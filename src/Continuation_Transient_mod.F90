
#include <petsc/finclude/petscksp.h>

module ContinuationTransient_mod

    use petscksp
    use SolutionVariables_mod, only: SolutionArraysType

    implicit none

    contains

    subroutine loopTransient()

        use omp_lib, only: omp_get_wtime
        use Tools_mod, only: writeElapsedTime
        use PhysicalParameters_mod, only: Total_time, Time_f
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Increment, Cvar2
        use Petsc_mod, only: writeMemoryUsage
        use Physical_mod, only: setCVar
        use ContinuationStability_mod, only: makeVar2Directories, readVar1Solution

        implicit none

        PetscReal :: tstart, tEnd

        tstart = omp_get_wtime()
        
        Increment = 0

        call setCVar('',Cvar2)
        if (Rank == 0) then
            call makeVar2Directories(Cvar2%name, Cvar2%p)
        end if

        call readVar1Solution()

        loop_continuation:do

            if (Total_time > Time_f) exit

            call continuationStepTransient()

            tEnd = omp_get_wtime()
            if (Rank == 0) call writeElapsedTime(tStart, tEnd, 'Total elapsed time')

            call writeMemoryUsage()

        end do loop_continuation

    end subroutine loopTransient

    !---------------------------------------------------------------------

    subroutine continuationStepTransient()

        use omp_lib, only: omp_get_wtime
        use Tools_mod, only: writeElapsedTime
        use ContinuationParameters_mod, only: adjust_dt
        use PhysicalParameters_mod, only: Total_time, Problem_Inflow, Problem_Main, Time0
        use FEMParameters_mod, only: FE_Inflow, FE_Main
        use MeshParameters_mod, only: Mesh_Inflow, Mesh_Main
        use BoundaryParameters_mod, only: Boundary_Inflow, Boundary_Main
        use GaussParameters_mod, only: GaussInt_Inflow, GaussInt_Main
        use MPIParameters_mod, only: Rank
        use NewtonRaphsonParameters_mod, only: LS_maxit
        use ContinuationVariables_mod, only: Increment, dt
        use ElementVariables_mod, only: NodeArrays_Inflow, NodeArrays_Main
        use SolutionVariables_mod, only: Sol_Inflow, Sol_Main
        use LinearSystemVariables_mod, only: LS_Inflow, LS_Main, LS_Main_proj, LS_Inflow_proj
        use NR_mod, only: loopNewtonRaphson, perturbSolution
        use Petsc_mod, only: reCreatePetscObjects
        use PostProcess_mod, only: postProcess
        use PostProcessStability_mod, only: readLeadingEigenvectorSolution

        implicit none

        PetscReal, save :: twrite0_Main = Time0, twrite0_Inflow = Time0
        PetscBool :: reCreate
        PetscReal :: tStart, tEnd
        PetscErrorCode :: ierr

        tStart = omp_get_wtime()

        Increment = Increment + 1

        Total_time = Total_time + dt

        if (Rank == 0) call printStepInfoTransient()

        if (Problem_Inflow%Ndim > 0) then
            call loopNewtonRaphson(Problem_Inflow, FE_Inflow, Mesh_Inflow, &
                Boundary_Inflow, GaussInt_Inflow, Sol_Inflow, NodeArrays_Inflow, LS_Inflow)

            call postProcess(Problem_Inflow, FE_Inflow, Mesh_Inflow, Boundary_Inflow, &
                GaussInt_Inflow, Sol_Inflow, NodeArrays_Inflow, LS_Inflow_proj, twrite0_Inflow)
        end if

        ! call readLeadingEigenvectorSolution(Problem_Main, FE_Main, Mesh_Main%Connectivity, Sol_Main)
        ! call perturbSolution(Sol_Main)

        reCreate = (LS_Main%name == 'Iterative')
        ! reCreate = reCreate .and. (Total_time > 20.0d0)
        reCreate = reCreate .and. (LS_Main%reCreate)
        reCreate = reCreate .or. (LS_Main%maxit >= LS_maxit/2)
        if (reCreate) then
            call reCreatePetscObjects(Problem_Main, FE_Main, Mesh_Main, Boundary_Main, LS_Main)
        end if

        call loopNewtonRaphson(Problem_Main, FE_Main, Mesh_Main, &
            Boundary_Main, GaussInt_Main, Sol_Main, NodeArrays_Main, LS_Main)

        call postProcess(Problem_Main, FE_Main, Mesh_Main, &
            Boundary_Main, GaussInt_Main, Sol_Main, NodeArrays_Main, LS_Main_proj, twrite0_Main)

        if (adjust_dt) then
            if (Rank == 0) call setTimeStep(Problem_Main%Neq, Sol_Main)
            call MPI_Bcast(dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
        end if

        if (Problem_Inflow%Ndim > 0) then
            call updateSolutionTransient(Sol_Inflow)
        end if
        call updateSolutionTransient(Sol_Main)

        tEnd = omp_get_wtime()
        if (Rank == 0) call writeElapsedTime(tStart, tEnd, 'Continuation step time')

    end subroutine continuationStepTransient

    !---------------------------------------------------------------------

    subroutine printStepInfoTransient()

        use ContinuationVariables_mod, only: Increment, Cvar1
        use SolutionVariables_mod, only: Sol_Main

        implicit none

        write(*,*)''
        write(*,'(a)')'======================================================================'
        write(*,60) Increment, trim(adjustl(Cvar1%name)), Cvar1%p, Sol_Main%EU(:)
        write(*,'(a)')'======================================================================'
        write(*,*)''

        60 format('STEP = ', i8, 5x, a, ' = ' f12.4, 5x, 'EU = ', *(f10.3))

    end subroutine printStepInfoTransient

    !---------------------------------------------------------------------

    Subroutine setTimeStep(Neq, Sol)

        use ContinuationVariables_mod, only: dt

        implicit none

        PetscInt, intent(in) :: Neq
        type(SolutionArraysType), intent(in) :: Sol

        PetscInt :: i
        PetscReal :: dif

        dif = 0.0d0
        do i = 1, Neq
            dif = dif + dot_product((Sol%TL(:,i)-Sol%TLp(:,i)), &
                                    (Sol%TL(:,i)-Sol%TLp(:,i)))
        end do
        dif = dif + dot_product((Sol%EU(:)-Sol%EUp(:)), &
                                (Sol%EU(:)-Sol%EUp(:)))
        dif = sqrt(dif)
        ! dif = dif/(dt**2)

        if (dif < 10.0d0) dt = 2.0d0*dt
        if (dif > 1000.0d0) dt = dt/2.0d0

        if (dt < 0.01d0) dt = 0.01d0
        if (dt > 0.05d0) dt = 0.05d0

        write(*,'(2(a,es12.4,3x))') 'Solution difference: ', dif, 'New time step: ', dt
        write(*, '(a)') '--------------------------------------------------------------'

    end subroutine setTimeStep

    !--------------------------------------------------------------------------

    subroutine updateSolutionTransient(Sol)

        use Tools_mod, only: extrapolationLagrange
        use PhysicalParameters_mod, only: Total_time, Time0
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: dt, Increment

        implicit none

        type(SolutionArraysType), intent(inout) :: Sol

        PetscReal :: Lb, Lo, L, max_dif
        PetscBool :: extrapolate
        PetscErrorCode :: ierr

        extrapolate = (Increment >= 3) .or. (Time0 > 0.0d0)

        if (Rank == 0) then
            max_dif = maxval( abs(Sol%TL(:,:)-Sol%TLo(:,:)) )
        end if

        call MPI_Bcast(max_dif, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)

        if ( max_dif < 1.0d-8 ) extrapolate = .false.

        if (Rank == 0) then
            write(*,'(2a,es12.4,a,l)') trim(Sol%name), ' Max solution difference: ', max_dif, &
                                    ', Extrapolate = ', extrapolate
            write(*, '(a)') '--------------------------------------------------------------'
        end if

        Lb = 0.0d0 ; Lo = 0.0d0 ; L = 1.0d0
        if (extrapolate) then
            call extrapolationLagrange(Total_time+dt, Total_time-(dt+dt), &
                Total_time-dt, Total_time, Lb, Lo, L)
        end if

        if (Rank == 0) then
            Sol%TLp(:,:) = Lb*Sol%TLb(:,:) + Lo*Sol%TLo(:,:) + L*Sol%TL(:,:)
            Sol%TLb(:,:) = Sol%TLo(:,:)
            Sol%TLo(:,:) = Sol%TL(:,:)
            Sol%TL(:,:) = Sol%TLp(:,:)
        end if

        Sol%TLp_Rank(:,:) = Lb*Sol%TLb_Rank(:,:) + Lo*Sol%TLo_Rank(:,:) &
                                + L*Sol%TL_Rank(:,:)
        Sol%TLb_Rank(:,:) = Sol%TLo_Rank(:,:)
        Sol%TLo_Rank(:,:) = Sol%TL_Rank(:,:)
        Sol%TL_Rank(:,:)  = Sol%TLp_Rank(:,:)

        Sol%EUp(:) = Lb*Sol%EUb(:) + Lo*Sol%EUo(:) + L*Sol%EU(:)
        Sol%EUb(:) = Sol%EUo(:)
        Sol%EUo(:) = Sol%EU(:)
        Sol%EU(:) = Sol%EUp(:)

    end subroutine updateSolutionTransient

end module ContinuationTransient_mod
