
#include <petsc/finclude/petscksp.h>

module Continuation_mod

    use petscksp
    use PhysicalParameters_mod, only: ProblemParameters
    use MeshParameters_mod, only: MeshParameters
    use BoundaryParameters_mod, only: BoundaryParameters
    use SolutionVariables_mod, only: SolutionArraysType

    implicit none

    contains

    subroutine loopContinuation(tStart)

        use omp_lib, only: omp_get_wtime
        use ContinuationParameters_mod, only: Continuation_Method
        use PhysicalParameters_mod, only: Problem_Main, Problem_Inflow
        use FEMParameters_mod, only: FE_Inflow, FE_Main
        use MeshParameters_mod, only: Mesh_Inflow, Mesh_Main
        use BoundaryParameters_mod, only: Boundary_Inflow, Boundary_Main
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Cvar1, Increment
        use SolutionVariables_mod, only: Sol_Inflow, Sol_Main
        use PostProcess_mod, only: readSolution, writeTecplotFile
        use ContinuationArclength_mod, only: loopArclength
        use ContinuationTransient_mod, only: loopTransient
        use ContinuationStability_mod, only: loopStability

        implicit none

        PetscReal, intent(in) :: tStart

        PetscReal :: tEnd
        PetscInt :: ierr
        character(100) :: fn, str
        character(len=:), allocatable :: val1, ival

        call initializeSolution(Problem_Inflow, Sol_Inflow)
        call initializeSolution(Problem_Main, Sol_Main)

        write(str,'(f12.4)') Cvar1%p
        val1 = trim(adjustl(str))

        if (Rank == 0) then

            write(str,'(f12.4)') real(Increment,8)
            ival = trim(adjustl(str))

            fn = 'TEC_'//trim(adjustl(Cvar1%name))//'_'//val1//'_iter_'//ival &
                //'_'//trim(adjustl(Problem_Main%name))//'.PLT'
            fn = trim(adjustl(fn))

            call writeTecplotFile(Problem_Main, FE_Main, Mesh_Main, &
                Sol_Main%TL, Sol_Main%TL_proj, fn, 0.0d0)

            tEnd = omp_get_wtime()
            write(*,*)
            write(*, '(a,f6.3,a)') 'Time prior to continuation: ', tEnd - tStart, ' s'
            write(*,*)

            deallocate(ival)

        end if

        deallocate(val1)

        call deallocateExtraArrays(Mesh_Inflow, Boundary_Inflow, Sol_Inflow)
        call deallocateExtraArrays(Mesh_Main, Boundary_Main, Sol_Main)

        call setGlobalPointers()
        
        call MPI_Barrier(MPI_COMM_WORLD,ierr)

        select case (Continuation_Method)
        case ('Stability')
            call loopStability()
        case ('Arclength')
            call loopArclength()
        case ('Transient')
            call loopTransient()
        case default
            write(*,'(a)') 'Wrong Continuation_Method in loopContinuation!'
            stop
        end select

    end subroutine loopContinuation

    ! ----------------------------------------------------------------------

    subroutine initializeSolution(Problem, Sol)

        use PhysicalParameters_mod, only: Stress_Reform
        use ContinuationParameters_mod, only: Continuation_Method, LN0
        use MPIParameters_mod, only: Rank
        use Tools_mod, only: getStressComponent

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(SolutionArraysType), intent(inout) :: Sol

        PetscInt :: Neq, Neq_s, istart, i, j, k, Nex
        PetscReal :: value
        PetscBool :: check

        value = 0.1d0
        if (Continuation_Method == 'Transient') value = 0.0d0

        check = (Rank == 0)

        if (check) Sol%TLo(:,:) = value

        Sol%TLo_Rank(:,:) = value

        if (Stress_Reform == 'SQRT') then
            Neq = Problem%Neq
            Neq_s = Problem%Neq_s
            istart = Neq - Neq_s + 1
            do i = istart, size(Sol%TLo_Rank,2)
                call getStressComponent(i,istart,j,k)
                if (j == k) then
                    Sol%TLo_Rank(:,i) = 1.0d0
                    if (check) Sol%TLo(:,i) = 1.0d0
                end if
            end do
        end if

        if (check) then
            Sol%TLb(:,:) = Sol%TLo(:,:)
            Sol%TL(:,:)  = Sol%TLo(:,:)
            Sol%TLp(:,:) = Sol%TL(:,:)
        end if

        Sol%TLb_Rank(:,:) = Sol%TLo_Rank(:,:)
        Sol%TL_Rank(:,:)  = Sol%TLo_Rank(:,:)
        Sol%TLp_Rank(:,:) = Sol%TL_Rank(:,:)

        Sol%EUo(:) = value
        check = (Continuation_Method == 'Arclength') .and. (Sol%name == 'Main')
        if (check) then
            Nex = size(Sol%EUo)
            Sol%EUo(Nex) = LN0
        end if
        Sol%EUb(:) = Sol%EUo(:)
        Sol%EU(:)  = Sol%EUo(:)
        Sol%EUp(:) = Sol%EU(:)

    end subroutine initializeSolution

    ! ----------------------------------------------------------------------

    subroutine deallocateExtraArrays(Mesh, Boundary, Sol)

        use MPIParameters_mod, only: Rank

        implicit none

        type(MeshParameters), intent(inout) :: Mesh
        type(BoundaryParameters), dimension(:), intent(inout) :: Boundary
        type(SolutionArraysType), intent(inout) :: Sol

        PetscInt :: ibnd

        do ibnd = 1, size(Boundary)
            deallocate(Boundary(ibnd)%xi)
        end do

        if (Rank == 0) return

        if (allocated(Mesh%Connectivity)) deallocate(Mesh%Connectivity)
        if (allocated(Mesh%Xi_Mesh)) deallocate(Mesh%Xi_Mesh)
        if (allocated(Sol%TL)) deallocate(Sol%TL)
        if (allocated(Sol%TLo)) deallocate(Sol%TLo)
        if (allocated(Sol%TLb)) deallocate(Sol%TLb)
        if (allocated(Sol%TLp)) deallocate(Sol%TLp)
        if (allocated(Sol%TL_proj)) deallocate(Sol%TL_proj)
        if (allocated(Sol%TL_b)) deallocate(Sol%TL_b)
        if (allocated(Sol%TL_d)) deallocate(Sol%TL_d)

    end subroutine deallocateExtraArrays

    ! ----------------------------------------------------------------------

    subroutine setGlobalPointers()

        use PhysicalParameters_mod, only: Stress_Reform, &
            Problem_Main, Problem_Inflow, Problem_Main_Stability
        use Residuals_mod, only: momentumMain, momentumMainSQRT, &
            stressesBulkMain, stressesBulkMainSQRT, momentumMain_p, stressMain_p
        use ResidualsStability_mod, only: momentumStabilityT, momentumStabilitySQRT, &
            stressStabilityT, stressStabilitySQRT, &
            momentumStability_p, stressStability_p
        use ElementCalculations_mod, only: setBasisPointers
        use ConstitutiveModels_mod, only: setConEQ_T
        use ConstitutiveModelsSQRT_mod, only: setConEQ_SQRT, setPointersSQRT
        use ConstitutiveModelsStability_mod, only: setConEQ_T_Stability
        use ConstitutiveModelsSQRTStability_mod, only: setConEQ_SQRT_Stability, &
            setPointersSQRT_Stability

        implicit none

        select case (Stress_Reform)
        case ('SQRT')
            call setConEQ_SQRT()
            call setConEQ_SQRT_Stability()
            call setPointersSQRT(Problem_Main%Ndim)
            call setPointersSQRT_Stability(Problem_Main_Stability%Ndim)
            momentumMain_p => momentumMainSQRT
            stressMain_p => stressesBulkMainSQRT
            momentumStability_p => momentumStabilitySQRT
            stressStability_p => stressStabilitySQRT
        case default
            call setConEQ_T()
            call setConEQ_T_Stability()
            momentumMain_p => momentumMain
            stressMain_p => stressesBulkMain
            momentumStability_p => momentumStabilityT
            stressStability_p => stressStabilityT
        end select

        call setBasisPointers(Problem_Main%Ndim, Problem_Inflow%Ndim)

    end subroutine setGlobalPointers

end module Continuation_mod
