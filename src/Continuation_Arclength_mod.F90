
#include <petsc/finclude/petscksp.h>

module ContinuationArclength_mod

    use petscksp
    use PhysicalParameters_mod, only: ProblemParameters
    use FEMParameters_mod, only: FEParameters
    use SolutionVariables_mod, only: SolutionArraysType
    use LinearSystemVariables_mod, only: LinearSystemType

    implicit none

    contains

    subroutine initializeVariablesArclength()

        use ContinuationParameters_mod, only: dS0
        use ContinuationVariables_mod, only: dSo, Increment

        implicit none

        Increment = 0

        dSo = dS0

    end subroutine initializeVariablesArclength

    !---------------------------------------------------------------------

    subroutine loopArclength()

        use omp_lib, only: omp_get_wtime
        use Tools_mod, only: writeElapsedTime
        use ContinuationParameters_mod, only: adjust_dt
        use PhysicalParameters_mod, only: dt0
        use MPIParameters_mod, only: Rank
        use Petsc_mod, only: writeMemoryUsage
        use ContinuationVariables_mod, only: Cvar2, dt
        use Physical_mod, only: setCVar
        use ContinuationStability_mod, only: makeVar2Directories, readVar1Solution

        implicit none

        PetscReal :: tstart, tEnd

        tstart = omp_get_wtime()

        ! if (adjust_dt) then
            ! dt = dt0
        ! end if

        call initializeVariablesArclength()

        call setCVar('',Cvar2)
        if (Rank == 0) then
            call makeVar2Directories(Cvar2%name, Cvar2%p)
        end if

        call readVar1Solution()

        loop_continuation:do

            call continuationStepArclength()

            tEnd = omp_get_wtime()
            if (Rank == 0) call writeElapsedTime(tStart, tEnd, 'Total elapsed time')

            call writeMemoryUsage()

        end do loop_continuation

    end subroutine loopArclength

    !---------------------------------------------------------------------

    subroutine continuationStepArclength()

        use omp_lib, only: omp_get_wtime
        use Tools_mod, only: writeElapsedTime
        use PhysicalParameters_mod, only: Problem_Inflow, Problem_Main
        use FEMParameters_mod, only: FE_Main, FE_Inflow
        use MeshParameters_mod, only: Mesh_Main, Mesh_Inflow
        use BoundaryParameters_mod, only: Boundary_Main, Boundary_Inflow
        use GaussParameters_mod, only: GaussInt_Main, GaussInt_Inflow
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Increment
        use ElementVariables_mod, only: NodeArrays_Main, NodeArrays_Inflow
        use SolutionVariables_mod, only: Sol_Inflow, Sol_Main
        use LinearSystemVariables_mod, only: LS_Main, LS_Main_proj, LS_Inflow_proj
        use NR_mod, only: loopNewtonRaphson
        use PostProcess_mod, only: postProcess
        use ContinuationStability_mod, only: updateBaseSolution

        implicit none

        PetscReal, save :: twrite0_Main = 0.0d0, twrite0_Inflow = 0.0d0
        PetscReal :: tstart, tEnd, dTLodS_Norm

        tstart = omp_get_wtime()

        Increment = Increment + 1

        if (Rank == 0) call printStepInfoArclength()

        call loopNewtonRaphson(Problem_Main, FE_Main, Mesh_Main, &
            Boundary_Main, GaussInt_Main, Sol_Main, NodeArrays_Main, LS_Main)

        call postProcess(Problem_Main, FE_Main, Mesh_Main, &
            Boundary_Main, GaussInt_Main, Sol_Main, NodeArrays_Main, LS_Main_proj, twrite0_Main)

        if (Increment > 1) then
            call setDirectionVectors(Problem_Main, FE_Main, &
                Mesh_Main%Connectivity_Rank, Sol_Main, LS_Main, dTLodS_Norm)
            call adjustContinuationStep(Problem_Main%Neq, Sol_Main)
        end if

        if (Problem_Inflow%Ndim > 0) then
            call postProcess(Problem_Inflow, FE_Inflow, Mesh_Inflow, &
                Boundary_Inflow, GaussInt_Inflow, Sol_Inflow, NodeArrays_Inflow, LS_Inflow_proj, twrite0_Inflow)
            call updateBaseSolution(Sol_Inflow)
        end if
        call updateSolutionArclength(dTLodS_Norm, Sol_Main)

        tEnd = omp_get_wtime()
        if (Rank == 0) call writeElapsedTime(tStart, tEnd, 'Continuation step time')

    end subroutine continuationStepArclength

    !---------------------------------------------------------------------

    subroutine printStepInfoArclength()

        use ContinuationVariables_mod, only: Increment, dSo, Cvar1
        use SolutionVariables_mod, only: Sol_Main

        implicit none

        write(*,*)''
        write(*,'(A)')'======================================================================'
        write(*,60) Increment, trim(adjustl(Cvar1%name)), Cvar1%p, dSo, Sol_Main%EU(:)
        write(*,'(A)')'======================================================================'
        write(*,*)''

        60 format('STEP = ',i5,3x,a,' = 'f10.3,3x,'dSo = ',f10.3,3x,'EU = ',*(f10.3))

    end subroutine printStepInfoArclength

    !----------------------------------------------------------------------

    subroutine setDirectionVectors(Problem, FE, Connectivity_Rank, Sol, LS, dTLodS_Norm)

        use Tools_mod, only: getRows, getRankNode
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Increment
        use SolutionVariables_mod, only: dTLodS, dEUodS, dTLodS_Rank

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:,:), intent(in) :: Connectivity_Rank
        type(SolutionArraysType), intent(in) :: Sol
        type(LinearSystemType), intent(inout) :: LS
        real(8), intent(out) :: dTLodS_Norm

        PetscInt :: i, inod, iel_rank, rnod, ieq, Nunknowns, Neq, Nex
        PetscReal :: dS_loc
        PetscScalar :: value
        PetscInt, dimension(1) :: gnod
        PetscInt, dimension(Problem%Neq) :: ieqs, irows_eqs
        Vec :: s_f
        PetscScalar, pointer :: s_f_pointer(:)
        VecScatter :: ctx
        PetscErrorCode :: ierr

        dTLodS_Norm = 0.0d0

        Nunknowns = Problem%Nunknowns
        Neq = Problem%Neq
        Nex = Problem%Nex

        if (Increment == 2) then

            if (Rank == 0) then

                dS_loc = 0.0d0
                do ieq = 1, Neq
                    dS_loc = dS_loc + dot_product( (Sol%TL(:,ieq)-Sol%TLo(:,ieq)), &
                                                    (Sol%TL(:,ieq)-Sol%TLo(:,ieq)) )
                end do

                dS_loc = dS_loc + dot_product( (Sol%EU-Sol%EUo),&
                                                (Sol%EU-Sol%EUo) )
                dS_loc = sqrt(dS_loc)
                
                dTLodS(:,:) = (Sol%TL(:,:)-Sol%TLo(:,:))/dS_loc

            end if

            !Send dS_loc to all processes
            call MPI_Bcast(dS_loc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)

            !Extra unknowns
            dEUodS(:) = (Sol%EU(:)-Sol%EUo(:))/dS_loc

            dTLodS_Rank(:,:) = (Sol%TL_rank(:,:)-Sol%TLo_rank(:,:))/dS_loc

            return

        end if

        ! Be_f = 0.0d0 ; Be_f(Nex) = 1.0d0
        
        ! call VecZeroEntries(LS%s_f, ierr)
        ! call VecZeroEntries(LS%b_f, ierr)

        value = 1.0d0
        call VecSetValues(LS%b_f, 1, Nunknowns+Nex-1, value, INSERT_VALUES, ierr)
        call VecAssemblyBegin(LS%b_f, ierr)
        call VecAssemblyEnd(LS%b_f, ierr)

        call MatAssemblyBegin(LS%A_f, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(LS%A_f, MAT_FINAL_ASSEMBLY, ierr)

        ! ! if (Flag_NR == 'MNR') then
        ! !     call KSPSetReusePreconditioner(LS%ksp, PETSC_TRUE, ierr)
        ! ! else
        !     call KSPSetReusePreconditioner(LS%ksp, PETSC_FALSE, ierr)
        ! ! end if

        call KSPSolve(LS%ksp, LS%b_f, LS%s_f, ierr)

        ! call KSPConvergedReasonView(LS%ksp, Petsc_VIEWER_STDOUT_WORLD, ierr)

        call VecScatterCreateToAll(LS%s_f, ctx, s_f, ierr)
        !INSERT_VALUES: any location non-scattered retains its old value
        !SCATTER_fORWARD: scatter from LS%s_f to s_f
        call VecScatterBegin(ctx, LS%s_f, s_f, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterEnd(ctx, LS%s_f, s_f, INSERT_VALUES, SCATTER_FORWARD, ierr)

        call VecGetArrayReadF90(s_f, s_f_pointer, ierr)

        !Norm 2 of LS%s_f
        call VecNorm(LS%s_f, NORM_2, dTLodS_Norm, ierr)

        if (Rank == 0) write(*,'(a, f12.6)') 'Norm = ', dTLodS_Norm

        !Normalize direction vectors
        if (Rank == 0) then
            i = 0
            do inod = 1, FE%Nnodes
                do ieq = 1, Neq
                    i = i + 1
                    dTLodS(inod,ieq) = s_f_pointer(i)/dTLodS_Norm
                end do
            end do
        end if

        dEUodS(:) = s_f_pointer(Nunknowns+1:Nunknowns+Nex)/dTLodS_Norm

        ieqs = [ (ieq, ieq = 1, Neq) ]
        
        do iel_rank = 1, FE%Nel_Rank
            do inod = 1, FE%Nbf
                gnod(1) = Connectivity_Rank(iel_rank,inod)
                rnod = getRankNode(iel_rank,inod, FE%Nbf)
                irows_eqs(:) = getRows(gnod(:),ieqs(:),Neq)
                dTLodS_Rank(rnod,:) = s_f_pointer(irows_eqs(:))/dTLodS_Norm
            end do
        end do

        call VecRestoreArrayReadF90(s_f, s_f_pointer, ierr)

        call VecScatterDestroy(ctx, ierr)
        call VecDestroy(s_f, ierr)

    end subroutine setDirectionVectors

    !----------------------------------------------------------------------

    subroutine adjustContinuationStep(Neq, Sol)

        use ContinuationParameters_mod, only: eps_pr, dS_max
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod,  only: Increment, dSo

        implicit none

        PetscInt, intent(in) :: Neq
        type(SolutionArraysType), intent(in) :: Sol

        PetscReal, dimension(Neq+1,3) :: erind
        PetscInt  :: i, idim_er
        PetscReal :: err_pr, trns
        PetscErrorCode :: ierr

        idim_er = Neq + 1

        if (Rank == 0) then
            
            do i = 1, Neq
                erind(i,1) = dot_product( Sol%TL(:,i)-Sol%TLp(:,i), Sol%TL(:,i)-Sol%TLp(:,i) )
                erind(i,2) = dot_product(Sol%TL(:,i), Sol%TL(:,i))
                erind(i,3) = sqrt(erind(i,1)/erind(i,2))
            end do

            erind(idim_er,1) = dot_product(Sol%EU(:)-Sol%EUp(:), Sol%EU(:)-Sol%EUp(:))
            erind(idim_er,2) = dot_product(Sol%EU(:), Sol%EU(:))
            erind(idim_er,3) = sqrt(erind(idim_er,1)/erind(idim_er,2))

            err_pr = maxval(erind(1:idim_er,3))
            trns = sqrt(eps_pr/err_pr)
            ! trns = sqrt(err_pr/eps_pr)
            ! trns = err_pr/eps_pr
            write(*,'(a, es12.4, a, f12.6)') 'Error in prediction:', err_pr, ', trns:', trns
            write(*,'(a)') '--------------------------------------------------------------'

        end if

        !Send trns to all processes
        call MPI_Bcast(trns, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)

        !MAX variation
        if (trns > 2.0d0) trns = 2.0d0

        ! !MIN variation
        ! if (trns < 0.5d0) trns = 0.5d0

        ! if (Increment <= 5) trns = 1.0d0

        dSo = min(dSo*trns, dS_max)
        ! dSo = dS0

    end subroutine adjustContinuationStep

    !--------------------------------------------------------------------------

    subroutine updateSolutionArclength(dTLodS_Norm, Sol)

        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: dSo, Increment
        use SolutionVariables_mod, only: dTLodS_Rank, dEUodS, dTLodS

        implicit none

        PetscReal, intent(in) :: dTLodS_Norm
        type(SolutionArraysType), intent(inout) :: Sol

        if (Rank == 0) then
            if (Increment == 1) then
                Sol%TLp(:,:) = Sol%TL(:,:)
            else
                Sol%TLp(:,:) = Sol%TL(:,:) + dTLodS(:,:)*dSo!*dTLodS_Norm
            end if
            Sol%TLb(:,:) = Sol%TLo(:,:)
            Sol%TLo(:,:) = Sol%TL(:,:)
            Sol%TL(:,:) = Sol%TLp(:,:)
        end if

        if (Increment == 1) then
            Sol%TLp_Rank(:,:) = Sol%TL_Rank(:,:)
            Sol%EUp(:) = Sol%EU(:)
        else
            Sol%TLp_Rank(:,:) = Sol%TL_Rank(:,:) + dTLodS_Rank(:,:)*dSo!*dTLodS_Norm
            Sol%EUp(:) = Sol%EU(:) + dEUodS(:)*dSo
        end if
        
        Sol%TLb_Rank(:,:) = Sol%TLo_Rank(:,:)
        Sol%TLo_Rank(:,:) = Sol%TL_Rank(:,:)
        Sol%TL_Rank(:,:) = Sol%TLp_Rank(:,:)

        Sol%EUb(:) = Sol%EUo(:)
        Sol%EUo(:) = Sol%EU(:)
        Sol%EU(:) = Sol%EUp(:)

    end subroutine updateSolutionArclength

end module ContinuationArclength_mod
