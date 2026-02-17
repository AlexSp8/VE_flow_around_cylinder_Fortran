
#include <petsc/finclude/petscksp.h>

module NR_mod

    use petscksp
    use PhysicalParameters_mod, only: ProblemParameters
    use FEMParameters_mod, only: FEParameters
    use MeshParameters_mod, only: MeshParameters
    use BoundaryParameters_mod, only: BoundaryParameters
    use GaussParameters_mod, only: GaussIntegration
    use ElementVariables_mod, only: NodeArrays
    use SolutionVariables_mod, only: SolutionArraysType
    use LinearSystemVariables_mod, only: LinearSystemType
    use NewtonRaphsonVariables_mod, only: NewtonType

    implicit none

    contains

    subroutine loopNewtonRaphson(Problem, FE, Mesh, Boundary, GaussInt, Sol, Elem, LS)

        use omp_lib, only: omp_get_wtime
        use Tools_mod, only: writeElapsedTime
        use ContinuationParameters_mod, only: Continuation_Method
        use MPIParameters_mod, only: Rank
        use NewtonRaphsonParameters_mod, only: Cor_Norm_eps, MNR_eps
        use ContinuationVariables_mod, only: Increment
        use Petsc_mod, only: initializeLSObjects, solveLinearSystem, &
            calculateMatrixDeterminant
        ! use PostProcess_mod, only: writeTecplotFile

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(inout) :: Sol
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS

        type(NewtonType) :: NR
        PetscInt, save :: iter0_Main = 0
        PetscReal, save :: cor_norm0_Main = 1.0d0
        PetscReal :: tStart, tEnd, tStart2, tStart3
        PetscBool :: check, output, emergency
        PetscErrorCode :: ierr
        PetscViewer :: viewer

        tStart = omp_get_wtime()

        NR%name = Sol%name

        output = (Rank == 0)

        if (output .and. NR%name == 'Main') then
            write(*, '(a)') '======================================================================'
            write(*, '(33x,a)') 'Main'
            write(*, '(a)') '======================================================================'
            open(20, file='Results/Base/Info/NR_data.dat', position='append')
            write(20,'(a,i0,a)') '------------ i = ', Increment, ' ---------------'
        end if

        if (output .and. NR%name == 'Inflow') then
            write(*, '(a)') '======================================================================'
            write(*, '(32x,a)') 'Inflow'
            write(*, '(a)') '======================================================================'
        end if

        NR%xF = 1.0d0
        call initializeNewtonVariables(NR)
        if (NR%name == 'Main') then
            NR%Flag_NR =  checkMNR(iter0_Main, cor_norm0_Main)
            if (LS%maxit == 0) then
                NR%Flag_NR = 'NRP'
            end if
        end if

        LOOP_NewtonRaphson:do

            check = NR%Cor_Norm_New_f < Cor_Norm_eps
            check = check .and. NR%Res_Norm < Cor_Norm_eps!/100.0d0
            if (check) then
                if (NR%name == 'Main') then
                    iter0_Main = NR%Iter_f
                end if
                exit
            end if

            tStart2 = omp_get_wtime()

            NR%Iter_f = NR%Iter_f + 1

            ! if (NR%Iter_f >= 10) then
            if (mod(NR%Iter_f,10) == 0) then
                NR%Flag_NR = 'NRP'
            end if

            call initializeLSObjects(NR%Flag_NR, LS)

            call assemblePetscObjects(output, Problem, FE, Mesh, Boundary, &
                GaussInt, Sol, NR, Elem, LS)

            call VecNormBegin(LS%b_f, NORM_2, NR%Res_Norm, ierr)
            call VecNormEnd(LS%b_f, NORM_2, NR%Res_Norm, ierr)

            tStart3 = omp_get_wtime()
            call solveLinearSystem(NR, LS)
            tEnd = omp_get_wtime()
            if (output) write(*,'(a,f8.3,a)') 'LS Solve time: ', tEnd - tStart3, ' s'

            check = (LS%name == 'Direct' .and. Sol%name == 'Main')
            check = check .and. (NR%Flag_NR == 'NRP')
            if (check) then

                call calculateMatrixDeterminant(LS%A_f_factored, output, NR%Iter_f)

                if (FE%Ndim == 3) then
                    call PetscLogDefaultBegin(ierr)
                    call PetscViewerASCIIOpen(PETSC_COMM_WORLD,&
                            'Results/Base/Info/log_info.dat',viewer,ierr)
                    call PetscLogView(viewer,ierr)
                    call PetscViewerDestroy(viewer, ierr)
                end if

            end if

            call VecNormBegin(LS%s_f, NORM_2, NR%Cor_Norm_New_f, ierr)
            call VecNormEnd(LS%s_f, NORM_2, NR%Cor_Norm_New_f, ierr)

            if (NR%name == 'Main') then

                if (output) write(20,'(i5,2x,a,2es12.4,i5)') NR%Iter_f, NR%Flag_NR, &
                    NR%Res_norm, NR%Cor_Norm_New_f, LS%iter

                if (NR%Iter_f == 1) cor_norm0_Main = NR%Cor_Norm_New_f

            end if

            call checkNewtonConvergence(output, NR, emergency)

            if (emergency) then
                call emergencyNewtonRaphson(output, NR, Sol)
                cycle LOOP_NewtonRaphson
            end if

            call updateNewton(Problem, FE, Mesh%Connectivity_Rank, NR%xF, LS%s_f, Sol)

            if (output) then
                write(*,'(a,*(f14.6))') 'EU = ', Sol%EU(:)
                ! call writeTecplotFile(Problem, FE, Mesh, &
                !     Sol%TL, Sol%TL_proj, fn, real(NR%Iter_f,8))
                ! pause
            end if

            tEnd = omp_get_wtime()
            if (output) call writeElapsedTime(tStart2, tEnd, 'Newton time')

        end do LOOP_NewtonRaphson

        if (output .and. NR%name == 'Main') close(20, status='keep')

        LS%reCreate = (cor_norm0_Main < MNR_eps/1.0d1)

        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        
        tEnd = omp_get_wtime()
        if (output) call writeElapsedTime(tStart, tEnd, 'Solution time')

    end subroutine loopNewtonRaphson

    !--------------------------------------------------------------------------

    subroutine loopNewtonRaphsonInflow(Problem, FE, Mesh, Boundary, GaussInt, Sol, Elem, LS)

        use omp_lib, only: omp_get_wtime
        use Tools_mod, only: writeElapsedTime
        use ContinuationParameters_mod, only: Continuation_Method
        use MPIParameters_mod, only: Rank
        use NewtonRaphsonParameters_mod, only: Cor_Norm_eps
        use Petsc_mod, only: initializeLSObjects, solveLinearSystem

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(inout) :: Sol
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS

        type(NewtonType) :: NR
        PetscReal :: tStart, tEnd, tStart2, tStart3
        PetscBool :: check, output, emergency
        PetscErrorCode :: ierr

        tStart = omp_get_wtime()

        NR%name = Sol%name

        output = (Rank == 0)
        check = (Continuation_Method == 'Arclength')
        if (check) output = .false.

        if (output) then
            write(*, '(a)') '======================================================================'
            write(*, '(32x,a)') 'Inflow'
            write(*, '(a)') '======================================================================'
        end if

        NR%xF = 1.0d0
        call initializeNewtonVariables(NR)

        LOOP_NewtonRaphson:do

            check = NR%Cor_Norm_New_f < Cor_Norm_eps
            check = check .and. NR%Res_Norm < Cor_Norm_eps!/100.0d0
            if (check) then
                exit
            end if

            tStart2 = omp_get_wtime()

            NR%Iter_f = NR%Iter_f + 1

            ! if (NR%Iter_f >= 10) then
            if (mod(NR%Iter_f,10) == 0) then
                NR%Flag_NR = 'NRP'
            end if

            call initializeLSObjects(NR%Flag_NR, LS)

            call assemblePetscObjects(output, Problem, FE, Mesh, Boundary, &
                GaussInt, Sol, NR, Elem, LS)

            call VecNormBegin(LS%b_f, NORM_2, NR%Res_Norm, ierr)
            call VecNormEnd(LS%b_f, NORM_2, NR%Res_Norm, ierr)

            tStart3 = omp_get_wtime()
            call solveLinearSystem(NR, LS)
            tEnd = omp_get_wtime()
            if (output) write(*,'(a,f8.3,a)') 'LS Solve time: ', tEnd - tStart3, ' s'

            call VecNormBegin(LS%s_f, NORM_2, NR%Cor_Norm_New_f, ierr)
            call VecNormEnd(LS%s_f, NORM_2, NR%Cor_Norm_New_f, ierr)

            call checkNewtonConvergence(output, NR, emergency)

            if (emergency) then
                call emergencyNewtonRaphson(output, NR, Sol)
                cycle LOOP_NewtonRaphson
            end if

            call updateNewton(Problem, FE, Mesh%Connectivity_Rank, NR%xF, LS%s_f, Sol)

            if (output) then
                write(*,'(a,*(f14.6))') 'EU = ', Sol%EU(:)
                ! call writeTecplotFile(Problem, FE, Mesh, &
                !     Sol%TL, Sol%TL_proj, fn, real(NR%Iter_f,8))
                ! pause
            end if

            tEnd = omp_get_wtime()
            if (output) call writeElapsedTime(tStart2, tEnd, 'Newton time')

        end do LOOP_NewtonRaphson

        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        
        tEnd = omp_get_wtime()
        if (output) call writeElapsedTime(tStart, tEnd, 'Solution time')

    end subroutine loopNewtonRaphsonInflow

    !--------------------------------------------------------------------------

    subroutine assemblePetscObjects(output, Problem, FE, Mesh, Boundary, &
        GaussInt, Sol, NR, Elem, LS)

        use omp_lib, only: omp_get_wtime
        use ContinuationParameters_mod, only: Continuation_Method
        use ElementCalculations_mod, only: basis, basisInflow_p, basisMain_p
        use Equations_mod, only: setEquationsBulk, setEquationsEdge, setEquationsFace
        use Extra_Equations_mod, only: setExtraEquationsBulk, setExtraEquationsEdge, &
            setExtraEquationsFace, arclength
        use Dirichlet_mod, only: clearEntries, applyDirichletBC, applyPeriodicBC, &
            applyDirichletBCFullyDeveloped
        use Petsc_mod, only: finalAssembly, flushAssembly

        use PhysicalParameters_mod, only: Problem_Inflow
        use FEMParameters_mod, only: FE_Inflow
        use MeshParameters_mod, only: Mesh_Inflow
        use BoundaryParameters_mod, only: Boundary_Inflow
        use GaussParameters_mod, only: GaussInt_Inflow
        use SolutionVariables_mod, only: Sol_inflow
        use ElementVariables_mod, only: NodeArrays_Inflow
        use LinearSystemVariables_mod, only: LS_Inflow

        implicit none
    
        PetscBool, intent(in) :: output
        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        type(NewtonType), intent(in) :: NR
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS

        PetscReal :: tStart, tEnd
        PetscInt :: ibnd, inex, i, ierr
        PetscBool :: check
        procedure(basis), pointer :: basis_p

        tStart = omp_get_wtime()

        if (Sol%name == 'Main') then
            basis_p => basisMain_p
        end if

        if (Sol%name == 'Inflow') then
            basis_p => basisInflow_p
        end if

        !Bulk
        call setEquationsBulk(Problem, FE, Mesh, Boundary, GaussInt, Sol, &
            basis_p, NR%Flag_NR, Elem, LS)
        call flushAssembly(NR%Flag_NR, LS)

        tEnd = omp_get_wtime()
        if (output) write(*,'(a,f8.3,a)', advance='no') 'Bulk EQs time: ', tEnd - tStart, ' s'

        do inex = 1, Problem%Nex_b
            call setExtraEquationsBulk(Problem, FE, Mesh, &
                GaussInt, Sol, basis_p, NR%Flag_NR, inex, Elem, LS)
        end do
        call flushAssembly(NR%Flag_NR, LS)

        !Boundary
        do ibnd = 1, size(Boundary)
            ! call finalAssembly(NR%Flag_NR, LS)
            ! call clearEntries(Problem%Neq, Boundary(ibnd), NR%Flag_NR, LS)
            call setEquationsEdge(Problem, FE, Mesh, Boundary(ibnd), Boundary, &
                GaussInt, Sol, basis_p, NR%Flag_NR, Elem, LS)
            call setEquationsFace(Problem, FE, Mesh, Boundary(ibnd), Boundary, &
                GaussInt, Sol, basis_p, NR%Flag_NR, Elem, LS)
            call setExtraEquationsEdge(Problem, FE, Mesh, Boundary(ibnd), &
                GaussInt, Sol, basis_p, NR%Flag_NR, Elem, LS)
            call setExtraEquationsFace(Problem, FE, Mesh, Boundary(ibnd), &
                GaussInt, Sol, basis_p, NR%Flag_NR, Elem, LS)
        end do
        call flushAssembly(NR%Flag_NR, LS)

        !Arclength
        check = (Continuation_Method == 'Arclength')
        check = check .and. (Sol%name == 'Main')
        if (check) then

            call arclength(Problem, FE%Nnodes, Sol, NR%Flag_NR, LS)
            call flushAssembly(NR%Flag_NR, LS)

            if (Problem_Inflow%Ndim > 0) then
                call loopNewtonRaphsonInflow(Problem_Inflow, FE_Inflow, Mesh_Inflow, &
                    Boundary_Inflow, GaussInt_Inflow, Sol_Inflow, NodeArrays_Inflow, LS_Inflow)
            end if

        end if

        !PBCs
        do i = 1, size(Problem%PeriodicBNDs,1)
            ibnd = Problem%PeriodicBNDs(i,1)
            call finalAssembly(NR%Flag_NR, LS)
            call applyPeriodicBC(Problem%Neq, Boundary(ibnd), Sol, NR%Flag_NR, LS)
        end do

        call finalAssembly(NR%Flag_NR, LS)

        !BCs
        check = (Sol%name == 'Main')
        check = check .and. ( len(LS_Inflow%name) > 0)
        if (check) then
            call applyDirichletBCFullyDeveloped(NR%Flag_NR, LS)
        end if

        do ibnd = 1, size(Boundary)
            call applyDirichletBC(Problem%Neq, Boundary(ibnd), Sol%TL_Rank, NR%Flag_NR, LS)
        end do
        
        call finalAssembly(NR%Flag_NR, LS)

        tEnd = omp_get_wtime()
        if (output) write(*, '(a,f8.3,a)') ' / Assembly time: ', tEnd - tStart, ' s'

    end subroutine assemblePetscObjects

    !--------------------------------------------------------------------------

    subroutine initializeNewtonVariables(NR)

        use NewtonRaphsonParameters_mod, only: Niter

        implicit none

        type(NewtonType), intent(inout) :: NR

        NR%Iter_f = 0

        NR%xNiter = Niter

        NR%Res_Norm = 1.0d0
        NR%Cor_Norm_Old_f = 1.0d0 ; NR%Cor_Norm_New_f = 1.0d0

        NR%Flag_NR = 'NRP'

    end subroutine initializeNewtonVariables

    !--------------------------------------------------------------------------

    function checkMNR(iter0, cor_norm0) result(Flag_NR)

        use NewtonRaphsonParameters_mod, only: MNR_eps
        use ContinuationVariables_mod, only: Increment

        implicit none

        PetscInt, intent(in) :: iter0
        PetscReal, intent(in) :: cor_norm0
        character(3) :: Flag_NR

        PetscBool :: check_MNR, check_NRP

        Flag_NR = 'NRP'

        check_MNR = (cor_norm0 < MNR_eps)
        if (check_MNR) Flag_NR = 'MNR'

        check_NRP = (mod(Increment,100) == 0 .and. cor_norm0 > MNR_eps/10.0d0)
        ! check_NRP = check_NRP .or. (iter0 >= 10)
        check_NRP = check_NRP .or. (mod(Increment,1000) == 0)

        if (check_NRP) Flag_NR = 'NRP'

    end function checkMNR

    !-----------------------------------------------------------------------

    subroutine checkNewtonConvergence(output, NR, emergency)

        use NewtonRaphsonParameters_mod, only: MNR_eps, Cor_Norm_div

        implicit none

        PetscBool, intent(in) :: output
        type(NewtonType), intent(inout) :: NR
        PetscBool, intent(out) :: emergency

        emergency = .false.

        if (NR%Cor_Norm_New_f /= NR%Cor_Norm_New_f) then
            if (output) write(*,'(a)') 'NaN ENCOUNTERED, GENERAL STOP!'
            stop
        end if

        KIND_OF_NEWTON_RAPHSON: select case (NR%Flag_NR)

        case ('NRP')

            if (output) then
                write(*,48) NR%Iter_f, NR%Flag_NR, NR%Res_Norm, NR%Cor_Norm_New_f
            end if

            if (NR%Cor_Norm_New_f < MNR_eps) then
                NR%Flag_NR = 'MNR'
            else
                NR%Flag_NR = 'NRP'
            end if

        case ('MNR')

            if (output) then
                write(*,48) NR%Iter_f, NR%Flag_NR, NR%Res_Norm, NR%Cor_Norm_New_f
            end if

            if (NR%Cor_Norm_New_f > NR%Cor_Norm_Old_f) then
                NR%Flag_NR = 'NRP'
            else
                NR%Flag_NR = 'MNR'
            end if

        case default

            if (output) write(*,'(a)')' INCORRECT CHOICE OF NEWTON RAPHSON FLAG '

        end select KIND_OF_NEWTON_RAPHSON

        ! !ONLY FULL NEWTON RAPHSON 
        ! NR%Flag_NR = 'NRP'

        if (NR%Cor_Norm_New_f > Cor_Norm_div) then
            if (output) then
                write(*,*)
                write(*,'(a)') 'TOO LARGE NORM. RESTARTING WITH RELAXATION!'
                write(*,'(a)') '--------------------------------------------------------------'
            end if
            emergency = .true.
            return
        end if

        if (NR%Iter_f > NR%xNiter) then
            if (output) then
                write(*,*)
                write(*,'(a)') 'TOO MANY ITERATIONS. RESTARTING WITH RELAXATION!'
                write(*,'(a)') '--------------------------------------------------------------'
            end if
            emergency = .true.
            return
        end if

        NR%Cor_Norm_Old_f = NR%Cor_Norm_New_f

        48 format(i5,'. ',a,': RES_NORM = ', es12.4, ', COR_NORM = ', es12.4)

    end subroutine checkNewtonConvergence

    !--------------------------------------------------------------------------

    subroutine emergencyNewtonRaphson(output, NR, Sol)

        use ContinuationParameters_mod, only: Continuation_Method
        use NewtonRaphsonParameters_mod, only: Niter
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: dSo, Cvar1

        implicit none

        PetscBool, intent(in) :: output
        type(NewtonType), intent(inout) :: NR
        type(SolutionArraysType), intent(inout) :: Sol

        PetscBool :: check        

        call initializeNewtonVariables(NR)

        NR%xNiter = 50*Niter
        NR%xF = 0.5d0*NR%xF
        if (NR%xF < 1.0d-2) then
            if (output) write(*,'(a)') 'VERY SMALL RELAXATION FACTOR, GENERAL STOP!'
            stop
        end if

        if (output) then
            write(*,60) NR%xF, trim(adjustl(Cvar1%name)), Cvar1%p
        end if

        if (Rank == 0) then
            if (NR%name == 'Projection') then
                Sol%TL(:,:) = 0.1d0
            else
                ! Sol%TL(:,:) = Sol%TLo(:,:)
                Sol%TL(:,:) = Sol%TLp(:,:)
            end if
        end if

        if (NR%name == 'Projection') then
            Sol%TL_Rank(:,:) = 0.1d0
        else
            ! Sol%TL_Rank(:,:) = Sol%TLo_Rank(:,:)
            ! Sol%EU(:) = Sol%EUo(:)
            Sol%TL_Rank(:,:) = Sol%TLp_Rank(:,:)
            Sol%EU(:) = Sol%EUp(:)
        end if

        check = (Continuation_Method == 'Arclength') .and. (NR%name == 'Main')
        if (check) dSo = 0.5d0*dSo

        60 format('xF = ', f8.4, ', ', a, ' = ', f12.4)

    end subroutine emergencyNewtonRaphson

    !---------------------------------------------------------------------

    subroutine perturbSolution(Sol)

        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Cvar1

        implicit none

        type(SolutionArraysType), intent(inout) :: Sol

        PetscReal :: factor

        factor = 0.5d0

        if (Rank == 0) then

            write(*,'(a,f8.4)') 'Perturbation: ', factor

            Sol%TL(:,:) = (1.0d0-factor)*(Sol%TL(:,:)) + factor*(Sol%TL_d(:,:))
            Sol%EU(:) = (1.0d0-factor)*(Sol%EU(:)) + factor*(Sol%EU_d(:))
            Sol%TLp(:,:) = Sol%TL(:,:)
            Sol%EUp(:) = Sol%EU(:)

        end if

        Sol%TL_Rank(:,:) = (1.0d0-factor)*(Sol%TL_Rank(:,:)) + factor*(Sol%TL_d_Rank(:,:))
        Sol%TLp_Rank(:,:) = Sol%TL_Rank(:,:)

    end subroutine perturbSolution

    !--------------------------------------------------------------------------

    subroutine updateNewton(Problem, FE, Connectivity_Rank, xF, LS_s_f, Sol)

        use Tools_mod, only: getRows, getRankNode
        use MPIParameters_mod, only: Rank

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:,:), intent(in) :: Connectivity_Rank
        PetscReal, intent(in) :: xF
        Vec, intent(in) :: LS_s_f
        type(SolutionArraysType), intent(inout) :: Sol

        PetscInt :: inod, ieq, irow, iel_rank, rnod
        PetscInt, dimension(Problem%Neq) :: ieqs
        PetscInt, dimension(Problem%Neq) :: irows_eqs
        PetscInt, dimension(1) :: gnod
        Vec :: s_f
        PetscScalar, pointer :: s_f_pointer(:)
        VecScatter :: ctx
        PetscErrorCode :: ierr

        call VecScatterCreateToAll(LS_s_f, ctx, s_f, ierr)
        !INSERT_VALUES: any location non-scattered retains its old value
        !SCATTER_fORWARD: scatter from LS_s_f to s_f
        call VecScatterBegin(ctx, LS_s_f, s_f, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterEnd(ctx, LS_s_f, s_f, INSERT_VALUES, SCATTER_FORWARD, ierr)

        !Here, we created the pointer to access the values of the solution vector s_f
        call VecGetArrayReadF90(s_f, s_f_pointer, ierr)

        if (Rank == 0) then
            irow = 0
            do inod = 1, FE%Nnodes
                do ieq = 1, Problem%Neq
                    irow = irow+1
                    Sol%TL(inod,ieq) = Sol%TL(inod,ieq) - xF*s_f_pointer(irow)
                end do
            end do
        end if

        ieqs = [ (ieq, ieq = 1, Problem%Neq) ]
        do iel_rank = 1, FE%Nel_Rank
            do inod = 1, FE%Nbf
                gnod(1) = Connectivity_Rank(iel_rank,inod)
                rnod = getRankNode(iel_rank,inod,FE%Nbf)
                irows_eqs(:) = getRows(gnod(:),ieqs(:),Problem%Neq)
                Sol%TL_Rank(rnod,:) = Sol%TL_Rank(rnod,:) - xF*s_f_pointer(irows_eqs(:))
            end do
        end do

        do ieq = 1, Problem%Nex
            Sol%EU(ieq) = Sol%EU(ieq) - xF*s_f_pointer(Problem%Nunknowns+ieq)
        end do

        call VecRestoreArrayReadF90(s_f, s_f_pointer, ierr)

        call VecScatterDestroy(ctx, ierr)
        call VecDestroy(s_f, ierr)

    end subroutine updateNewton

end module NR_mod
