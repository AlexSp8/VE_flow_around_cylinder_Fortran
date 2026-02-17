
#include <petsc/finclude/petscksp.h>

module Petsc_mod

#include <slepc/finclude/slepceps.h>
    use slepceps
    use petscksp
    use PhysicalParameters_mod, only: ProblemParameters
    use FEMParameters_mod, only: FEParameters
    use MeshParameters_mod, only: MeshParameters
    use BoundaryParameters_mod, only: BoundaryParameters
    use NewtonRaphsonVariables_mod, only: NewtonType
    use LinearSystemVariables_mod, only: LinearSystemType
    use SlepcVariables_mod, only: EVPType

    implicit none

    contains

    subroutine initializePetsc(continuation_method)

        use MPIParameters_mod, only: Rank, NRanks
        use Tools_mod, only: makeDirectories

        implicit none

        character(*), intent(in) :: continuation_method

        PetscErrorCode :: ierr

        ! if (continuation_method == 'Stability') then
            call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
        ! else
        !     call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
        ! end if

        call MPI_Comm_Rank(PETSC_COMM_WORLD, Rank, ierr)
        call MPI_Comm_size(PETSC_COMM_WORLD, NRanks, ierr)

        if (Rank == 0) then
            call makeDirectories(continuation_method)
        end if

        call PetscMemorySetGetMaximumUsage(ierr)

    end subroutine initializePetsc

    !---------------------------------------------------------------------
    
    subroutine createPetscObjects(Problem, FE, Mesh, Boundary, LS)

        use ContinuationParameters_mod, only: Continuation_Method
        use NewtonRaphsonParameters_mod, only: Iterative_Solver, Monitor_Convergence, &
            Solver_Name, PC_Name

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(LinearSystemType), intent(inout) :: LS

        PetscInt ::  N
        PetscBool :: check
        PetscViewer :: viewer
        PetscViewerAndFormat :: vf
        PetscErrorCode :: ierr

        N = Problem%Nunknowns+Problem%Nex
        if (N == 0) return
        
        call KSPCreate(PETSC_COMM_WORLD, LS%ksp, ierr)

        call createPetscVector(LS%b_f, N)
        call createPetscVector(LS%s_f, N)
        call createPetscMatrix(LS%A_f, N)

        call KSPSetOperators(LS%ksp, LS%A_f, LS%A_f, ierr)

        LS%maxit = 0
        LS%iter = 0
        LS%reCreate = .false.

        check = (FE%Ndim == 3)
        check = check .and. (Iterative_Solver)
        ! check = check .and. (Continuation_Method == 'Transient')
        if (check) then

            LS%name = 'Iterative'

            !Select the iterative solver from Krylov space: 
            !KSPGMRES (default), KSPFGMRES (Flexible)
            !KSPCG, KSPBCGS, KSPCGS, KSPBICG
            !KSPRICHARDSON, KSPCHEBYSHEV, KSPTCQMP, KSPTFQMR, KSPCR, KSPLSQR
            !Flexible methods allow the preconditioner to be non-linear
            call setKSPSolver(Solver_Name, LS)

            !PCJACOBI, PCBJACOBI, PCSOR
            !PCICC: Incomplete Cholesky, PCCHOLESKY: Cholesky
            !PCILU: Incomplete LU
            !PCASM: Additive Schwarz, PCGASM: Generalized Addititve Schwarz
            !PCGAMG: Algebraic Multigrid
            !PCSHELL: Shell for user-defined
            call setKSPPreconditioner(PC_Name, Problem, LS)

            !Selecting the norm for convergence test
            !KSP_NORM_UNPRECONDITIONED: 12 norm of the true b - Ax residual
            call KSPSetNormType(LS%ksp,KSP_NORM_UNPRECONDITIONED,ierr)

            if (Monitor_Convergence) then
                call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'Results/Base/Info/GMRES.dat',viewer,ierr) 
                call PetscViewerAndFormatCreate(viewer,PETSC_VIEWER_DEFAULT,vf,ierr)
                !Monitor the convergence history of the iterative solver
                call KSPMonitorSet(LS%ksp,KSPMonitorTrueResidual,vf,PetscViewerAndFormatDestroy,ierr)
            end if

            call KSPSetFromOptions(LS%ksp, ierr)

            call preAllocatePetsc(Problem, FE, Mesh, Boundary, LS%A_f)

            return

        end if

        LS%name = 'Direct'

        !Direct solver: KSPPREONLY, KSPNONE
        call KSPSetType(LS%ksp, KSPPREONLY, ierr)

        call setKSPPreconditioner('PCLU', Problem, LS)

        !Set direct solver: MATSOLVERPASTIX, MATSOLVERMUMPS, MATSOLVERMKL_CPARDISO
        call PCFactorSetMatSolverType(LS%pc, MATSOLVERMUMPS, ierr)

        !For factorization: to compute det(A)
        call PCFactorSetUpMatSolverType(LS%pc,ierr)
        call PCFactorGetMatrix(LS%pc, LS%A_f_factored, ierr)
        !A_f, ICNTL index, value
        call MatMumpsSetIcntl(LS%A_f_factored,33,1,ierr) !For the determinant
        ! call MatMumpsSetIcntl(LS%A_f_factored,24,1,ierr) !For null pivots
        ! call MatMumpsSetIcntl(LS%A_f_factored,11,1,ierr) !For the condition numbers

        call PetscLogDefaultBegin(ierr)
        call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'Results/Base/Info/log_info.dat',viewer,ierr)
        call PetscLogView(viewer,ierr)
        call PetscViewerDestroy(viewer, ierr)

        call KSPSetFromOptions(LS%ksp, ierr)

        call preAllocatePetsc(Problem, FE, Mesh, Boundary, LS%A_f)

    end subroutine createPetscObjects

    !---------------------------------------------------------------------

    subroutine createPetscObjectsProjection(Problem, Ndim, LS)

        use ContinuationParameters_mod, only: Continuation_Method
        use NewtonRaphsonParameters_mod, only: Iterative_Solver, &
            Solver_Name, PC_Name

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Ndim
        type(LinearSystemType), intent(inout) :: LS

        PetscInt ::  N
        PetscBool :: check
        PetscErrorCode :: ierr

        N = Problem%Nunknowns_proj
        if (N == 0) return

        call KSPCreate(PETSC_COMM_WORLD, LS%ksp, ierr)

        call createPetscVector(LS%b_f, N)
        call createPetscVector(LS%s_f, N)
        call createPetscMatrix(LS%A_f, N)

        call KSPSetOperators(LS%ksp, LS%A_f, LS%A_f, ierr)

        LS%maxit = 0
        LS%iter = 0
        LS%reCreate = .false.

        check = (Ndim == 3)
        check = check .and. (Iterative_Solver)
        ! check = check .and. (Continuation_Method == 'Transient')
        if (check) then
            LS%name = 'Iterative'
            call setKSPSolver(Solver_Name, LS)
            call setKSPPreconditioner(PC_Name, Problem, LS)
            call KSPSetNormType(LS%ksp,KSP_NORM_UNPRECONDITIONED,ierr)
            call KSPSetFromOptions(LS%ksp, ierr)
            return
        end if

        LS%name = 'Direct'
        call KSPSetType(LS%ksp, KSPPREONLY, ierr)
        call setKSPPreconditioner('PCLU', Problem, LS)
        call PCFactorSetMatSolverType(LS%pc, MATSOLVERMUMPS, ierr)

        call KSPSetFromOptions(LS%ksp, ierr)

    end subroutine createPetscObjectsProjection

    !---------------------------------------------------------------------
    
    subroutine reCreatePetscObjects(Problem, FE, Mesh, Boundary, LS)

        use MPIParameters_mod, only: Rank

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(LinearSystemType), intent(inout) :: LS

        PetscInt ::  N
        PetscErrorCode :: ierr
        PetscViewer :: viewer

        if (Rank == 0) then
            write(*,'(a)') 'Changing LS from Iterative to Direct!'
        end if
        
        call destroyLinearSystemObjects(LS)

        N = Problem%Nunknowns+Problem%Nex
        call KSPCreate(PETSC_COMM_WORLD, LS%ksp, ierr)
        call createPetscVector(LS%b_f, N)
        call createPetscVector(LS%s_f, N)
        call createPetscMatrix(LS%A_f, N)
        call KSPSetOperators(LS%ksp, LS%A_f, LS%A_f, ierr)

        LS%name = 'Direct'

        LS%maxit = 0
        LS%iter = 0
        LS%reCreate = .false.
        
        !Direct solver: KSPPREONLY, KSPNONE
        call KSPSetType(LS%ksp, KSPPREONLY, ierr)

        call setKSPPreconditioner('PCLU', Problem, LS)

        !Set direct solver: MATSOLVERPASTIX, MATSOLVERMUMPS, MATSOLVERMKL_CPARDISO
        call PCFactorSetMatSolverType(LS%pc, MATSOLVERMUMPS, ierr)
        !For factorization: to compute det(A)
        call PCFactorSetUpMatSolverType(LS%pc,ierr)
        call PCFactorGetMatrix(LS%pc, LS%A_f_factored, ierr)
        call MatMumpsSetIcntl(LS%A_f_factored, 33, 1,ierr) ! preparing for the determinant
        call MatMumpsSetIcntl(LS%A_f_factored, 24, 1,ierr) ! preparing for null pivots

        call PetscLogDefaultBegin(ierr)
        call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'Results/Base/Info/log_info.dat',viewer,ierr)
        call PetscLogView(viewer,ierr)
        call PetscViewerDestroy(viewer, ierr)

        call KSPSetFromOptions(LS%ksp, ierr)

        call preAllocatePetsc(Problem, FE, Mesh, Boundary, LS%A_f)

    end subroutine reCreatePetscObjects

    ! ----------------------------------------------------------------------

    subroutine createPetscVector(petsc_Vec, size_Vec)

        implicit none

        PetscInt, intent(in)  :: size_Vec
        Vec, intent(out) :: petsc_Vec

        PetscErrorCode :: ierr

        call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, size_Vec, petsc_Vec, ierr)
        call VecSetFromOptions(petsc_Vec, ierr)
        call VecZeroEntries(petsc_Vec, ierr)

    end subroutine createPetscVector

    ! ----------------------------------------------------------------------

    subroutine createPetscMatrix(petsc_Mat, size_Mat)

        use MPIParameters_mod, only: Rank, NRanks

        implicit none

        PetscInt, intent(in) :: size_Mat
        Mat, intent(inout) :: petsc_Mat

        PetscInt :: grow1, growf, nrows, ncolumns
        PetscInt :: localDiagonalNonZeros, localOffDiagonalNonZeros
        PetscErrorCode :: ierr

        call MatCreate(PETSC_COMM_WORLD, petsc_Mat, ierr)
        call MatSetType(petsc_Mat, MATMPIAIJ, ierr)
        call MatSetSizes(petsc_Mat, PETSC_DECIDE, PETSC_DECIDE, size_Mat, size_Mat, ierr)

        ! nrows = size_Mat
        ! ncolumns = size_Mat
        ! localDiagonalNonZeros = (size_Mat)/NRanks
        ! localOffDiagonalNonZeros = (size_Mat)/NRanks
        ! call MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE,PETSC_DECIDE, nrows,ncolumns, &
        !                 localDiagonalNonZeros, PETSC_NULL_INTEGER, &
        !                 localOffDiagonalNonZeros,PETSC_NULL_INTEGER, petsc_Mat, ierr)

        call MatSetFromOptions(petsc_Mat, ierr)

        call MatSetUp(petsc_Mat, ierr)

        ! call MatGetOwnershipRange(petsc_Mat, grow1, growf, ierr)
        ! print*, Rank, size_Mat, grow1, growf

        call MatZeroEntries(petsc_Mat, ierr)

        ! !Force diagonal entries to be allocated (PETSC_TRUE)
        ! call MatSetOption(petsc_Mat, MAT_FORCE_DIAGONAL_ENTRIES, PETSC_TRUE, ierr)

        ! !Ignore new non-zero insertions (PETSC_FALSE)
        ! call MatSetOption(petsc_Mat, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE, ierr)

        !Any new entry in the structure will produce an error (PETSC_TRUE)
        call MatSetOption(petsc_Mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE, ierr)

        !Any new entry that has not been pre-allocated will produce an error (PETSC_TRUE)
        call MatSetOption(petsc_Mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE, ierr)

        !Zero entries (when MatZeroRows is called) are kept in the non-zero structure (PETSC_TRUE)
        call MatSetOption(petsc_Mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr)

    end subroutine createPetscMatrix

    !--------------------------------------------------------------------------

    subroutine setKSPSolver(solver_name, LS)

        implicit none

        character(*), intent(in) :: solver_name
        type(LinearSystemType), intent(inout) :: LS

        PetscInt :: maxit
        PetscErrorCode :: ierr

        select case (solver_name)

        case ('KSPGMRES')

            call KSPSetType(LS%ksp,KSPGMRES,ierr)
            maxit = 1000
            call KSPGMRESSetRestart(LS%ksp,maxit,ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_gmres_restart","1000",ierr)

            ! !Pre-allocate all needed vectors
            ! call KSPGMRESSetPreAllocateVectors(LS%ksp,ierr)
            ! ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_gmres_preallocate",PETSC_NULL_OPTIONS,ierr)

            !Orthogonalization of Krylov vectors method
            !KSPGMRESClassicalGramSchmidtOrthogonalization (default)
            ! call KSPGMRESSetOrthogonalization(LS%ksp,KSPGMRESClassicalGramSchmidtOrthogonalization,ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_gmres_classicalgramschmidt",ierr)
            !Modified Gram-Schmidt
            ! call KSPGMRESSetOrthogonalization(LS%ksp,KSPGMRESModifiedGramSchmidtOrthogonalization,ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_gmres_modifiedgramschmidt",ierr)

            !Iterative refinement (not good for orthogonalization stability)
            !KSP_GMRES_CGS_REFINE_NEVER (default)
            !KSP_GMRES_CGS_REFINE_IFNEEDED: 1-step refinement if orthogonality is poor
            !KSP_GMRES_CGS_REFINE_ALWAYS: 2-step refinement
            call KSPGMRESSetCGSRefinementType(LS%ksp,KSP_GMRES_CGS_REFINE_IFNEEDED,ierr)

        case ('KSPFGMRES')

            call KSPSetType(LS%ksp,KSPFGMRES,ierr)

            maxit = 1000
            call KSPGMRESSetRestart(LS%ksp,maxit,ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_gmres_restart","1000",ierr)

            ! !Pre-allocate all needed vectors
            ! call KSPGMRESSetPreAllocateVectors(LS%ksp,ierr)
            ! ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_gmres_preallocate",PETSC_NULL_OPTIONS,ierr)

            !Orthogonalization of Krylov vectors method
            !KSPGMRESClassicalGramSchmidtOrthogonalization (default)
            ! call KSPGMRESSetOrthogonalization(LS%ksp,KSPGMRESClassicalGramSchmidtOrthogonalization,ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_gmres_classicalgramschmidt",ierr)
            !Modified Gram-Schmidt
            ! call KSPGMRESSetOrthogonalization(LS%ksp,KSPGMRESModifiedGramSchmidtOrthogonalization,ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-ksp_gmres_modifiedgramschmidt",ierr)

            !Iterative refinement (not good for orthogonalization stability)
            !KSP_GMRES_CGS_REFINE_NEVER (default)
            !KSP_GMRES_CGS_REFINE_IFNEEDED: 1-step refinement if orthogonality is poor
            !KSP_GMRES_CGS_REFINE_ALWAYS: 2-step refinement
            call KSPGMRESSetCGSRefinementType(LS%ksp,KSP_GMRES_CGS_REFINE_IFNEEDED,ierr)

        case ('KSPPREONLY')

            call KSPSetType(LS%ksp, KSPPREONLY, ierr)

        case default
            write(*,'(a)') 'Not supported solver_name in setKSPSolver!'
            stop
        end select

    end subroutine setKSPSolver

    !--------------------------------------------------------------------------

    subroutine setKSPPreconditioner(pc_name, Problem, LS)

        use MPIParameters_mod, only: Rank

        implicit none

        character(*), intent(in) :: pc_name
        type(ProblemParameters), intent(in) :: Problem
        type(LinearSystemType), intent(inout) :: LS

        !Index Sets for each block
        IS :: IS_velocity, IS_stress
        PetscInt, dimension(:), allocatable :: irows_velocity, irows_stress
        PetscErrorCode :: ierr

        !Set Preconditioner
        call KSPGetPC(LS%ksp,LS%pc,ierr)

        select case (pc_name)

        case ('PCFIELDSPLIT')

            call PCSetType(LS%pc,PCFIELDSPLIT,ierr)

            call setPetscISFieldSplit(Problem, LS, irows_velocity, irows_stress)

            call ISCreateGeneral(MPI_COMM_WORLD, size(irows_velocity), irows_velocity, &
                                PETSC_COPY_VALUES, IS_velocity, ierr)
            call ISCreateGeneral(MPI_COMM_WORLD, size(irows_stress), irows_stress, &
                                PETSC_COPY_VALUES, IS_stress, ierr)

            deallocate(irows_velocity)
            deallocate(irows_stress)

            call PCFieldSplitSetIS(LS%pc,"0",IS_velocity,ierr)
            call PCFieldSplitSetIS(LS%pc,"1",IS_stress,ierr)

            !Type of relaxation/factorization splitting
            !additive: addition of all preconditioners
            !multiplicative (default): new preconditioner computation 
            !                          and sequential application to the residual
            !symmetric_multiplicative: multiplicative for symmetric matrices
            !schur: Schur complement of the matrix from 2 blocks
            !gkb: generalized Golub-Kahan bidiagonalization preconditioner
            ! call PCFieldSplitSetType(LS%pc,PC_COMPOSITE_SCHUR,ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-pc_fieldsplit_type","schur",ierr)

            !The operator to construct the Schur complement
            !a11 (default): from A11
            !self: from the symbolic representation of the Schur complement matrix
            !selfp: from an explicitly-assembled approximation 
            !       Sp=A11−A10*diag(A00)^(-1)*A01
            !user: from user provided matrix
            !full: from the exact Schur complement of the matrix (expensive)
            ! call PCFieldSplitSetSchurPre(LS%pc,PC_FIELDSPLIT_SCHUR_PRE_A11,PETSC_NULL_MAT,ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-pc_fieldsplit_schur_precondition","a11",ierr)

            !Which blocks of the approximate block factorization to retain in the preconditioner
            !diag: D
            !lower: L D
            !upper: D U
            !full (default): L (D U)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-pc_fieldsplit_schur_fact_type","lower",ierr)

            !Block 0 options
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_ksp_type","gmres",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_ksp_max_it","10",ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_ksp_ksp_rtol","1e-4",ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_ksp_ksp_atol","1e-4",ierr)

            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_pc_type","asm",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_pc_asm_type","interpolate",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_pc_asm_overlap","1",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_pc_asm_local_type","additive",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_sub_ksp_type","preonly",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_sub_pc_type","ilu",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_sub_pc_factor_levels","1",ierr)

            ! !Block 1 options
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_1_ksp_type","preonly",ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_1_pc_type","bjacobi",ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_1_sub_ksp_type","preonly",ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_1_sub_pc_type","ilu",ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_1_sub_pc_factor_levels","1",ierr)

            ! !Check
            ! !Block 0 options
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_ksp_type","gmres",ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_ksp_max_it","10",ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_ksp_ksp_rtol","1e-4",ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_ksp_ksp_atol","1e-4",ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_pc_type","asm",ierr)
            ! !basic, interpolate, restrict, none
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_pc_asm_type","interpolate",ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_pc_asm_overlap","1",ierr)
            ! !additive, multiplicative
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_pc_asm_local_type","additive",ierr)

            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_sub_ksp_type","preonly",ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_sub_pc_type","ilu",ierr)
            ! call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_0_sub_pc_factor_levels","1",ierr)

            !Block 1 options
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_1_ksp_type","preonly",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_1_pc_type","asm",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_1_pc_asm_type","interpolate",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_1_pc_asm_overlap","2",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_1_pc_asm_local_type","additive",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_1_sub_ksp_type","preonly",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_1_sub_pc_type","ilu",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-fieldsplit_1_sub_pc_factor_levels","1",ierr)

        case ('PCASM')

            call PCSetType(LS%pc,PCASM,ierr)

            !ASM options
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-pc_asm_type","interpolate",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-pc_asm_local_type","additive",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-sub_ksp_type","preonly",ierr)  
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-sub_pc_type","ilu",ierr)
            call PetscOptionsSetValue(PETSC_NULL_OPTIONS,"-sub_pc_factor_levels","1",ierr)

            !Overlap to compute constructing subdomains
            call PCASMSetOverlap(LS%pc,1,ierr)

            !PC_ASM_BASIC: full restriction and interpolation operators
            !PC_ASM_RESTRICT: full restriction, ignores off-process interpolation
            !PC_ASM_INTERPOLATE: limited restriction with full interpolation
            !PC_ASM_NONE: ignores off-process values for both restriction and interpolation
            call PCASMSetType(LS%pc, PC_ASM_INTERPOLATE, ierr)

        case ('PCLU')

            call PCSetType(LS%pc, PCLU, ierr)

        case default
            write(*,'(a)') 'Not supported pc_name in setKSPPreconditioner!'
            stop
        end select

    end subroutine setKSPPreconditioner

    !--------------------------------------------------------------------------

    subroutine setPetscISFieldSplit(Problem, LS, irows_velocity, irows_stress)

        use MPIParameters_mod, only: Rank, NRanks

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(LinearSystemType), intent(in) :: LS
        PetscInt, dimension(:), allocatable, intent(out) :: irows_velocity, irows_stress

        PetscInt :: istart, iend, isize, irow, i_Velocity, i_stress, inex
        PetscErrorCode :: ierr

        !Range of rows owned by MPI process (global indices)
        call MatGetOwnershipRange(LS%A_f, istart, iend, ierr)

        iend = iend-1

        isize = iend - istart
        !Extra unknowns added to velocity block
        allocate(irows_velocity(0:isize+Problem%Nex)) ; irows_velocity(:) = -1
        allocate(irows_stress(0:isize)) ; irows_stress(:) = -1

        !Substract Nex (owned by the last process)
        if (Rank == NRanks-1) then
            iend = iend - Problem%Nex
        end if

        i_Velocity = 0 ;  i_stress = 0
        do irow = istart, iend
            !Using 0-based index
            if (irow < Problem%Neq) then
                if (irow >= Problem%Neq_f+1) then
                    irows_stress(i_stress) = irow
                    i_stress = i_stress + 1
                else
                    irows_velocity(i_Velocity) = irow 
                    i_Velocity = i_Velocity + 1
                end if
            else
                if (mod(irow,Problem%Neq) >= Problem%Neq_f+1) then
                    irows_stress(i_stress) = irow
                    i_stress = i_stress + 1
                else
                    irows_velocity(i_Velocity) = irow 
                    i_Velocity = i_Velocity + 1
                end if
            end if

        end do

        if (rank == NRanks-1) then
            do inex = 1, Problem%Nex
                irows_velocity(i_Velocity - 1 + inex) = Problem%Nunknowns - Problem%Nex + inex
                print*, i_Velocity -1 + inex
            end do
        end if

        irows_velocity = pack(irows_velocity, irows_velocity /= -1)
        irows_stress = pack(irows_stress, irows_stress /= -1)

    end subroutine setPetscISFieldSplit

    !--------------------------------------------------------------------------

    subroutine preAllocatePetsc(Problem, FE, Mesh, Boundary, A_f)

        use Tools_mod, only: getRows, indexInArray
        use Storage_mod, only: storeJacobian, storeJacobian_Ac, storeJacobian_Ar, &
            storeJacobian_Ah, storeJacobianPeriodic
        use ContinuationParameters_mod, only: Continuation_Method

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        Mat, intent(inout) :: A_f

        PetscInt, dimension(1), parameter :: ieqs = [1]
        PetscScalar, parameter :: value = 0.0d0
        PetscInt, dimension(1) :: gnod_periodic, grows_periodic
        PetscInt :: iel_rank, inod, i, j, ibnd, gnod, gnod_index, jcol
        PetscInt :: Nrows, Mcols, inex, irow, Nex, ieq, Neq
        PetscErrorCode :: ierr
        PetscBool :: check
        PetscInt, dimension(FE%Nbf) :: gnodes, grows
        PetscInt, dimension(Problem%Nunknowns+Problem%Nex) :: petsc_cols
        PetscScalar, dimension(Problem%Nunknowns+Problem%Nex) :: values

        PetscInt, dimension(Problem%Neq) :: irows, jeqs
        PetscInt, dimension(1) :: jcol_p

        PetscScalar, dimension(FE%Nbf, Problem%Neq, FE%Nbf, Problem%Neq) :: TEMP_A_f
        PetscScalar, dimension(FE%Nbf, Problem%Neq) :: TEMP_Ac_f
        PetscScalar, dimension(FE%Nbf, Problem%Neq) :: TEMP_Ar_f
        PetscScalar, dimension(Problem%Nex) :: TEMP_Ah_f

        Neq = Problem%Neq
        Nex = Problem%Nex

        TEMP_A_f(:,:,:,:) = 0.0d0
        TEMP_Ac_f(:,:) = 0.0d0
        TEMP_Ar_f(:,:) = 0.0d0
        TEMP_Ah_f(:) = 0.0d0

        loop_elements:do iel_rank = 1, FE%Nel_Rank

            gnodes(:) = Mesh%Connectivity_Rank(iel_rank,:)

            grows(:) = getRows(gnodes(:),ieqs(:),Problem%Neq)

            call storeJacobian(TEMP_A_f, grows, A_f)

            call MatGetSize(A_f, Nrows, Mcols, ierr)
            
            do inex = 1, Nex
                jcol = Mcols - Nex + inex -1
                call storeJacobian_Ac(TEMP_Ac_f, grows, A_f, jcol)
                irow = Nrows - Nex + inex -1
                call storeJacobian_Ar(grows, TEMP_Ar_f, A_f, irow)
                call storeJacobian_Ah(TEMP_Ah_f, A_f, inex, irow)
            end do

            do i = 1, size(Problem%PeriodicBNDs,1)
                ibnd = Problem%PeriodicBNDs(i,1)
                do inod = 1, size(gnodes)
                    gnod = gnodes(inod)
                    gnod_index = indexInArray(Boundary(ibnd)%gnodes_total, gnod)
                    if (gnod_index /= 0) then
                        gnod_periodic(1) = Boundary(ibnd)%gnodes_periodic(gnod_index)
                        grows_periodic(:) = getRows(gnod_periodic,ieqs,Problem%Neq)
                        call storeJacobianPeriodic(TEMP_A_f(:,:,inod:inod,:), &
                                            grows, grows_periodic, A_f)
                        do inex = 1, Nex
                            jcol = Mcols - Nex + inex -1
                            call storeJacobian_Ac(TEMP_Ac_f(inod:inod,:), grows_periodic, A_f, jcol)
                        end do
                    end if
                end do
            end do

        end do loop_elements

        !PBCs
        jeqs = [ (ieq, ieq = 1, size(jeqs)) ]
        do i = 1, size(Problem%PeriodicBNDs,1)
            ibnd = Problem%PeriodicBNDs(i,1)
            do inod = 1, size(Boundary(ibnd)%gnodes_total)

                irows(:) = getRows(Boundary(ibnd)%gnodes_total(inod:inod),jeqs,Neq) - 1

                do j = 1, size(irows)
                    jcol_p(:) = getRows(Boundary(ibnd)%gnodes_periodic(inod:inod),jeqs(j:j),Neq) - 1
                    !Insert +1 and -1
                    call MatSetValues(A_f, 1, irows(j), 1, irows(j), value, ADD_VALUES, ierr)
                    call MatSetValues(A_f, 1, irows(j), 1, jcol_p(1), value, ADD_VALUES, ierr)
                end do

            end do
        end do

        check = (Continuation_Method == 'Arclength')
        check = check .and. (Problem%name == 'Main')
        if (check) then

            call MatGetSize(A_f, Nrows, Mcols, ierr)
            irow = Nrows - 1

            petsc_cols = [ (i, i = 0, size(petsc_cols)-1) ]
            call MatSetValues(A_f, size(petsc_cols), petsc_cols, 1, irow, values, ADD_VALUES, ierr)
            call MatSetValues(A_f, 1, irow, size(petsc_cols), petsc_cols, values, ADD_VALUES, ierr)

        end if

        call MatAssemblyBegin(A_f, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(A_f, MAT_FINAL_ASSEMBLY, ierr)

    end subroutine preAllocatePetsc

    !--------------------------------------------------------------------------

    subroutine initializeLSObjects(Flag_NR, LS)

        implicit none

        character(*), intent(in) :: Flag_NR
        type(LinearSystemType), intent(inout) :: LS

        PetscErrorCode :: ierr

        call VecZeroEntries(LS%b_f, ierr)
        call VecZeroEntries(LS%s_f, ierr)

        if (Flag_NR == 'NRP') then
            call MatZeroEntries(LS%A_f, ierr)
        end if

    end subroutine initializeLSObjects

    !-----------------------------------------------------------------------

    subroutine flushAssembly(Flag_NR, LS)

        implicit none

        character(*), intent(in) :: Flag_NR
        type(LinearSystemType), intent(inout) :: LS

        PetscErrorCode :: ierr

        call VecAssemblyBegin(LS%b_f,ierr)
        if (Flag_NR == 'NRP') then
            call MatAssemblyBegin(LS%A_f, MAT_FLUSH_ASSEMBLY, ierr)
            call MatAssemblyEnd(LS%A_f, MAT_FLUSH_ASSEMBLY, ierr)
        end if
        call VecAssemblyEnd(LS%b_f, ierr)

    end subroutine flushAssembly

    !-----------------------------------------------------------------------

    subroutine finalAssembly(Flag_NR, LS)

        implicit none

        character(*), intent(in) :: Flag_NR
        type(LinearSystemType), intent(inout) :: LS

        PetscErrorCode :: ierr

        call VecAssemblyBegin(LS%b_f,ierr)
        if (Flag_NR == 'NRP') then
            call MatAssemblyBegin(LS%A_f, MAT_FINAL_ASSEMBLY, ierr)
            call MatAssemblyEnd(LS%A_f, MAT_FINAL_ASSEMBLY, ierr)
        end if
        call VecAssemblyEnd(LS%b_f, ierr)

    end subroutine finalAssembly

    !-----------------------------------------------------------------------

    subroutine solveLinearSystem(NR, LS)

        use NewtonRaphsonParameters_mod, only: LS_maxit, LS_eps
        
        implicit none

        type(NewtonType), intent(in) :: NR
        type(LinearSystemType), intent(inout) :: LS

        PetscInt  :: ncomp
        PetscReal :: rtol, atol, dtol
        PetscInt, parameter :: n = 10
        PetscReal, dimension(n) :: Kr, Ki
        PetscErrorCode :: ierr

        if (NR%Flag_NR == 'MNR') then
            call KSPSetReusePreconditioner(LS%ksp, PETSC_TRUE, ierr)
        else
            call KSPSetReusePreconditioner(LS%ksp, PETSC_FALSE, ierr)
        end if

        ! ! We scale the matrix AND the residual vector
        ! ! ATTENTION: The matrix and the residual vector are changed!
        ! call KSPSetDiagonalScale(LS%ksp, PETSC_TRUE,ierr)

        rtol = 1.0d-16
        atol = max( min(1.0d-4, NR%Res_norm*1.0d-4), LS_eps )
        ! atol = 1.0d-12

        if (NR%Res_norm < LS_eps) atol = NR%Res_Norm*1.0d-2
        
        ! if (NR%Iter_f == 1) atol = 1.0d-5
        ! if (NR%Iter_f == 2) atol = 1.0d-6

        dtol = 1.0d8

        !Setting the tolerance
        !rtol: relative convergence tolerance, relative decrease in the residual norm
        !atol: absolute convergence tolerance absolute size of the residual norm
        !dtol: divergence tolerance before KSPConvergedDefault() concludes that the method is diverging
        !maxit: maximum number of iterations
        call KSPSetTolerances(LS%ksp, rtol, atol, dtol, LS_maxit, ierr)

        ! !Eigenvalue approximation of the preconditioner (called before KSPSetUp and KSPSolve)
        ! call KSPSetComputeEigenvalues(LS%ksp, PETSC_TRUE, ierr)

        call KSPSolve(LS%ksp, LS%b_f, LS%s_f, ierr)

        call KSPGetIterationNumber(LS%ksp, LS%iter, ierr)

        LS%maxit = max(LS%maxit,LS%iter)

        ! !n: size of the the Kr, Ki arrays
        ! !Kr, Ki: real and complex part of the eigenvalues
        ! !ncomp: computed eigenvalues
        ! call KSPComputeEigenvalues(LS%ksp, n, Kr, Ki, ncomp, ierr)

        if (NR%name == 'Main') then
            call KSPConvergedReasonView(LS%ksp, PETSC_VIEWER_STDOUT_WORLD, ierr)
        end if

        ! call KSPGetSolution(LS%ksp, LS%s_f)
        ! call KSPGetRhs(LS%ksp, LS%b_f)

    end subroutine solveLinearSystem

    !---------------------------------------------------------------------

    subroutine calculateMatrixDeterminant(A_f_factored, output, iter_f)

        use ContinuationVariables_mod, only: Cvar1

        implicit none

        Mat, intent(in) :: A_f_factored
        PetscBool, intent(in) :: output
        PetscInt, intent(in) :: iter_f

        PetscInt :: c, ierr
        PetscReal :: a, b, detA_sign, detA_signo, cond1, cond2
        PetscScalar :: detA
        PetscScalar, save :: detAo = 0.0d0
        PetscBool :: detA_sign_change

        !RINFOG
        !1: floating-point operations for the elimination process after analysis
        !2: floating-point operations for the assembly process after factorization
        !3: floating-point operations for the elimination process after factorization
        !4-8: For ICNTL(11) = 1 or 2:
        !   4: infinite norm of A, ||A||oo
        !   5: infinite norm of x solution
        !   6: scaled residual ||Ax-b|| / (||A||oo * ||x||oo)
        !   7: ω1 backward error: max( ((|b-Ax|)i) / (|b|+|A||x|)i )
        !   8: ω2 backward error: max( ((|b-Ax|)i) / ((|A||x|)i+(|A||x|)oo )
        !9-11: For ICNTL(11) = 1:
        !   9: ω1*cond1+ω2*cond2
        !  10: cond1: ||A||*(||A||^-1)
        !  11: cond2: σmax(A)/σmin(A), σ singular values of A (sqrt of non-negative eigenvalues)
        !12: real part of det, a
        !13: imaginary part of det, b
        !14: effective floating-point operations for the elimination process after factorization
        !15: disk space in MB for out-of-core factorization
        !16: disk space in MB for out-of-core execution
        !17: file sizes in MB
        !18: MUMPS structures sizes in MB
        !19: smallest absolute pivot (including null pivots)
        !20: smallest absolute pivot (excluding null pivots)
        !21: largest absolute pivot
        !24-40: not used

        !INFOG
        !1: 0 (successful MUMPS call), negative (error), positive (warning)
        !2: information about the error/warning
        !8: structural symmetric % (0: fully unsymmetric, 100: symmetric)
        !15: number of iterative refinements
        !34: exponent of determinant, c

        ! call MatMumpsGetRinfog(A_f_factored,10,cond1,ierr)
        ! call MatMumpsGetRinfog(A_f_factored,11,cond2,ierr)

        !(a+ib)x2^c, a = RINFOG(12), b = RINFOG(13), c = RINFOG(34)
        call MatMumpsGetRinfog(A_f_factored,12,a,ierr)
        call MatMumpsGetRinfog(A_f_factored,13,b,ierr)
        ! call MatMumpsGetInfog(A_f_factored,34,c,ierr)

        detA = (a+b*PETSC_i)!*(2.0d0**c)
        ! call MatMumpsGetNullPivots(A_f_factored, numb_null_pivots, numb_null_pivots_arrays,ierr)

        if (output) then

            ! write(*,'(2(a,f16.8),a,i0)') "a+bi = ", a, " + ", b, " i / c = ", c
            write(*,'(2(a,es16.8),a)'), "detA = ", a, " + ", b, " i"
            ! write(*,'(2(a,es16.8),a)'), "detA = ", PetscRealPart(detA), &
            !                             " + ", PetscImaginaryPart(detA), " i"

            ! write(*,'(a,es16.8)'), "cond1 = ", cond1
            ! write(*,'(a,es16.8)'), "cond2 = ", cond2

            detA_sign = sign(1.0d0,PetscRealPart(detA))
            detA_signo = sign(1.0d0,PetscRealPart(detAo))
            detA_sign_change = (detA_sign*detA_signo < 0.0d0)
            if (detA_sign_change) then
                write(*,'(3a,f12.4)') 'detA < 0 at ', Cvar1%name, ' = ', Cvar1%p
            end if

            write(*, '(a)') '--------------------------------------------------------------'

            open(2,file='Results/Base/Info/detA.dat',position='append')
            write(2,'(f12.4,i5,2(es16.4,a))') Cvar1%p, iter_f, a, " + ", b, " i"
            close(2,status='keep')

            ! print*, "NUMBER OF NULL PIVOTS = ", numb_null_pivots

        end if

        detAo = detA

    end subroutine calculateMatrixDeterminant

    !---------------------------------------------------------------------
    
    subroutine createSlepcObjects(Problem_Stability, FE, Mesh, Boundary, EVP)

        implicit none

        type(ProblemParameters), intent(in) :: Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(EVPType), intent(inout) :: EVP

        PetscInt :: N
        PetscErrorCode :: ierr
        PetscViewer :: viewer

        N = Problem_Stability%Nunknowns + Problem_Stability%Nex
        if (N == 0) return

        call KSPCreate(PETSC_COMM_WORLD, EVP%ksp, ierr)
        call EPSCreate(PETSC_COMM_WORLD, EVP%eps, ierr)
        
        call createPetscMatrix(EVP%A_f, N)
        call createPetscMatrix(EVP%B_f, N)

        call MatCreateVecs(EVP%A_f, EVP%Xr, PETSC_NULL_VEC, ierr)

        call PetscLogDefaultBegin(ierr)
        call PetscViewerASCIIOpen(PETSC_COMM_WORLD,'Results/Base/Info/log_info.dat',viewer,ierr)
        call PetscLogView(viewer,ierr)
        call PetscViewerDestroy(viewer, ierr)

        call preAllocatePetsc(Problem_Stability, FE, Mesh, Boundary, EVP%A_f)
        call preAllocatePetsc(Problem_Stability, FE, Mesh, Boundary, EVP%B_f)

    end subroutine createSlepcObjects

    !-----------------------------------------------------------------------

    subroutine flushAssemblyEVPObjects(A_f, B_f)

        implicit none

        Mat, intent(inout) :: A_f, B_f

        PetscErrorCode :: ierr

        call MatAssemblyBegin(A_f, MAT_FLUSH_ASSEMBLY, ierr)
        call MatAssemblyBegin(B_f, MAT_FLUSH_ASSEMBLY, ierr)
        
        call MatAssemblyEnd(A_f, MAT_FLUSH_ASSEMBLY, ierr)
        call MatAssemblyEnd(B_f, MAT_FLUSH_ASSEMBLY, ierr)

    end subroutine flushAssemblyEVPObjects
    
    !-----------------------------------------------------------------------

    subroutine finalAssemblyEVPObjects(A_f, B_f)

        implicit none

        Mat, intent(inout) :: A_f, B_f

        PetscErrorCode :: ierr

        call MatAssemblyBegin(A_f, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyBegin(B_f, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(A_f, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd(B_f, MAT_FINAL_ASSEMBLY, ierr)

    end subroutine finalAssemblyEVPObjects

    ! ----------------------------------------------------------------------

    subroutine writeMemoryUsage()

        use MPIParameters_mod, only: Rank

        implicit none

        PetscErrorCode :: ierr
        PetscLogDouble :: mem, malloc_mem, max_memory, total_malloc, total_mem

        call PetscMemoryGetCurrentUsage(mem, ierr)
        total_mem = mem!-Mem0
        call MPI_Allreduce(MPI_IN_PLACE, total_mem, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ierr)

        call PetscMallocGetCurrentUsage(malloc_mem, ierr)
        total_malloc = malloc_mem!-Malloc_Mem0
        call MPI_Allreduce(MPI_IN_PLACE, total_malloc, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ierr)

        call PetscMemoryGetMaximumUsage(max_memory, ierr)
        call MPI_Allreduce(MPI_IN_PLACE, max_memory, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ierr)

        if (Rank == 0) then
            open(8, file='Results/Base/Info/Memory_usage.dat',position='append')
            write(8,'(3es12.4)') total_malloc/(1024**3), total_mem/(1024**3), max_memory/(1024**3)
            close(8, status='keep')
            write(*,'(a,es12.4,a)') "Max Memory used: ", max_memory/(1024**3), " GB"
        end if

    end subroutine writeMemoryUsage

    !---------------------------------------------------------------------
    
    subroutine destroyPetscObjects(continuation_method)

        use SlepcVariables_mod, only: EVP_Inflow, EVP_Main
        use LinearSystemVariables_mod, only: LS_Inflow, LS_Main

        implicit none

        character(*), intent(in) :: continuation_method

        call destroyLinearSystemObjects(LS_Inflow)

        call destroyLinearSystemObjects(LS_Main)

        if (continuation_method == 'Stability') then
            call destroyEVPSystemObjects(EVP_Inflow)
            call destroyEVPSystemObjects(EVP_Main)
        end if

    end subroutine destroyPetscObjects

    !---------------------------------------------------------------------
    
    subroutine destroyLinearSystemObjects(LS)

        implicit none

        type(LinearSystemType), intent(out) :: LS

        PetscErrorCode :: ierr

        call KSPDestroy(LS%ksp, ierr)
        call PCDestroy(LS%pc, ierr)

        call MatDestroy(LS%A_f, ierr)
        call MatDestroy(LS%A_f_factored, ierr)

        call VecDestroy(LS%b_f, ierr)
        call VecDestroy(LS%s_f, ierr)

    end subroutine destroyLinearSystemObjects

    ! ----------------------------------------------------------------------

    subroutine destroyEVPSystemObjects(EVP)

        use SlepcVariables_mod, only: EVPType

        implicit none

        type(EVPType), intent(out) :: EVP

        PetscErrorCode :: ierr
        
        call KSPDestroy(EVP%ksp, ierr)
        call EPSDestroy(EVP%eps, ierr)

        call MatDestroy(EVP%A_f, ierr)
        call MatDestroy(EVP%B_f, ierr)

        call VecDestroy(EVP%Xr, ierr)

    end subroutine destroyEVPSystemObjects

    ! ----------------------------------------------------------------------

    subroutine finalizePetsc(continuation_method)

        implicit none

        character(*), intent(in) :: continuation_method
        
        PetscErrorCode :: ierr

        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        
        call destroyPetscObjects(continuation_method)
        
        ! if (continuation_method == 'Stability') then
            call SlepcFinalize(ierr)
        ! else
            call PetscFinalize(ierr)
        ! end if

    end subroutine finalizePetsc

end module Petsc_mod
