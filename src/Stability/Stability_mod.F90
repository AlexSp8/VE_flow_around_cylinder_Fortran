
#include <slepc/finclude/slepceps.h>

module Stability_mod

    use slepceps
    use PhysicalParameters_mod, only: ProblemParameters
    use FEMParameters_mod, only: FEParameters
    use MeshParameters_mod, only: MeshParameters
    use BoundaryParameters_mod, only: BoundaryParameters
    use GaussParameters_mod, only: GaussIntegration
    use SolutionVariables_mod, only: SolutionArraysType
    use ElementVariables_mod, only: NodeArrays
    use SlepcVariables_mod, only: EVPType

    implicit none

    contains

    subroutine stabilityProblem(Problem, Problem_Stability, FE, Mesh, &
        Boundary, GaussInt, Sol, Elem, EVP)

        use omp_lib, only: omp_get_wtime
        use Tools_mod, only: writeElapsedTime, updateEMA, interpolationLinear
        use ContinuationParameters_mod, only: Cvar2_Names
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Cvar1

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(inout) :: Sol
        type(NodeArrays), intent(inout) :: Elem
        type(EVPType), intent(inout) :: EVP

        !smoothing factor (0 < a < 1): a = 1 only current value
        PetscReal :: tStart, tEnd, Kr_Lead, a
        PetscErrorCode :: ierr

        tstart = omp_get_wtime()

        if (Rank == 0) then
            write(*, '(a)') '======================================================================'
            write(*, '(30x,a)') 'Stability'
            write(*, '(a)') '======================================================================'
        end if

        EVP%Kr_Leado = PetscRealPart(EVP%K_Lead)

        if (Rank == 0) then
            Sol%TL_d(:,:) = 0.0d0
        end if
        Sol%TL_d_rank(:,:) = 0.0d0
        Sol%EU_d(:) = 0.0d0

        call assembleSlepcObjects(Problem, Problem_Stability, FE, Mesh, &
            Boundary, GaussInt, Sol, Elem, EVP)

        call loopShifts(Problem_Stability, FE, Mesh, Sol, EVP)

        Kr_Lead = PetscRealPart(EVP%K_Lead)

        a = 0.5d0
        ! if (size(Cvar2_Names) > 0) a = 0.5d0
        EVP%target_real = updateEMA(Kr_Lead, a, EVP%target_real)
        ! EVP%target_real = Kr_Lead

        EVP%critical = (Kr_Lead*EVP%Kr_Leado < 0.0d0)
        if (EVP%critical) then
            EVP%xCr = interpolationLinear(Cvar1%p-Cvar1%dvar, Cvar1%p, &
                EVP%Kr_Leado, Kr_Lead, 0.0d0)
        end if

        if (size(Cvar2_Names) > 0) then
            call checkCriticalConditions(EVP)
        end if

        call postProcessStability(Problem, Problem_Stability, FE, Mesh, &
            GaussInt, Sol, Elem, Kr_Lead)

        call MPI_Barrier(MPI_COMM_WORLD,ierr)

        if (Rank == 0) then
            call writeStabilityInfo(Problem_Stability%name, EVP)
            tEnd = omp_get_wtime()
            call writeElapsedTime(tStart, tEnd, 'Stability solution time')
        end if

    end subroutine stabilityProblem

    !--------------------------------------------------------------------------

    subroutine assembleSlepcObjects(Problem, Problem_Stability, FE, Mesh, &
        Boundary, GaussInt, Sol, Elem, EVP)

        use omp_lib, only: omp_get_wtime
        use MPIParameters_mod, only: Rank
        use EquationsStability_mod, only: setEquationsBulkStability, &
            setEquationsEdgeStability, setEquationsFaceStability
        use Extra_EquationsStability_mod, only: setExtraEquationsBulkStability, &
            setExtraEquationsEdgeStability, setExtraEquationsFaceStability
        use Dirichlet_mod, only: clearEntriesStability, &
            applyDirichletBCStability, applyPeriodicBCStability
        use Petsc_mod, only: finalAssemblyEVPObjects, flushAssemblyEVPObjects
        use ElementCalculations_mod, only: basis, basisInflow_p, basisMain_p

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(inout) :: Sol
        type(NodeArrays), intent(inout) :: Elem
        type(EVPType), intent(inout) :: EVP

        PetscInt :: ibnd, inex, i
        PetscReal :: t_start, t_end, norm_A, norm_B
        PetscReal, save :: norm_Ao = 0.0d0, norm_Bo = 0.0d0
        PetscErrorCode :: ierr
        procedure(basis), pointer :: basis_p

        t_start = omp_get_wtime()

        if (Sol%name == 'Main') then
            basis_p => basisMain_p
        end if

        if (Sol%name == 'Inflow') then
            basis_p => basisInflow_p
        end if

        call MatZeroEntries(EVP%A_f, ierr)
        call MatZeroEntries(EVP%B_f, ierr)

        !Bulk EQs
        call setEquationsBulkStability(Problem, Problem_Stability, FE, Mesh, &
            Boundary, GaussInt, Sol, basis_p, Elem, EVP)
        call flushAssemblyEVPObjects(EVP%A_f, EVP%B_f)

        do inex = 1, Problem_Stability%Nex_b
            call setExtraEquationsBulkStability(Problem, Problem_Stability, FE, Mesh, &
                GaussInt, Sol, basis_p, Elem, inex, EVP)
        end do
        call flushAssemblyEVPObjects(EVP%A_f, EVP%B_f)

        !Boundary EQs
        do ibnd = 1, size(Boundary)
            call finalAssemblyEVPObjects(EVP%A_f, EVP%B_f)
            call clearEntriesStability(Problem_Stability%Neq, &
                Boundary(ibnd), EVP%A_f, EVP%B_f)
            call setEquationsEdgeStability(Problem, Problem_Stability, FE, Mesh, &
                Boundary(ibnd), Boundary, GaussInt, Sol, basis_p, Elem, EVP)
            call setEquationsFaceStability(Problem, Problem_Stability, FE, Mesh, &
                Boundary(ibnd), Boundary, GaussInt, Sol, basis_p, Elem, EVP)
            call setExtraEquationsEdgeStability(Problem, Problem_Stability, FE, Mesh, &
                Boundary(ibnd), GaussInt, Sol, basis_p, Elem, EVP)
            call setExtraEquationsFaceStability(Problem, Problem_Stability, FE, Mesh, &
                Boundary(ibnd), GaussInt, Sol, basis_p, Elem, EVP)
        end do
        call flushAssemblyEVPObjects(EVP%A_f, EVP%B_f)
        
        !PBCs
        do i = 1, size(Problem_Stability%PeriodicBNDs,1)
            ibnd = Problem_Stability%PeriodicBNDs(i,1)
            call finalAssemblyEVPObjects(EVP%A_f, EVP%B_f)
            call applyPeriodicBCStability(Problem_Stability%Neq, Boundary(ibnd), EVP)
        end do

        !Dirichlet
        call finalAssemblyEVPObjects(EVP%A_f, EVP%B_f)
        do ibnd = 1, size(Boundary)
            call applyDirichletBCStability(Problem_Stability%Neq, Boundary(ibnd), &
                EVP%A_f, EVP%B_f)
        end do
        
        call finalAssemblyEVPObjects(EVP%A_f, EVP%B_f)

        ! call MatView(EVP%A_f, PETSC_VIEWER_DRAW_WORLD, ierr)
        ! call MatView(EVP%B_f, PETSC_VIEWER_STDOUT_WORLD, ierr)

        ! call MatNorm(EVP%A_f, NORM_FROBENIUS, norm_A, ierr)
        ! call MatNorm(EVP%B_f, NORM_FROBENIUS, norm_B, ierr)

        t_end = omp_get_wtime()
        if (Rank == 0) then
        !     write(*, '(a,2es16.8)') 'Norm A = ', norm_A, abs(norm_A-norm_Ao)
        !     write(*, '(a,2es16.8)') 'Norm B = ', norm_B, abs(norm_B-norm_Bo)
            write(*, '(a,f8.3,a)') 'EVP Assembly time: ', t_end-t_start, ' s'
        end if
        ! norm_Ao = norm_A
        ! norm_Bo = norm_B

    end subroutine assembleSlepcObjects

    !--------------------------------------------------------------------------

    subroutine loopShifts(Problem_Stability, FE, Mesh, Sol, EVP)

        use omp_lib, only: omp_get_wtime
        use MPIParameters_mod, only: Rank
        use SlepcParameters_mod, only: Nshifts
        use PostProcessStability_mod, only: writeEigenpairs, storeEigenvector

        implicit none

        type(ProblemParameters), intent(in) :: Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(SolutionArraysType), intent(inout) :: Sol
        type(EVPType), intent(inout) :: EVP

        PetscInt :: ishift, nconv, ilead
        PetscReal :: t_start, t_end
        PetscReal :: Kr_Lead, Kr_Leado, Ki_maxo, term
        PetscBool :: check

        EVP%K_Lead = -1.0d8 + PETSC_i*0.0d0
        EVP%Ki_max = 0.0d0
        ! EVP%LwN_Lead = 0.0d0
        do ishift = 1, Nshifts

            t_start = omp_get_wtime()

            call solveEVP(ishift, EVP%Ki_max, EVP)

            t_end = omp_get_wtime()
            if (Rank == 0) write(*, '(a,f8.3,a)') 'EVP Solution time: ', t_end-t_start, ' s'

            nconv = getEVPSolutionInfo(EVP%eps)

            Kr_Leado = PetscRealPart(EVP%K_Lead)
            Ki_maxo = EVP%Ki_max

            !Write eigenvalues, find max real, imaginary parts
            call writeEigenpairs(Problem_Stability, FE, Mesh, EVP, nconv, ishift, ilead)

            Kr_Lead = PetscRealPart(EVP%K_Lead)

            check = (Kr_Lead - Kr_Leado > 1.0d-6)
            if (check) then
                call storeEigenvector(Problem_Stability%Neq, FE, &
                    Mesh%Connectivity_Rank, Sol%TL_d, Sol%TL_d_rank, ilead, EVP)
            end if

            term = abs((EVP%Ki_max-Ki_maxo)/(EVP%Ki_max+1.0d-8))
            check = (term < 1.0d-3)
            if (check) then
                if (Rank == 0) then
                    write(*,*)
                    write(*,2) term, EVP%Ki_max
                    write(*,*)
                end if
                EVP%Ki_max = 2.0d0*EVP%Ki_max
                ! exit
            end if

        end do

        2 format('Minimal shift change: ', es12.4, 5x, 'Ki_max: ', f12.8)

    end subroutine loopShifts
    
    !--------------------------------------------------------------------------

    subroutine solveEVP(ishift, Ki_Lead, EVP)

        use MPIParameters_mod, only: Rank
        use SlepcParameters_mod, only: maxit, nev, EVP_tol

        implicit none

        PetscInt, intent(in) :: ishift
        PetscReal, intent(in) :: Ki_Lead
        type(EVPType), intent(inout) :: EVP

        PetscInt :: ncv, mpd
        PetscReal :: a_int, b_int
        PetscScalar :: shift, anti_shift, target_value
        PetscReal, parameter :: rtol = 1.0d-10, abs_tol = 1.0d-10, div_tol = 1.0d8
        PetscBool, parameter :: true_res = .false.
        PetscErrorCode :: ierr
        PetscViewer :: viewer
        ST :: st

        !EPSPOWER, EPSSUBSPACE, EPSARNOLDI, EPSLANCZOS(not recommended), EPSKRYLOVSCHUR(default)
        !EPSGD, EPSJD, EPSRQCG, EPSLOBPCG, EPSCISS, EPSLYAPII
        !External packages: EPSLAPACK, EPSARPACK, EPSPRIMME, EPSEVSL, EPSTRLAN, EPSBLOPEX
        !EPSSCALAPACK, EPSELPA, EPSELEMENTAL, EPSFEAST
        call EPSSetType(EVP%eps, EPSKRYLOVSCHUR, ierr)

        !Set operators A and B. Generalized eigenvalue problem (Ax = λBx)
        call EPSSetOperators(EVP%eps, EVP%A_f, EVP%B_f, ierr)

        ! !If the user provided initial guesses or constraints, pass them here
        ! call EPSSetInitialSpace(EVP%eps,nini,Iv,ierr)
        ! call EPSSetDeflationSpace(EVP%eps,ncon,Cv,ierr)

        !Set type of EVP problem:
        !EPS_HEP: Hermitian: Ax = λx, A = (AT)*
        !EPS_NHEP(default): non-Hermitian: Ax = λx
        !EPS_GHEP: generalized Hermitian: Ax = λBx, A = (AT)*, B = (BT)*
        !EPS_GHIEP: generalized Hermitian indefinite
        !EPS_GNHEP: generalized non-Hermitian: Ax = λBx
        !EPS_PGNHEP: generalized non-Hermitian with positive semi-definite B
        call EPSSetProblemType(EVP%eps, EPS_GNHEP, ierr)

        !----------------------------------------------------------------------------
        !Set Spectral Transformation (default is no ST to be performed)
        !Extract ST context
        call EPSGetST(EVP%eps, st, ierr)

        !Set ST type
        !STSHIFT: shift: (A-σΙ)x = θx, θ = λ-σ (shift eigenvalues by σ, not altering eigenvectors)
        !STSINVERT: shift+invert (for eigenvalues near σ): ((A-σΙ)**(-1))x = θx
        !           θ = 1/(λ-σ) (shift eigenvalues by σ, not altering eigenvectors)
        !STCAYLEY: (for eigenvalues near σ): ((A-σΙ)**(-1) +(A+vI))x = θx, 
        !           θ = (λ+v)/(λ-σ), v: anti-shift
        !STPRECOND: for preconditioned iterative solvers
        !STFILTER: polynomial filtering: p(A)x = θx
        !STSHELL
        call STSetType(st, STSINVERT, ierr)

        ! !For shift (usually not called, but EPSSetTarget is used to control the shift)
        ! if (ishift == 1) then
        !     shift = 0.0d0
        ! else
        !     shift = CMPLX(0.0d0,ki_max)
        ! end if
        ! call STSetShift(st,shift,ierr)

        ! !For anti-shift
        ! anti_shift = 1.0d0
        ! call STCayleySetAntiShift(st,anti_shift,ierr)

        !Set which eigenpairs are of interest
        !EPS_LARGEST_MAGNITUDE: largest |λ|, EPS_SMALLEST_MAGNITUDE: smallest |λ|
        !EPS_LARGEST_REAL: largest Re(λ), EPS_SMALLEST_REAL: smallest Re(λ)
        !EPS_LARGEST_IMAGINARY: largest Im(λ), EPS_SMALLEST_IMAGINARY: smallest Im(λ)

        ! !EPS_ALL:Compute all eigenvalues in [a,b]
        ! a_int = -100.0d0 ; b_int = 1.0d0
        ! call EPSSetInterval(EVP%eps,a_int,b_int,ierr)

        !For EPSSetTarget
        !EPS_TARGET_MAGNITUDE: smallest |λ-τ|
        !EPS_TARGET_REAL: smallest |Re(λ-τ)|, EPS_TARGET_IMAGINARY: smallest |Im(λ-τ)|
        target_value = EVP%target_real + 1.0d0*(Ki_Lead)*PETSC_i
        call EPSSetTarget(EVP%eps, target_value, ierr)

        if (Rank == 0) then
            write(*,*)
            write(*,'(i3, a)', advance='no') ishift, ". Shift value: "
            write(*,*) target_value
            write(*,*)
        end if

        ! call EPSSetWhichEigenpairs(EVP%eps,EPS_TARGET_REAL,ierr)

        !----------------------------------------------------------------------------
        !ncv = 2*nev ; mpd = nev
        !ncv: the maximum dimension of the subspace to be used by the solver (number of column vectors)
        !ncv should always be between nev and (nev+mpd), typically ncv=nev+mpd
        !mpd: the maximum dimension allowed for the projected problem
        !The parameters ncv and mpd are intimately related (set one of them at most). 
        !Normal usage: (a) in cases where nev is small, the user sets ncv (a reasonable default is 2*nev)
        !and (b) in cases where nev is large, the user sets mpd
        !When computing all eigenvalues in an interval, see EPSSetInterval(), 
        !these parameters lose relevance, and tuning must be done with EPSKrylovSchurSetdimensions()
        call EPSSetdimensions(EVP%eps, nev, PETSC_DEFAULT_INTEGER, PETSC_DEFAULT_INTEGER, ierr)

        !EPS_CONV_ABS: absolute residual: ||r||
        !EPS_CONV_REL: relative residual to eigenvalue : ||r||/|λ| (default, not for λ->0)
        !EPS_CONV_NORM: relative residual to matrix norms: ||r||/(||A|| + |λ|*||B||)
        !EPS_CONV_USER: user-defined
        call EPSSetConvergenceTest(EVP%eps, EPS_CONV_ABS, ierr)

        call EPSSetTolerances(EVP%eps, EVP_tol, maxit, ierr)

        !Tolerance based on true residual, not based on shifted problem
        call EPSSetTrueResidual(EVP%eps, true_res, ierr)

        ! !CANNOT USE WITH ST
        ! !EPS_RITZ: standard projection (default)
        ! !EPS_HARMONIC: harmonic projection: to stabilize the eigenvalues close to the target first
        ! call EPSSetExtraction(EVP%eps,EPS_HARMONIC,ierr)

        ! !Set KSP for operator to vector operation
        ! call STGetKSP(st,EVP%ksp,ierr)
        ! call KSPSetType(EVP%ksp,KSPPREONLY,ierr)
        ! !Preconditioner
        ! call KSPGetPC(EVP%ksp,pc,ierr)
        ! call PCSetType(pc,PCLU,ierr)
        ! !MATSOLVERPASTIX, MATSOLVERMUMPS, MATSOLVERMKL_CPARDISO
        ! call PCFactorSetMatSolverType(pc,MATSOLVERMUMPS,ierr)

        ! !Set tolerances, usually stricter compared to EVP problem
        ! call KSPSetTolerances(EVP%ksp, rtol, abs_tol, div_tol, maxit, ierr)

        ! !Balancing to reduce the matrix norm
        ! !EPS_BALANCE_NONE, EPS_BALANCE_ONESIDE, EPS_BALANCE_TWOSIDE, or EPS_BALANCE_USER
        ! maxit = 100
        ! !cutoof = PETSC_DEFAULT_REAL
        ! call EPSSetBalance(EVP%eps, EPS_BALANCE_ONESIDE, maxit, PETSC_DEFAULT_REAL, ierr)

        !This is the last command to set up all the previously stated parameters
        call KSPSetFromOptions(EVP%ksp,ierr)

        !Set solver parameters at runtime
        call EPSSetFromOptions(EVP%eps,ierr)

        call EPSSolve(EVP%eps,ierr)

        call PetscLogDefaultBegin(ierr)
        call PetscViewerASCIIOpen(PETSC_COMM_WORLD, &
            'Results/Base/Info/log_info.dat',viewer,ierr)
        call PetscLogView(viewer,ierr)
        call PetscViewerDestroy(viewer, ierr)

    end subroutine solveEVP

    !----------------------------------------------------------------------

    function getEVPSolutionInfo(eps) result(nconv)

        use SlepcParameters_mod, only: nev

        implicit none

        EPS, intent(in) :: eps
        PetscInt :: nconv

        PetscBool :: terse
        character(len=100) :: str
        PetscErrorCode :: ierr

        write(str,'(a,g,a)') "Number of requested eigenvalues: ", nev, "\n"
        call PetscPrintf(PETSC_COMM_WORLD, str, ierr)

        call EPSGetConverged(eps, nconv, ierr)
        write(str,'(a,g,a)') "Number of computed eigenvalues: ", nconv, "\n"
        call PetscPrintf(PETSC_COMM_WORLD, str, ierr)

        call PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-terse", terse, ierr)
        if (terse) then
            call EPSErrorView(eps, EPS_ERROR_RELATIVE, PETSC_VIEWER_STDOUT_WORLD, ierr)
        else
            call PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO_DETAIL, ierr)
            call EPSConvergedReasonView(eps, PETSC_VIEWER_STDOUT_WORLD, ierr)
            ! call EPSErrorView(eps, EPS_ERROR_RELATIVE, PETSC_VIEWER_STDOUT_WORLD, ierr)
            call EPSErrorView(eps, EPS_ERROR_ABSOLUTE, PETSC_VIEWER_STDOUT_WORLD, ierr)
            call PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD, ierr)
        end if

    end function getEVPSolutionInfo

    !--------------------------------------------------------------------------

    subroutine checkCriticalConditions(EVP)

        use ContinuationVariables_mod, only: Cvar1
        use Tools_mod, only: interpolationLinear

        implicit none

        type(EVPType), intent(inout) :: EVP
        PetscBool :: critical

        PetscReal, parameter :: eps = 1.0d-8
        PetscInt, save :: Npos = 0
        PetscReal :: Kr_Lead, Kr_Leado, Cvar1_o
        PetscBool :: check

        EVP%critical = .false.

        Kr_Lead = PetscRealPart(EVP%K_Lead)
        Kr_Leado = EVP%Kr_Leado

        check = (Kr_Lead*Kr_Leado < 0.0d0)
        check = check .or. ( Kr_Lead > 0.0d0 .and. Kr_Leado > 0.0d0 )
        if (check) Npos = Npos+1
        if (Npos == 1) then
            call startingSolution(EVP%K_Lead, 'write')
        end if

        if (checkVar2Change(Kr_Lead, Kr_Leado)) return
        if (checkQChange(Kr_Lead, Kr_Leado)) return

        EVP%critical = checkSignChange(Kr_Lead, Kr_Leado)

        if (EVP%critical) then

            Cvar1_o = Cvar1%p - Cvar1%dvar
            check = (Kr_Lead < 0.0d0) .and. (Kr_Leado > 0.0d0)
            if (check) then
                Cvar1_o = Cvar1%p + Cvar1%dvar
            end if

            EVP%xCr = interpolationLinear(Cvar1_o, Cvar1%p, Kr_Leado, Kr_Lead, 0.0d0)

            ! write(*,'(*(es12.4))') Cvar1_o, Cvar1%p, Kr_Leado, Kr_Lead

        end if

    end subroutine checkCriticalConditions

    !--------------------------------------------------------------------------

    subroutine postProcessStability(Problem, Problem_Stability, FE, Mesh, &
        GaussInt, Sol, Elem, Kr_Lead)

        use MPIParameters_mod, only: Rank
        use PostProcessStability_mod, only: writeFinalSolution, energyAnalysis, &
            writeLeadingEigenvectorSolution

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        type(NodeArrays), intent(inout) :: Elem
        PetscReal, intent(in) :: Kr_Lead

        PetscBool :: check

        ! check = (Rank == 0) .and. (Kr_Lead > 0.0d0)
        ! check = check .and. (Problem%Ndim == Problem_Stability%Ndim)
        ! if (check) then
        !     call writeFinalSolution(Problem_Stability, FE, Mesh, Sol, Kr_Lead)
        !     call writeLeadingEigenvectorSolution(FE, Connectivity, Sol, Cvar1%p)
        ! end if

        call energyAnalysis(Problem, Problem_Stability, FE, Mesh, &
                GaussInt, Sol, Elem, Kr_Lead)

    end subroutine postProcessStability

    !--------------------------------------------------------------------------

    subroutine writeStabilityInfo(Problem_Stability_name, EVP)

        use Tools_mod, only: sortFile
        use PhysicalParameters_mod, only: QbN
        use ContinuationVariables_mod, only: Cvar1, Cvar2

        implicit none

        character(*), intent(in) :: Problem_Stability_name
        type(EVPType), intent(inout) :: EVP

        PetscInt :: Ncols_sort
        PetscReal :: Kr_Lead, Ki_Lead
        character(200) :: fn, folder
        character(len=:), allocatable :: val2, Q_val, val1

        Kr_Lead = PetscRealPart(EVP%K_Lead)
        Ki_Lead = PetscImaginaryPart(EVP%K_Lead)

        write(*,*)
        write(*,'(a)', advance='no') 'Leading eigenpair:'
        write(*,'(3f10.4,2f12.6,a)') Cvar2%p, QbN, Cvar1%p, Kr_Lead, Ki_Lead, ' i'
        write(*,*)

        write(fn,'(f12.4)') Cvar2%p
        val2 = trim(adjustl(fn))

        write(fn,'(f12.4)') QbN
        Q_val = trim(adjustl(fn))
        
        write(fn,'(f12.4)') Cvar1%p
        val1 = trim(adjustl(fn))

        folder = 'Results/Stability/'//trim(adjustl(Cvar2%name))//'/' &
            //trim(adjustl(Cvar2%name))//'_'//trim(adjustl(val2))
        fn = trim(adjustl(folder))//'/Kr_LwN_'//trim(adjustl(Cvar1%name))//'_'//val1//'.dat'
        open(20, file=fn, position='append')
        write(20,'(f12.4,3f16.8)') QbN, Kr_Lead, Ki_Lead, EVP%LwN_Lead
        close(20)

        if (EVP%critical) then

            fn = trim(adjustl(folder))//'/Critical_conditions.dat'
            open(20, file=fn, position='append')
            write(20,'(f12.4,2f16.8)') QbN, EVP%xCr, EVP%LwN_Lead
            close(20)

            write(*,2) trim(adjustl(Cvar2%name)), Cvar2%p, QbN, trim(adjustl(Cvar1%name)), EVP%xCr
            2 format ('Critical conditions at ',a,' = ',f8.4,', QbN = ',f8.4,', ',a,' = ',f8.4)
            write(*,*)

        end if

        folder = trim(adjustl(folder))//'/QbN_'//trim(adjustl(Q_val))//'/Eigenvalues/'
        fn = trim(adjustl(folder))//'Leading_eigenvalues.dat'
        open(16, file=fn, position='append')
        write(16,'(f12.4,4f16.8)') Cvar1%p, Kr_Lead, Ki_Lead, EVP%LwN_Lead, EVP%target_real
        close(16)

        ! Ncols_sort = 5
        ! call sortFile(fn, Ncols_sort)

        deallocate(val2) ; deallocate(Q_val) ; deallocate(val1)

    end subroutine writeStabilityInfo

    !--------------------------------------------------------------------------

    subroutine startingSolution(K_Lead, read_write)

        use MPIParameters_mod, only: Rank
        use FEMParameters_mod, only: FE_Inflow, FE_Main
        use MeshParameters_mod, only: Mesh_Inflow, Mesh_Main
        use SolutionVariables_mod, only: Sol_Inflow, Sol_Main
        use ContinuationVariables_mod, only: Cvar1
        use PostProcess_mod, only: writeSolution, readSolution

        implicit none

        PetscScalar, intent(inout) :: K_Lead
        character(*), intent(in) :: read_write

        PetscReal, save :: Cvar1_Start = 0.0d0
        PetscScalar, save :: K_Lead_Start = 0.0d0
        PetscReal :: iter_dump
        character(50) :: fn, str
        character(len=:), allocatable :: val1

        if (read_write == 'write') then

            Cvar1_Start = Cvar1%p
            K_Lead_Start = K_Lead

            if (Rank /= 0) return

            ! write(*,'(a)') 'Writing starting solution!'

            write(str,'(f12.4)') Cvar1%p
            val1 = trim(adjustl(str))

            fn = 'SOL_'//val1//'_Inflow_Start'
            call writeSolution(FE_Inflow, Mesh_Inflow%Connectivity, Sol_Inflow, fn)

            fn = 'SOL_'//val1//'_Main_Start'
            call writeSolution(FE_Main, Mesh_Main%Connectivity, Sol_Main, fn)

        end if

        if (read_write == 'read') then

            Cvar1%p = Cvar1_Start
            K_Lead = K_Lead_Start

            ! if (Rank == 0) then
            !     write(*,'(a)') 'Reading starting solution!'
            ! end if

            write(str,'(f12.4)') Cvar1%p
            val1 = trim(adjustl(str))

            fn = 'SOL_'//val1//'_Inflow_Start'
            call readSolution(FE_Inflow, Mesh_Inflow%Connectivity, iter_dump, &
                                fn, Sol_Inflow)

            fn = 'SOL_'//val1//'_Main_Start'
            call readSolution(FE_Main, Mesh_Main%Connectivity, iter_dump, fn, Sol_Main)

        end if

    end subroutine startingSolution

    !--------------------------------------------------------------------------

    subroutine previousSolution(K_Lead, read_write)

        use MPIParameters_mod, only: Rank
        use FEMParameters_mod, only: FE_Inflow, FE_Main
        use MeshParameters_mod, only: Mesh_Inflow, Mesh_Main
        use SolutionVariables_mod, only: Sol_Inflow, Sol_Main
        use ContinuationVariables_mod, only: Cvar1
        use PostProcess_mod, only: writeSolution, readSolution

        implicit none

        PetscScalar, intent(inout) :: K_Lead
        character(*), intent(in) :: read_write

        PetscReal, save :: Cvar1_Previous = 0.0d0
        PetscScalar, save :: K_Lead_Previous = 0.0d0
        PetscReal :: iter_dump
        character(50) :: fn, str, cmd
        character(len=:), allocatable :: val1

        if (read_write == 'write') then

            Cvar1_Previous = Cvar1%p
            K_Lead_Previous = K_Lead

            if (Rank /= 0) return

            cmd = 'rm -f *Previous.DTA*'
            call execute_command_line(cmd)

            ! write(*,'(a)') 'Writing previous solution!'

            write(str,'(f12.4)') Cvar1%p
            val1 = trim(adjustl(str))

            fn = 'SOL_'//val1//'_Inflow_Previous'
            call writeSolution(FE_Inflow, Mesh_Inflow%Connectivity, Sol_Inflow, fn)

            fn = 'SOL_'//val1//'_Main_Previous'
            call writeSolution(FE_Main, Mesh_Main%Connectivity, Sol_Main, fn)

        end if

        if (read_write == 'read') then

            Cvar1%p = Cvar1_Previous
            K_Lead = K_Lead_Previous

            ! if (Rank == 0) then
            !     write(*,'(a)') 'Reading previous solution!'
            ! end if

            write(str,'(f12.4)') Cvar1%p
            val1 = trim(adjustl(str))

            fn = 'SOL_'//val1//'_Inflow_Previous'
            call readSolution(FE_Inflow, Mesh_Inflow%Connectivity, iter_dump, &
                                fn, Sol_Inflow)

            fn = 'SOL_'//val1//'_Main_Previous'
            call readSolution(FE_Main, Mesh_Main%Connectivity, iter_dump, fn, Sol_Main)

        end if

    end subroutine previousSolution

    !--------------------------------------------------------------------------

    function checkVar2Change(Kr_Lead, Kr_Leado)

        use ContinuationVariables_mod, only: Cvar1, Cvar2
        use MPIParameters_mod, only: Rank

        implicit none

        PetscReal, intent(in) :: Kr_Lead, Kr_Leado
        PetscBool :: checkVar2Change

        PetscReal, parameter :: eps = 1.0d-8
        PetscReal :: dvar2
        PetscBool :: check_ok, check_var1

        checkVar2Change = .false.

        Cvar1%dvar = abs(Cvar1%dvar)

        dvar2 = abs(Cvar2%p-Cvar2%var_o)

        Cvar2%var_o = Cvar2%p

        check_ok = ( (Kr_Lead < 0.0d0) .and. (Kr_Leado > 0.0d0) .and. (dvar2 > eps) )
        if (check_ok) then

            if (Rank == 0) then
                write(*,'(4a)') 'Eigenvalue < 0 when ', trim(adjustl(Cvar2%name)), &
                    ' changed. Increasing ', trim(adjustl(Cvar1%name))
            end if

            checkVar2Change = .true.
            return

        end if

        check_var1 = ( (Kr_Lead > 0.0d0) .and. (Kr_Leado < 0.0d0) .and. (dvar2 > eps) )
        if (check_var1) then

            if (Rank == 0) then
                write(*,'(4a)') 'Eigenvalue > 0 when ', trim(adjustl(Cvar2%name)), &
                    ' changed. Decreasing ', trim(adjustl(Cvar1%name))
            end if

            Cvar1%dvar = -Cvar1%dvar

            checkVar2Change = .true.
            return

        end if

    end function checkVar2Change

    !--------------------------------------------------------------------------

    function checkQChange(Kr_Lead, Kr_Leado)

        use ContinuationVariables_mod, only: Cvar1
        use PhysicalParameters_mod, only: QbN0, QbN
        use MPIParameters_mod, only: Rank

        implicit none

        PetscReal, intent(in) :: Kr_Lead, Kr_Leado
        PetscBool :: checkQChange

        PetscReal, parameter :: eps = 1.0d-8
        PetscReal, save :: QbNo = QbN0
        PetscReal :: dQ
        PetscBool :: check_ok, check_var1

        checkQChange = .false.

        Cvar1%dvar = abs(Cvar1%dvar)

        dQ = abs(QbN-QbNo)

        QbNo = QbN

        check_ok = ( (Kr_Lead < 0.0d0) .and. (Kr_Leado > 0.0d0) .and. (dQ > eps) )
        if (check_ok) then

            if (Rank == 0) then
                write(*,'(2a)') 'Eigenvalue < 0 when QbN changed. Increasing ', &
                    trim(adjustl(Cvar1%name))
            end if

            checkQChange = .true.
            return

        end if

        check_var1 = ( (Kr_Lead > 0.0d0) .and. (Kr_Leado < 0.0d0) .and. (dQ > eps) )
        if (check_var1) then

            if (Rank == 0) then
                write(*,'(2a)') 'Eigenvalue > 0 when QbN changed. Decreasing ', &
                    trim(adjustl(Cvar1%name))
            end if

            Cvar1%dvar = -Cvar1%dvar

            checkQChange = .true.
            return

        end if

    end function checkQChange

    !--------------------------------------------------------------------------

    function checkSignChange(Kr_Lead, Kr_Leado)

        use ContinuationVariables_mod, only: Cvar1
        use MPIParameters_mod, only: Rank

        implicit none

        PetscReal, intent(in) :: Kr_Lead, Kr_Leado
        PetscBool :: checkSignChange

        PetscBool :: check_ok, check_var1

        checkSignChange = .true.

        Cvar1%dvar = abs(Cvar1%dvar)

        check_ok = (Kr_Lead < 0.0d0) .and. (Kr_Leado < 0.0d0)
        if (check_ok) then

            if (Rank == 0) then
                write(*,'(2a)') 'Consecutive Eigenvalues < 0. Increasing ', &
                    trim(adjustl(Cvar1%name))
            end if

            checkSignChange = .false.
            return

        end if

        check_var1 = (Kr_Lead > 0.0d0) .and. (Kr_Leado > 0.0d0)
        if (check_var1) then

            if (Rank == 0) then
                write(*,'(2a)') 'Consecutive Eigenvalues > 0. Decreasing ', &
                    trim(adjustl(Cvar1%name))
            end if

            Cvar1%dvar = -Cvar1%dvar

            checkSignChange = .false.
            return

        end if

    end function checkSignChange

end module Stability_mod
