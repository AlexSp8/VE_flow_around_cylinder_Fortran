
#include <petsc/finclude/petscksp.h>

module PostProcessStability_mod

    use slepceps
    use PhysicalParameters_mod, only: ProblemParameters
    use FEMParameters_mod, only: FEParameters
    use MeshParameters_mod, only: MeshParameters
    use BoundaryParameters_mod, only: BoundaryParameters
    use GaussParameters_mod, only: GaussIntegration
    use SolutionVariables_mod, only: SolutionArraysType
    use ElementVariables_mod, only: NodeArrays
    use SlepcVariables_mod, only: EVPType

    contains

    subroutine writeEigenpairs(Problem_Stability, FE, Mesh, &
        EVP, nconv, ishift, ilead)

        use PhysicalParameters_mod, only: QbN, Nbd_periodic
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Cvar1, Cvar2

        implicit none

        type(ProblemParameters), intent(in) :: Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(EVPType), intent(inout) :: EVP
        PetscInt, intent(in)  :: nconv, ishift
        PetscInt, intent(out) :: ilead

        PetscInt :: i
        PetscReal :: eigen_err, Kr, Ki, Kr_Lead
        PetscScalar :: K_eig
        PetscErrorCode :: ierr
        character(200) :: fn, folder
        character(len=:), allocatable :: val2, Q_val, val1

        if (Rank == 0) then

            write(fn,'(f12.4)') Cvar2%p
            val2 = trim(adjustl(fn))

            write(fn,'(f12.4)') QbN
            Q_val = trim(adjustl(fn))

            folder = 'Results/Stability/'//trim(adjustl(Cvar2%name))//'/' &
                //trim(adjustl(Cvar2%name))//'_'//trim(adjustl(val2)) &
                //'/QbN_'//trim(adjustl(Q_val))//'/Eigenvalues/'

            write(fn,'(f10.4)') Cvar1%p
            val1 = trim(adjustl(fn))

            fn = trim(adjustl(folder))//trim(adjustl(Cvar1%name))//'_'//val1//'.dat'

            open(2, file=fn, position='append')

            deallocate(val1) ; deallocate(val2) ; deallocate(Q_val)

        end if

        ilead = 0
        do i = 0, nconv-1

            call EPSGetEigenvalue(EVP%eps, i, K_eig, PETSC_NULL_SCALAR, ierr)

            call EPSComputeError(EVP%eps, i, EPS_ERROR_ABSOLUTE, eigen_err, ierr)

            if (abs(eigen_err) > 1.0d-6) cycle

            Kr = PetscRealPart(K_eig)
            Ki = PetscImaginaryPart(K_eig)

            if (Rank == 0) then
                write(2,'(2f24.12)') Kr, Ki
            end if

            if (Kr > 0.0d0) then
                call plotEigenvector(Problem_Stability, FE, Mesh, ishift, i, EVP)
            end if

            Kr_Lead = PetscRealPart(EVP%K_Lead)
            if (Kr > Kr_Lead) then
                EVP%K_Lead = K_eig
                ilead = i
            end if

            if (abs(Ki) > EVP%Ki_max) then
                EVP%Ki_max = abs(Ki)
            end if

        end do

        if (Rank == 0) then
            write(2,*)
            close(2)
        end if

    end subroutine writeEigenpairs

    !-----------------------------------------------------------------------

    subroutine plotEigenvector(Problem_Stability, FE, Mesh, ishift, i, EVP)

        use PhysicalParameters_mod, only: QbN, Nbd_periodic
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Cvar1, Increment, Cvar2
        use PostProcess_mod, only: writeTecplotFile

        implicit none

        type(ProblemParameters), intent(in) :: Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        PetscInt, intent(in) :: ishift, i
        type(EVPType), intent(inout) :: EVP

        PetscInt :: K_eig_ID
        PetscScalar :: K_eig
        PetscReal :: Kr, LwN, Kr_Lead
        PetscReal, dimension(FE%Nnodes,Problem_Stability%Neq) :: TL_dum
        PetscScalar, dimension(FE%Nnodes,Problem_Stability%Neq) :: TL_d
        PetscScalar, dimension(FE%Nnodes_Rank,Problem_Stability%Neq) :: TL_d_rank
        PetscErrorCode :: ierr
        character(200) :: fn, folder
        character(len=:), allocatable :: val2, Q_val, val1, K_val

        call EPSGetEigenvalue(EVP%eps, i, K_eig, PETSC_NULL_SCALAR, ierr)

        call storeEigenvector(Problem_Stability%Neq, FE, &
                    Mesh%Connectivity_Rank, TL_d, TL_d_rank, i, EVP)

        if (Rank /= 0) return

        Kr = PetscRealPart(K_eig)

        write(fn,'(f12.4)') Cvar2%p
        val2 = trim(adjustl(fn))

        write(fn,'(f12.4)') QbN
        Q_val = trim(adjustl(fn))

        write(fn,'(f12.4)') Cvar1%p
        val1 = trim(adjustl(fn))

        folder = 'Results/Stability/'//trim(adjustl(Cvar2%name))//'/' &
            //trim(adjustl(Cvar2%name))//'_'//trim(adjustl(val2)) &
            //'/QbN_'//trim(adjustl(Q_val))//'/Eigenvectors/' &
            //trim(adjustl(Cvar1%name))//'_'//val1//'/'

        call execute_command_line('mkdir -p '//folder)
        call execute_command_line('cp -r ConvertToBinary '//folder)

        write(fn,'(f16.8)') Kr
        K_val = trim(adjustl(fn))

        fn = trim(adjustl(folder))//'Kr_'//K_val//'.PLT'

        deallocate(val1) ; deallocate(val2) ; deallocate(Q_val) ; deallocate(K_val)

        TL_dum(:,:) = 0.0d0

        call writeTecplotFile(Problem_Stability, FE, Mesh, PetscRealPart(TL_d), &
            TL_dum, fn, real(Increment,8))

        if (Nbd_periodic > 0) then
            K_eig_ID = 1000*ishift+i+1
            LwN = longwaveDisturbance(Problem_Stability, Mesh%Xi_Mesh, TL_d, Kr, K_eig_ID)
            Kr_Lead = PetscRealPart(EVP%K_Lead)
            if (Kr > Kr_Lead) then
                EVP%LwN_Lead = LwN
            end if
        end if

    end subroutine plotEigenvector

    !-----------------------------------------------------------------------

    subroutine storeEigenvector(Neq, FE, Connectivity_Rank, TL_d, TL_d_rank, iEig, EVP)

        use MPIParameters_mod, only: Rank
        use Tools_mod, only: getRows, getRankNode

        implicit none

        PetscInt, intent(in) :: Neq
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:,:), intent(in) :: Connectivity_Rank
        PetscScalar, dimension(:,:), intent(out) :: TL_d, TL_d_rank
        PetscInt, intent(in) :: iEig
        type(EVPType), intent(in) :: EVP

        PetscInt :: ieq, inod, irow, iel_rank, rnod
        PetscInt, dimension(Neq) :: ieqs, irows_eqs
        PetscInt, dimension(1) :: gnod
        PetscScalar, pointer :: Xr_scat_p(:)
        Vec :: Xr_scat
        VecScatter :: ctx
        PetscErrorCode :: ierr

        call EPSGetEigenvector(EVP%eps, iEig, EVP%Xr, PETSC_NULL_VEC, ierr)

        call VecScatterCreateToAll(EVP%Xr, ctx, Xr_scat, ierr)
        call VecScatterBegin(ctx, EVP%Xr, Xr_scat, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterEnd(ctx, EVP%Xr, Xr_scat, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecGetArrayReadF90(Xr_scat, Xr_scat_p, ierr)

        if (Rank == 0) then
            irow = 0
            do inod = 1, FE%Nnodes
                do ieq = 1, Neq
                    irow = irow+1
                    ! TL_d(inod,ieq) = PetscRealPart(Xr_scat_p(irow))
                    TL_d(inod,ieq) = Xr_scat_p(irow)
                end do
            end do
        end if

        ieqs = [ (ieq, ieq = 1, Neq) ]
        do iel_rank = 1, FE%Nel_Rank
            do inod = 1, FE%Nbf
                gnod(1) = Connectivity_Rank(iel_rank,inod)
                rnod = getRankNode(iel_rank,inod,FE%Nbf)
                irows_eqs(:) = getRows(gnod(:),ieqs(:),Neq)
                TL_d_rank(rnod,:) = Xr_scat_p(irows_eqs(:))
            end do
        end do

        call VecRestoreArrayReadF90(Xr_scat, Xr_scat_p, ierr)
        call VecScatterDestroy(ctx, ierr)
        call VecDestroy(Xr_scat, ierr)

    end subroutine storeEigenvector

    !--------------------------------------------------------------------------

    function longwaveDisturbance(Problem_Stability, Xi_Mesh, TL_d, Kr, K_eig_ID) result(LwN)

        use petscksp, only: PETSC_i
        use PhysicalParameters_mod, only: QbN, Length, Height, Width, Periodicity, PI
        use ContinuationVariables_mod, only: Cvar1, Cvar2
        use Tools_mod, only: arrayClosestValue, sortFile

        implicit none

        type(ProblemParameters), intent(in) :: Problem_Stability
        PetscReal, dimension(:,:), intent(in) :: Xi_Mesh
        PetscScalar, dimension(:,:), intent(in) :: TL_d
        PetscReal, intent(in) :: Kr
        PetscInt, intent(in) :: K_eig_ID
        PetscReal :: LwN

        PetscInt :: i, j, k, inod, Neq_f, Ncells, Nnodes, Ncols_sort, p_int
        PetscReal :: Lc, x, p
        PetscReal, parameter, dimension(1) :: x_target = [1.1d0]
        PetscReal, parameter :: eps = 1.0d-8
        PetscReal, dimension(:), allocatable :: x_nodes
        PetscReal, dimension(3) :: xi
        PetscScalar, dimension(3) :: U_f, U_f0
        PetscReal, dimension(size(x_target)) :: x_close
        PetscBool :: check
        character(150) :: fn, folder
        character(len=:), allocatable :: val2, Q_val, val1, K_val

        write(fn,'(f12.4)') Cvar2%p
        val2 = trim(adjustl(fn))

        write(fn,'(f12.4)') QbN
        Q_val = trim(adjustl(fn))

        write(fn,'(f12.4)') Cvar1%p
        val1 = trim(adjustl(fn))

        folder = 'Results/Stability/'//trim(adjustl(Cvar2%name))//'/' &
            //trim(adjustl(Cvar2%name))//'_'//trim(adjustl(val2)) &
            //'/QbN_'//trim(adjustl(Q_val))//'/Line/' &
            //trim(adjustl(Cvar1%name))//'_'//val1//'/'

        call execute_command_line('mkdir -p '//folder)

        ! write(fn,'(f16.8)') Kr
        write(fn,'(i8)') K_eig_ID
        K_val = trim(adjustl(fn))

        fn = trim(adjustl(folder))//'Kr_'//K_val//'.dat'

        open(12, file=fn, position='append')

        deallocate(val1) ; deallocate(val2) ; deallocate(Q_val) ; deallocate(K_val)

        Neq_f = Problem_Stability%Neq_f
        Nnodes = size(Xi_Mesh,1)

        x_nodes = pack(Xi_Mesh(:,1), Xi_Mesh(:,2) > 0.0d0)
        do i = 1, size(x_target)
            x_close(i) = arrayClosestValue(x_nodes,x_target(i))
        end do
        deallocate(x_nodes)

        ! LwN = 1.0d0/(QbN+1.0d-12)
        ! if ( abs(QbN) < eps ) LwN = 1.0d0-eps
        ! Ncells = ceiling(LwN)

        Ncells = 1
        if (abs(QbN) > eps) then
            do
                p = Ncells*QbN
                p_int = nint(p)
                check = ( abs(p-p_int) < eps)
                if (check) then
                    exit
                end if
                Ncells = Ncells+1
            end do
        end if

        do i = 1, size(Periodicity)

            select case (Periodicity(i))
            case ('X')
                Lc = Length
            case ('Y')
                Lc = Height
            case ('Z')
                Lc = Width
            case default
                Lc = 0.0d0
            end select
            Lc = 2.0d0*Lc

            do k = 1, size(x_close)
                do inod = 1, Nnodes

                    xi(:) = Xi_Mesh(inod,:)

                    check = ( abs(xi(1)-x_close(k)) < eps ) .and. ( xi(2) > 0.0d0 )
                    if ( .not. check ) cycle

                    U_f0(:) = 0.0d0
                    U_f0(1:Neq_f) = TL_d(inod,1:Neq_f)

                    do j = 1, Ncells

                        x = xi(3) + (j-1)*Lc

                        U_f(:) = U_f0(:)*exp(2.0d0*PI*QbN*PETSC_i*(j-1))
                        write(12,'(*(es16.8))') x, PetscRealPart(U_f)

                    end do

                end do
            end do

        end do

        close(12)

        Ncols_sort = 4
        call sortFile(fn, Ncols_sort)

        LwN = wavelengthDataStability(Problem_Stability%name, folder, fn, Kr, K_eig_ID)

    end function longwaveDisturbance

    !--------------------------------------------------------------------------

    function wavelengthDataStability(name, folder, fn, Kr, K_eig_ID) result(LwN)

        use ContinuationVariables_mod, only: Cvar1, Cvar2
        use PhysicalParameters_mod, only: QbN, PI, Width
        use Tools_mod, only: interpolationLinear, sortFile

        implicit none

        character(*), intent(in) :: name, fn, folder
        PetscReal, intent(in) :: Kr
        PetscInt, intent(in) :: K_eig_ID
        PetscReal :: LwN

        PetscReal, parameter :: eps = 1.0d-6
        PetscInt :: ioerr, Nchanges, np, Ncols_sort
        PetscReal :: xo, xi, x0, x1, x2, LwN_max, Lc
        PetscReal :: U_max, U_min, U_mean, alpha, phi, C
        PetscBool :: check
        PetscReal, dimension(:), allocatable :: x
        PetscReal, dimension(:,:), allocatable :: U
        PetscReal, dimension(3) :: U_f, U_fo
        PetscReal, dimension(4) :: Xparam
        character(100) :: fn1

        Lc = 2.0d0*Width

        ! LwN_max = Lc/(QbN+1.0d-12)
        ! if ( abs(QbN) < eps ) LwN_max = Lc

        open(20, file=fn, action='read', iostat=ioerr, position='rewind')

        np = 0
        Nchanges = 0
        U_fo(:) = 0.0d0
        xo = 0.0d0
        x1 = 0.0d0
        x2 = 0.0d0
        U_max = -huge(1.0d0)
        U_min = huge(1.0d0)
        U_mean = 0.0d0

        read(20,*,iostat=ioerr) x0, U_f
        rewind(20)
        do

            read(20,*,iostat=ioerr) xi, U_f(1), U_f(2), U_f(3)

            check = (ioerr /= 0)
            ! check = check .or. (xi > x0+LwN_max+eps)
            if (check) exit

            if (U_f(3) > U_max) U_max = U_f(3)
            if (U_f(3) < U_min) U_min = U_f(3)

            U_mean = U_mean + U_f(3)

            np = np+1

            if (U_fo(3)*U_f(3) < 0.0d0) then
                Nchanges = Nchanges+1
                if (Nchanges == 1) then
                    x1 = interpolationLinear(xo, xi, U_fo(3), U_f(3), 0.0d0)
                end if
                if (Nchanges == 2) then
                    x2 = interpolationLinear(xo, xi, U_fo(3), U_f(3), 0.0d0)
                end if
            end if

            xo = xi
            U_fo(:) = U_f(:)

            LwN_max = xi-x0

        end do

        if (Nchanges == 1) then
            x2 = x0+LwN_max
        end if

        U_mean = U_mean/max(1,np)

        rewind(20)

        allocate(x(np)) ; allocate(U(np,3))

        do np = 1, size(x)
            read(20,*) x(np), U(np,1), U(np,2), U(np,3)
        end do

        close(20, status='keep')

        LwN = 2.0d0*(x2-x1)
        ! LwN = Lc/(Nwaves+QbN*Lc)
        phi = PI/LwN - 2.0d0*PI*x1/LwN
        alpha = (U_max-U_min)/2.0d0
        C = U_mean

        Xparam(:) = [LwN, phi, alpha, C]

        call fitSine(x, U(:,3), Xparam)

        deallocate(x)
        deallocate(U)

        LwN = Xparam(1)

        fn1 = trim(adjustl(folder))//'Wave_'//trim(adjustl(name))//'.dat'
        open(15, file=fn1, position='append')
        write(15,'(es16.8,i8,*(es16.8))') Kr, K_eig_ID, LwN_max/LwN, LwN, phi/PI, alpha, C
        close(15, status='keep')

        Ncols_sort = 7
        call sortFile(fn1, Ncols_sort)

    end function wavelengthDataStability

    !--------------------------------------------------------------------------

    subroutine fitSine(x, y, Xparam)

        use Tools_mod, only: gaussElimination

        implicit none

        PetscReal, dimension(:), intent(in) :: x, y
        PetscReal, dimension(:), intent(inout) :: Xparam

        PetscInt, parameter :: max_iter = 1000
        PetscReal, parameter :: pi = 4.0d0*atan(1.0d0), tol = 1.0d-10
        PetscInt :: iter, i, N, M
        PetscReal :: JTR_norm, R_norm, Cor_norm, relax
        PetscReal :: A, L, phi, C, theta, yfit
        PetscBool :: check
        PetscReal, dimension(size(y)) :: Res
        PetscReal, dimension(size(y),size(Xparam)) :: Jac
        PetscReal, dimension(size(Xparam),size(y)) :: JT
        PetscReal, dimension(size(Xparam)) :: JTR, dx
        PetscReal, dimension(size(Xparam),size(Xparam)) :: JTJ

        ! open(24,file='sine.dat',position='append')
        ! write(24,'(a)') '---------------------------------------------------'
        ! write(24,'(*(f16.10))') Xparam

        ! open(23,file='regression.dat', position='append')

        N = size(Xparam)
        M = size(x)

        relax = 1.0d0

        do iter = 1, max_iter

            L = Xparam(1)
            phi = Xparam(2)
            A = Xparam(3)
            C = Xparam(4)

            do i = 1, size(x)

                theta = 2.0d0*pi*x(i)/L + phi

                yfit = A*sin(theta) + C

                Res(i) = yfit - y(i)

                Jac(i,1) = -A*cos(theta)*(2.0d0*pi*x(i)/(L**2)) !d/dL
                Jac(i,2) = A*cos(theta)                         !d/dphi
                Jac(i,3) = sin(theta)                           !d/dA
                Jac(i,4) = 1.0d0                                !d/dC

            end do

            JT(:,:) = transpose(Jac)
            JTR(:) = matmul(JT,Res)

            JTR_norm = norm2(JTR)
            R_norm = norm2(Res)

            check = (JTR_norm < tol) .and. (R_norm < tol)
            if (check) then
                ! write(*,'(a,i5)')'Non linear regression converged at i = ', iter
                ! write(*,'(a,es12.4)')'JTR norm:        ', JTR_norm
                ! write(*,'(a,es12.4)')'R norm:          ', R_norm
                exit
            end if

            JTJ(:,:) = matmul(JT,Jac)
            do i = 1, N
                JTJ(i,i) = JTJ(i,i) + JTR_norm
            end do

            dx(:) = gaussElimination(JTJ, -JTR)

            Cor_norm = norm2(dx)

            ! write(23,'(i5,*(es20.10))') iter, JTR_norm, R_norm, Cor_norm, L

            check = (JTR_norm < tol) .and. (Cor_norm < tol)
            if (check) then
                ! write(*,'(a,i5)')'Non linear regression converged at i = ', iter
                ! write(*,'(a,es12.4)')'JTR norm:        ', JTR_norm
                ! write(*,'(a,es12.4)')'Cor norm:        ', Cor_norm
                exit
            end if

            Xparam(:) = Xparam(:) + relax*dx(:)

        end do

        ! write(23,'(a)') '---------------------------------------------------'
        ! close(23)

        ! write(24,'(*(f16.10))') Xparam
        ! write(24,'(a)') '---------------------------------------------------'

        ! do i = 1, size(x)

        !     L = Xparam(1)
        !     phi = Xparam(2)
        !     A = Xparam(3)
        !     C = Xparam(4)

        !     yfit = A*sin(2.0d0*pi*x(i)/L + phi) + C

        !     write(24,'(3f16.10)') x(i), yfit, y(i)

        ! end do

        ! write(24,*)
        ! close(24)

    end subroutine fitSine

    ! ----------------------------------------------------------------------

    subroutine writeFinalSolution(Problem_Stability, FE, Mesh, Sol, Kr_Lead)

        use PhysicalParameters_mod, only: QbN
        use ContinuationVariables_mod, only: Cvar1, Cvar2
        use PostProcess_mod, only: writeTecplotFile

        implicit none

        type(ProblemParameters), intent(in) :: Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(SolutionArraysType), intent(in) :: Sol
        PetscReal, intent(in) :: Kr_Lead

        PetscInt :: i, ieq, k, factor, n
        PetscReal :: eig_norm
        PetscReal, dimension(11) :: iter_time
        character(200) :: str, fn, folder
        character(len=:), allocatable :: val2, Q_val, val1, t_val
        PetscReal, dimension(FE%Nnodes,Problem_Stability%Neq) :: TL, TL_p, TL_dum

        iter_time(:) = [(i-1, i = 1, size(iter_time))]

        ! eig_norm = 0.0d0
        ! do ieq = 1, Problem_Main_Stability%Neq
        !     eig_norm = eig_norm + dot_product(Sol%TL_d(:,ieq),Sol%TL_d(:,ieq))
        ! end do
        ! eig_norm = sqrt(eig_norm)

        factor = 1/Kr_Lead
        n = 0
        do while(factor > 0)
            factor = factor/10
            n = n+1
        end do
        factor = 10**(n-1)
        iter_time(:) = factor*iter_time(:)


        write(str,'(f12.4)') Cvar2%p
        val2 = trim(adjustl(str))

        write(str,'(f12.4)') QbN
        Q_val = trim(adjustl(str))

        write(str,'(f12.4)') Cvar1%p
        val1 = trim(adjustl(str))

        folder = 'Results/Stability/'//trim(adjustl(Cvar2%name))//'/' &
            //trim(adjustl(Cvar2%name))//'_'//trim(adjustl(val2)) &
            //'/QbN_'//trim(adjustl(Q_val))//'/Final_solution/' &
            //trim(adjustl(Cvar1%name))//'_'//val1//'/'

        call execute_command_line('mkdir -p '//folder)
        call execute_command_line('cp -r ConvertToBinary '//folder)

        deallocate(val1) ; deallocate(val2) ; deallocate(Q_val)

        TL_dum(:,:) = 0.0d0

        do k = 1, size(iter_time)

            TL_p(:,:) = Sol%TL_d(:,:)*exp(Kr_Lead*iter_time(k))!/eig_norm
            TL_p(:,:) = PetscRealPart(TL_p)

            TL(:,:) = Sol%TL_b(:,:) + TL_p(:,:)

            write(str,'(f12.4)') iter_time(k)
            t_val = trim(adjustl(str))

            fn = trim(adjustl(folder))//'Time_'//t_val//'.PLT'

            deallocate(t_val)

            call writeTecplotFile(Problem_Stability, FE, Mesh, &
                TL, TL_dum, fn, iter_time(k))

        end do

    end subroutine writeFinalSolution

    !--------------------------------------------------------------------------

    subroutine energyAnalysis(Problem, Problem_Stability, FE, Mesh, &
        GaussInt, Sol, Elem, Kr_Lead)

        use PhysicalParameters_mod, only: BvN, WiN, ReN, Stress_Reform, QbN
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Cvar2
        use ElementVariables_mod, only: GaussPointQuantities
        use ElementCalculations_mod, only: copyToLocal, setElementBasisFunctions, &
            basis, basisInflow_p, basisMain_p
        use ElementCalculationsStability_mod, only: setFlowQuantitiesStability
        use Tools_mod, only: doubleDotProduct

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        type(NodeArrays), intent(inout) :: Elem
        PetscReal, intent(in) :: Kr_Lead

        type(GaussPointQuantities) :: GsPt_b, GsPt_d
        PetscInt  :: ig, ieq, iel_rank, Nbf, j
        PetscReal :: phi_in, phi_pr, phi_vis, phi_ps1, phi_pu1, phi_pu2, phi_ps2, phi_rel
        PetscReal :: dKEdt, dEvedt, dEedt, dEvdt, dVDdt, phi_jump, dEdt, balance
        PetscReal :: term, eig_norm, wet
        PetscReal, dimension(3,3,3) :: dTdXi_d
        PetscReal, dimension(3,3) :: term2, Se_d, GU_d, GUT_d, G_dot_d, T_ve_d
        PetscReal, dimension(3,3) :: T_ve_b
        PetscBool :: check
        PetscErrorCode :: ierr
        procedure(basis), pointer :: basis_p
        character(200) :: folder, fn
        character(len=:), allocatable :: val2, Q_val

        select case (Sol%name)
        case ('Inflow')
            basis_p => basisInflow_p
        case ('Main')
            basis_p => basisMain_p
        case default
            write(*,'(a)') 'Wrong Sol%name in energyAnalysis!'
        end select

        if (Rank == 0) then
            eig_norm = 0.0d0
            do ieq = 1, Problem_Stability%Neq
                eig_norm = eig_norm + dot_product(Sol%TL_d(:,ieq),Sol%TL_d(:,ieq))
            end do
            eig_norm = sqrt(eig_norm)
            ! print*, eig_norm
        end if
        call MPI_Bcast(eig_norm, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)

        Nbf = FE%Nbf
        
        phi_in = 0.0d0 ; phi_pr = 0.0d0 ; phi_ps1 = 0.0d0 ; phi_pu1 = 0.0d0
        phi_ps2 = 0.0d0 ; phi_pu2 = 0.0d0 ; phi_rel = 0.0d0 ; phi_vis = 0.0d0
        dKEdt = 0.0d0 ; dEvedt = 0.0d0 ; dEedt = 0.0d0 ; dEvdt = 0.0d0 ; phi_jump = 0.0d0

        do iel_rank = 1, FE%Nel_Rank

            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)

            Elem%TEMP_TL_d(:,:) = Elem%TEMP_TL_d(:,:)/eig_norm

            check = checkBoundary(Elem%Xi_loc(:,1),Elem%Xi_loc(:,2))
            if (.not. check) cycle

            ! write(*,'(*(es12.4))') Xi_loc(:,1), Xi_loc(:,2)

            LOOP_GAUSS: do ig = 1, GaussInt%NGauss_b

                call setElementBasisFunctions(ig, GaussInt, Elem)

                call basis_p(Problem%idir, Nbf, Elem, GsPt_b)

                call setFlowQuantitiesStability(Problem, Problem_Stability, &
                    Nbf, Elem, GsPt_b, GsPt_d)

                wet = (GaussInt%Weights_b(ig))*abs(GsPt_b%detJ)

                GU_d(:,:) = PetscRealPart(GsPt_d%GU(:,:))
                GUT_d(:,:) = PetscRealPart(GsPt_d%GUT(:,:))
                G_dot_d(:,:) = GU_d(:,:) + GUT_d(:,:)
                dTdXi_d(:,:,:) = PetscRealPart(GsPt_d%dTdXi(:,:,:))

                T_ve_b(:,:) = GsPt_b%T_ve(:,:)
                T_ve_d(:,:) = PetscRealPart(GsPt_d%T_ve(:,:))

                if (Stress_Reform == 'SQRT') then
                    call calculateViscoelasticFunction(Problem%Ndim, &
                        Problem_Stability%Ndim, GsPt_b, GsPt_d)
                end if
            
                !dKE/dt
                dKEdt = dKEdt + ReN*Kr_Lead*dot_product(GsPt_d%U_f,GsPt_d%U_f)*wet

                !dEp/dt
                term = doubleDotProduct(T_ve_d,GU_d)
                dEvedt = dEvedt - WiN*Kr_Lead*term*wet/(GsPt_b%f(2))

                !inertia
                phi_in = phi_in - ReN*dot_product(GsPt_d%Convection,GsPt_d%U_f)*wet

                !pressure
                phi_pr = phi_pr - dot_product(GsPt_d%dPdXi,GsPt_d%U_f)*wet

                !ps1
                term2(:,:) = 0.0d0
                do j = 1, 3
                    term2(:,:) = term2(:,:) + (GsPt_b%U_f(j))*dTdXi_d(:,:,j)
                end do
                term = doubleDotProduct(term2,GU_d)
                phi_ps1 = phi_ps1 + WiN*term*wet/(GsPt_b%f(2))

                !pu1
                term2(:,:) = 0.0d0
                do j = 1, 3
                    term2(:,:) = term2(:,:) + (GsPt_d%U_f(j))*(GsPt_b%dTdXi(:,:,j))
                end do
                term = doubleDotProduct(term2,GU_d)
                phi_pu1 = phi_pu1 + WiN*term*wet/(GsPt_b%f(2))

                !ps2
                term2(:,:) = matmul(GsPt_b%GUT,T_ve_d) + matmul(T_ve_d,GsPt_b%GU)
                term = doubleDotProduct(term2,GU_d)
                phi_ps2 = phi_ps2 - WiN*term*wet/(GsPt_b%f(2))

                !pu2
                term2(:,:) = matmul(GUT_d,T_ve_b) + matmul(T_ve_b,GU_d)
                term = doubleDotProduct(term2,GU_d)
                phi_pu2 = phi_pu2 - WiN*term*wet/(GsPt_b%f(2))

                !relaxation
                term = doubleDotProduct(T_ve_b,GU_d)
                phi_rel = phi_rel + (GsPt_d%f(2))*term*wet/(GsPt_b%f(2))

                !viscous
                term = doubleDotProduct(G_dot_d,GU_d)
                term = term*( BvN + ( (GsPt_b%f(1))*(1.0d0-BvN)/(GsPt_b%f(2)) ) )
                phi_vis = phi_vis - term*wet
                term = doubleDotProduct(GsPt_b%G_dot,GU_d)
                term = term*( (GsPt_d%f(1))*(1.0d0-BvN)/(GsPt_b%f(2)) )
                phi_vis = phi_vis - term*wet

                !elastic: dEe/dt
                Se_d(:,:) = T_ve_d(:,:) - (1.0d0-BvN)*G_dot_d(:,:)
                term = doubleDotProduct(Se_d,GU_d)
                dEedt = dEedt - WiN*Kr_Lead*term*wet/(GsPt_b%f(2))

                !viscous: dEv/dt
                term2(:,:) = Kr_Lead*G_dot_d(:,:)
                ! term2(:,:) = Kr_Lead*GUT_d(:,:)
                term = doubleDotProduct(term2,GU_d)
                dEvdt = dEvdt - (1.0d0-BvN)*WiN*term*wet/(GsPt_b%f(2))

                !jump: dEv,G/dt

            end do LOOP_GAUSS

        end do

        call MPI_Allreduce(MPI_IN_PLACE,dKEdt,   1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,dEvedt,  1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,phi_in,  1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,phi_pr,  1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,phi_ps1, 1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,phi_pu1, 1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,phi_ps2, 1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,phi_pu2, 1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,phi_rel, 1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,phi_vis, 1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,dEedt,   1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,dEvdt,   1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,phi_jump,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

        ! dEvedt = -dKEdt + (phi_in+phi_pr+phi_ps1+phi_pu1+phi_ps2+phi_pu2+phi_rel+phi_vis)
        ! balance = dKEdt+dEvedt-(phi_in+phi_pr+phi_ps1+phi_pu1+phi_ps2+phi_pu2+phi_rel+phi_vis)

        dVDdt = -dEvdt
        dEdt = dKEdt+dEedt+dEvdt+phi_jump
        balance = dEdt-(phi_in+phi_pr+phi_ps1+phi_pu1+phi_ps2+phi_pu2+phi_rel+phi_vis)

        if (Rank /= 0) return

        write(fn,'(f12.4)') Cvar2%p
        val2 = trim(adjustl(fn))

        write(fn,'(f12.4)') QbN
        Q_val = trim(adjustl(fn))

        folder = 'Results/Stability/'//trim(adjustl(Cvar2%name))//'/'
        folder = trim(adjustl(folder))//trim(adjustl(Cvar2%name))//'_'//trim(adjustl(val2))
        folder = trim(adjustl(folder))//'/QbN_'//trim(adjustl(Q_val))//'/Stability_data/'

        fn = trim(adjustl(folder))//'energy_data.dat'
        open(16, file=fn, position='append')
        write(16,'(f12.4,*(es14.4))') WiN, dEvedt, phi_ps1, phi_pu1, phi_ps2, phi_pu2, phi_rel
        close(16)

        fn = trim(adjustl(folder))//'energy_data_2.dat'
        open(18, file=fn, position='append')
        write(18,'(*(es14.4))') phi_vis, dEvdt, dEedt, dKEdt, phi_in, phi_pr, phi_jump
        close(18)

        fn = trim(adjustl(folder))//'energy_data_3.dat'
        open(20, file=fn, position='append')
        write(20,'(*(es14.4))') dEdt, dVDdt, balance, dEdt-dEvedt
        close(20)

        deallocate(val2) ; deallocate(Q_val)

    end subroutine energyAnalysis

    !--------------------------------------------------------------------------

    subroutine calculateViscoelasticFunction(Ndim, Ndim_S, GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: Model_Name, WiN, BvN, VeN
        use ElementVariables_mod, only: GaussPointQuantities
        use Tools_mod, only: traceTensor
        use ConstitutiveModels_mod, only: setFENEMaterialDerivative
        use ConstitutiveModelsStability_mod, only: setFENEMaterialDerivativeStability

        implicit none

        PetscInt, intent(in) :: Ndim, Ndim_S
        type(GaussPointQuantities), intent(inout) :: GsPt_b, GsPt_d

        PetscReal :: trT_b, MDfDt_b
        PetscReal :: trT_d, MDfDt_d

        trT_b = traceTensor(GsPt_b%T_ve)
        trT_d = traceTensor(PetscRealPart(GsPt_d%T_ve))

        select case (Model_Name)

        case ('L-PTT')

            GsPt_b%f_T = 1.0d0 + (VeN*WiN*trT_b)/(1.0d0-BvN)
            GsPt_b%f(1) = 1.0d0
            GsPt_b%f(2) = GsPt_b%f_T
            GsPt_b%f(3) = 0.0d0

            GsPt_d%f_T = VeN*WiN*trT_d/(1.0d0-BvN)
            GsPt_d%f(1) = 0.0d0
            GsPt_d%f(2) = GsPt_d%f_T
            GsPt_d%f(3) = 0.0d0

        case ('m-L-PTT')

            GsPt_b%f_T = 1.0d0 + (VeN*WiN*trT_b)/(1.0d0-BvN)
            GsPt_b%f(1) = GsPt_b%f_T
            GsPt_b%f(2) = GsPt_b%f_T
            GsPt_b%f(3) = 0.0d0

            GsPt_d%f_T = VeN*WiN*trT_d/(1.0d0-BvN)
            GsPt_d%f(1) = GsPt_d%f_T
            GsPt_d%f(2) = GsPt_d%f_T
            GsPt_d%f(3) = 0.0d0

        case ('e-PTT')

            GsPt_b%f_T = exp(VeN*WiN*trT_b/(1.0d0-BvN))
            GsPt_b%f(1) = 1.0d0
            GsPt_b%f(2) = GsPt_b%f_T
            GsPt_b%f(3) = 0.0d0

            GsPt_d%f_T = (GsPt_b%f_T)*VeN*WiN*trT_d/(1.0d0-BvN)
            GsPt_d%f(1) = 0.0d0
            GsPt_d%f(2) = GsPt_d%f_T
            GsPt_d%f(3) = 0.0d0

        case ('m-e-PTT')

            GsPt_b%f_T = exp(VeN*WiN*trT_b/(1.0d0-BvN))
            GsPt_b%f(1) = GsPt_b%f_T
            GsPt_b%f(2) = GsPt_b%f_T
            GsPt_b%f(3) = 0.0d0

            GsPt_d%f_T = (GsPt_b%f_T)*VeN*WiN*trT_d/(1.0d0-BvN)
            GsPt_d%f(1) = GsPt_d%f_T
            GsPt_d%f(2) = GsPt_d%f_T
            GsPt_d%f(3) = 0.0d0

        case ('FENE-CR')

            GsPt_b%f_T = 1.0d0 + ((WiN*trT_b/(1.0d0-BvN)))/(VeN-Ndim)
            GsPt_b%f(1) = GsPt_b%f_T
            MDfDt_b = setFENEMaterialDerivative(GsPt_b)
            GsPt_b%f(2) = GsPt_b%f_T - (WiN*MDfDt_b/GsPt_b%f_T)
            GsPt_b%f(3) = 0.0d0

            GsPt_d%f_T = (WiN*trT_d/(1.0d0-BvN))/(VeN - Ndim_S)
            GsPt_d%f(1) = GsPt_d%f_T
            MDfDt_d = setFENEMaterialDerivativeStability(GsPt_b, GsPt_d)
            GsPt_d%f(2) = GsPt_d%f_T - WiN*MDfDt_d/(GsPt_b%f_T)
            GsPt_d%f(3) = 0.0d0

        case ('FENE-P')

            GsPt_b%f_T = 1.0d0 + ((WiN*trT_b/(1.0d0-BvN)))/VeN
            GsPt_b%f(1) = 1.0d0
            MDfDt_b = setFENEMaterialDerivative(GsPt_b)
            GsPt_b%f(2) = GsPt_b%f_T - (WiN*MDfDt_b/GsPt_b%f_T)
            GsPt_b%f(3) = 0.0d0

            GsPt_d%f_T = (WiN*trT_d/(1.0d0-BvN))/VeN
            GsPt_d%f(1) = 0.0d0
            MDfDt_d = setFENEMaterialDerivativeStability(GsPt_b, GsPt_d)
            GsPt_d%f(2) = GsPt_d%f_T - WiN*MDfDt_d/(GsPt_b%f_T)
            GsPt_d%f(3) = 0.0d0

        case ('GIESEKUS')
            write(*,'(a)') 'energyAnalysisSQRT not defined for Giesekus!'
        case default
            write(*,'(a)') 'Wrong Model_Name in energyAnalysisSQRT!'
        end select

    end subroutine calculateViscoelasticFunction

    !--------------------------------------------------------------------------

    logical function checkBoundary(x,y)

        implicit none

        PetscReal, dimension(:), intent(in) :: x, y

        checkBoundary = .true.!all( x(:) >= abs(y(:)) .and. x(:) <= 10.0d0 )

    end function checkBoundary

    ! ----------------------------------------------------------------------

    subroutine writeLeadingEigenvectorSolution(FE, Connectivity, Sol, value)

        implicit none

        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:,:), intent(in) :: Connectivity
        type(SolutionArraysType), intent(in) :: Sol
        PetscReal, intent(in) :: value

        PetscInt :: inod, ieq, ierror, iel, gnod
        character(len=50) :: str, fn
        character(len=:), allocatable :: str1

        write(str,'(f12.4)') value
        str1 = trim(adjustl(str))

        fn = 'Eigenvectors/SOL_Eigenvector_'//str1//'.DTA'
        fn = trim(adjustl(fn))

        deallocate(str1)

        open(12, file=fn, status='unknown', action='write', iostat=ierror, &
            position='rewind', form='unformatted')

        do ieq = 1, size(Sol%EU_d)
            write(12) Sol%EU_d(ieq)
        end do

        do iel = 1, FE%Nel
            do inod = 1, FE%Nbf
                gnod = Connectivity(iel,inod)
                do ieq = 1, size(Sol%TL_d,2)
                    write(12) Sol%TL_d(gnod,ieq)
                end do
            end do
        end do

        close(12, status='keep')

    end subroutine writeLeadingEigenvectorSolution

    !-----------------------------------------------------------------------

    subroutine readLeadingEigenvectorSolution(Problem, FE, Connectivity, Sol)

        use MPIParameters_mod, only: Rank, NRanks
        use Tools_mod, only: allocateArray, getRankNode

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:,:), intent(in) :: Connectivity
        type(SolutionArraysType), intent(inout) :: Sol

        PetscInt :: inod, ieq, ierr, iel_rank, iel, iel_start, rnod, gnod, Nel_Rank0
        PetscReal, dimension(FE%Nnodes, Problem%Neq) :: TL_d
        PetscReal, dimension(FE%Nnodes_Rank, Problem%Neq) :: TL_d_rank
        PetscReal, dimension(Problem%Nex) :: EU_d

        if (allocated(Sol%EU_d)) deallocate(Sol%EU_d)
        call allocateArray(Sol%EU_d, Problem%Nex)

        if (allocated(Sol%TL_d_rank)) deallocate(Sol%TL_d_rank)
        call allocateArray(Sol%TL_d_rank, FE%Nnodes_Rank, Problem%Neq)

        open(13, file='Disturbance.DTA', status='unknown', action='read', iostat=ierr, &
            position='rewind', form='unformatted')

        if (ierr /= 0) then
            return
        end if

        if (Rank == 0) then

            if (allocated(Sol%TL_d)) deallocate(Sol%TL_d)
            call allocateArray(Sol%TL_d, FE%Nnodes, Problem%Neq)

            do ieq = 1, size(Sol%EU_d)
                read(13) EU_d(ieq)
            end do

            do iel = 1, FE%Nel
                do inod = 1, FE%Nbf
                    gnod = Connectivity(iel,inod)
                    do ieq = 1, size(Sol%TL_d,2)
                        read(13) TL_d(gnod,ieq)
                    end do
                end do
            end do

            rewind(13)

            Sol%TL_d(:,:) = TL_d(:,:)

        end if

        do ieq = 1, size(Sol%EU_d)
            read(13) EU_d(ieq)
        end do

        Nel_Rank0 = FE%Nel/NRanks

        !Skip
        iel_start = Rank*Nel_Rank0+1
        do iel = 1, iel_start-1
            do inod = 1, FE%Nbf
                do ieq = 1, size(Sol%TL_d_rank,2)
                    read(13)
                end do
            end do            
        end do

        do iel_rank = 1, FE%Nel_Rank
            do inod = 1, FE%Nbf
                rnod = getRankNode(iel_rank,inod,FE%Nbf)
                do ieq = 1, size(Sol%TL_d_rank,2)
                    read(13) TL_d_rank(rnod,ieq)
                end do
            end do
        end do
        close(13, status='keep')

        Sol%TL_d_Rank(:,:) = TL_d_Rank(:,:)
        Sol%EU_d(:) = EU_d(:)

        call MPI_Barrier(MPI_COMM_WORLD,ierr)

    end subroutine readLeadingEigenvectorSolution

end module PostProcessStability_mod
