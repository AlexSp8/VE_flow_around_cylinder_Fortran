
#include <petsc/finclude/petscksp.h>

module PostProcess_mod

    use petscksp
    use PhysicalParameters_mod, only: ProblemParameters
    use FEMParameters_mod, only: FEParameters
    use MeshParameters_mod, only: MeshParameters
    use BoundaryParameters_mod, only: BoundaryParameters
    use GaussParameters_mod, only: GaussIntegration
    use SolutionVariables_mod, only: SolutionArraysType
    use LinearSystemVariables_mod, only: LinearSystemType
    use ElementVariables_mod, only: GaussPointQuantities, NodeArrays

    contains

    subroutine postProcess(Problem, FE, Mesh, Boundary, &
        GaussInt, Sol, Elem, LS_proj, twrite0)

        use ContinuationParameters_mod, only: Dx_write, Continuation_Method
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Increment, Cvar1, dSo, Cvar2

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(inout) :: Sol
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS_proj
        PetscReal, intent(inout) :: twrite0

        PetscReal, parameter :: eps = 1.0d-8
        PetscInt :: Nex
        PetscReal :: dl
        PetscBool :: check
        character(200) :: fn, str, folder, cmd
        character(len=:), allocatable :: val2, val1

        check = (Continuation_Method == 'Arclength')
        check = check .and. (Rank == 0)
        check = check .and. (Sol%name == 'Main')
        if (check) then
            Nex = size(Sol%EU)
            dl = Sol%EU(Nex)
            fn = 'Results/Base/Info/Arclength_data.dat'
            open(20, file=fn,position='append')
            write(20,'(i5, f24.12,2f16.6)') Increment, dl, Cvar1%p, dSo
            close(20, status='keep')
        end if

        check = (abs(Cvar1%p - twrite0) > Dx_write - eps)
        if (Cvar1%p > 10.0d0 + eps) then
            check = ( abs(Cvar1%p - twrite0) > 10.0d0*Dx_write - eps )
        end if
        if (Cvar1%p > 1000.0d0 + eps) then
            check = ( abs(Cvar1%p - twrite0) > 50.0d0*Dx_write - eps )
        end if
        if (Cvar1%p > 5000.0d0 + eps) then
            check = ( abs(Cvar1%p - twrite0) > 100.0d0*Dx_write - eps )
        end if
        if (Cvar1%p > 10000.0d0 + eps) then
            check = ( abs(Cvar1%p - twrite0) > 500.0d0*Dx_write - eps )
        end if

        if (check) then

            twrite0 = Cvar1%p

            if (Problem%Neq_proj > 0) then
                call projectQuantities(Problem, FE, Mesh, GaussInt, Sol, Elem, LS_proj)
            end if

            call calculatePostProcessQuantities(Problem, FE, Mesh, &
                Boundary, GaussInt, Sol, Elem)
            
            check = (Problem%Ndim > 0)
            check = check .and. (Rank == 0)
            if (.not. check) return

            write(str,'(f12.4)') Cvar2%p
            val2 = trim(adjustl(str))

            write(str,'(f12.4)') Cvar1%p
            val1 = trim(adjustl(str))

            !Tec
            folder = 'Results/Base/'//trim(adjustl(Cvar2%name))//'/' &
                //trim(adjustl(Cvar2%name))//'_'//trim(adjustl(val2))//'/'

            cmd = 'cp -r ConvertToBinary '//trim(adjustl(folder))//'Tec/'
            call execute_command_line(cmd)

            fn = trim(adjustl(folder))//'Tec/'//trim(adjustl(Problem%name))//'_' &
                    //trim(adjustl(Cvar1%name))//'_'//val1//'.PLT'
            fn = trim(adjustl(fn))

            call writeTecplotFile(Problem, FE, Mesh, Sol%TL, Sol%TL_proj, &
                fn, real(Increment,8))

            !Sol
            fn = trim(adjustl(folder))//'Sol/'//trim(adjustl(Problem%name))//'_' &
                    //trim(adjustl(Cvar1%name))//'_'//val1
            fn = trim(adjustl(fn))

            call writeSolution(FE, Mesh%Connectivity, Sol, fn)

            deallocate(val1) ; deallocate(val2)

        end if

    end subroutine postProcess

    !--------------------------------------------------------------------------

    subroutine writeTecplotFile(Problem, FE, Mesh, TL, TL_proj, fn, iter1)

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        PetscReal, dimension(:,:), intent(in) :: TL, TL_proj
        character(len=100), intent(in) :: fn
        PetscReal, intent(in) :: iter1

        PetscInt :: ierror, inod, iel, idim, Nel_v
        character(len=20) :: zonetype
        character(len=200) :: variables
        PetscBool :: check
        PetscInt, dimension(:,:), allocatable :: Connectivity_v

        open(14, file=fn, status='unknown', action='write', iostat=ierror, position='rewind')

        if (ierror /= 0) then
            write(*,'(a)')'Cannot open file in writeTecplotFile!'
            stop
        end if

        variables = setVariablesTecplot(Problem,FE%Ndim)

        zonetype = setZonetypeTecplot(FE%name)

        Connectivity_v = setConnectivityTecplot(FE, Mesh%Connectivity)

        Nel_v = size(Connectivity_v,1)
        
        write(14,'(a)')'TITLE="FLOW"'
        write(14,'(2a)') 'VARIABLES=', adjustl(trim(variables))
        write(14,'(a,es12.4,2(a,i0),2a)') 'ZONE T="iter=', iter1, '", DATAPACKING=BLOCK, N=', &
                    FE%Nnodes, ', E=', Nel_v, ', ZONETYPE=', zonetype
        write(14,*)

        do idim = 1, 3
            check = any (Problem%idir == idim)
            if (check) then
                write(14,'(f32.15)') (Mesh%Xi_Mesh(inod,idim), inod = 1, size(Mesh%Xi_Mesh,1))
            end if
        end do

        call writeTecplotUnknowns(Problem, TL, TL_proj)

        do iel = 1, Nel_v
            write(14,*)(Connectivity_v(iel,inod), inod = 1, size(Connectivity_v,2))
        end do

        deallocate(Connectivity_v)

        close(14, status='keep')

    end subroutine writeTecplotFile

    !-----------------------------------------------------------------------

    function setZonetypeTecplot(name) result(zonetype)

        implicit none

        character(*), intent(in) :: name
        character(20) :: zonetype

        select case(name)
        case('Line')
            zonetype = 'FELINESEG'
        case('Triangle')
            zonetype = 'FETRIANGLE'
        case('Quadrangle')
            zonetype = 'FEQUADRILATERAL'
        case('Tetrahedron')
            zonetype = 'FETETRAHEDRON'
        case('Hexahedron')
            zonetype = 'FEBRICK'
        case default
            write(*,'(a)') 'Wrong name in setZonetypeTecplot!'
            stop
        end select

    end function setZonetypeTecplot

    !-----------------------------------------------------------------------

    function setVariablesTecplot(Problem, Ndim) result(variables)

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Ndim
        character(200) :: variables

        PetscBool :: check

        select case (Ndim)
        case(0)
            variables = ''
        case(1)
            variables = '"X"'
        case(2)
            variables = '"X", "Y"'
        case(3)
            variables = '"X", "Y", "Z"'
        case default
            write(*,'(a)') 'Wrong Ndim in setVariablesTecplot!'
            stop
        end select

        select case (Problem%Neq_f)
        case(1)
            variables = trim(adjustl(variables))//', "Ux"'
        case(2)
            variables = trim(adjustl(variables))//', "Ux", "Uy"'
        case(3)
            variables = trim(adjustl(variables))//', "Ux", "Uy", "Uz"'
        case default
            write(*,'(a)') 'Wrong Problem%Neq_f in setVariablesTecplot!'
            stop
        end select

        check = (Problem%name == 'Inflow_1D') .or. (Problem%name == 'Inflow_2D')
        if (.not. check) then
            variables = trim(adjustl(variables))//', "P"'
        end if

        select case (Problem%Neq_s)
        case(3)
            variables = trim(adjustl(variables))//', "Txx", "Txy", "Tyy"'
        case(6)
            variables = trim(adjustl(variables))//', "Txx", "Txy", "Tyy", "Txz", "Tyz", "Tzz"'
        case default
            write(*,'(a)') 'Wrong Problem%Neq_s in setVariablesTecplot!'
            stop
        end select

        variables = trim(adjustl(variables))//', "detC"'

        if (Problem%Neq_proj > 0) then
            select case (Problem%Neq_proj)
            case(1)
                variables = trim(adjustl(variables))//', "Gxx"'
            case(4)
                variables = trim(adjustl(variables))//', "Gxx", "Gxy", "Gyx", "Gyy"'
            case(9)
                variables = trim(adjustl(variables))//', "Gxx", "Gxy", "Gyx", "Gyy"'// &
                    ', "Gxz", "Gzx", "Gyz", "Gzy", "Gzz"'
            case default
                write(*,'(a)') 'Wrong Problem%Neq_proj in setVariablesTecplot!'
            end select
            variables = trim(adjustl(variables))//', "Mc_local", "Flow_type"'
        end if

    end function setVariablesTecplot

    !-----------------------------------------------------------------------

    subroutine writeTecplotUnknowns(Problem, TL, TL_proj)

        use PhysicalParameters_mod, only: Stress_Reform

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscReal, dimension(:,:), intent(in) :: TL, TL_proj

        PetscInt :: inod, ieq, Nnodes_Total
        PetscReal, dimension(size(TL,1),Problem%Neq_s) :: T_Tensor
        PetscReal, dimension(size(TL,1)) :: Mc_local, FT_parameter, detC
        PetscBool :: check
        ! PetscReal, dimension(size(TL,1))  :: theta, Ur, Uth

        Nnodes_Total = size(TL,1)

        select case (Stress_Reform)

        case ('SQRT')

            !Ui
            do ieq = 1, Problem%Neq_f
                write(14,'(f32.15)') (TL(inod,ieq), inod = 1, size(TL,1))
            end do

            !P
            check = (Problem%name == 'Inflow_1D') .or. (Problem%name == 'Inflow_2D')
            if (.not. check) then
                ieq = Problem%Neq_f+1
                write(14,'(f32.15)') (TL(inod,ieq), inod = 1, size(TL,1))
            end if

            !Tij
            check = (Problem%name == 'Cylinder_Stability_2D')
            check = check .or. (Problem%name == 'Cylinder_Stability_2D_3D')
            check = check .or. (Problem%name == 'Cylinder_Stability_3D')
            if (check) then
                T_Tensor = calculateStressTensor_SQRT_Stability(Problem, TL)
            else
                T_Tensor = calculateStressTensor_SQRT(Problem, TL)
            end if
            do ieq = 1, size(T_Tensor,2)
                write(14,'(f32.15)') (T_Tensor(inod,ieq), inod = 1, size(T_Tensor,1))
            end do

            check = (Problem%name == 'Cylinder_Stability_2D')
            check = check .or. (Problem%name == 'Cylinder_Stability_2D_3D')
            check = check .or. (Problem%name == 'Cylinder_Stability_3D')
            if (check) then
                detC = calculateDetC_SQRT_Stability(Problem, TL)                
            else
                detC = calculateDetC_SQRT(Problem, TL)
            end if
            do inod = 1, size(detC)
                write(14,'(f32.15)') detC(inod)
            end do

        case default

            do ieq = 1, Problem%Neq
                write(14,'(f32.15)') (TL(inod,ieq), inod = 1, size(TL,1))
            end do

            check = (Problem%name == 'Cylinder_Stability_2D')
            check = check .or. (Problem%name == 'Cylinder_Stability_2D_3D')
            check = check .or. (Problem%name == 'Cylinder_Stability_3D')
            if (check) then
                detC = calculateDetC_T_Stability(Problem, TL)
            else
                detC = calculateDetC_T(Problem, TL)
            end if
            do inod = 1, size(detC)
                write(14,'(f32.15)') detC(inod)
            end do

        end select

        if (Problem%Neq_proj > 0) then

            do ieq = 1, Problem%Neq_proj
                write(14,"(f32.15)") (TL_proj(inod,ieq), inod = 1, size(TL_proj,1))
            end do

            Mc_local = calculateLocalM(Problem, TL, TL_proj)
            write(14,"(f32.15)") ( Mc_local(inod), inod = 1, size(Mc_local))

            FT_parameter = calculateFlowTypeParameter(TL_proj)
            write(14,"(f32.15)") ( FT_parameter(inod), inod = 1, size(FT_parameter))

        end if

        ! !Ur, Uth
        ! theta(:) = atan2(Y_Mesh,X_Mesh)
        ! Ur(:) =  TL(:,1)*cos(theta(:)) + TL(:,2)*sin(theta(:))
        ! Uth(:) = -TL(:,1)*sin(theta(:)) + TL(:,2)*cos(theta(:))

        ! write(14,'(f32.15)') ( Ur(inod), inod = 1, size(Ur))
        ! write(14,'(f32.15)') (Uth(inod), inod = 1, size(Uth))

    end subroutine writeTecplotUnknowns

    !-----------------------------------------------------------------------

    function calculateStressTensor_SQRT(Problem, TL) result(T_Tensor)

        use Tools_mod, only: getStressComponent

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscReal, dimension(:,:), intent(in) :: TL
        PetscReal, dimension(size(TL,1),Problem%Neq_s) :: T_Tensor

        PetscInt :: inod, istart, Neq, Neq_s, i, j, k, i2, Nnodes_Total
        PetscReal, dimension(3,3) :: T_Tensor_inod

        Nnodes_Total = size(TL,1)
        Neq = Problem%Neq
        Neq_s = Problem%Neq_s

        istart = Neq - Neq_s + 1
        do inod = 1, Nnodes_Total

            T_Tensor_inod(:,:) = nodalStressTensorSQRT(Problem, TL(inod,:))
            
            i2 = 0
            do i = istart, size(TL,2)
                call getStressComponent(i,istart,j,k)
                i2 = i2+1
                T_Tensor(inod,i2) = T_Tensor_inod(j,k)
            end do
            
        end do

    end function calculateStressTensor_SQRT

    !-----------------------------------------------------------------------

    function calculateStressTensor_SQRT_Stability(Problem, TL) result(T_p)

        use Tools_mod, only: getStressComponent
        use SolutionVariables_mod, only: Sol_Main

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscReal, dimension(:,:), intent(in) :: TL
        PetscReal, dimension(size(TL,1),Problem%Neq_s) :: T_p

        PetscInt :: inod, istart, Neq, Neq_s, i, j, k, i2, Nnodes_Total
        PetscReal, dimension(3,3) :: T_p_inod

        Nnodes_Total = size(TL,1)
        Neq = Problem%Neq
        Neq_s = Problem%Neq_s

        istart = Neq - Neq_s + 1
        
        do inod = 1, Nnodes_Total

            T_p_inod(:,:) = nodalStressTensorSQRT_Stability(Problem, TL(inod,:), Sol_Main%TL(inod,:))

            i2 = 0
            do i = istart, size(TL,2)
                call getStressComponent(i,istart,j,k)
                i2 = i2+1
                T_p(inod,i2) = T_p_inod(j,k)
            end do
            
        end do

    end function calculateStressTensor_SQRT_Stability

    !-----------------------------------------------------------------------

    function calculateDetC_SQRT(Problem, TL) result(detC)

        use Tools_mod, only: getStressComponent, determinant2D, determinant3D
        use ContinuationVariables_mod, only: Cvar1

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscReal, dimension(:,:), intent(in) :: TL
        PetscReal, dimension(size(TL,1)) :: detC

        PetscInt :: inod, istart, Neq, Neq_s, i, j, k, Nnodes_Total
        PetscReal, dimension(3,3) :: S_Tensor, C_Tensor
        PetscReal, dimension(2,2) :: C_Tensor_2D
        PetscBool :: check

        Nnodes_Total = size(TL,1)
        Neq = Problem%Neq
        Neq_s = Problem%Neq_s

        istart = Neq - Neq_s + 1
        do inod = 1, Nnodes_Total

            S_Tensor(:,:) = 0.0d0
            do i = istart, size(TL,2)
                call getStressComponent(i,istart,j,k)
                S_Tensor(j,k) = TL(inod,i)
                S_Tensor(k,j) = S_Tensor(j,k)
            end do
            
            C_Tensor(:,:) = matmul(S_Tensor,S_Tensor)

            select case (Problem%Ndim)
            case(1)
                detC(inod) = C_Tensor(1,1)
            case (2)
                C_Tensor_2D = C_Tensor(1:2,1:2)
                detC(inod) = determinant2D(C_Tensor_2D)
            case (3)
                detC(inod) = determinant3D(C_Tensor)
            case default
                write(*,'(a)') 'Wrong Problem%Ndim in calculateDetC_SQRT!'
                stop
            end select

            if (detC(inod) < 0.0d0) then
                write(*,'(a,es16.8,a,i8,3a,f12.4)') 'detC = ', detC(inod), &
                    ' at inod = ', inod, ' and ', Cvar1%name, ' = ', Cvar1%p
            end if

        end do

    end function calculateDetC_SQRT

    !-----------------------------------------------------------------------

    function calculateDetC_SQRT_Stability(Problem, TL) result(detC_p)

        use Tools_mod, only: getStressComponent, determinant2D, determinant3D
        use PhysicalParameters_mod, only: PM => Problem_Main
        use ContinuationVariables_mod, only: Cvar1
        use SolutionVariables_mod, only: Sol_Main

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscReal, dimension(:,:), intent(in) :: TL
        PetscReal, dimension(size(TL,1)) :: detC_p

        PetscInt :: inod, istart, Neq, Neq_s, i, j, k, Neq_M, Neq_s_M, Nnodes_Total
        PetscReal, dimension(3,3) :: S_p, C_p, S_b
        PetscReal, dimension(2,2) :: C_p_2D
        PetscBool :: check

        Nnodes_Total = size(TL,1)
        Neq = Problem%Neq
        Neq_s = Problem%Neq_s
        Neq_M = PM%Neq
        Neq_s_M = PM%Neq_s

        do inod = 1, Nnodes_Total

            istart = Neq_M - Neq_s_M + 1
            S_b(:,:) = 0.0d0
            do i = istart, size(Sol_Main%TL,2)
                call getStressComponent(i,istart,j,k)
                S_b(j,k) = Sol_Main%TL(inod,i)
                S_b(k,j) = S_b(j,k)
            end do

            istart = Neq - Neq_s + 1
            S_p(:,:) = 0.0d0
            do i = istart, size(TL,2)
                call getStressComponent(i,istart,j,k)
                S_p(j,k) = TL(inod,i)
                S_p(k,j) = S_p(j,k)
            end do

            C_p(:,:) = matmul(S_b,S_p) + matmul(S_p,S_b)

            select case (Problem%Ndim)
            case(1)
                detC_p(inod) = C_p(1,1)
            case (2)
                C_p_2D = C_p(1:2,1:2)
                detC_p(inod) = determinant2D(C_p_2D)
            case (3)
                detC_p(inod) = determinant3D(C_p)
            case default
                write(*,'(a)') 'Wrong Problem%Ndim in calculateDetC_SQRT_Stability!'
                stop
            end select

            ! if (detC_p(inod) < 0.0d0) then
            !     write(*,'(a,es16.8,a,i8,3a,f12.4)') 'detC_p = ', detC_p(inod), &
            !         ' at inod = ', inod, ' and ', Cvar1%name, ' = ', Cvar1%p
            ! end if

        end do

    end function calculateDetC_SQRT_Stability

    !-----------------------------------------------------------------------

    function calculateDetC_T(Problem, TL) result(detC)

        use Tools_mod, only: getStressComponent, traceTensor, determinant2D, determinant3D
        use PhysicalParameters_mod, only: BvN, WiN, Model_Name, VeN
        use ContinuationVariables_mod, only: Cvar1

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscReal, dimension(:,:), intent(in) :: TL
        PetscReal, dimension(size(TL,1)) :: detC

        PetscInt :: inod, istart, Neq, Neq_s, i, j, k, N, Nnodes_Total
        PetscReal :: trT, g1, g2
        PetscReal, dimension(3,3) :: T_Tensor, C_Tensor
        PetscReal, dimension(2,2) :: C_Tensor_2D
        PetscReal, dimension(3,3), parameter :: &
            I_tensor = reshape( [1,0,0,0,1,0,0,0,1], shape = shape(I_tensor) )
        PetscBool :: check

        Nnodes_Total = size(TL,1)
        Neq = Problem%Neq
        Neq_s = Problem%Neq_s
        N = Problem%Ndim

        istart = Neq - Neq_s + 1
        do inod = 1, Nnodes_Total

            T_Tensor = 0.0d0
            do i = istart, size(TL,2)
                call getStressComponent(i,istart,j,k)
                T_Tensor(j,k) = TL(inod,i)
                T_Tensor(k,j) = T_Tensor(j,k)
            end do

            trT = traceTensor(T_Tensor)

            select case (Model_Name)
            case ('m-L-PTT')
                g1 = 1.0d0 + (VeN*WiN*trT)/(1.0d0-BvN)
                g2 = g1
            case ('m-e-PTT')
                g1 = exp(VeN*WiN*trT/(1.0d0-BvN))
                g2 = g1
            case ('FENE-CR')
                g1 = 1.0d0 + ((WiN*trT/(1.0d0-BvN)))/(VeN**2-N)
                g2 = g1
            case ('FENE-P')
                g1 = 1.0d0
                g2 = 1.0d0 + ((WiN*trT/(1.0d0-BvN)))/(VeN**2)
            case default
                g1 = 1.0d0
                g2 = 1.0d0
            end select

            C_Tensor(:,:) = (WiN/(1.0d0-BvN))*T_Tensor(:,:) + g1*I_tensor(:,:)
            C_Tensor(:,:) = (1.0d0/g2)*C_Tensor(:,:)
            
            select case (Problem%Ndim)
            case(1)
                detC(inod) = C_Tensor(1,1)
            case (2)
                C_Tensor_2D = C_Tensor(1:2,1:2)
                detC(inod) = determinant2D(C_Tensor_2D)
            case (3)
                detC(inod) = determinant3D(C_Tensor)
            case default
                write(*,'(a)') 'Wrong Problem%Ndim in calculateDetC_T!'
                stop
            end select

            check = (Problem%name == 'Cylinder_Stability_2D')
            check = check .or. (Problem%name == 'Cylinder_Stability_2D_3D')
            check = check .or. (Problem%name == 'Cylinder_Stability_3D')
            check = (.not. check)
            check = check .and. (detC(inod) < 0.0d0)
            if (check) then
                write(*,'(a,es16.8,a,i8,3a,f12.4)') 'detC = ', detC(inod), &
                    ' at inod = ', inod, ' and ', Cvar1%name, ' = ', Cvar1%p
            end if
            
        end do

    end function calculateDetC_T

    !-----------------------------------------------------------------------

    function calculateDetC_T_Stability(Problem, TL) result(detC)

        use Tools_mod, only: getStressComponent, traceTensor, determinant2D, determinant3D
        use PhysicalParameters_mod, only: BvN, WiN, Model_Name, VeN, PM => Problem_Main
        use ContinuationVariables_mod, only: Cvar1
        use SolutionVariables_mod, only: Sol_Main

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscReal, dimension(:,:), intent(in) :: TL
        PetscReal, dimension(size(TL,1)) :: detC

        PetscInt :: inod, istart, Neq, Neq_s, i, j, k, N, Neq_M, Neq_s_M, Nnodes_Total
        PetscReal :: trT_b, trT_p
        PetscReal, dimension(2) :: g_b, g_p
        PetscReal, dimension(3,3) :: T_p, C_p, T_b, C_b
        PetscReal, dimension(2,2) :: C_p_2D
        PetscReal, dimension(3,3), parameter :: &
            I_tensor = reshape( [1,0,0,0,1,0,0,0,1], shape = shape(I_tensor) )
        PetscBool :: check

        Nnodes_Total = size(TL,1)
        Neq = Problem%Neq
        Neq_s = Problem%Neq_s
        N = Problem%Ndim
        Neq_M = PM%Neq
        Neq_s_M = PM%Neq_s

        do inod = 1, Nnodes_Total

            istart = Neq_M - Neq_s_M + 1
            T_b(:,:) = 0.0d0
            do i = istart, size(Sol_Main%TL,2)
                call getStressComponent(i,istart,j,k)
                T_b(j,k) = Sol_Main%TL(inod,i)
                T_b(k,j) = T_b(j,k)
            end do

            istart = Neq - Neq_s + 1
            T_p(:,:) = 0.0d0
            do i = istart, size(TL,2)
                call getStressComponent(i,istart,j,k)
                T_p(j,k) = TL(inod,i)
                T_p(k,j) = T_p(j,k)
            end do

            select case (Model_Name)
            case ('m-L-PTT')

                trT_b = traceTensor(T_b)
                g_b(1) = 1.0d0 + (VeN*WiN*trT_b)/(1.0d0-BvN)
                g_b(2) = g_b(1)

                trT_p = traceTensor(T_p)
                g_p(1) = VeN*WiN*trT_p/(1.0d0-BvN)
                g_p(2) = g_p(1)

            case ('m-e-PTT')

                trT_b = traceTensor(T_b)
                g_b(1) = exp(VeN*WiN*trT_b/(1.0d0-BvN))
                g_b(2) = g_b(1)

                trT_p = traceTensor(T_p)
                g_p(1) = g_b(1)*VeN*WiN*trT_p/(1.0d0-BvN)
                g_p(2) = g_p(1)

            case ('FENE-P')

                trT_b = traceTensor(T_b)
                g_b(1) = 1.0d0
                g_b(2) = 1.0d0 + (WiN*trT_b/(1.0d0-BvN))/(VeN**2)

                trT_p = traceTensor(T_p)
                g_p(1) = 0.0d0
                g_p(2) =(WiN*trT_p/(1.0d0-BvN))/(VeN**2)

            case default
                g_b(1) = 1.0d0
                g_b(2) = 1.0d0
                g_p(1) = 0.0d0
                g_p(2) = 0.0d0
            end select

            C_p(:,:) = WiN*(T_p(:,:)/(1.0d0-BvN)) + g_p(1)*I_tensor(:,:) - g_p(2)*C_b(:,:)
            C_p(:,:) = (1.0d0/g_b(2))*C_p(:,:)
            
            select case (Problem%Ndim)
            case(1)
                detC(inod) = C_p(1,1)
            case (2)
                C_p_2D = C_p(1:2,1:2)
                detC(inod) = determinant2D(C_p_2D)
            case (3)
                detC(inod) = determinant3D(C_p)
            case default
                write(*,'(a)') 'Wrong Problem%Ndim in calculateDetC_T_Stability!'
                stop
            end select

            check = (Problem%name == 'Cylinder_Stability_2D')
            check = check .or. (Problem%name == 'Cylinder_Stability_2D_3D')
            check = check .or. (Problem%name == 'Cylinder_Stability_3D')
            check = (.not. check)
            check = check .and. (detC(inod) < 0.0d0)
            if (check) then
                write(*,'(a,es16.8,a,i8,3a,f12.4)') 'detC = ', detC(inod), &
                    ' at inod = ', inod, ' and ', Cvar1%name, ' = ', Cvar1%p
            end if
            
        end do

    end function calculateDetC_T_Stability

    !-----------------------------------------------------------------------

    function calculateLocalM(Problem, TL, TL_proj) result(Mc)

        use Tools_mod, only: vectorMagnitude, crossProduct3D, &
            normFrobenius, traceTensor, tensorMagnitude, &
            getStressComponent, getTensorComponent
        use PhysicalParameters_mod, only: WiN, BvN, Model_Name

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscReal, dimension(:,:), intent(in) :: TL, TL_proj
        PetscReal, dimension(size(TL,1)) :: Mc

        PetscInt :: inod, i, j, k, istart, N, Nnodes_Total
        PetscReal, dimension(3,3) :: S, GU, GUT, G_dot
        PetscReal, dimension(3) :: U_f, Convection, cross_vector
        PetscReal :: R_curv, trace_T, f_T, U_mag, Tss, T_mag, t_rel0, t_rel_eff, G_dot_mag

        Nnodes_Total = size(TL,1)
        N = Problem%Ndim

        istart = Problem%Neq - Problem%Neq_s + 1

        U_f(:) = 0.0d0 ; GU(:,:) = 0.0d0 ; S(:,:) = 0.0d0
        do inod = 1, Nnodes_Total

            U_f(1:N) = TL(inod,1:N)

            do i = istart, size(TL,2)
                call getStressComponent(i,istart,j,k)
                S(j,k) = TL(inod,i)
            end do

            do i = 1, size(TL_proj,2)
                call getTensorComponent(i,j,k)
                GU(j,k) = TL_proj(inod,i)
            end do

            U_mag = vectorMagnitude(U_f)

            Convection(:) = matmul(U_f,GU)
            cross_vector(:) = crossProduct3D(U_f,Convection)

            R_curv = (U_mag**3)/(vectorMagnitude(cross_vector) + 1.0d-12)

            Tss = dot_product(matmul(U_f,S),U_f)/(U_mag**2 + 1.0d-12)
            ! Tss = S(1,1)

            T_mag = normFrobenius(S)

            trace_T = traceTensor(S)

            ! select case (Model_Name)
            ! case ('FENE-CR')
            !     f_T = 1.0d0 + ((WiN*trace_T/(1.0d0-BvN)))/(VeN**2-N)
            ! case ('FENE-P')
            !     f_T = 1.0d0 + ((WiN*trace_T/(1.0d0-BvN)))/(VeN**2)
            ! case ('L-PTT', 'm-L-PTT')
            !     f_T = 1.0d0 + (VeN*WiN*trace_T)/(1.0d0-BvN)
            ! case ('e-PTT', 'm-e-PTT')
            !     f_T = exp(VeN*WiN*trace_T/(1.0d0-BvN))
            ! case default
            !     write(*,'(A)') 'Wrong Model_Name in calculateLocalM!'
            !     stop
            ! end select

            t_rel0 = WiN
            t_rel_eff = t_rel0!/f_T

            Mc(inod) = sqrt( (t_rel_eff*U_mag/(R_curv+1.0d-12))*abs(Tss/(T_mag+1.0d-12)) )

            ! GUT(:,:) = transpose(GU)
            ! G_dot(:,:) = GU(:,:) + GUT(:,:)
            ! G_dot_mag = tensorMagnitude(G_dot)
            ! Mc(inod) = sqrt( (t_rel_eff*U_mag/(R_curv+1.0d-12))*abs(Tss/(G_dot_mag+1.0d-12)) )

        end do

    end function calculateLocalM

    !-----------------------------------------------------------------------

    function calculateFlowTypeParameter(TL_proj) result(FT)

        use Tools_mod, only: tensorMagnitude, getTensorComponent

        implicit none

        PetscReal, dimension(:,:), intent(in) :: TL_proj
        PetscReal, dimension(size(TL_proj,1)) :: FT

        PetscReal, dimension(3,3) :: GU, GUT, D_tensor, Omega_tensor
        PetscReal :: D_mag, Omega_mag
        PetscInt :: inod, i, j ,k, Nnodes_Total

        Nnodes_Total = size(TL_proj,1)
        GU(:,:) = 0.0d0 ; GUT(:,:) = 0.0d0
        do inod = 1, Nnodes_Total

            do i = 1, size(TL_proj,2)
                call getTensorComponent(i,j,k)
                GU(j,k) = TL_proj(inod,i)
            end do

            GUT(:,:) = transpose(GU)

            D_tensor(:,:) = 0.5d0*(GU(:,:) + GUT(:,:))
            Omega_tensor(:,:) = 0.5d0*(GU(:,:) - GUT(:,:))
            ! Omega_tensor(:,:) = 0.5d0*(GUT(:,:) - GU(:,:))

            D_mag = tensorMagnitude(D_tensor)
            Omega_mag = tensorMagnitude(Omega_tensor)

            FT(inod) = (D_mag - Omega_mag)/(D_mag + Omega_mag + 1.0d-12)

        end do

    end function calculateFlowTypeParameter

    !-----------------------------------------------------------------------

    function setConnectivityTecplot(FE, Connectivity) result(Con)

        use Tools_mod, only: allocateArray

        implicit none

        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:,:), intent(in) :: Connectivity
        PetscInt, dimension(:,:), allocatable :: Con

        PetscInt :: Nel, Nnodes
        PetscBool :: check_quad

        if (FE%order == 'Linear') then
            Nel = size(Connectivity,1)
            Nnodes = size(Connectivity,2)
            call allocateArray(Con, Nel, Nnodes)
            Con(:,:) = Connectivity(:,:)
            return
        end if

        check_quad = (FE%order == 'Quadratic' .or. FE%order == 'Serendipity')
        if (.not. check_quad) then
            write(*,'(a)') 'Wrong FE%order in setConnectivityTecplot!'
            return
        end if

        select case(FE%name)
        case('Line')
            Con = setConnectivityQuadraticLine(FE%Nel, Connectivity)
        case('Triangle')
            Con = setConnectivityQuadraticTriangle(FE%Nel, Connectivity)
        case('Quadrangle')
            Con = setConnectivityQuadraticQuadrangle(FE%Nel, Connectivity)
        case('Tetrahedron')
            Con = setConnectivityQuadraticTetrahedron(FE%Nel, Connectivity)
        case('Hexahedron')
            Con = setConnectivityQuadraticHexahedron(FE%Nel, Connectivity)
        case default
            write(*,'(a)') 'Wrong FE%name in setConnectivityTecplot!'
            stop
        end select

    end function setConnectivityTecplot

    !-----------------------------------------------------------------------

    function setConnectivityQuadraticLine(Nel, Connectivity) result(Con)

        use Tools_mod, only: allocateArray

        implicit none

        PetscInt, intent(in) :: Nel
        PetscInt, dimension(:,:), intent(in) :: Connectivity
        PetscInt, allocatable, dimension(:,:) :: Con

        PetscInt, parameter :: Nbf_v = 2, Nel_mult = 2
        PetscInt :: Nel_v, iel

        Nel_v = Nel_mult*Nel

        call allocateArray(Con, Nel_v, Nbf_v)

        do iel = 1, Nel
            !Line 1
            Con(Nel_mult*iel-1,1) = Connectivity(iel,1)
            Con(Nel_mult*iel-1,2) = Connectivity(iel,3)

            !Line 2
            Con(Nel_mult*iel,1) = Connectivity(iel,3)
            Con(Nel_mult*iel,2) = Connectivity(iel,2)
        end do

    end function setConnectivityQuadraticLine

    !-----------------------------------------------------------------------

    function setConnectivityQuadraticTriangle(Nel, Connectivity) result(Con)

        use Tools_mod, only: allocateArray

        implicit none

        PetscInt, intent(in) :: Nel
        PetscInt, dimension(:,:), intent(in) :: Connectivity
        PetscInt, allocatable, dimension(:,:) :: Con

        PetscInt, parameter :: Nbf_v = 3, Nel_mult = 4
        PetscInt :: Nel_v, iel

        Nel_v = Nel_mult*Nel

        call allocateArray(Con, Nel_v, Nbf_v)

        do iel = 1, Nel
            !Triangle 1
            Con(Nel_mult*iel-3,1) = Connectivity(iel,5)
            Con(Nel_mult*iel-3,2) = Connectivity(iel,3)
            Con(Nel_mult*iel-3,3) = Connectivity(iel,6)

            !Triangle 2
            Con(Nel_mult*iel-2,1) = Connectivity(iel,4)
            Con(Nel_mult*iel-2,2) = Connectivity(iel,2)
            Con(Nel_mult*iel-2,3) = Connectivity(iel,5)

            !Triangle 3
            Con(Nel_mult*iel-1,1) = Connectivity(iel,4)
            Con(Nel_mult*iel-1,2) = Connectivity(iel,5)
            Con(Nel_mult*iel-1,3) = Connectivity(iel,6)

            !Triangle 4
            Con(Nel_mult*iel,1) = Connectivity(iel,1)
            Con(Nel_mult*iel,2) = Connectivity(iel,4)
            Con(Nel_mult*iel,3) = Connectivity(iel,6)
        end do

    end function setConnectivityQuadraticTriangle

    !-----------------------------------------------------------------------

    function setConnectivityQuadraticQuadrangle(Nel, Connectivity) result(Con)

        use Tools_mod, only: allocateArray

        implicit none

        PetscInt, intent(in) :: Nel
        PetscInt, dimension(:,:), intent(in) :: Connectivity
        PetscInt, allocatable, dimension(:,:) :: Con

        PetscInt, parameter :: Nbf_v = 4, Nel_mult = 4
        PetscInt :: Nel_v, iel

        Nel_v = Nel_mult*Nel

        call allocateArray(Con, Nel_v, Nbf_v)

        do iel = 1, Nel

            !Quadrangle 1
            Con(Nel_mult*iel-3,1) = Connectivity(iel,8)
            Con(Nel_mult*iel-3,2) = Connectivity(iel,9)
            Con(Nel_mult*iel-3,3) = Connectivity(iel,7)
            Con(Nel_mult*iel-3,4) = Connectivity(iel,4)

            !Quadrangle 2
            Con(Nel_mult*iel-2,1) = Connectivity(iel,9)
            Con(Nel_mult*iel-2,2) = Connectivity(iel,6)
            Con(Nel_mult*iel-2,3) = Connectivity(iel,3)
            Con(Nel_mult*iel-2,4) = Connectivity(iel,7)

            !Quadrangle 3
            Con(Nel_mult*iel-1,1) = Connectivity(iel,5)
            Con(Nel_mult*iel-1,2) = Connectivity(iel,2)
            Con(Nel_mult*iel-1,3) = Connectivity(iel,6)
            Con(Nel_mult*iel-1,4) = Connectivity(iel,9)

            !Quadrangle 4
            Con(Nel_mult*iel,1) = Connectivity(iel,1)
            Con(Nel_mult*iel,2) = Connectivity(iel,5)
            Con(Nel_mult*iel,3) = Connectivity(iel,9)
            Con(Nel_mult*iel,4) = Connectivity(iel,8)

        end do

    end function setConnectivityQuadraticQuadrangle

    !-----------------------------------------------------------------------

    function setConnectivityQuadraticTetrahedron(Nel, Connectivity) result(Con)

        use Tools_mod, only: allocateArray

        implicit none

        PetscInt, intent(in) :: Nel
        PetscInt, dimension(:,:), intent(in) :: Connectivity
        PetscInt, allocatable, dimension(:,:) :: Con

        PetscInt, parameter :: Nbf_v = 4, Nel_mult = 8
        PetscInt :: Nel_v, iel

        Nel_v = Nel_mult*Nel

        call allocateArray(Con, Nel_v, Nbf_v)

        do iel = 1, Nel

            !Tetrahedron 1
            Con(Nel_mult*iel-7,1) = Connectivity(iel,1)
            Con(Nel_mult*iel-7,2) = Connectivity(iel,5)
            Con(Nel_mult*iel-7,3) = Connectivity(iel,7)
            Con(Nel_mult*iel-7,4) = Connectivity(iel,8)

            !Tetrahedron 2
            Con(Nel_mult*iel-6,1) = Connectivity(iel,5)
            Con(Nel_mult*iel-6,2) = Connectivity(iel,2)
            Con(Nel_mult*iel-6,3) = Connectivity(iel,6)
            Con(Nel_mult*iel-6,4) = Connectivity(iel,9)

            !Tetrahedron 3
            Con(Nel_mult*iel-5,1) = Connectivity(iel,3)
            Con(Nel_mult*iel-5,2) = Connectivity(iel,6)
            Con(Nel_mult*iel-5,3) = Connectivity(iel,7)
            Con(Nel_mult*iel-5,4) = Connectivity(iel,10)

            !Tetrahedron 4
            Con(Nel_mult*iel-4,1) = Connectivity(iel,4)
            Con(Nel_mult*iel-4,2) = Connectivity(iel,8)
            Con(Nel_mult*iel-4,3) = Connectivity(iel,9)
            Con(Nel_mult*iel-4,4) = Connectivity(iel,10)

            !Tetrahedron 5
            Con(Nel_mult*iel-3,1) = Connectivity(iel,9)
            Con(Nel_mult*iel-3,2) = Connectivity(iel,5)
            Con(Nel_mult*iel-3,3) = Connectivity(iel,7)
            Con(Nel_mult*iel-3,4) = Connectivity(iel,8)

            !Tetrahedron 6
            Con(Nel_mult*iel-2,1) = Connectivity(iel,5)
            Con(Nel_mult*iel-2,2) = Connectivity(iel,7)
            Con(Nel_mult*iel-2,3) = Connectivity(iel,6)
            Con(Nel_mult*iel-2,4) = Connectivity(iel,9)

            !Tetrahedron 7
            Con(Nel_mult*iel-1,1) = Connectivity(iel,7)
            Con(Nel_mult*iel-1,2) = Connectivity(iel,6)
            Con(Nel_mult*iel-1,3) = Connectivity(iel,9)
            Con(Nel_mult*iel-1,4) = Connectivity(iel,10)

            !Tetrahedron 8
            Con(Nel_mult*iel,1) = Connectivity(iel,8)
            Con(Nel_mult*iel,2) = Connectivity(iel,9)
            Con(Nel_mult*iel,3) = Connectivity(iel,10)
            Con(Nel_mult*iel,4) = Connectivity(iel,7)

        end do

    end function setConnectivityQuadraticTetrahedron

    !-----------------------------------------------------------------------

    function setConnectivityQuadraticHexahedron(Nel, Connectivity) result(Con)

        use Tools_mod, only: allocateArray

        implicit none

        PetscInt, intent(in) :: Nel
        PetscInt, dimension(:,:), intent(in) :: Connectivity
        PetscInt, allocatable, dimension(:,:) :: Con

        PetscInt, parameter :: Nbf_v = 8, Nel_mult = 8
        PetscInt :: Nel_v, iel

        Nel_v = Nel_mult*Nel

        call allocateArray(Con, Nel_v, Nbf_v)

        do iel = 1, Nel

            !Hexahedron 1
            Con(Nel_mult*iel-7,1) = Connectivity(iel,1)
            Con(Nel_mult*iel-7,2) = Connectivity(iel,9)
            ! Con(Nel_mult*iel-7,3) = Connectivity(iel,25)
            Con(Nel_mult*iel-7,3) = Connectivity(iel,21)  !Salome
            Con(Nel_mult*iel-7,4) = Connectivity(iel,12)
            Con(Nel_mult*iel-7,5) = Connectivity(iel,17)
            ! Con(Nel_mult*iel-7,6) = Connectivity(iel,21)
            Con(Nel_mult*iel-7,6) = Connectivity(iel,22)  !Salome
            Con(Nel_mult*iel-7,7) = Connectivity(iel,27)
            ! Con(Nel_mult*iel-7,8) = Connectivity(iel,24)
            Con(Nel_mult*iel-7,8) = Connectivity(iel,25)  !Salome

            !Hexahedron 2
            Con(Nel_mult*iel-6,1) = Connectivity(iel,9)
            Con(Nel_mult*iel-6,2) = Connectivity(iel,2)
            Con(Nel_mult*iel-6,3) = Connectivity(iel,10)
            ! Con(Nel_mult*iel-6,4) = Connectivity(iel,25)
            Con(Nel_mult*iel-6,4) = Connectivity(iel,21)  !Salome
            ! Con(Nel_mult*iel-6,5) = Connectivity(iel,21)
            Con(Nel_mult*iel-6,5) = Connectivity(iel,22)  !Salome
            Con(Nel_mult*iel-6,6) = Connectivity(iel,18)
            ! Con(Nel_mult*iel-6,7) = Connectivity(iel,22)
            Con(Nel_mult*iel-6,7) = Connectivity(iel,23)  !Salome
            Con(Nel_mult*iel-6,8) = Connectivity(iel,27)

            !Hexahedron 3
            ! Con(Nel_mult*iel-5,1) = Connectivity(iel,25)
            Con(Nel_mult*iel-5,1) = Connectivity(iel,21)  !Salome
            Con(Nel_mult*iel-5,2) = Connectivity(iel,10)
            Con(Nel_mult*iel-5,3) = Connectivity(iel,3)
            Con(Nel_mult*iel-5,4) = Connectivity(iel,11)
            Con(Nel_mult*iel-5,5) = Connectivity(iel,27)
            ! Con(Nel_mult*iel-5,6) = Connectivity(iel,22)
            Con(Nel_mult*iel-5,6) = Connectivity(iel,23)  !Salome
            Con(Nel_mult*iel-5,7) = Connectivity(iel,19)
            ! Con(Nel_mult*iel-5,8) = Connectivity(iel,23)
            Con(Nel_mult*iel-5,8) = Connectivity(iel,24)  !Salome

            !Hexahedron 4
            Con(Nel_mult*iel-4,1) = Connectivity(iel,12)
            ! Con(Nel_mult*iel-4,2) = Connectivity(iel,25)
            Con(Nel_mult*iel-4,2) = Connectivity(iel,21)  !Salome
            Con(Nel_mult*iel-4,3) = Connectivity(iel,11)
            Con(Nel_mult*iel-4,4) = Connectivity(iel,4)
            ! Con(Nel_mult*iel-4,5) = Connectivity(iel,24)
            Con(Nel_mult*iel-4,5) = Connectivity(iel,25)  !Salome
            Con(Nel_mult*iel-4,6) = Connectivity(iel,27)
            ! Con(Nel_mult*iel-4,7) = Connectivity(iel,23)
            Con(Nel_mult*iel-4,7) = Connectivity(iel,24)  !Salome
            Con(Nel_mult*iel-4,8) = Connectivity(iel,20)

            !Hexahedron 5
            Con(Nel_mult*iel-3,1) = Connectivity(iel,17)
            ! Con(Nel_mult*iel-3,2) = Connectivity(iel,21)
            Con(Nel_mult*iel-3,2) = Connectivity(iel,22)  !Salome
            Con(Nel_mult*iel-3,3) = Connectivity(iel,27)
            ! Con(Nel_mult*iel-3,4) = Connectivity(iel,24)
            Con(Nel_mult*iel-3,4) = Connectivity(iel,25)  !Salome
            Con(Nel_mult*iel-3,5) = Connectivity(iel,5)
            Con(Nel_mult*iel-3,6) = Connectivity(iel,13)
            Con(Nel_mult*iel-3,7) = Connectivity(iel,26)
            Con(Nel_mult*iel-3,8) = Connectivity(iel,16)

            !Hexahedron 6
            ! Con(Nel_mult*iel-2,1) = Connectivity(iel,21)
            Con(Nel_mult*iel-2,1) = Connectivity(iel,22)  !Salome
            Con(Nel_mult*iel-2,2) = Connectivity(iel,18)
            ! Con(Nel_mult*iel-2,3) = Connectivity(iel,22)
            Con(Nel_mult*iel-2,3) = Connectivity(iel,23)  !Salome
            Con(Nel_mult*iel-2,4) = Connectivity(iel,27)
            Con(Nel_mult*iel-2,5) = Connectivity(iel,13)
            Con(Nel_mult*iel-2,6) = Connectivity(iel,6)
            Con(Nel_mult*iel-2,7) = Connectivity(iel,14)
            Con(Nel_mult*iel-2,8) = Connectivity(iel,26)

            !Hexahedron 7
            Con(Nel_mult*iel-1,1) = Connectivity(iel,27)
            ! Con(Nel_mult*iel-1,2) = Connectivity(iel,22)
            Con(Nel_mult*iel-1,2) = Connectivity(iel,23)  !Salome
            Con(Nel_mult*iel-1,3) = Connectivity(iel,19)
            ! Con(Nel_mult*iel-1,4) = Connectivity(iel,23)
            Con(Nel_mult*iel-1,4) = Connectivity(iel,24)  !Salome
            Con(Nel_mult*iel-1,5) = Connectivity(iel,26)
            Con(Nel_mult*iel-1,6) = Connectivity(iel,14)
            Con(Nel_mult*iel-1,7) = Connectivity(iel,7)
            Con(Nel_mult*iel-1,8) = Connectivity(iel,15)

            !Hexahedron 8
            ! Con(Nel_mult*iel,1) = Connectivity(iel,24)
            Con(Nel_mult*iel,1) = Connectivity(iel,25)  !Salome
            Con(Nel_mult*iel,2) = Connectivity(iel,27)
            ! Con(Nel_mult*iel,3) = Connectivity(iel,23)
            Con(Nel_mult*iel,3) = Connectivity(iel,24)  !Salome
            Con(Nel_mult*iel,4) = Connectivity(iel,20)
            Con(Nel_mult*iel,5) = Connectivity(iel,16)
            Con(Nel_mult*iel,6) = Connectivity(iel,26)
            Con(Nel_mult*iel,7) = Connectivity(iel,15)
            Con(Nel_mult*iel,8) = Connectivity(iel,8)

        end do

    end function setConnectivityQuadraticHexahedron

    !-----------------------------------------------------------------------

    subroutine writeSolution(FE, Connectivity, Sol, fn)

        use ContinuationParameters_mod, only: Continuation_Method

        implicit none

        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:,:), intent(in) :: Connectivity
        type(SolutionArraysType), intent(in) :: Sol
        character(*), intent(in) :: fn

        PetscInt :: inod, ieq, ierror, iel, gnod
        character(len=150) :: fn_full

        fn_full = trim(adjustl(fn))//'.DTA'
        fn_full = trim(adjustl(fn_full))

        open(12, file=fn_full, status='unknown', action='write', iostat=ierror, &
            position='rewind', form='unformatted')

        ! do inod = 1, size(Sol%TL,1)
        !     do ieq = 1, size(Sol%TL,2)
        !         write(12) Sol%TL(inod,ieq)
        !     end do
        ! end do

        do iel = 1, FE%Nel
            do inod = 1, FE%Nbf
                gnod = Connectivity(iel,inod)
                do ieq = 1, size(Sol%TL,2)
                    write(12) Sol%TL(gnod,ieq)
                end do
            end do
        end do

        do ieq = 1, size(Sol%EU)
            write(12) Sol%EU(ieq)
        end do

        close(12, status='keep')

        if (Continuation_Method == 'Transient') then
            !Previous
            fn_full = trim(adjustl(fn))//'_o.DTA'
            fn_full = trim(adjustl(fn_full))

            open(12, file=fn_full, status='unknown', action='write', iostat=ierror, &
                position='rewind', form='unformatted')

            do iel = 1, FE%Nel
                do inod = 1, FE%Nbf
                    gnod = Connectivity(iel,inod)
                    do ieq = 1, size(Sol%TLo,2)
                        write(12) Sol%TLo(gnod,ieq)
                    end do
                end do
            end do

            do ieq = 1, size(Sol%EUo)
                write(12) Sol%EUo(ieq)
            end do

            close(12, status='keep')

            !Back
            fn_full = trim(adjustl(fn))//'_b.DTA'
            fn_full = trim(adjustl(fn_full))

            open(12, file=fn_full, status='unknown', action='write', iostat=ierror, &
                position='rewind', form='unformatted')

            do iel = 1, FE%Nel
                do inod = 1, FE%Nbf
                    gnod = Connectivity(iel,inod)
                    do ieq = 1, size(Sol%TLb,2)
                        write(12) Sol%TLb(gnod,ieq)
                    end do
                end do
            end do

            do ieq = 1, size(Sol%EUb)
                write(12) Sol%EUb(ieq)
            end do

            close(12, status='keep')
        end if

    end subroutine writeSolution

    !-----------------------------------------------------------------------

    subroutine readSolution(FE, Connectivity, iter, fn, Sol)

        use Tools_mod, only: extrapolationLagrange, getRankNode
        use ContinuationParameters_mod, only: Continuation_Method
        use PhysicalParameters_mod, only: dt0
        use MPIParameters_mod, only: Rank, NRanks

        implicit none

        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:,:), intent(in) :: Connectivity
        PetscReal, intent(in) :: iter
        character(*), intent(in) :: fn
        type(SolutionArraysType), intent(inout) :: Sol

        PetscInt :: inod, ieq, ierr, iel, gnod
        PetscInt :: Nel_Rank0, iel_rank, iel_start, rnod
        PetscReal :: Lb, Lo, L
        PetscBool :: check
        character(len=150) :: fn_full

        fn_full = trim(adjustl(fn))//'.DTA'
        fn_full = trim(adjustl(fn_full))

        open(13, file=fn_full, status='unknown', action='read', iostat=ierr, &
            position='rewind', form='unformatted')

        if (ierr /= 0) return

        if (Rank == 0) then

            do iel = 1, FE%Nel
                do inod = 1, FE%Nbf
                    gnod = Connectivity(iel,inod)
                    do ieq = 1, size(Sol%TL,2)
                        read(13) Sol%TL(gnod,ieq)
                    end do
                end do
            end do

            do ieq = 1, size(Sol%EU)
                read(13, iostat=ierr) Sol%EU(ieq)
                if (ierr /= 0) exit
            end do

            rewind(13)

            Sol%TLb(:,:) = Sol%TL(:,:)
            Sol%TLo(:,:) = Sol%TL(:,:)
            Sol%TLp(:,:) = Sol%TL(:,:)

        end if

        call MPI_Bcast(Sol%EU, size(Sol%EU), MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)

        !For each rank
        Nel_Rank0  = FE%Nel/NRanks

        !Skip
        iel_start = Rank*Nel_Rank0+1
        do iel = 1, iel_start-1
            do inod = 1, FE%Nbf
                do ieq = 1, size(Sol%TL_Rank,2)
                    read(13)
                end do
            end do
        end do

        do iel_rank = 1, FE%Nel_Rank
            do inod = 1, FE%Nbf
                rnod = getRankNode(iel_rank,inod,FE%Nbf)
                do ieq = 1, size(Sol%TL_Rank,2)
                    read(13) Sol%TL_Rank(rnod,ieq)
                end do
            end do
        end do

        Sol%TLb_Rank(:,:) = Sol%TL_Rank(:,:)
        Sol%TLo_Rank(:,:) = Sol%TL_Rank(:,:)
        Sol%TLp_Rank(:,:) = Sol%TL_Rank(:,:)

        Sol%EUb(:) = Sol%EU(:)
        Sol%EUo(:) = Sol%EU(:)
        Sol%EUp(:) = Sol%EU(:)

        close(13, status='keep')

        check = (Continuation_Method == 'Transient')
        if (check) then
            call readSolutionPrevious(FE, Connectivity, fn, Sol)
            call readSolutionBack(FE, Connectivity, fn, Sol)
            call extrapolationLagrange(iter+dt0, iter-(dt0+dt0), iter-dt0, iter, Lb, Lo, L)

            if (Rank == 0) then
                Sol%TLp(:,:) = Lb*Sol%TLb(:,:) + Lo*Sol%TLo(:,:) + L*Sol%TL(:,:)
                Sol%TLb(:,:) = Sol%TLo(:,:)
                Sol%TLo(:,:) = Sol%TL(:,:)
                Sol%TL(:,:) = Sol%TLp(:,:)
            end if

            Sol%TLp_Rank(:,:) = Lb*Sol%TLb_Rank(:,:) + Lo*Sol%TLo_Rank(:,:) + L*Sol%TL_Rank(:,:)
            Sol%TLb_Rank(:,:) = Sol%TLo_Rank(:,:)
            Sol%TLo_Rank(:,:) = Sol%TL_Rank(:,:)
            Sol%TL_Rank(:,:)  = Sol%TLp_Rank(:,:)

            Sol%EUp(:) = Lb*Sol%EUb(:) + Lo*Sol%EUo(:) + L*Sol%EU(:)
            Sol%EUb(:) = Sol%EUo(:) ; Sol%EUo(:) = Sol%EU(:) ; Sol%EU(:) = Sol%EUp(:)
        end if
        
        call MPI_Barrier(MPI_COMM_WORLD,ierr)

    end subroutine readSolution

    !-----------------------------------------------------------------------

    subroutine readSolutionPrevious(FE, Connectivity, fn, Sol)

        use Tools_mod, only: getRankNode
        use MPIParameters_mod, only: Rank, NRanks

        implicit none

        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:,:), intent(in) :: Connectivity
        character(*), intent(in) :: fn
        type(SolutionArraysType), intent(inout) :: Sol

        PetscInt :: inod, ieq, ierr, iel, gnod, Nel_Rank0, iel_rank, iel_start, rnod
        character(150) :: fn_full
        PetscBool :: check

        fn_full = trim(adjustl(fn))//'.DTA'
        fn_full = trim(adjustl(fn_full))

        open(13, file=fn_full, status='unknown', action='read', iostat=ierr, &
            position='rewind', form='unformatted')

        if (ierr /= 0) return

        if (Rank == 0) then

            do iel = 1, FE%Nel
                do inod = 1, FE%Nbf
                    gnod = Connectivity(iel,inod)
                    do ieq = 1, size(Sol%TLo,2)
                        read(13) Sol%TLo(gnod,ieq)
                    end do
                end do
            end do

            do ieq = 1, size(Sol%EUo)
                read(13, iostat=ierr) Sol%EUo(ieq)
                if (ierr /= 0) exit
            end do

            rewind(13)

        end if

        call MPI_Bcast(Sol%EUo, size(Sol%EUo), MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)

        !For each rank
        Nel_Rank0  = FE%Nel/NRanks

        !Skip
        iel_start = Rank*Nel_Rank0+1
        do iel = 1, iel_start-1
            do inod = 1, FE%Nbf
                do ieq = 1, size(Sol%TLo_Rank,2)
                    read(13)
                end do
            end do
        end do

        do iel_rank = 1, FE%Nel_Rank
            do inod = 1, FE%Nbf
                rnod = getRankNode(iel_rank,inod,FE%Nbf)
                do ieq = 1, size(Sol%TLo_Rank,2)
                    read(13) Sol%TLo_Rank(rnod,ieq)
                end do
            end do
        end do

        close(13, status='keep')

    end subroutine readSolutionPrevious

    !-----------------------------------------------------------------------

    subroutine readSolutionBack(FE, Connectivity, fn, Sol)

        use Tools_mod, only: getRankNode
        use MPIParameters_mod, only: Rank, NRanks

        implicit none

        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:,:), intent(in) :: Connectivity
        character(*), intent(in) :: fn
        type(SolutionArraysType), intent(inout) :: Sol

        PetscInt :: inod, ieq, ierr, iel, gnod, Nel_Rank0, iel_rank, iel_start, rnod
        character(150) :: fn_full
        PetscBool :: check

        fn_full = trim(adjustl(fn))//'.DTA'
        fn_full = trim(adjustl(fn_full))

        open(13, file=fn_full, status='unknown', action='read', iostat=ierr, &
            position='rewind', form='unformatted')

        if (ierr /= 0) return

        if (Rank == 0) then

            do iel = 1, FE%Nel
                do inod = 1, FE%Nbf
                    gnod = Connectivity(iel,inod)
                    do ieq = 1, size(Sol%TLb,2)
                        read(13) Sol%TLb(gnod,ieq)
                    end do
                end do
            end do

            do ieq = 1, size(Sol%EUb)
                read(13, iostat=ierr) Sol%EUb(ieq)
                if (ierr /= 0) exit
            end do

            rewind(13)

        end if

        call MPI_Bcast(Sol%EUb, size(Sol%EUb), MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)

        !For each rank
        Nel_Rank0  = FE%Nel/NRanks

        !Skip
        iel_start = Rank*Nel_Rank0+1
        do iel = 1, iel_start-1
            do inod = 1, FE%Nbf
                do ieq = 1, size(Sol%TLb_Rank,2)
                    read(13)
                end do
            end do
        end do

        do iel_rank = 1, FE%Nel_Rank
            do inod = 1, FE%Nbf
                rnod = getRankNode(iel_rank,inod,FE%Nbf)
                do ieq = 1, size(Sol%TLb_Rank,2)
                    read(13) Sol%TLb_Rank(rnod,ieq)
                end do
            end do
        end do

        close(13, status='keep')

    end subroutine readSolutionBack

    !--------------------------------------------------------------------------

    subroutine calculatePostProcessQuantities(Problem, FE, Mesh, Boundary, GaussInt, Sol, Elem)

        use MPIParameters_mod, only: Rank
        use ElementCalculations_mod, only: basis, basisMain_p, basisInflow_p

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        type(NodeArrays), intent(inout) :: Elem

        PetscInt :: ibnd, ieq, ieq_ed, ieq_fc
        PetscReal :: term
        procedure(basis), pointer :: basis_p

        if (Sol%name == 'Main') then
            basis_p => basisMain_p
        end if

        if (Sol%name == 'Inflow') then
            basis_p => basisInflow_p
        end if

        ! if (Sol%name == 'Main') then
            call calculateBulkQuantities(Problem, FE, Mesh, GaussInt, Sol, basis_p, Elem)
        ! end if

        do ibnd = 1, size(Boundary)
            do ieq = 1, size(Boundary(ibnd)%epp_ed)
                ieq_ed = Boundary(ibnd)%epp_ed(ieq)
                call calculateEdgeQuantities(Problem, FE, Mesh, Boundary(ibnd), &
                    GaussInt, Sol, basis_p, ieq_ed, Elem, term)
            end do
            do ieq = 1, size(Boundary(ibnd)%epp_fc)
                ieq_fc = Boundary(ibnd)%epp_fc(ieq)
                call calculateFaceQuantities(Problem, FE, Mesh, Boundary(ibnd), &
                    GaussInt, Sol, basis_p, ieq_fc, Elem, term)
            end do
        end do

        if (Rank == 0) then
            call calculateNodalQuantities(Problem, FE%Nnodes, Mesh%Xi_Mesh, Sol%TL)
            ! call compareStrongWeakFormConvergence(Problem, FE, Mesh, Boundary, &
            !     GaussInt, basis_p, Sol, Elem)
            write(*,'(a)') '--------------------------------------------------------------'
        end if

    end subroutine calculatePostProcessQuantities

    !--------------------------------------------------------------------------

    subroutine calculateBulkQuantities(Problem, FE, Mesh, GaussInt, Sol, basis_p, Elem)

        use MPIParameters_mod,  only: Rank
        use ContinuationVariables_mod, only: Cvar1, Cvar2
        use ElementCalculations_mod, only: copyToLocal, basis

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem

        PetscInt :: ierr, iel_rank, Nbf
        PetscReal :: u1, u2, u1_total, u2_total, Ap_2D, abs_Uz, Ap_3D
        character(200) :: name, fn, folder, str
        character(len=:), allocatable :: val2

        Nbf = FE%Nbf

        u1_total = 0.0d0 ; u2_total = 0.0d0
        Ap_2D = 0.0d0 ; Ap_3D = 0.0d0
        loop_elements:do iel_rank = 1, FE%Nel_Rank

            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)

            call setQuantitiesBulk(Problem, Nbf, GaussInt, basis_p, Elem, u1, u2, abs_Uz)

            u1_total = u1_total + u1
            u2_total = u2_total + u2
            Ap_3D = Ap_3D + abs_Uz

        end do loop_elements

        call MPI_Allreduce(MPI_IN_PLACE,u1_total,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,u2_total,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call MPI_Allreduce(MPI_IN_PLACE,Ap_3D,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

        Ap_2D = (u1_total-u2_total)/(u1_total+u2_total+1.0d-8)

        if (Rank == 0) then

            write(str,'(f12.4)') Cvar2%p
            val2 = trim(adjustl(str))

            name = Cvar2%name
            folder = 'Results/Base/'//trim(adjustl(name))//'/' &
                //trim(adjustl(name))//'_'//trim(adjustl(val2))//'/Flow_data/'

            name = Problem%name
            fn = trim(adjustl(folder))//'Ap_'//trim(adjustl(name))//'.dat'
            open(24, file=fn, iostat=ierr, position='append')
            write(24,'(f16.8,2es24.12)') Cvar1%p, Ap_2D, Ap_3D
            close(24, status='keep')

            write(*,'(a,es16.8,10x)', advance='no') 'Ap_2D = ', Ap_2D
            write(*,'(a,es16.8)') 'Ap_3D = ', Ap_3D

        end if

    end subroutine calculateBulkQuantities

    !--------------------------------------------------------------------------

    subroutine setQuantitiesBulk(Problem, Nbf, GaussInt, basis_p, Elem, u1, u2, abs_Uz)

        use ElementCalculations_mod, only: setElementBasisFunctions, &
            setFlowQuantities, basis

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nbf
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscReal, intent(out) :: u1, u2, abs_Uz

        type(GaussPointQuantities) :: GsPt
        PetscInt :: ig
        PetscBool :: check
        PetscReal, dimension(Nbf) :: y

        u1 = 0.0d0 ; u2 = 0.0d0 ; abs_Uz = 0.0d0
        loop_gauss:do ig = 1, GaussInt%NGauss_b

            call setElementBasisFunctions(ig, GaussInt, Elem)
            
            call basis_p(Problem%idir, Nbf, Elem, GsPt)

            call setFlowQuantities(Problem, Nbf, Elem, GsPt)

            GsPt%WET = GaussInt%Weights_b(ig)*abs(GsPt%detJ)

            y(:) = Elem%Xi_loc(:,2)
            check = ( all(y(:) > 0.0d0) )
            if (check) then
                u1 = u1 + GsPt%U_f(1)*(GsPt%WET)
            end if

            check = ( all(y(:) < 0.0d0) )
            if (check) then
                u2 = u2 + GsPt%U_f(1)*(GsPt%WET)
            end if

            abs_Uz = abs_Uz + abs(GsPt%U_f(3))*(GsPt%WET)
        
        end do loop_gauss

    end subroutine setQuantitiesBulk

    !--------------------------------------------------------------------------

    subroutine calculateEdgeQuantities(Problem, FE, Mesh, bound, GaussInt, Sol, &
        basis_p, ieq_ed, Elem, term)

        use PhysicalParameters_mod, only: ReN
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Cvar1, Cvar2
        use ElementCalculations_mod, only: copyToLocal, basis, &
            setElementBasisFunctionsEdge, setFlowQuantities, setEdgeQuantities

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        PetscInt, intent(in) :: ieq_ed
        type(NodeArrays), intent(inout) :: Elem
        PetscReal, intent(out) :: term

        PetscInt :: iel_rank, ied, iel, Nbf, ig
        PetscReal :: term_el, drag_coef
        PetscBool :: check
        character(100) :: name, fn, str, folder
        character(:), allocatable :: val2
        type(GaussPointQuantities) :: GsPt
        PetscErrorCode :: ierr

        Nbf = FE%Nbf
        term = 0.0d0
        loop_elements:do iel = 1, size(bound%edge_elements_rank)

            iel_rank = bound%edge_elements_rank(iel)

            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)

            ied = bound%edges(iel)

            term_el = 0.0d0
            loop_gauss:do ig = 1, GaussInt%NGauss_1D

                call setElementBasisFunctionsEdge(ig, ied, GaussInt, Elem)

                call basis_p(Problem%idir, Nbf, Elem, GsPt)

                call setFlowQuantities(Problem, Nbf, Elem, GsPt)

                call setEdgeQuantities(ied, Problem%idir, FE%name, GsPt)

                GsPt%WET = (GaussInt%Weights_1D(ig))*(GsPt%dL)

                select case (ieq_ed)
                case (0)
                    return
                case (1)
                    term_el = term_el + dot_product(-GsPt%normal,GsPt%P_tot(:,1))*(GsPt%WET)
                case default
                    write(*,'(a)') "Wrong ieq_ed choice in calculateEdgeQuantities!"
                    stop
                end select

            end do loop_gauss

            term = term + term_el

        end do loop_elements

        call MPI_Allreduce(MPI_IN_PLACE,term,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

        check = (Rank == 0) .and. (ieq_ed == 1)
        if (check) then

            drag_coef = 2.0d0*term/(2.0d0*ReN + 1.0d-8)
            drag_coef = drag_coef/2.0d0

            write(str,'(f12.4)') Cvar2%p
            val2 = trim(adjustl(str))

            name = Cvar2%name
            folder = 'Results/Base/'//trim(adjustl(name))//'/' &
                //trim(adjustl(name))//'_'//trim(adjustl(val2))//'/Flow_data/'

            name = Problem%name
            fn = trim(adjustl(folder))//'Fd_'//trim(adjustl(name))//'.dat'
            open(24, file=fn, iostat=ierr, position='append')
            write(24,'(f16.8,es24.12)') Cvar1%p, term
            close(24, status='keep')

            write(*,'(a,f15.8)') 'Fd = ', term

        end if

    end subroutine calculateEdgeQuantities

    !--------------------------------------------------------------------------

    subroutine calculateFaceQuantities(Problem, FE, Mesh, bound, GaussInt, Sol, &
        basis_p, ieq_fc, Elem, term)

        use PhysicalParameters_mod, only: ReN
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Cvar1, Cvar2
        use ElementCalculations_mod, only: copyToLocal, basis, &
            setElementBasisFunctionsFace, setFlowQuantities, setFaceQuantities

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        PetscInt, intent(in) :: ieq_fc
        type(NodeArrays), intent(inout) :: Elem
        PetscReal, intent(out) :: term

        PetscInt :: iel_rank, ifc, iel, Nbf, ig
        PetscReal :: term_el, drag_coef
        PetscBool :: check
        character(100) :: name, fn, folder, str
        character(:), allocatable :: val2
        type(GaussPointQuantities) :: GsPt
        PetscErrorCode :: ierr

        Nbf = FE%Nbf
        term = 0.0d0
        loop_elements:do iel = 1, size(bound%face_elements_rank)

            iel_rank = bound%face_elements_rank(iel)

            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)

            ifc = bound%faces(iel)

            term_el = 0.0d0
            loop_gauss:do ig = 1, GaussInt%NGauss_2D

                call setElementBasisFunctionsFace(ig, ifc, GaussInt, Elem)

                call basis_p(Problem%idir, Nbf, Elem, GsPt)

                call setFlowQuantities(Problem, Nbf, Elem, GsPt)

                call setFaceQuantities(ifc, FE%name, GsPt)

                GsPt%WET = (GaussInt%Weights_2D(ig))*(GsPt%dL)

                select case (ieq_fc)
                case (0)
                    return
                case (1)
                    term_el = term_el + dot_product(-GsPt%normal,GsPt%P_tot(:,1))*(GsPt%WET)
                case default
                    write(*,'(a)') "Wrong ieq_fc choice in calculateFaceQuantities!"
                    stop
                end select

            end do loop_gauss

            term = term + term_el

        end do loop_elements

        call MPI_Allreduce(MPI_IN_PLACE,term,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

        check = (Rank == 0) .and. (ieq_fc == 1)
        if (check) then

            drag_coef = 2.0d0*term/(2.0d0*ReN + 1.0d-8)
            drag_coef = drag_coef/2.0d0

            write(str,'(f12.4)') Cvar2%p
            val2 = trim(adjustl(str))

            name = Cvar2%name
            folder = 'Results/Base/'//trim(adjustl(name))//'/' &
                //trim(adjustl(name))//'_'//trim(adjustl(val2))//'/Flow_data/'

            name = Problem%name
            fn = trim(adjustl(folder))//'Fd_'//trim(adjustl(name))//'.dat'
            open(24, file=fn, iostat=ierr, position='append')
            write(24,'(f16.8,es24.12)') Cvar1%p, term
            close(24, status='keep')

            write(*,'(a,f15.8)') 'Fd = ', term

        end if

    end subroutine calculateFaceQuantities

    !--------------------------------------------------------------------------

    subroutine calculateNodalQuantities(Problem, Nnodes, Xi_Mesh, TL)

        use Tools_mod, only: arrayClosestValue, sortFile, getStressComponent
        use PhysicalParameters_mod, only: r => Radius, PI, Stress_Reform, Length, QbN
        use MeshParameters_mod, only: eps => Mesh_Tolerance
        use ContinuationVariables_mod, only: Cvar1, Cvar2

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nnodes
        PetscReal, dimension(:,:), intent(in) :: Xi_Mesh, TL

        PetscInt, parameter :: Nbd = 3, Ntargets = 1
        PetscReal, dimension(Ntargets), parameter :: x_target = [1.1d0]
        PetscReal, dimension(Ntargets) :: x_close
        PetscReal, dimension(:), allocatable :: x_nodes
        PetscInt :: inod, ibnd, Ncols_sort, Neq_f, i
        PetscReal :: theta, L_arc, P
        PetscReal, dimension(3) :: xi, U_f
        PetscReal, dimension(3,3) :: T_Tensor
        PetscBool :: check
        character(200) :: str, fn1, fn2, folder
        character(len=:), allocatable :: val2, val1

        x_nodes = pack(Xi_Mesh(:,1), Xi_Mesh(:,2) > 0.0d0)
        do i = 1, size(x_target)
            x_close(i) = arrayClosestValue(x_nodes,x_target(i))
        end do
        deallocate(x_nodes)

        check = all( abs(x_close-Xi_Mesh(1,1)) < eps )
        if (check) then
            ! x_close(:) = +huge(1.0d0)
            return
        end if

        write(str,'(f12.4)') Cvar2%p
        val2 = trim(adjustl(str))

        write(str,'(f10.4)') Cvar1%p
        val1 = trim(adjustl(str))

        folder = 'Results/Base/'//trim(adjustl(Cvar2%name))//'/' &
            //trim(adjustl(Cvar2%name))//'_'//trim(adjustl(val2))//'/'

        fn1 = trim(adjustl(folder))//'Txx/'//trim(adjustl(Problem%name))//'_'
        fn1 = trim(adjustl(fn1))//trim(adjustl(Cvar1%name))//'_'//val1//'.dat'
        open(24, file=fn1, action='write', position='rewind')

        fn2 = trim(adjustl(folder))//'Line/'//trim(adjustl(Problem%name))//'_'
        fn2 = trim(adjustl(fn2))//trim(adjustl(Cvar1%name))//'_'//val1//'.dat'
        open(20, file=fn2, action='write', position='rewind')

        deallocate(val2) ; deallocate(val1)

        Neq_f = Problem%Neq_f

        do inod = 1, Nnodes

            U_f(:) = 0.0d0
            U_f(1:Neq_f) = TL(inod,1:Neq_f)

            P = TL(inod,Neq_f+1)

            select case (Stress_Reform)
            case ('SQRT')
                T_Tensor(:,:) = nodalStressTensorSQRT(Problem, TL(inod,:))
            case default
                T_Tensor(:,:) = nodalStressTensorT(Problem, TL(inod,:))
            end select

            xi(:) = Xi_Mesh(inod,:)
            do ibnd = 1, Nbd

                select case (ibnd)
                case (1)
                    if ( distanceToExtraBND(ibnd, xi, x_close(1)) ) then
                        theta = atan2(xi(2),xi(1))
                        L_arc = theta*r
                        if (abs(xi(1)+r) < eps) L_arc = 0.0d0
                        if (abs(xi(1)-r) < eps) L_arc = PI*r
                        write(24,'(*(es16.8))') L_arc, T_Tensor(1,1)
                    end if
                case (2)
                    if ( distanceToExtraBND(ibnd, xi, x_close(1)) ) then
                        L_arc = PI*r+(xi(1)-r)
                        write(24,'(*(es16.8))') L_arc, T_Tensor(1,1)
                    end if
                case (3)
                    do i = 1, size(x_close)
                        if ( distanceToExtraBND(ibnd, xi, x_close(i)) ) then  
                            write(20,'(*(es16.8))') xi(3), U_f, P, &
                                T_Tensor(1,1), T_Tensor(1,2), T_Tensor(2,2), &
                                T_Tensor(1,3), T_Tensor(2,3), T_Tensor(3,3)
                        end if
                    end do
                case default

                end select

            end do

        end do

        close(24, status='keep')
        close(20, status='keep')

        Ncols_sort = 2
        call sortFile(fn1, Ncols_sort)

        Ncols_sort = 11
        call sortFile(fn2, Ncols_sort)
        call wavelengthData(Problem%name, folder, fn2)

    end subroutine calculateNodalQuantities

    !--------------------------------------------------------------------------

    function nodalStressTensorSQRT(Problem, TL_inod) result(T_Tensor)

        use Tools_mod, only: getStressComponent, traceTensor
        use PhysicalParameters_mod, only: BvN, WiN, VeN, Model_Name

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscReal, dimension(:), intent(in) :: TL_inod
        PetscReal, dimension(3,3) :: T_tensor

        PetscReal, dimension(3,3), parameter :: &
            I_tensor = reshape( [1,0,0,0,1,0,0,0,1], shape = shape(I_tensor) )
        PetscInt :: i, j, k, istart
        PetscInt :: Neq, Neq_s, N
        PetscReal ::  g1, g2, trC
        PetscReal, dimension(3,3) :: S_Tensor, C_Tensor

        Neq = Problem%Neq
        Neq_s = Problem%Neq_s
        N = Problem%Ndim

        istart = Neq - Neq_s + 1

        S_Tensor(:,:) = 0.0d0
        do i = istart, size(TL_inod,1)
            call getStressComponent(i,istart,j,k)
            S_Tensor(j,k) = TL_inod(i)
            S_Tensor(k,j) = S_Tensor(j,k)
        end do
        C_Tensor(:,:) = matmul(S_Tensor,S_Tensor)

        trC = traceTensor(C_Tensor)
        select case (Model_Name)
        case ('m-L-PTT')
            g1 = 1.0d0+VeN*(trC-N)
            g2 = g1
        case ('m-e-PTT')
            g1 = exp(VeN*(trC-N))
            g2 = g1
        case ('FENE-CR')
            g1 = (VeN-N)/(VeN-trC)
            g2 = g1
        case ('FENE-P')
            g1 = 1.0d0
            g2 = VeN/(VeN-trC)
        case default
            g1 = 1.0d0
            g2 = 1.0d0
        end select

        T_Tensor(:,:) = g2*C_Tensor(:,:)-g1*I_Tensor(:,:)
        T_Tensor(:,:) = ((1.0d0-BvN)/WiN)*T_Tensor(:,:)

    end function nodalStressTensorSQRT

    !--------------------------------------------------------------------------

    function nodalStressTensorSQRT_Stability(Problem, TL_inod, TLb_inod) result(T_p)

        use Tools_mod, only: getStressComponent, traceTensor
        use PhysicalParameters_mod, only: BvN, WiN, VeN, Model_Name, PM => Problem_Main

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscReal, dimension(:), intent(in) :: TL_inod, TLb_inod
        PetscReal, dimension(3,3) :: T_p

        PetscReal, dimension(3,3), parameter :: &
            I_tensor = reshape( [1,0,0,0,1,0,0,0,1], shape = shape(I_tensor) )
        PetscInt :: i, j, k, istart
        PetscInt :: Neq, Neq_s, N, Neq_M, Neq_s_M
        PetscReal :: trace_C_p, trace_C_b
        PetscReal, dimension(2) :: g_b, g_p
        PetscReal, dimension(3,3) :: S_p, C_p, S_b, C_b

        Neq = Problem%Neq
        Neq_s = Problem%Neq_s
        N = Problem%Ndim
        Neq_M = PM%Neq
        Neq_s_M = PM%Neq_s

        istart = Neq_M - Neq_s_M + 1
        S_b(:,:) = 0.0d0
        do i = istart, size(TLb_inod)
            call getStressComponent(i,istart,j,k)
            S_b(j,k) = TLb_inod(i)
            S_b(k,j) = S_b(j,k)
        end do

        C_b(:,:) = matmul(S_b,S_b)

        istart = Neq - Neq_s + 1
        S_p(:,:) = 0.0d0
        do i = istart, size(TL_inod,1)
            call getStressComponent(i,istart,j,k)
            S_p(j,k) = TL_inod(i)
            S_p(k,j) = S_p(j,k)
        end do

        C_p(:,:) = matmul(S_b,S_p) + matmul(S_p,S_b)

        select case (Model_Name)
        case ('FENE-CR')
            
            trace_C_b = traceTensor(C_b)
            g_b(1) = (VeN - N)/(VeN-trace_C_b)
            g_b(2) = g_b(1)

            trace_C_p = traceTensor(C_p)
            g_p(1) = g_b(1)*trace_C_p/(VeN-trace_C_b)
            g_p(2) = g_p(1)

        case ('FENE-P')
        
            trace_C_b = traceTensor(C_b)
            g_b(1) = 1.0d0
            g_b(2) = (VeN - N)/(VeN-trace_C_b)

            trace_C_p = traceTensor(C_p)
            g_p(1) = 0.0d0
            g_p(2) = g_b(1)*trace_C_p/(VeN-trace_C_b)

        case default
            g_b(1) = 1.0d0
            g_b(2) = 1.0d0
            g_p(1) = 0.0d0
            g_p(2) = 0.0d0
        end select

        T_p(:,:) = g_b(2)*C_p(:,:)+g_p(2)*C_b(:,:)-g_p(1)*I_Tensor(:,:)
        T_p(:,:) = ((1.0d0-BvN)/WiN)*T_p(:,:)

    end function nodalStressTensorSQRT_Stability

    !--------------------------------------------------------------------------

    function nodalStressTensorT(Problem, TL_inod) result(T_Tensor)

        use Tools_mod, only: getStressComponent

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscReal, dimension(:), intent(in) :: TL_inod
        PetscReal, dimension(3,3) :: T_tensor

        PetscInt :: i, j, k, istart
        PetscInt :: Neq, Neq_s

        Neq = Problem%Neq
        Neq_s = Problem%Neq_s
        istart = Neq-Neq_s+1
        do i = istart, size(TL_inod,1)
            call getStressComponent(i,istart,j,k)
            T_Tensor(j,k) = TL_inod(i)
            T_Tensor(k,j) = T_Tensor(j,k)
        end do

    end function nodalStressTensorT

    !--------------------------------------------------------------------------

    subroutine wavelengthData(Problem_name, folder, fn)

        use ContinuationVariables_mod, only: Cvar1

        implicit none

        character(*), intent(in) :: Problem_name, folder, fn

        PetscInt :: ioerr, nchanges, n
        PetscReal, dimension(3) :: U_f, U_f_mean
        PetscReal :: Uzo, xi
        character(100) :: fn1

        open(20, file=fn, action='read', iostat=ioerr, position='rewind')

        nchanges = 0
        Uzo = 0.0d0
        U_f_mean(:) = 0.0d0
        n = 0
        do
            read(20,*,iostat=ioerr) xi, U_f(1), U_f(2), U_f(3)
            if (ioerr /= 0) exit

            U_f_mean(:) = U_f_mean(:) + U_f(:)
            n = n+1

            if (Uzo*U_f(3) < 0.0d0) then
                nchanges = nchanges+1
            end if

            Uzo = U_f(3)

        end do
        close(20, status='keep')

        U_f_mean(:) = U_f_mean(:)/max(1,n)

        fn1 = trim(adjustl(folder))//'Line/Wave_'//trim(adjustl(Problem_name))//'.dat'

        open(15, file=fn1, action='write', iostat=ioerr, position='append')
        write(15,'(4es20.8,f6.1)') Cvar1%p, U_f_mean, nchanges/2.0d0
        close(15, status='keep')

    end subroutine wavelengthData

    !--------------------------------------------------------------------------

    function distanceToExtraBND(ibnd, x, x_close) result(check)

        use PhysicalParameters_mod, only: r => Radius
        use MeshParameters_mod, only: eps => Mesh_Tolerance

        implicit none

        PetscInt, intent(in) :: ibnd
        PetscReal, dimension(:), intent(in) :: x
        PetscReal, intent(in) :: x_close
        PetscBool :: check

        select case (ibnd)
        case (1)
            check = (x(2) >= 0.0d0) &
                .and. (abs(x(1)**2+x(2)**2-r**2) < eps)
        case (2)
            check = (abs(x(2)-0.0d0) < eps) &
                .and. (x(1) >= r) .and. (x(1) <= 8.0d0*r)
        case (3)
            check = (abs(x(1)-x_close) < eps) .and. (x(2) > 0.0d0) !&
                    ! (abs(x(2)-0.0d0) < eps)
        case default
            write(*,'(a)') 'Wrong ibnd in distanceToExtraBND!'
            stop
        end select

    end function distanceToExtraBND

    !--------------------------------------------------------------------------

    subroutine compareStrongWeakFormConvergence(Problem, FE, Mesh, Boundary, &
        GaussInt, basis_p, Sol, Elem)

        use PhysicalParameters_mod, only: Problem_Inflow
        use FEMParameters_mod, only: FE_Inflow
        use ContinuationVariables_mod, only: Increment, Cvar2, Cvar1
        use SolutionVariables_mod, only: Sol_Inflow
        use ElementCalculations_mod, only: copyToLocalGlobal, basis
        use Equations_mod, only: residualBulk
        use Dirichlet_mod, only: Inflow_Nodes_Main

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem

        PetscInt :: iel, Nel, Nbf, inod, ieq, irow, gnod, ibnd, jeq
        PetscReal :: res_norm_weak, res_norm_strong
        PetscInt, dimension(Problem%Neq) :: node_weak, node_strong
        PetscReal, dimension(Problem%Neq) :: Res_temp, max_res_weak, max_res_strong
        PetscScalar, dimension(FE%Nbf, Problem%Neq) :: Res_weak_el, Res_strong_el
        PetscReal, dimension(FE%Nnodes, Problem%Neq) :: Res_weak, Res_strong
        character(100) :: str, fn, folder
        character(len=:), allocatable :: val

        write(str,'(f12.4)') Cvar2%p
        val = trim(adjustl(str))

        folder = 'Results/Base/'//trim(adjustl(Cvar2%name))//'/'
        folder = trim(adjustl(folder))//trim(adjustl(Cvar2%name)) &
            //'_'//trim(adjustl(val))//'/Convergence/'
        folder = trim(adjustl(folder))//trim(adjustl(Problem%name))//'/'

        deallocate(val)

        call execute_command_line('mkdir -p '//folder)

        write(str,'(f12.4)') Cvar1%p

        fn = trim(adjustl(folder))//'WeakForm_'// &
            trim(adjustl(Cvar1%name))//'_'//trim(adjustl(str))//'.dat'
        open(2, file=fn,position='rewind')

        fn = trim(adjustl(folder))//'StrongForm_'// &
            trim(adjustl(Cvar1%name))//'_'//trim(adjustl(str))//'.dat'
        open(3, file=fn,position='rewind')

        Nel = FE%Nel
        Nbf = FE%Nbf

        !Bulk
        Res_weak(:,:) = 0.0d0
        Res_strong(:,:) = 0.0d0
        do iel = 1, Nel

            call copyToLocalGlobal(Nbf, Mesh, Sol, iel, Elem)

            call residualBulk(Problem, Nbf, GaussInt, &
                    Sol%name, basis_p, Elem, Res_weak_el)
            call residualBulk(Problem, Nbf, GaussInt, &
                    'Strong', basis_p, Elem, Res_strong_el)

            do inod = 1, Nbf

                gnod = Mesh%Connectivity(iel,inod)

                Res_temp(:) = PetscRealPart(Res_weak_el(inod,:))
                Res_weak(gnod,:) = Res_weak(gnod,:) + Res_temp(:)

                Res_temp(:) = PetscRealPart(Res_strong_el(inod,:))
                Res_strong(gnod,:) = Res_strong(gnod,:) + Res_temp(:)

            end do

        end do

        !Residual = 0 at boundaries
        do ibnd = 1, size(Boundary)
            do inod = 1, size(Boundary(ibnd)%gnodes_total)
                gnod = Boundary(ibnd)%gnodes_total(inod)
                Res_weak(gnod,:) = 0.0d0
                Res_strong(gnod,:) = 0.0d0
            end do
        end do

        ! !Residual at boundaries
        ! bnd_loop:do ibnd = 1, size(Boundary)

        !     !Inflow
        !     if (Boundary(ibnd)%name == 'Inlet') then
        !         do inod = 1, FE_Inflow%Nnodes
        !             gnod = Inflow_Nodes_Main(inod)
        !             do ieq = 1, Problem_Inflow%Neq
        !                 jeq = ieq
        !                 !Skip P
        !                 if (ieq > Problem_Inflow%Neq_f) jeq = ieq+1
        !                 Res_weak(gnod,jeq) = Sol%TL(gnod,jeq) - Sol_Inflow%TL(inod,ieq)
        !             end do
        !         end do
        !         cycle bnd_loop
        !     end if

        !     !Dirichlet
        !     do inod = 1, size(Boundary(ibnd)%gnodes_total)
        !         gnod = Boundary(ibnd)%gnodes_total(inod)
        !         do ieq = 1, size(Boundary(ibnd)%values_Dirichlet)
        !             jeq = Boundary(ibnd)%eq_Dirichlet(ieq)
        !             Res_weak(gnod,jeq) = Sol%TL(gnod,jeq) &
        !                     - Boundary(ibnd)%values_Dirichlet(ieq)
        !         end do
        !     end do

        ! end do bnd_loop

        do ieq = 1, Problem%Neq
            max_res_weak(ieq) = maxval( abs(Res_weak(:,ieq)) )
            node_weak(ieq) = maxloc( abs(Res_weak(:,ieq)), dim=1 )
            max_res_strong(ieq) = maxval( abs(Res_strong(:,ieq)) )
            node_strong(ieq) = maxloc( abs(Res_strong(:,ieq)), dim=1 )

            write(2,'(i6,(*(es12.4)))') node_weak(ieq), max_res_weak(ieq)
            write(3,'(i6,(*(es12.4)))') node_strong(ieq), max_res_strong(ieq)
        end do

        res_norm_weak = 0.0d0
        res_norm_strong = 0.0d0
        do inod = 1, FE%Nnodes

            ! !Actual residual (Dirichlet etc included)
            ! do ieq = 1, Problem%Neq
            !     irow = (inod-1)*Problem%Neq+ieq
            !     Res_weak(inod,ieq) = PetscRealPart(B_f_pointer(irow))
            ! end do
            
            Res_temp(:) = Res_weak(inod,:)
            write(2,'(i6,(*(es12.4)))') inod, Res_temp
            res_norm_weak = res_norm_weak + dot_product(Res_temp,Res_temp)

            Res_temp(:) = Res_strong(inod,:)
            write(3,'(i6,(*(es12.4)))') inod, Res_temp
            res_norm_strong = res_norm_strong + dot_product(Res_temp,Res_temp)

        end do
        res_norm_weak = sqrt(res_norm_weak)
        res_norm_strong = sqrt(res_norm_strong)

        write(*,'(a,es20.12)') 'Res_norm_weak = ', res_norm_weak
        write(*,'(a,es20.12)') 'Res_norm_strong = ', res_norm_strong

        close(2)
        close(3)

    end subroutine compareStrongWeakFormConvergence
    
    !--------------------------------------------------------------------------

    subroutine projectQuantities(Problem, FE, Mesh, GaussInt, Sol, Elem, LS)

        use MPIParameters_mod,  only: Rank
        use PhysicalParameters_mod, only: ProblemParameters
        use FEMParameters_mod, only: FEParameters
        use MeshParameters_mod, only: MeshParameters
        use GaussParameters_mod, only: GaussIntegration
        use SolutionVariables_mod, only: SolutionArraysType
        use ElementVariables_mod, only: NodeArrays
        use LinearSystemVariables_mod, only: LinearSystemType
        use Projection_mod, only: loopNewtonRaphsonProjection

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(inout) :: Sol
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS

        if (Rank == 0) then
            Sol%TL_proj(:,:) = 0.1d0
        end if

        Sol%TL_proj_Rank(:,:) = 0.1d0

        call loopNewtonRaphsonProjection(Problem, FE, Mesh, GaussInt, Sol, Elem, LS)

    end subroutine projectQuantities

end module PostProcess_mod
