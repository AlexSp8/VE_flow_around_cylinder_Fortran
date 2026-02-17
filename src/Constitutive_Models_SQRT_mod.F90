
#include <petsc/finclude/petscksp.h>

module ConstitutiveModelsSQRT_mod

    use petscksp
    use Tools_mod, only: traceTensor
    use ElementVariables_mod, only: GaussPointQuantities

    implicit none

    interface omegaNominator2D
        module procedure omegaNominator2DReal
        module procedure omegaNominator2DRealComplex
        module procedure omegaNominator2DComplexReal
    end interface omegaNominator2D

    interface omegaDenominator2D
        module procedure omegaDenominator2DReal
        module procedure omegaDenominator2DComplex
    end interface omegaDenominator2D

    interface WTensor3D_RHS
        module procedure WTensor3D_RHS_Real
        module procedure WTensor3D_RHS_RealComplex
        module procedure WTensor3D_RHS_ComplexReal
    end interface

    interface WTensor3D_LHS
        module procedure WTensor3D_LHS_Real
        module procedure WTensor3D_LHS_Complex
    end interface

    interface
        subroutine conEQ(GsPt)
            use ElementVariables_mod, only: GaussPointQuantities
            implicit none
            type(GaussPointQuantities), intent(inout) :: GsPt
        end subroutine conEQ

        subroutine W_Tensor(S, GUT, W)
            implicit none
            PetscReal, dimension(:,:), intent(in) :: S, GUT
            PetscReal, dimension(:,:), intent(out) :: W
        end subroutine W_Tensor

        subroutine SI_Tensor(S,SI)
            implicit none
            PetscReal, dimension(:,:), intent(in) :: S
            PetscReal, dimension(:,:), intent(out) :: SI
        end subroutine SI_Tensor
    end interface

    procedure(conEQ), pointer :: conEQ_p
    procedure(W_Tensor), pointer :: WTensor_p
    procedure(SI_Tensor), pointer :: inverseTensor_p

    contains

    subroutine setConEQ_SQRT()

        use PhysicalParameters_mod, only: Model_Name

        implicit none

        select case (Model_Name)
        case ('FENE-CR')
            conEQ_p => conEQ_FENE_CR_SQRT
        case ('FENE-P')
            conEQ_p => conEQ_FENE_P_SQRT
        case ('L-PTT')
            conEQ_p => conEQ_L_PTT_SQRT
        case ('m-L-PTT')
            conEQ_p => conEQ_m_L_PTT_SQRT
        case ('e-PTT')
            conEQ_p => conEQ_e_PTT_SQRT
        case ('m-e-PTT')
            conEQ_p => conEQ_m_e_PTT_SQRT
        case ('GIESEKUS')
            conEQ_p => conEQ_Giesekus_SQRT
        case default
            write(*,'(a)') 'Wrong Model_Name in setConEQ_SQRT!'
            stop
        end select

    end subroutine setConEQ_SQRT

    !---------------------------------------------------------------------

    subroutine setPointersSQRT(Ndim)

        use Tools_mod, only: inverseTensor2D, inverseTensor3D

        implicit none

        PetscInt, intent(in) :: Ndim

        select case (Ndim)
        case (0)
            WTensor_p => NULL()
            inverseTensor_p => NULL()
        case (1)
            WTensor_p => WTensor2D
            inverseTensor_p => inverseTensor2D
        case (2)
            WTensor_p => WTensor2D
            inverseTensor_p => inverseTensor2D
        case (3)
            WTensor_p => WTensor3D
            inverseTensor_p => inverseTensor3D
        case default
            write(*,'(a)') 'Wrong Ndim in setPointersSQRT!'
            stop
        end select

    end subroutine setPointersSQRT

    !---------------------------------------------------------------------

    subroutine strongFormConEQ_SQRT(GsPt)

        use PhysicalParameters_mod, only: WiN

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscReal :: trC
        PetscReal, dimension(3,3) :: SI, UCD, S3, term_S, S_real
        
        call setCTensorSQRT(GsPt)

        call conEQ_p(GsPt)

        call setStressTensorSQRT(GsPt)

        UCD(:,:) = setUpperConvectedDerivativeSQRT(GsPt)

        S_real(:,:) = PetscRealPart(GsPt%S)
        SI(:,:) = setSInverse(S_real)

        S3(:,:) = matmul( GsPt%S,matmul(GsPt%S,GsPt%S) )

        term_S(:,:) = (GsPt%f(3))*S3(:,:)+(GsPt%f(2))*GsPt%S(:,:)-(GsPt%f(1))*SI(:,:)

        GsPt%Con_EQ(:,:) = WiN*UCD(:,:)+term_S(:,:)/2.0d0
        GsPt%Con_EQ_time(:,:) = 0.0d0

    end subroutine strongFormConEQ_SQRT

    !---------------------------------------------------------------------

    subroutine conEQ_L_PTT_SQRT(GsPt)

        use PhysicalParameters_mod, only: ePTT => VeN, PM => Problem_Main

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt :: N
        PetscReal :: trC

        N = PM%Ndim

        trC = traceTensor(GsPt%C_Tensor)

        GsPt%f_T = 1.0d0 + ePTT*(trC - N)

        GsPt%f(1) = GsPt%f_T
        GsPt%f(2) = GsPt%f(1)
        GsPt%f(3) = 0.0d0

        GsPt%g(1) = 1.0d0
        GsPt%g(2) = 1.0d0

        GsPt%dgdxi(1,:) = 0.0d0
        GsPt%dgdxi(2,:) = 0.0d0

    end subroutine conEQ_L_PTT_SQRT

    !---------------------------------------------------------------------

    subroutine conEQ_m_L_PTT_SQRT(GsPt)

        use PhysicalParameters_mod, only: ePTT => VeN, PM => Problem_Main

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt :: i, N
        PetscReal :: trC

        N = PM%Ndim

        trC = traceTensor(GsPt%C_Tensor)

        GsPt%f_T = 1.0d0 + ePTT*(trC - N)

        GsPt%f(1) = GsPt%f_T
        GsPt%f(2) = GsPt%f(1)
        GsPt%f(3) = 0.0d0

        GsPt%g(1) = GsPt%f_T
        GsPt%g(2) = GsPt%f_T

        GsPt%dgdxi(1,:) = 0.0d0
        do i = 1, 3
            GsPt%dgdxi(1,:) = GsPt%dgdxi(1,:) + GsPt%dCdXi(i,i,:)
        end do
        GsPt%dgdxi(1,:) = ePTT*(GsPt%dgdxi(1,:))
        GsPt%dgdxi(2,:) = GsPt%dgdxi(1,:)

    end subroutine conEQ_m_L_PTT_SQRT

    !---------------------------------------------------------------------

    subroutine conEQ_e_PTT_SQRT(GsPt)

        use PhysicalParameters_mod, only: ePTT => VeN, PM => Problem_Main

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt :: N
        PetscReal :: trC

        N = PM%Ndim

        trC = traceTensor(GsPt%C_Tensor)

        GsPt%f_T = exp(ePTT*(trC-N))

        GsPt%f(1) = GsPt%f_T
        GsPt%f(2) = GsPt%f(1)
        GsPt%f(3) = 0.0d0

        GsPt%g(1) = 1.0d0
        GsPt%g(2) = 1.0d0

        GsPt%dgdxi(1,:) = 0.0d0
        GsPt%dgdxi(2,:) = 0.0d0

    end subroutine conEQ_e_PTT_SQRT

    !---------------------------------------------------------------------

    subroutine conEQ_m_e_PTT_SQRT(GsPt)

        use PhysicalParameters_mod, only: ePTT => VeN, PM => Problem_Main

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt :: i, N
        PetscReal :: trC

        N = PM%Ndim

        trC = traceTensor(GsPt%C_Tensor)

        GsPt%f_T = exp(ePTT*(trC-N))

        GsPt%f(1) = GsPt%f_T
        GsPt%f(2) = GsPt%f(1)
        GsPt%f(3) = 0.0d0

        GsPt%g(1) = GsPt%f_T
        GsPt%g(2) = GsPt%f_T

        GsPt%dgdxi(1,:) = 0.0d0
        do i = 1, 3
            GsPt%dgdxi(1,:) = GsPt%dgdxi(1,:) + GsPt%dCdXi(i,i,:)
        end do
        GsPt%dgdxi(1,:) = ePTT*(GsPt%f_T)*(GsPt%dgdxi(1,:))
        GsPt%dgdxi(2,:) = GsPt%dgdxi(1,:)

    end subroutine conEQ_m_e_PTT_SQRT

    !---------------------------------------------------------------------

    subroutine conEQ_Giesekus_SQRT(GsPt)

        use PhysicalParameters_mod, only: aG => VeN

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        GsPt%f_T = 1.0d0-2.0d0*aG

        GsPt%f(1) = 1.0d0-aG
        GsPt%f(2) = 1.0d0-2.0d0*aG
        GsPt%f(3) = aG

        GsPt%g(1) = 1.0d0
        GsPt%g(2) = 1.0d0

        GsPt%dgdxi(1,:) = 0.0d0
        GsPt%dgdxi(2,:) = 0.0d0

    end subroutine conEQ_Giesekus_SQRT

    !---------------------------------------------------------------------

    subroutine conEQ_FENE_CR_SQRT(GsPt)

        use PhysicalParameters_mod, only: Lamda2 => VeN, PM => Problem_Main

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt :: i, N
        PetscReal :: trC

        N = PM%Ndim

        trC = traceTensor(GsPt%C_Tensor)

        GsPt%f_T = (Lamda2 - N)/(Lamda2 - trC)

        GsPt%f(1) = GsPt%f_T
        GsPt%f(2) = GsPt%f(1)
        GsPt%f(3) = 0.0d0

        GsPt%g(1) = GsPt%f_T
        GsPt%g(2) = GsPt%f_T

        GsPt%dgdxi(1,:) = 0.0d0
        do i = 1, 3
            GsPt%dgdxi(1,:) = GsPt%dgdxi(1,:) + GsPt%dCdXi(i,i,:)
        end do
        GsPt%dgdxi(1,:) = (GsPt%f_T)*(GsPt%dgdxi(1,:))/(Lamda2-trC)
        GsPt%dgdxi(2,:) = GsPt%dgdxi(1,:)

    end subroutine conEQ_FENE_CR_SQRT

    !---------------------------------------------------------------------

    subroutine conEQ_FENE_P_SQRT(GsPt)

        use PhysicalParameters_mod, only: Lamda2 => VeN, PM => Problem_Main

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt  :: i, N
        PetscReal :: trC

        N = PM%Ndim
        
        trC = traceTensor(GsPt%C_Tensor)

        GsPt%f_T = (Lamda2 - N)/(Lamda2 - trC)

        GsPt%f(1) = 1.0d0
        GsPt%f(2) = GsPt%f_T
        GsPt%f(3) = 0.0d0

        GsPt%g(1) = 1.0d0
        GsPt%g(2) = GsPt%f_T

        GsPt%dgdxi(1,:) = 0.0d0
        GsPt%dgdxi(2,:) = 0.0d0
        do i = 1, 3
            GsPt%dgdxi(2,:) = GsPt%dgdxi(2,:) + GsPt%dCdXi(i,i,:)
        end do
        GsPt%dgdxi(2,:) = (GsPt%f_T)*(GsPt%dgdxi(2,:))/(Lamda2-trC)

    end subroutine conEQ_FENE_P_SQRT

    !---------------------------------------------------------------------

    subroutine setCTensorSQRT(GsPt)

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt :: k

        GsPt%C_Tensor(:,:) = matmul(GsPt%S,GsPt%S)

        do k = 1, 3
            GsPt%dCdXi(:,:,k) = matmul(GsPt%dSdXi(:,:,k),GsPt%S) + matmul(GsPt%S,GsPt%dSdXi(:,:,k))
        end do

    end subroutine setCTensorSQRT

    !---------------------------------------------------------------------

    function setUpperConvectedDerivativeSQRT(GsPt) result(UCD)

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt
        PetscReal, dimension(3,3) :: UCD

        PetscInt :: k
        PetscReal, dimension(3,3) :: W

        UCD(:,:) = GsPt%dSdt(:,:)
        do k = 1, 3
            UCD(:,:) = UCD(:,:) + (GsPt%U_f(k))*(GsPt%dSdXi(:,:,k))
        end do

        W(:,:) = setOmegaTensor(GsPt)

        UCD(:,:) = UCD(:,:) - matmul(GsPt%S,GsPt%GU_tr) - matmul(W,GsPt%S)

    end function setUpperConvectedDerivativeSQRT

    !---------------------------------------------------------------------

    function setOmegaTensor(GsPt) result(W)

        use PhysicalParameters_mod, only: PM => Problem_Main

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt
        PetscReal, dimension(3,3) :: W

        PetscInt :: N
        PetscReal, dimension(PM%Ndim,PM%Ndim) :: TEMP_S, TEMP_GUT, TEMP_W

        N = PM%Ndim

        TEMP_S(:,:) = GsPt%S(1:N,1:N)
        TEMP_GUT(:,:) = GsPt%GUT_tr(1:N,1:N)
        ! TEMP_GUT(:,:) = GsPt%GUT(1:N,1:N)

        call WTensor_p(TEMP_S, TEMP_GUT, TEMP_W)
        ! select case (N)
        ! case (0)
        ! case (1)
        !     call WTensor2D(TEMP_S, TEMP_GUT, TEMP_W)
        ! case (2)
        !     call WTensor2D(TEMP_S, TEMP_GUT, TEMP_W)
        ! case (3)
        !     call WTensor3D(TEMP_S, TEMP_GUT, TEMP_W)
        ! case default
        !     write(*,'(a)') 'Wrong N in setOmegaTensor!'
        !     stop
        ! end select

        W(:,:) = 0.0d0
        W(1:N,1:N) = TEMP_W(:,:)

    end function setOmegaTensor

    !---------------------------------------------------------------------

    function setSInverse(S) result(SI)

        use PhysicalParameters_mod, only: PM => Problem_Main

        implicit none

        PetscReal, dimension(3,3), intent(in) :: S
        PetscReal, dimension(3,3) :: SI

        PetscInt :: N
        PetscReal, dimension(PM%Ndim,PM%Ndim) :: TEMP_S, TEMP_SI

        N = PM%Ndim

        TEMP_S(:,:) = S(1:N,1:N)

        SI(:,:) = 0.0d0
        call inverseTensor_p(TEMP_S,TEMP_SI)
        SI(1:N,1:N) = TEMP_SI(:,:)

    end function setSInverse
    
    !---------------------------------------------------------------------

    subroutine setStressTensorSQRT(GsPt)

        use PhysicalParameters_mod, only: WiN, BvN

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt :: k

        PetscReal, dimension(3,3), parameter :: &
            I_tensor = reshape( [1,0,0,0,1,0,0,0,1], shape = shape(I_tensor) )

        GsPt%T_ve(:,:) = (GsPt%g(2))*(GsPt%C_Tensor(:,:)) - (GsPt%g(1))*I_tensor(:,:)
        GsPt%T_ve(:,:) = ((1.0d0-BvN)/WiN)*(GsPt%T_ve(:,:))

        GsPt%P_tot(:,:) = GsPt%T_ve(:,:) + BvN*(GsPt%G_dot_tr(:,:)) - (GsPt%P)*I_tensor(:,:)

        do k = 1, 3
            GsPt%dTdXi(:,:,k) = (GsPt%dgdxi(2,k))*(GsPt%C_Tensor(:,:)) &
                              + (GsPt%g(2))*(GsPt%dCdXi(:,:,k)) &
                              - (GsPt%dgdxi(1,k))*I_Tensor(:,:)
        end do

        GsPt%dTdXi(:,:,:) = ((1.0d0-BvN)/WiN)*(GsPt%dTdXi(:,:,:))

    end subroutine setStressTensorSQRT

    !----------------------------------------------------------------------- 

    subroutine WTensor2D(S, GUT, W)

        implicit none

        PetscReal, dimension(:,:), intent(in) :: S, GUT
        PetscReal, dimension(:,:), intent(out) :: W

        PetscReal :: n, d

        W(:,:) = 0.0d0

        n = omegaNominator2D(S,GUT)
        d = omegaDenominator2D(S)

        W(1,2) = n/d
        W(2,1) = -W(1,2)

    end subroutine WTensor2D

    !----------------------------------------------------------------------- 

    function omegaNominator2DReal(S, GUT) result(n)

        implicit none

        PetscReal, dimension(2,2), intent(in) :: S, GUT
        PetscReal :: n

        n = -S(1,1)*GUT(2,1) + S(1,2)*(GUT(1,1)-GUT(2,2)) + S(2,2)*GUT(1,2)

    end function omegaNominator2DReal

    !----------------------------------------------------------------------- 

    function omegaNominator2DRealComplex(S, GUT) result(n)

        implicit none

        PetscReal, dimension(2,2), intent(in) :: S
        PetscScalar, dimension(2,2), intent(in) :: GUT
        PetscScalar :: n

        n = -S(1,1)*GUT(2,1) + S(1,2)*(GUT(1,1)-GUT(2,2)) + S(2,2)*GUT(1,2)

    end function omegaNominator2DRealComplex

    !----------------------------------------------------------------------- 

    function omegaNominator2DComplexReal(S, GUT) result(n)

        implicit none

        PetscScalar, dimension(2,2), intent(in) :: S
        PetscReal, dimension(2,2), intent(in) :: GUT
        PetscScalar :: n

        n = -S(1,1)*GUT(2,1) + S(1,2)*(GUT(1,1)-GUT(2,2)) + S(2,2)*GUT(1,2)

    end function omegaNominator2DComplexReal

    !----------------------------------------------------------------------- 

    function omegaDenominator2DReal(S) result(d)

        implicit none

        PetscReal, dimension(2,2), intent(in) :: S
        PetscReal :: d

        d = S(1,1)+S(2,2)

    end function omegaDenominator2DReal

    !----------------------------------------------------------------------- 

    function omegaDenominator2DComplex(S) result(d)

        implicit none

        PetscScalar, dimension(2,2), intent(in) :: S
        PetscScalar :: d

        d = S(1,1)+S(2,2)

    end function omegaDenominator2DComplex

    !----------------------------------------------------------------------- 

    subroutine WTensor3D(S, GUT, W)

        PetscReal, dimension(:,:), intent(in) :: S, GUT
        PetscReal, dimension(:,:), intent(out) :: W

        PetscReal :: w12, w13, w23, d, T1, T2, T3, B3, B2, B1, w1, w2, w3

        T1 = S(2,2)+S(3,3)
        T2 = S(1,1)+S(3,3)
        T3 = S(1,1)+S(2,2)

        B3 = S(1,2)
        B2 = S(1,3)
        B1 = S(2,3)

        d = T1*(T2*T3-B1*B1) - B2*(B2*T2+B1*B3) - B3*(B2*B1+B3*T3)

        w1 = (S(1,2)*GUT(1,1)-S(1,1)*GUT(2,1)) + (S(2,2)*GUT(1,2)-S(1,2)*GUT(2,2)) &
            + (S(2,3)*GUT(1,3)-S(1,3)*GUT(2,3))
        w2 = (S(1,3)*GUT(1,1)-S(1,1)*GUT(3,1)) + (S(3,3)*GUT(1,3)-S(1,3)*GUT(3,3)) &
            + (S(2,3)*GUT(1,2)-S(1,2)*GUT(3,2))
        w3 = (S(1,3)*GUT(2,1)-S(1,2)*GUT(3,1)) + (S(2,3)*GUT(2,2)-S(2,2)*GUT(3,2)) &
            + (S(3,3)*GUT(2,3)-S(2,3)*GUT(3,3))

        w12 = (  (T1*T2-B3*B3)*w1 - (B1*T1+B3*B2)*w2 + (B2*T2+B1*B3)*w3)/d
        w13 = (- (B1*T1+B3*B2)*w1 + (T1*T3-B2*B2)*w2 - (B2*B1+B3*T3)*w3)/d
        w23 = (  (B2*T2+B1*B3)*w1 - (B2*B1+B3*T3)*w2 + (T2*T3-B1*B1)*w3)/d

        W(:,:) = 0.0d0

        W(1,2) =  w12
        W(2,1) = -w12

        W(1,3) =  w13
        W(3,1) = -w13

        W(2,3) =  w23
        W(3,2) = -w23

    end subroutine WTensor3D

    ! !----------------------------------------------------------------------- 

    ! subroutine WTensor3D(S, GUT, W)

    !     use Tools_mod, only: gaussElimination

    !     implicit none

    !     PetscReal, dimension(:,:), intent(in) :: S, GUT
    !     PetscReal, dimension(:,:), intent(out) :: W

    !     PetscReal, dimension(3) :: b, x
    !     PetscReal, dimension(3,3) :: A

    !     A(:,:) = WTensor3D_LHS(S)
    !     b(:) = WTensor3D_RHS(S,GUT)

    !     x(:) = gaussElimination(A,b)

    !     W(:,:) = 0.0d0
    !     W(1,2) = x(1) ; W(2,1) = -W(1,2)
    !     W(1,3) = x(2) ; W(3,1) = -W(1,3)
    !     W(2,3) = x(3) ; W(3,2) = -W(2,3)

    ! end subroutine WTensor3D

    !-----------------------------------------------------------------------

    function WTensor3D_RHS_Real(S,GUT) result(b)

        implicit none

        PetscReal, dimension(3,3), intent(in) :: S, GUT
        PetscReal, dimension(3) :: b

        b(1) = S(1,2)*GUT(1,1) - S(1,1)*GUT(2,1) &
             + S(2,2)*GUT(1,2) - S(1,2)*GUT(2,2) &
             + S(2,3)*GUT(1,3) - S(1,3)*GUT(2,3)

        b(2) = S(1,3)*GUT(1,1)-S(1,1)*GUT(3,1) &
             + S(3,3)*GUT(1,3)-S(1,3)*GUT(3,3) &
             + S(2,3)*GUT(1,2)-S(1,2)*GUT(3,2)
                
        b(3) = S(1,3)*GUT(2,1)-S(1,2)*GUT(3,1) &
             + S(2,3)*GUT(2,2)-S(2,2)*GUT(3,2) &
             + S(3,3)*GUT(2,3)-S(2,3)*GUT(3,3)

    end function WTensor3D_RHS_Real

    !-----------------------------------------------------------------------

    function WTensor3D_RHS_RealComplex(S,GUT) result(b)

        implicit none

        PetscReal, dimension(3,3), intent(in) :: S
        PetscScalar, dimension(3,3), intent(in) :: GUT
        PetscScalar, dimension(3) :: b

        b(1) = S(1,2)*GUT(1,1) - S(1,1)*GUT(2,1) &
             + S(2,2)*GUT(1,2) - S(1,2)*GUT(2,2) &
             + S(2,3)*GUT(1,3) - S(1,3)*GUT(2,3)

        b(2) = S(1,3)*GUT(1,1)-S(1,1)*GUT(3,1) &
             + S(3,3)*GUT(1,3)-S(1,3)*GUT(3,3) &
             + S(2,3)*GUT(1,2)-S(1,2)*GUT(3,2)
                
        b(3) = S(1,3)*GUT(2,1)-S(1,2)*GUT(3,1) &
             + S(2,3)*GUT(2,2)-S(2,2)*GUT(3,2) &
             + S(3,3)*GUT(2,3)-S(2,3)*GUT(3,3)

    end function WTensor3D_RHS_RealComplex

    !-----------------------------------------------------------------------

    function WTensor3D_RHS_ComplexReal(S,GUT) result(b)

        implicit none

        PetscScalar, dimension(3,3), intent(in) :: S
        PetscReal, dimension(3,3), intent(in) :: GUT
        PetscScalar, dimension(3) :: b

        b(1) = S(1,2)*GUT(1,1) - S(1,1)*GUT(2,1) &
             + S(2,2)*GUT(1,2) - S(1,2)*GUT(2,2) &
             + S(2,3)*GUT(1,3) - S(1,3)*GUT(2,3)

        b(2) = S(1,3)*GUT(1,1)-S(1,1)*GUT(3,1) &
             + S(3,3)*GUT(1,3)-S(1,3)*GUT(3,3) &
             + S(2,3)*GUT(1,2)-S(1,2)*GUT(3,2)
                
        b(3) = S(1,3)*GUT(2,1)-S(1,2)*GUT(3,1) &
             + S(2,3)*GUT(2,2)-S(2,2)*GUT(3,2) &
             + S(3,3)*GUT(2,3)-S(2,3)*GUT(3,3)

    end function WTensor3D_RHS_ComplexReal

    !-----------------------------------------------------------------------

    function WTensor3D_LHS_Real(S) result(A)

        implicit none

        PetscReal, dimension(3,3), intent(in) :: S
        PetscReal, dimension(3,3) :: A

        A(1,1) =  S(1,1)+S(2,2) ; A(1,2) = S(2,3)        ; A(1,3) = -S(1,3)
        A(2,1) =  S(2,3)        ; A(2,2) = S(1,1)+S(3,3) ; A(2,3) =  S(1,2)
        A(3,1) = -S(1,3)        ; A(3,2) = S(1,2)        ; A(3,3) =  S(2,2)+S(3,3)

    end function WTensor3D_LHS_Real

    !-----------------------------------------------------------------------

    function WTensor3D_LHS_Complex(S) result(A)

        implicit none

        PetscScalar, dimension(3,3), intent(in) :: S
        PetscScalar, dimension(3,3) :: A

        A(1,1) =  S(1,1)+S(2,2) ; A(1,2) = S(2,3)        ; A(1,3) = -S(1,3)
        A(2,1) =  S(2,3)        ; A(2,2) = S(1,1)+S(3,3) ; A(2,3) =  S(1,2)
        A(3,1) = -S(1,3)        ; A(3,2) = S(1,2)        ; A(3,3) =  S(2,2)+S(3,3)

    end function WTensor3D_LHS_Complex
    
end module ConstitutiveModelsSQRT_mod
