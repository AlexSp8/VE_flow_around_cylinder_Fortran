
#include <petsc/finclude/petscksp.h>

module ConstitutiveModelsSQRTStability_mod

    use petscksp
    use Tools_mod, only: traceTensor
    use ElementVariables_mod, only: GaussPointQuantities

    implicit none

    interface
        subroutine conEQ_Stability(GsPt_b, GsPt_d)
            use ElementVariables_mod, only: GaussPointQuantities
            implicit none
        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d
        end subroutine conEQ_Stability

        subroutine W_Tensor_Stability(S_b, S_d, GUT_b, GUT_d, W_b, W_d)
            implicit none
            PetscReal, dimension(:,:), intent(in) :: S_b, GUT_b, W_b
            PetscScalar, dimension(:,:), intent(in) :: S_d, GUT_d
            PetscScalar, dimension(:,:), intent(out) :: W_d
        end subroutine W_Tensor_Stability
    end interface

    procedure(conEQ_Stability), pointer :: conEQ_Stability_p
    procedure(W_Tensor_Stability), pointer :: WTensorStability_p

    contains

    subroutine setConEQ_SQRT_Stability()

        use PhysicalParameters_mod, only: Model_Name

        implicit none

        select case (Model_Name)
        case ('FENE-CR')
            conEQ_Stability_p => conEQ_FENE_CR_Stability_SQRT
        case ('FENE-P')
            conEQ_Stability_p => conEQ_FENE_P_Stability_SQRT
        case ('L-PTT')
            conEQ_Stability_p => conEQ_L_PTT_Stability_SQRT
        case ('m-L-PTT')
            conEQ_Stability_p => conEQ_m_L_PTT_Stability_SQRT
        case ('e-PTT')
            conEQ_Stability_p => conEQ_e_PTT_Stability_SQRT
        case ('m-e-PTT')
            conEQ_Stability_p => conEQ_m_e_PTT_Stability_SQRT
        case ('GIESEKUS')
            conEQ_Stability_p => conEQ_Giesekus_Stability_SQRT
        case default
            write(*,'(a)') 'Wrong Model_Name in setConEQ_SQRT_Stability!'
            stop
        end select

    end subroutine setConEQ_SQRT_Stability

    !---------------------------------------------------------------------

    subroutine setPointersSQRT_Stability(Ndim)

        implicit none

        PetscInt, intent(in) :: Ndim

        select case (Ndim)
        case (0)
            WTensorStability_p => NULL()
        case (1)
            WTensorStability_p => WTensor2DStability
        case (2)
            WTensorStability_p => WTensor2DStability
        case (3)
            WTensorStability_p => WTensor3DStability
        case default
            write(*,'(a)') 'Wrong Ndim in setPointersSQRT!'
            stop
        end select

    end subroutine setPointersSQRT_Stability

    !---------------------------------------------------------------------

    subroutine strongFormConEQ_Stability_SQRT(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: WiN
        use ConstitutiveModelsSQRT_mod, only: setSInverse

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscReal, dimension(3,3) :: SI_b, S3_b, S_b_real
        PetscScalar, dimension(3,3) :: UCD_d, SI_d, S3_d, term_S

        call setCTensorStabilitySQRT(GsPt_b, GsPt_d)

        call conEQ_Stability_p(GsPt_b, GsPt_d)

        call setStressTensorStabilitySQRT(GsPt_b, GsPt_d)

        UCD_d(:,:) = setUCDStabilitySQRT(GsPt_b, GsPt_d)

        S_b_real(:,:) = PetscRealPart(GsPt_b%S)
        SI_b(:,:) = setSInverse(S_b_real)
        SI_d(:,:) = setSInverseStability(SI_b, GsPt_d%S)

        S3_b(:,:) = matmul( GsPt_b%S,matmul(GsPt_b%S,GsPt_b%S) )
        S3_d(:,:) = matmul( matmul(GsPt_b%S,GsPt_b%S),GsPt_d%S ) &
                  + matmul( matmul(GsPt_b%S,GsPt_d%S),GsPt_b%S ) &
                  + matmul( GsPt_d%S,matmul(GsPt_b%S,GsPt_b%S) )

        term_S(:,:) = (GsPt_b%f(3))*S3_d(:,:) + (GsPt_d%f(3))*S3_b(:,:) &
                    + (GsPt_d%f(2))*(GsPt_b%S(:,:)) + (GsPt_b%f(2))*(GsPt_d%S(:,:)) &
                    - (GsPt_b%f(1))*SI_d(:,:) - (GsPt_d%f(1))*SI_b(:,:)

        GsPt_d%Con_EQ(:,:) = WiN*UCD_d(:,:) + term_S(:,:)/2.0d0

        GsPt_d%Con_EQ_time(:,:) = WiN*(GsPt_d%S(:,:))

    end subroutine strongFormConEQ_Stability_SQRT

    !---------------------------------------------------------------------

    subroutine setCTensorStabilitySQRT(GsPt_b, GsPt_d)

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscInt :: k

        GsPt_d%C_Tensor(:,:) = matmul(GsPt_b%S,GsPt_d%S)+matmul(GsPt_d%S,GsPt_b%S)

        do k = 1, 3
            GsPt_d%dCdXi(:,:,k) = matmul(GsPt_b%dSdXi(:,:,k),GsPt_d%S) &
                                + matmul(GsPt_b%S,GsPt_d%dSdXi(:,:,k)) &
                                + matmul(GsPt_d%dSdXi(:,:,k),GsPt_b%S) &
                                + matmul(GsPt_d%S,GsPt_b%dSdXi(:,:,k))
        end do

    end subroutine setCTensorStabilitySQRT

    !---------------------------------------------------------------------

    subroutine conEQ_L_PTT_Stability_SQRT(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: ePTT => VeN

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscScalar :: trC_d

        trC_d = traceTensor(GsPt_d%C_Tensor)

        GsPt_d%f_T = ePTT*trC_d

        GsPt_d%f(1) = GsPt_d%f_T
        GsPt_d%f(2) = GsPt_d%f(1)
        GsPt_d%f(3) = 0.0d0

        GsPt_d%g(1) = 0.0d0
        GsPt_d%g(2) = 0.0d0
        GsPt_d%dgdxi(1,:) = 0.0d0
        GsPt_d%dgdxi(2,:) = 0.0d0

    end subroutine conEQ_L_PTT_Stability_SQRT

    !---------------------------------------------------------------------

    subroutine conEQ_m_L_PTT_Stability_SQRT(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: ePTT => VeN

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscInt :: i
        PetscScalar :: trC_d

        trC_d = traceTensor(GsPt_d%C_Tensor)

        GsPt_d%f_T = ePTT*trC_d

        GsPt_d%f(1) = GsPt_d%f_T
        GsPt_d%f(2) = GsPt_d%f(1)
        GsPt_d%f(3) = 0.0d0

        GsPt_d%g(1) = GsPt_d%f_T
        GsPt_d%g(2) = GsPt_d%f_T

        GsPt_d%dgdxi(1,:) = 0.0d0
        do i = 1, 3
            GsPt_d%dgdxi(1,:) = GsPt_d%dgdxi(1,:) + GsPt_d%dCdXi(i,i,:)
        end do
        GsPt_d%dgdxi(1,:) = ePTT*(GsPt_d%dgdxi(1,:))
        GsPt_d%dgdxi(2,:) = GsPt_d%dgdxi(1,:)

    end subroutine conEQ_m_L_PTT_Stability_SQRT

    !---------------------------------------------------------------------

    subroutine conEQ_e_PTT_Stability_SQRT(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: ePTT => VeN

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscScalar :: trC_d

        trC_d = traceTensor(GsPt_d%C_Tensor)

        GsPt_d%f_T = (GsPt_b%f_T)*ePTT*trC_d

        GsPt_d%f(1) = GsPt_d%f_T
        GsPt_d%f(2) = GsPt_d%f(1)
        GsPt_d%f(3) = 0.0d0

        GsPt_d%g(1) = 0.0d0
        GsPt_d%g(2) = 0.0d0
        GsPt_d%dgdxi(1,:) = 0.0d0
        GsPt_d%dgdxi(2,:) = 0.0d0

    end subroutine conEQ_e_PTT_Stability_SQRT

    !---------------------------------------------------------------------

    subroutine conEQ_m_e_PTT_Stability_SQRT(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: ePTT => VeN

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscInt :: i
        PetscScalar :: trC_d

        trC_d = traceTensor(GsPt_d%C_Tensor)

        GsPt_d%f_T = (GsPt_b%f_T)*ePTT*trC_d

        GsPt_d%f(1) = GsPt_d%f_T
        GsPt_d%f(2) = GsPt_d%f(1)
        GsPt_d%f(3) = 0.0d0

        GsPt_d%g(1) = GsPt_d%f_T
        GsPt_d%g(2) = GsPt_d%f_T
        
        GsPt_d%dgdxi(1,:) = 0.0d0
        do i = 1, 3
            GsPt_d%dgdxi(1,:) = GsPt_d%dgdxi(1,:) + GsPt_d%dCdXi(i,i,:)
        end do
        GsPt_d%dgdxi(1,:) = (GsPt_b%f_T)*ePTT*(GsPt_d%dgdxi(1,:))
        GsPt_d%dgdxi(1,:) = GsPt_d%dgdxi(1,:) + (GsPt_b%dgdxi(1,:))*(GsPt_d%f_T)
        GsPt_d%dgdxi(2,:) = GsPt_d%dgdxi(1,:)

    end subroutine conEQ_m_e_PTT_Stability_SQRT

    !---------------------------------------------------------------------

    subroutine conEQ_Giesekus_Stability_SQRT(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: aG => VeN

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        GsPt_d%f_T = 0.0d0

        GsPt_d%f(1) = 0.0d0
        GsPt_d%f(2) = 0.0d0
        GsPt_d%f(3) = 0.0d0

        GsPt_d%g(1) = 0.0d0
        GsPt_d%g(2) = 0.0d0

        GsPt_d%dgdxi(1,:) = 0.0d0
        GsPt_d%dgdxi(2,:) = 0.0d0

    end subroutine conEQ_Giesekus_Stability_SQRT

    !---------------------------------------------------------------------

    subroutine conEQ_FENE_CR_Stability_SQRT(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: Lamda2 => VeN

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscInt :: i
        PetscReal :: trC_b
        PetscScalar :: trC_d

        trC_d = traceTensor(GsPt_d%C_Tensor)
        trC_b = traceTensor(GsPt_b%C_Tensor)

        GsPt_d%f_T = (GsPt_b%f_T)*trC_d/(Lamda2 - trC_b)

        GsPt_d%f(1) = GsPt_d%f_T
        GsPt_d%f(2) = GsPt_d%f(1)
        GsPt_d%f(3) = 0.0d0

        GsPt_d%g(1) = GsPt_d%f_T
        GsPt_d%g(2) = GsPt_d%f_T

        GsPt_d%dgdxi(1,:) = 0.0d0
        do i = 1, 3
            GsPt_d%dgdxi(1,:) = GsPt_d%dgdxi(1,:) + GsPt_d%dCdXi(i,i,:)
        end do
        GsPt_d%dgdxi(1,:) = (GsPt_b%f_T)*(GsPt_d%dgdxi(1,:)) + (GsPt_b%dgdxi(1,:))*trC_d
        GsPt_d%dgdxi(1,:) = (GsPt_d%dgdxi(1,:))/(Lamda2-trC_b)
        GsPt_d%dgdxi(1,:) = (GsPt_d%dgdxi(1,:)) + (GsPt_b%dgdxi(1,:))*(GsPt_d%f_T)/(GsPt_b%f_T)
        GsPt_d%dgdxi(2,:) = GsPt_d%dgdxi(1,:)

    end subroutine conEQ_FENE_CR_Stability_SQRT

    !---------------------------------------------------------------------

    subroutine conEQ_FENE_P_Stability_SQRT(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: Lamda2 => VeN

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscInt :: i
        PetscReal :: trC_b
        PetscScalar :: trC_d

        trC_d = traceTensor(GsPt_d%C_Tensor)
        trC_b = traceTensor(GsPt_b%C_Tensor)

        GsPt_d%f_T = (GsPt_b%f_T)*trC_d/(Lamda2 - trC_b)

        GsPt_d%f(1) = 0.0d0
        GsPt_d%f(2) = GsPt_d%f(1)
        GsPt_d%f(3) = 0.0d0

        GsPt_d%g(1) = 0.0d0
        GsPt_d%g(2) = GsPt_d%f_T

        GsPt_d%dgdxi(1,:) = 0.0d0
        do i = 1, 3
            GsPt_d%dgdxi(2,:) = GsPt_d%dgdxi(2,:) + GsPt_d%dCdXi(i,i,:)
        end do
        GsPt_d%dgdxi(2,:) = (GsPt_b%f_T)*(GsPt_d%dgdxi(2,:)) + (GsPt_b%dgdxi(2,:))*trC_d
        GsPt_d%dgdxi(2,:) = (GsPt_d%dgdxi(2,:))/(Lamda2-trC_b)
        GsPt_d%dgdxi(2,:) = (GsPt_d%dgdxi(2,:)) + (GsPt_b%dgdxi(2,:))*(GsPt_d%f_T)/(GsPt_b%f_T)

    end subroutine conEQ_FENE_P_Stability_SQRT

    !---------------------------------------------------------------------

    function setUCDStabilitySQRT(GsPt_b, GsPt_d) result(UCD_d)

        use PhysicalParameters_mod, only: PS => Problem_Main_Stability
        use ConstitutiveModelsSQRT_mod, only: setOmegaTensor

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d
        PetscScalar, dimension(3,3) :: UCD_d

        PetscInt :: k, N
        PetscReal, dimension(3,3) :: W_b
        PetscScalar, dimension(3,3) :: W_d
        PetscReal, dimension(PS%Ndim,PS%Ndim) :: TEMP_S_b, TEMP_GUT_b
        PetscScalar, dimension(PS%Ndim,PS%Ndim) :: TEMP_S_d, TEMP_GUT_d, TEMP_W_d

        UCD_d(:,:) = -matmul(GsPt_b%S,GsPt_d%GU)-matmul(GsPt_d%S,GsPt_b%GU)
        do k = 1, 3
            UCD_d(:,:) = UCD_d(:,:) + (GsPt_b%U_f(k))*(GsPt_d%dSdXi(:,:,k)) &
                                    + (GsPt_d%U_f(k))*(GsPt_b%dSdXi(:,:,k))
        end do

        W_b(:,:) = setOmegaTensor(GsPt_b)

        N = PS%Ndim
        TEMP_S_b = GsPt_b%S(1:N,1:N)
        TEMP_S_d = GsPt_d%S(1:N,1:N)
        TEMP_GUT_b = GsPt_b%GUT(1:N,1:N)
        TEMP_GUT_d = GsPt_d%GUT(1:N,1:N)

        call WTensorStability_p(TEMP_S_b, TEMP_S_d, TEMP_GUT_b, TEMP_GUT_d, W_b, TEMP_W_d)
        W_d(:,:) = 0.0d0
        W_d(1:N,1:N) = TEMP_W_d(:,:)

        UCD_d(:,:) = UCD_d(:,:) - matmul(W_b,GsPt_d%S)- matmul(W_d,GsPt_b%S)

    end function setUCDStabilitySQRT

    !---------------------------------------------------------------------

    function setSInverseStability(SI_b, S_d) result(SI_d)

        implicit none

        PetscReal, dimension(3,3), intent(in) :: SI_b
        PetscScalar, dimension(3,3), intent(in) :: S_d
        PetscScalar, dimension(3,3) :: SI_d

        SI_d(:,:) = -matmul(SI_b,matmul(S_d,SI_b))

    end function setSInverseStability

    !---------------------------------------------------------------------

    subroutine setStressTensorStabilitySQRT(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: WiN, BvN

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscInt :: k
        PetscReal, dimension(3,3), parameter :: &
            I_tensor = reshape( [1,0,0,0,1,0,0,0,1], shape = shape(I_tensor) )

        GsPt_d%T_ve(:,:) = (GsPt_b%g(2))*(GsPt_d%C_Tensor(:,:)) &
                         + (GsPt_d%g(2))*(GsPt_b%C_Tensor(:,:)) - (GsPt_d%g(1))*I_tensor(:,:)
        GsPt_d%T_ve(:,:) = (1.0d0-BvN)*(GsPt_d%T_ve(:,:))/WiN

        GsPt_d%P_tot(:,:) = GsPt_d%T_ve(:,:) + BvN*(GsPt_d%G_dot(:,:)) - (GsPt_d%P)*I_tensor(:,:)

        do k = 1, 3
            GsPt_d%dTdXi(:,:,k) = (GsPt_b%dgdxi(2,k))*(GsPt_d%C_Tensor(:,:)) &
                                + (GsPt_d%dgdxi(2,k))*(GsPt_b%C_Tensor(:,:)) &
                                + (GsPt_b%g(2))*(GsPt_d%dCdXi(:,:,k)) &
                                + (GsPt_d%g(2))*(GsPt_b%dCdXi(:,:,k)) &
                                - (GsPt_d%dgdxi(1,k))*I_Tensor(:,:)
        end do
        GsPt_d%dTdXi(:,:,:) = (1.0d0-BvN)*(GsPt_d%dTdXi(:,:,:))/WiN

    end subroutine setStressTensorStabilitySQRT

    !----------------------------------------------------------------------- 

    subroutine WTensor2DStability(S_b, S_d, GUT_b, GUT_d, W_b, W_d)

        use ConstitutiveModelsSQRT_mod, only: omegaNominator2D, omegaDenominator2D
        
        implicit none

        PetscReal, dimension(:,:), intent(in) :: S_b, GUT_b, W_b
        PetscScalar, dimension(:,:), intent(in) :: S_d, GUT_d
        PetscScalar, dimension(:,:), intent(out) :: W_d

        PetscReal :: n_b, d_b
        PetscScalar :: n_d, d_d

        n_b = omegaNominator2D(S_b,GUT_b)
        d_b = omegaDenominator2D(S_b)

        n_d = omegaNominator2D(S_b,GUT_d) + omegaNominator2D(S_d,GUT_b)
        d_d = omegaDenominator2D(S_d)

        W_d(:,:) = 0.0d0
        W_d(1,2) = (n_d*d_b-n_b*d_d)/(d_b**2)
        W_d(2,1) = -W_d(1,2)

    end subroutine WTensor2DStability

    !-----------------------------------------------------------------------

    subroutine WTensor3DStability(S_b, S_d, GUT_b, GUT_d, W_b, W_d)

        use ConstitutiveModelsSQRT_mod, only: WTensor3D_RHS, WTensor3D_LHS
        use Tools_mod, only: gaussElimination

        implicit none

        PetscReal, dimension(:,:), intent(in) :: S_b, GUT_b, W_b
        PetscScalar, dimension(:,:), intent(in) :: S_d, GUT_d
        PetscScalar, dimension(:,:), intent(out) :: W_d

        PetscReal, dimension(3,3) :: A
        PetscScalar, dimension(3) :: x, b

        A(:,:) = WTensor3D_LHS(S_b)

        b(:) = WTensor3D_RHS(S_b,GUT_d) + WTensor3D_RHS(S_d,GUT_b)

        b(1) = b(1) - W_b(1,2)*(S_d(1,1)+S_d(2,2))
        b(2) = b(2) - W_b(1,2)*S_d(2,3)
        b(3) = b(3) - W_b(1,2)*(-S_d(1,3))

        x(:) = gaussElimination(A,b)

        W_d(:,:) = 0.0d0
        W_d(1,2) = x(1) ; W_d(2,1) = -W_d(1,2)
        W_d(1,3) = x(2) ; W_d(3,1) = -W_d(1,3)
        W_d(2,3) = x(3) ; W_d(3,2) = -W_d(2,3)

    end subroutine WTensor3DStability

end module ConstitutiveModelsSQRTStability_mod
