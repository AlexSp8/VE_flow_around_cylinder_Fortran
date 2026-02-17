
#include <petsc/finclude/petscksp.h>

module ConstitutiveModelsStability_mod

    use petscksp
    use Tools_mod, only: traceTensor
    use PhysicalParameters_mod, only: WiN, BvN
    use ElementVariables_mod, only: GaussPointQuantities

    implicit none

    interface
        subroutine conEQ_Stability(GsPt_b, GsPt_d)
            use ElementVariables_mod, only: GaussPointQuantities
            implicit none
        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d
        end subroutine conEQ_Stability
    end interface

    procedure(conEQ_Stability), pointer :: conEQ_Stability_p

    contains

    subroutine setConEQ_T_Stability()

        use PhysicalParameters_mod, only: Model_Name

        implicit none

        select case (Model_Name)
        case ('FENE-CR')
            conEQ_Stability_p => conEQ_FENE_CR_Stability
        case ('FENE-P')
            conEQ_Stability_p => conEQ_FENE_P_Stability
        case ('L-PTT')
            conEQ_Stability_p => conEQ_L_PTT_Stability
        case ('m-L-PTT')
            conEQ_Stability_p => conEQ_m_L_PTT_Stability
        case ('e-PTT')
            conEQ_Stability_p => conEQ_e_PTT_Stability
        case ('m-e-PTT')
            conEQ_Stability_p => conEQ_m_e_PTT_Stability
        case ('GIESEKUS')
            conEQ_Stability_p => conEQ_Giesekus_Stability
        case default
            write(*,'(a)') 'Wrong Model_Name in setConEQ_T_Stability!'
            stop
        end select

    end subroutine setConEQ_T_Stability

    !---------------------------------------------------------------------

    subroutine strongFormConEQ_Stability(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: Model_Name

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscScalar :: term_I, MDfDt_d
        PetscScalar, dimension(3,3) :: UCD_d, term_T, term_G, term_dt
        PetscReal, dimension(3,3), parameter :: &
            I_tensor = reshape( [1,0,0,0,1,0,0,0,1], shape = shape(I_tensor) )

        GsPt_d%T_ve(:,:) = GsPt_d%S(:,:)
        
        GsPt_d%dTdXi(:,:,:) = GsPt_d%dSdXi(:,:,:)

        call conEQ_Stability_p(GsPt_b, GsPt_d)

        call setStressTensorStability(GsPt_b,GsPt_d)

        UCD_d(:,:) = setUCDStability(GsPt_b, GsPt_d)

        term_T(:,:) = (GsPt_b%f(3))*( matmul(GsPt_b%S,GsPt_d%S)+matmul(GsPt_d%S,GsPt_b%S) ) &
                    + (GsPt_d%f(3))*matmul(GsPt_b%S,GsPt_b%S) &
                    + (GsPt_d%f(2))*GsPt_b%S(:,:) + (GsPt_b%f(2))*GsPt_d%S(:,:)

        term_I = 0.0d0
        if (Model_Name == 'FENE-P') then
            MDfDt_d = setFENEMaterialDerivativeStability(GsPt_b, GsPt_d)
            term_I = MDfDt_d/(GsPt_b%f_T)
        end if

        term_G(:,:) = (GsPt_b%f(1))*(GsPt_d%G_dot(:,:)) + (GsPt_d%f(1))*(GsPt_b%G_dot(:,:))

        GsPt_d%Con_EQ(:,:) = WiN*UCD_d(:,:) + term_T(:,:) &
                            - (1.0d0-BvN)*( term_G(:,:)+term_I*I_tensor(:,:) )

        term_I = 0.0d0
        term_T(:,:) = WiN*(GsPt_d%S(:,:))
        if (Model_Name == 'FENE-CR') then
            term_T(:,:) = term_T(:,:) - WiN*(GsPt_d%f_T)*GsPt_b%S(:,:)/(GsPt_b%f_T)
        end if
        if (Model_Name == 'FENE-P') then
            term_I = -(1.0d0-BvN)*(GsPt_d%f_T)/(GsPt_b%f_T)
            term_T(:,:) = term_T(:,:) - WiN*(GsPt_d%f_T)*GsPt_b%S(:,:)/(GsPt_b%f_T)
        end if
        GsPt_d%Con_EQ_time(:,:) = term_T(:,:) + term_I*I_tensor(:,:)

    end subroutine strongFormConEQ_Stability

    !----------------------------------------------------------------------

    subroutine conEQ_L_PTT_Stability(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: ePTT => VeN

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscScalar :: trT_d

        trT_d = traceTensor(GsPt_d%S)

        GsPt_d%f_T = ePTT*WiN*trT_d/(1.0d0-BvN)

        GsPt_d%f(1) = 0.0d0
        GsPt_d%f(2) = GsPt_d%f_T
        GsPt_d%f(3) = 0.0d0

        GsPt_d%g(1) = 0.0d0
        GsPt_d%g(2) = 0.0d0

    end subroutine conEQ_L_PTT_Stability

    !----------------------------------------------------------------------

    subroutine conEQ_m_L_PTT_Stability(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: ePTT => VeN

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscScalar :: trT_d

        trT_d = traceTensor(GsPt_d%S)

        GsPt_d%f_T = ePTT*WiN*trT_d/(1.0d0-BvN)

        GsPt_d%f(1) = GsPt_d%f_T
        GsPt_d%f(2) = GsPt_d%f(1)
        GsPt_d%f(3) = 0.0d0

        GsPt_d%g(1) = GsPt_d%f_T
        GsPt_d%g(2) = GsPt_d%g(1)

    end subroutine conEQ_m_L_PTT_Stability

    !----------------------------------------------------------------------

    subroutine conEQ_e_PTT_Stability(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: ePTT => VeN

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscScalar :: trT_d

        trT_d = traceTensor(GsPt_d%S)

        GsPt_d%f_T = (GsPt_b%f_T)*ePTT*WiN*trT_d/(1.0d0-BvN)

        GsPt_d%f(1) = 0.0d0
        GsPt_d%f(2) = GsPt_d%f_T
        GsPt_d%f(3) = 0.0d0

        GsPt_d%g(1) = 0.0d0
        GsPt_d%g(2) = 0.0d0

    end subroutine conEQ_e_PTT_Stability

    !----------------------------------------------------------------------

    subroutine conEQ_m_e_PTT_Stability(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: ePTT => VeN

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscScalar :: trT_d

        trT_d = traceTensor(GsPt_d%S)

        GsPt_d%f_T = (GsPt_b%f_T)*ePTT*WiN*trT_d/(1.0d0-BvN)

        GsPt_d%f(1) = GsPt_d%f_T
        GsPt_d%f(2) = GsPt_d%f(1)
        GsPt_d%f(3) = 0.0d0

        GsPt_d%g(1) = GsPt_d%f_T
        GsPt_d%g(2) = GsPt_d%g(1)

    end subroutine conEQ_m_e_PTT_Stability

    !---------------------------------------------------------------------

    subroutine conEQ_Giesekus_Stability(GsPt_b, GsPt_d)

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

    end subroutine conEQ_Giesekus_Stability

    !---------------------------------------------------------------------

    subroutine conEQ_FENE_CR_Stability(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: Lamda2 => VeN, PS => Problem_Main

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscInt :: N
        PetscScalar :: trT_d, MDfDt_d

        N = PS%Ndim

        trT_d = traceTensor(GsPt_d%S)

        GsPt_d%f_T = (WiN*trT_d/(1.0d0-BvN))/(Lamda2 - N)

        GsPt_d%f(1) = GsPt_d%f_T
        MDfDt_d = setFENEMaterialDerivativeStability(GsPt_b, GsPt_d)
        GsPt_d%f(2) = GsPt_d%f_T - WiN*MDfDt_d/(GsPt_b%f_T)
        GsPt_d%f(3) = 0.0d0

        GsPt_d%g(1) = GsPt_d%f_T
        GsPt_d%g(2) = GsPt_d%g(1)

    end subroutine conEQ_FENE_CR_Stability

    !----------------------------------------------------------------------

    subroutine conEQ_FENE_P_Stability(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: Lamda2 => VeN

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscScalar :: trT_d, MDfDt_d

        trT_d = traceTensor(GsPt_d%S)

        GsPt_d%f_T = (WiN*trT_d/(1.0d0-BvN))/Lamda2

        GsPt_d%f(1) = 0.0d0
        MDfDt_d = setFENEMaterialDerivativeStability(GsPt_b, GsPt_d)
        GsPt_d%f(2) = GsPt_d%f_T - WiN*MDfDt_d/(GsPt_b%f_T)
        GsPt_d%f(3) = 0.0d0

        GsPt_d%g(1) = 0.0d0
        GsPt_d%g(2) = GsPt_d%f_T

    end subroutine conEQ_FENE_P_Stability

    !---------------------------------------------------------------------

    function setFENEMaterialDerivativeStability(GsPt_b, GsPt_d) result(MDfDt_d)

        use PhysicalParameters_mod, only: Lamda2 => VeN, Model_Name, PS => Problem_Main

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d
        PetscScalar :: MDfDt_d

        PetscInt  :: i, N
        PetscScalar :: term
        PetscReal, dimension(3) :: dfdXi_b
        PetscScalar, dimension(3) :: dfdXi_d

        dfdXi_b(:) = 0.0d0 ; dfdXi_d(:) = 0.0d0
        do i = 1, 3
            dfdXi_b(:) = dfdXi_b(:) + GsPt_b%dTdXi(i,i,:)
            dfdXi_d(:) = dfdXi_d(:) + GsPt_d%dTdXi(i,i,:)
        end do
        if (Model_Name == 'FENE-P') then
            dfdXi_b(:) = (WiN/(1.0d0-BvN))*dfdXi_b(:)/Lamda2
            dfdXi_d(:) = (WiN/(1.0d0-BvN))*dfdXi_d(:)/Lamda2
        else
            N = PS%Ndim
            dfdXi_b(:) = (WiN/(1.0d0-BvN))*dfdXi_b(:)/(Lamda2-N)
            dfdXi_d(:) = (WiN/(1.0d0-BvN))*dfdXi_d(:)/(Lamda2-N)
        end if

        MDfDt_d = 0.0d0
        do i = 1, 3
            term = ( (GsPt_b%f_T)*dfdXi_d(i)-(GsPt_d%f_T)*dfdXi_b(i) )/(GsPt_b%f_T)
            MDfDt_d = MDfDt_d + (GsPt_d%U_f(i))*dfdXi_b(i) + (GsPt_b%U_f(i))*term
        end do

    end function setFENEMaterialDerivativeStability

    !----------------------------------------------------------------------

    function setUCDStability(GsPt_b, GsPt_d) result(UCD_d)

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d
        PetscScalar, dimension(3,3) :: UCD_d

        PetscInt :: k

        UCD_d(:,:) = 0.0d0
        do k = 1, 3
            UCD_d(:,:) = UCD_d(:,:) + (GsPt_b%U_f(k))*(GsPt_d%dTdXi(:,:,k)) &
                                    + (GsPt_d%U_f(k))*(GsPt_b%dTdXi(:,:,k))
        end do

        UCD_d(:,:) = UCD_d(:,:) - matmul(GsPt_b%GUT,GsPt_d%S) - matmul(GsPt_d%GUT,GsPt_b%S) &
                                - matmul(GsPt_b%S,GsPt_d%GU)  - matmul(GsPt_d%S,GsPt_b%GU)

    end function setUCDStability

    !---------------------------------------------------------------------

    subroutine setStressTensorStability(GsPt_b,GsPt_d)

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscReal, dimension(3,3), parameter :: &
            I_tensor = reshape( [1,0,0,0,1,0,0,0,1], shape = shape(I_tensor) )

        GsPt_d%P_tot(:,:) = GsPt_d%T_ve(:,:) + BvN*(GsPt_d%G_dot(:,:)) - (GsPt_d%P)*I_tensor(:,:)

        GsPt_d%C_Tensor(:,:) = WiN*(GsPt_d%T_ve(:,:)/(1.0d0-BvN)) + GsPt_d%g(1)*I_tensor(:,:) &
                            - (GsPt_d%g(2))*(GsPt_b%C_Tensor(:,:))
        GsPt_d%C_Tensor(:,:) = (1.0d0/(GsPt_b%g(2)))*(GsPt_d%C_Tensor(:,:))

    end subroutine setStressTensorStability

end module ConstitutiveModelsStability_mod
