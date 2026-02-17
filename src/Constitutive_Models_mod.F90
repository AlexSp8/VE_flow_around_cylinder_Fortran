
#include <petsc/finclude/petscksp.h>

module ConstitutiveModels_mod

    use petscksp
    use Tools_mod, only: traceTensor
    use PhysicalParameters_mod, only: WiN, BvN
    use ElementVariables_mod, only: GaussPointQuantities

    implicit none
    
    interface
        subroutine conEQ(GsPt)
            use ElementVariables_mod, only: GaussPointQuantities
            implicit none
            type(GaussPointQuantities), intent(inout) :: GsPt
        end subroutine conEQ
    end interface

    procedure(conEQ), pointer :: conEQ_p

    contains

    subroutine setConEQ_T()

        use PhysicalParameters_mod, only: Model_Name

        implicit none

        select case (Model_Name)
        case ('FENE-CR')
            conEQ_p => conEQ_FENE_CR
        case ('FENE-P')
            conEQ_p => conEQ_FENE_P
        case ('L-PTT')
            conEQ_p => conEQ_L_PTT
        case ('m-L-PTT')
            conEQ_p => conEQ_m_L_PTT
        case ('e-PTT')
            conEQ_p => conEQ_e_PTT
        case ('m-e-PTT')
            conEQ_p => conEQ_m_e_PTT
        case ('GIESEKUS')
            conEQ_p => conEQ_Giesekus
        case default
            write(*,'(a)') 'Wrong Model_Name in setConEQ_T!'
            stop
        end select

    end subroutine setConEQ_T

    !---------------------------------------------------------------------

    subroutine strongFormConEQ(GsPt)

        use PhysicalParameters_mod, only: Model_Name

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscReal :: term_I, MDfDt
        PetscReal, dimension(3,3) :: UCD, Tve2
        PetscReal, dimension(3,3), parameter :: &
            I_tensor = reshape( [1,0,0,0,1,0,0,0,1], shape = shape(I_tensor) )

        GsPt%T_ve(:,:) = GsPt%S(:,:)

        GsPt%dTdXi(:,:,:) = GsPt%dSdXi(:,:,:)
        
        call conEQ_p(GsPt)

        call setStressTensor(GsPt)

        UCD(:,:) = setUpperConvectedDerivativeT(GsPt)
        Tve2(:,:) = matmul(GsPt%S,GsPt%S)

        term_I = 0.0d0
        if (Model_Name == 'FENE-P') then
            MDfDt = setFENEMaterialDerivative(GsPt)
            term_I = MDfDt/(GsPt%f_T)
        end if

        GsPt%Con_EQ(:,:) = WiN*UCD(:,:) + (GsPt%f(3))*Tve2(:,:)
        GsPt%Con_EQ(:,:) = GsPt%Con_EQ(:,:) + (GsPt%f(2))*GsPt%S(:,:)
        GsPt%Con_EQ(:,:) = GsPt%Con_EQ(:,:) &
                         - GsPt%f(1)*(1.0d0-BvN)*( GsPt%G_dot_tr(:,:)+term_I*I_tensor(:,:) )
        
        GsPt%Con_EQ_time(:,:) = 0.0d0

    end subroutine strongFormConEQ

    !---------------------------------------------------------------------

    subroutine conEQ_L_PTT(GsPt)

        use PhysicalParameters_mod, only: ePTT => VeN

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscReal :: trT

        trT = traceTensor(GsPt%S)

        GsPt%f_T = 1.0d0 + (ePTT*WiN*trT)/(1.0d0-BvN)

        GsPt%f(1) = 1.0d0
        GsPt%f(2) = GsPt%f_T
        GsPt%f(3) = 0.0d0

        GsPt%g(1) = 1.0d0
        GsPt%g(2) = 1.0d0

    end subroutine conEQ_L_PTT

    !---------------------------------------------------------------------

    subroutine conEQ_m_L_PTT(GsPt)

        use PhysicalParameters_mod, only: ePTT => VeN

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscReal :: trT

        trT = traceTensor(GsPt%S)

        GsPt%f_T = 1.0d0 + (ePTT*WiN*trT)/(1.0d0-BvN)

        GsPt%f(1) = GsPt%f_T
        GsPt%f(2) = GsPt%f(1)
        GsPt%f(3) = 0.0d0

        GsPt%g(1) = GsPt%f_T
        GsPt%g(2) = GsPt%g(1)

    end subroutine conEQ_m_L_PTT

    !---------------------------------------------------------------------

    subroutine conEQ_e_PTT(GsPt)

        use PhysicalParameters_mod, only: ePTT => VeN

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscReal :: trT

        trT = traceTensor(GsPt%S)

        GsPt%f_T = exp(ePTT*WiN*trT/(1.0d0-BvN))

        GsPt%f(1) = 1.0d0
        GsPt%f(2) = GsPt%f_T
        GsPt%f(3) = 0.0d0

        GsPt%g(1) = 1.0d0
        GsPt%g(2) = 1.0d0

    end subroutine conEQ_e_PTT

    !---------------------------------------------------------------------

    subroutine conEQ_m_e_PTT(GsPt)

        use PhysicalParameters_mod, only: ePTT => VeN

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscReal :: trT

        trT = traceTensor(GsPt%S)

        GsPt%f_T = exp(ePTT*WiN*trT/(1.0d0-BvN))

        GsPt%f(1) = GsPt%f_T
        GsPt%f(2) = GsPt%f(1)
        GsPt%f(3) = 0.0d0

        GsPt%g(1) = GsPt%f_T
        GsPt%g(2) = GsPt%g(1)

    end subroutine conEQ_m_e_PTT

    !---------------------------------------------------------------------

    subroutine conEQ_Giesekus(GsPt)

        use PhysicalParameters_mod, only: aG => VeN

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        GsPt%f_T = aG*WiN/(1.0d0-BvN)

        GsPt%f(1) = 1.0d0
        GsPt%f(2) = 1.0d0
        GsPt%f(3) = GsPt%f_T

        GsPt%g(1) = 1.0d0
        GsPt%g(2) = 1.0d0

    end subroutine conEQ_Giesekus

    !---------------------------------------------------------------------

    subroutine conEQ_FENE_CR(GsPt)

        use PhysicalParameters_mod, only: Lamda2 => VeN, PM => Problem_Main

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt :: N
        PetscReal :: trT, MDfDt

        N = PM%Ndim

        trT = traceTensor(GsPt%S)

        GsPt%f_T = 1.0d0 + ((WiN*trT/(1.0d0-BvN)))/(Lamda2-N)

        GsPt%f(1) = GsPt%f_T
        MDfDt = setFENEMaterialDerivative(GsPt)
        GsPt%f(2) = GsPt%f_T - (WiN*MDfDt/GsPt%f_T)
        GsPt%f(3) = 0.0d0

        GsPt%g(1) = GsPt%f_T
        GsPt%g(2) = GsPt%g(1)

    end subroutine conEQ_FENE_CR

    !---------------------------------------------------------------------

    subroutine conEQ_FENE_P(GsPt)

        use PhysicalParameters_mod, only: Lamda2 => VeN

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscReal :: trT, MDfDt

        trT = traceTensor(GsPt%S)

        GsPt%f_T = 1.0d0 + (WiN*trT/(1.0d0-BvN))/Lamda2

        GsPt%f(1) = 1.0d0
        MDfDt = setFENEMaterialDerivative(GsPt)
        GsPt%f(2) = GsPt%f_T - (WiN*MDfDt/GsPt%f_T)
        GsPt%f(3) = 0.0d0

        GsPt%g(1) = 1.0d0
        GsPt%g(2) = GsPt%f_T

    end subroutine conEQ_FENE_P

    !---------------------------------------------------------------------

    function setFENEMaterialDerivative(GsPt) result(MDfDt)

        use ContinuationParameters_mod, only: Continuation_Method
        use PhysicalParameters_mod, only: Lamda2 => VeN, Model_Name, PM => Problem_Main
        use ContinuationVariables_mod, only: dt

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt
        PetscReal :: MDfDt

        PetscInt  :: i, N
        PetscReal :: f_To, trTo
        PetscReal, dimension(3) :: dfdXi

        trTo = traceTensor(GsPt%So)

        dfdXi(:) = 0.0d0
        do i = 1, 3
            dfdXi(:) = dfdXi(:) + GsPt%dTdXi(i,i,:)
        end do
        if (Model_Name == 'FENE-P') then
            dfdXi(:) = (WiN/(1.0d0-BvN))*dfdXi(:)/Lamda2
            f_To = 1.0d0 + ((WiN*trTo/(1.0d0-BvN)))/Lamda2
        else
            N = PM%Ndim
            dfdXi(:) = (WiN/(1.0d0-BvN))*dfdXi(:)/(Lamda2-N)
            f_To = 1.0d0 + ((WiN*trTo/(1.0d0-BvN)))/(Lamda2-N)           
        end if

        MDfDt = 0.0d0
        if (Continuation_Method == 'Transient') then
            MDfDt = (GsPt%f_T - f_To)/dt
        end if

        do i = 1, 3
            MDfDt = MDfDt + GsPt%U_f(i)*dfdXi(i)
        end do

    end function setFENEMaterialDerivative

    !---------------------------------------------------------------------

    function setUpperConvectedDerivativeT(GsPt) result(UCD)

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt
        PetscReal, dimension(3,3) :: UCD

        PetscInt :: k

        UCD(:,:) = GsPt%dSdt(:,:)
        do k = 1, 3
            UCD(:,:) = UCD(:,:) + GsPt%U_f(k)*GsPt%dSdXi(:,:,k)
        end do

        UCD(:,:) = UCD(:,:) - matmul(GsPt%GUT_tr,GsPt%S) - matmul(GsPt%S,GsPt%GU_tr)

    end function setUpperConvectedDerivativeT

    !---------------------------------------------------------------------

    subroutine setStressTensor(GsPt)

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscReal, dimension(3,3), parameter :: &
            I_tensor = reshape( [1,0,0,0,1,0,0,0,1], shape = shape(I_tensor) )

        GsPt%P_tot(:,:) = GsPt%T_ve(:,:) + BvN*GsPt%G_dot_tr(:,:) - GsPt%P*I_tensor(:,:)

        GsPt%C_Tensor(:,:) = WiN*(GsPt%T_ve(:,:)/(1.0d0-BvN)) + GsPt%g(1)*I_tensor(:,:)
        GsPt%C_Tensor(:,:) = (1.0d0/GsPt%g(2))*GsPt%C_Tensor(:,:)

    end subroutine setStressTensor

end module ConstitutiveModels_mod
