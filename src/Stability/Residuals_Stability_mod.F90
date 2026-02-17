
#include <slepc/finclude/slepceps.h>

module ResidualsStability_mod

    use PhysicalParameters_mod, only: Problem_Main, Problem_Main_Stability
    use ElementVariables_mod, only: GaussPointQuantities, Elem => NodeArrays_Main

    implicit none

    interface
        subroutine momentumStability(inod, GsPt_b, GsPt_d, TERM_RES_a, TERM_RES_b)
            use ElementVariables_mod, only: GaussPointQuantities
            implicit none
            PetscInt, intent(in) :: inod
            type(GaussPointQuantities), intent(in) :: GsPt_b, GsPt_d
            PetscScalar, dimension(:), intent(inout) :: TERM_RES_a, TERM_RES_b
        end subroutine momentumStability

        subroutine stressStability(inod, GsPt_b, GsPt_d, TERM_RES_a, TERM_RES_b)
            use ElementVariables_mod, only: GaussPointQuantities
            implicit none
            PetscInt, intent(in) :: inod
            type(GaussPointQuantities), intent(in) :: GsPt_b, GsPt_d
            PetscScalar, dimension(:), intent(inout) :: TERM_RES_a, TERM_RES_b
        end subroutine stressStability
    end interface

    procedure(momentumStability), pointer :: momentumStability_p
    procedure(stressStability), pointer :: stressStability_p

    contains

    subroutine momentumStabilityT(inod, GsPt_b, GsPt_d, TERM_RES_a, TERM_RES_b)

        use PhysicalParameters_mod, only: ReN
        use Tools_mod, only: doubleDotProduct
        use Residuals_mod, only: LSMEMomentum, LSICMomentum, LSCEMomentumT

        implicit none

        PetscInt, intent(in) :: inod
        type(GaussPointQuantities), intent(in) :: GsPt_b, GsPt_d
        PetscScalar, dimension(:), intent(inout) :: TERM_RES_a, TERM_RES_b

        PetscInt  :: i, Neq_f_S, Neq_f_M
        PetscReal :: bifn
        PetscReal, dimension(3) :: dFdXi
        PetscScalar :: lsme_term, lsic_term, lsce_term
        PetscScalar, dimension(3,3) :: lsce

        Neq_f_S = Problem_Main_Stability%Neq_f

        bifn = Elem%bfn_p(inod)

        dFdXi(:) = Elem%dFdXi(:,inod)

        do i = 1, Neq_f_S

            lsme_term = LSMEMomentum(dFdXi,GsPt_b)

            lsic_term = LSICMomentum(dFdXi(i),GsPt_b%tlsic)
            lsic_term = lsic_term*(GsPt_d%Continuity)

            lsce(:,:) = LSCEMomentumT(i,dFdXi,GsPt_b%tlsce)
            lsce_term = doubleDotProduct(lsce,GsPt_d%Con_EQ)

            TERM_RES_a(i) = ReN*(GsPt_d%Convection(i))*bifn &
                          + dot_product( GsPt_d%P_tot(i,:),dFdXi ) &
                          !Stabilization terms
                          + lsme_term*(GsPt_d%Momentum(i)) + lsic_term + lsce_term

            lsce_term = doubleDotProduct(lsce,GsPt_d%Con_EQ_time)

            TERM_RES_b(i) = GsPt_d%Momentum_time(i)*bifn &
                          + lsme_term*(GsPt_d%Momentum_time(i)) + lsce_term

        end do

        Neq_f_M = Problem_Main%Neq_f
        if (Neq_f_S > Neq_f_M) then
            call momentumStability3D(inod, GsPt_d, TERM_RES_a)
        end if

    end subroutine momentumStabilityT

    !----------------------------------------------------------------------

    subroutine momentumStability3D(inod, GsPt_d, TERM_RES_a)

        use petscksp, only: PETSC_i
        use PhysicalParameters_mod, only: BvN, KwN
        use ContinuationParameters_mod, only: Stability_Im

        PetscInt, intent(in) :: inod
        type(GaussPointQuantities), intent(in) :: GsPt_d
        PetscScalar, dimension(:), intent(inout) :: TERM_RES_a

        PetscInt  :: i, N
        PetscReal :: bifn
        PetscScalar :: term

        N = Problem_Main_Stability%Ndim

        bifn = Elem%bfn_p(inod)

        do i = 1, N
            term = bifn*KwN
            if (Stability_Im) then
                term = term*PETSC_i
            else
                if (i == N) then
                    term = -term
                end if
            end if
            TERM_RES_a(i) = TERM_RES_a(i) - term*(GsPt_d%P_tot(i,N))
        end do                

        ! - bifn*PETSC_i*KwN*( GsPt_d%T_ve(N,i) + BvN*(PETSC_i*KwN)*GsPt_d%U_f(i) )

    end subroutine momentumStability3D

    !----------------------------------------------------------------------

    subroutine continuityStability(inod, GsPt_b, GsPt_d, TERM_RES_a, TERM_RES_b)

        use Residuals_mod, only: LSMEContinuity

        implicit none

        PetscInt, intent(in) :: inod
        type(GaussPointQuantities), intent(in) :: GsPt_b, GsPt_d
        PetscScalar, dimension(:), intent(inout) :: TERM_RES_a, TERM_RES_b

        PetscInt :: Neq_f_S
        PetscReal :: bifn
        PetscReal, dimension(3) :: dFdXi
        PetscScalar :: lsme_term

        Neq_f_S = Problem_Main_Stability%Neq_f

        bifn = Elem%bfn_p(inod)

        dFdXi(:) = Elem%dFdXi(:,inod)

        lsme_term = LSMEContinuity(dFdXi,GsPt_b%tlsme,GsPt_d%Momentum)
        TERM_RES_a(Neq_f_S+1) = (GsPt_d%Continuity)*bifn + lsme_term

        lsme_term = LSMEContinuity(dFdXi,GsPt_b%tlsme,GsPt_d%Momentum_time)
        TERM_RES_b(Neq_f_S+1) = lsme_term

    end subroutine continuityStability

    !----------------------------------------------------------------------

    subroutine stressStabilityT(inod, GsPt_b, GsPt_d, TERM_RES_a, TERM_RES_b)

        use Tools_mod, only: getStressComponent
        use PhysicalParameters_mod, only: DfN
        use Residuals_mod, only: LSMEStressT, LSCEStressT, DCSStressT

        implicit none

        PetscInt, intent(in) :: inod
        type(GaussPointQuantities), intent(in) :: GsPt_b, GsPt_d
        PetscScalar, dimension(:), intent(inout) :: TERM_RES_a, TERM_RES_b

        PetscInt  :: i, j, k, istart, Neq_s
        PetscReal :: bifn
        PetscReal, dimension(3) :: dFdXi
        PetscScalar :: lsme_term, lsce_term, dcs_term

        Neq_s = Problem_Main_Stability%Neq_s

        bifn = Elem%bfn_p(inod)

        dFdXi(:) = Elem%dFdXi(:,inod)

        istart = size(TERM_RES_a) - Neq_s + 1
        do i = istart, size(TERM_RES_a)

            call getStressComponent(i,istart,j,k)

            lsme_term = LSMEStressT(j,k,dFdXi,GsPt_b%tlsme,GsPt_d%Momentum)

            lsce_term = LSCEStressT(bifn,dFdXi,GsPt_b)

            dcs_term = DCSStressT(j,k,dFdXi,GsPt_b)

            TERM_RES_a(i) = (GsPt_d%Con_EQ(j,k))*bifn &
                            + DfN*dot_product(GsPt_d%dTdXi(j,k,:),dFdXi) &
                            !Stabilization terms
                            + lsme_term + lsce_term*(GsPt_d%Con_EQ(j,k)) &
                            + dcs_term*(GsPt_d%mag_con_EQ)
            
            lsme_term = LSMEStressT(j,k,dFdXi,GsPt_b%tlsme,GsPt_d%Momentum_time)

            TERM_RES_b(i) = (GsPt_d%Con_EQ_time(j,k))*bifn &
                            !Stabilization terms
                            + lsme_term + lsce_term*(GsPt_d%Con_EQ_time(j,k)) &
                            + dcs_term*(GsPt_d%mag_con_EQ_time)

        end do

    end subroutine stressStabilityT

    !----------------------------------------------------------------------

    subroutine momentumStabilitySQRT(inod, GsPt_b, GsPt_d, TERM_RES_a, TERM_RES_b)

        use PhysicalParameters_mod, only: ReN
        use Tools_mod, only: doubleDotProduct
        use Residuals_mod, only: LSMEMomentum, LSICMomentum, LSCEMomentumSQRT

        implicit none

        PetscInt, intent(in) :: inod
        type(GaussPointQuantities), intent(in) :: GsPt_b, GsPt_d
        PetscScalar, dimension(:), intent(inout) :: TERM_RES_a, TERM_RES_b

        PetscInt  :: i, Neq_f_S, Neq_f_M
        PetscReal :: bifn
        PetscReal, dimension(3) :: dFdXi
        PetscScalar :: lsme_term, lsic_term, lsce_term
        PetscScalar, dimension(3,3) :: lsce

        Neq_f_S = Problem_Main_Stability%Neq_f
        Neq_f_M = Problem_Main%Neq_f

        bifn = Elem%bfn_p(inod)

        dFdXi(:) = Elem%dFdXi(:,inod)

        do i = 1, Neq_f_S

            lsme_term = LSMEMomentum(dFdXi,GsPt_b)

            lsic_term = LSICMomentum(dFdXi(i),GsPt_b%tlsic)
            lsic_term = lsic_term*(GsPt_d%Continuity)

            lsce(:,:) = LSCEMomentumSQRT(i,dFdXi,GsPt_b)
            lsce_term = doubleDotProduct(lsce,GsPt_d%Con_EQ)

            TERM_RES_a(i) = ReN*(GsPt_d%Convection(i))*bifn &
                          + dot_product( GsPt_d%P_tot(i,:),dFdXi ) &
                          !Stabilization terms
                          + lsme_term*(GsPt_d%Momentum(i)) + lsic_term + lsce_term

            lsce_term = doubleDotProduct(lsce,GsPt_d%Con_EQ_time)

            TERM_RES_b(i) = GsPt_d%Momentum_time(i)*bifn &
                          + lsme_term*(GsPt_d%Momentum_time(i)) + lsce_term

        end do

        if (Neq_f_S > Neq_f_M) then
            call momentumStability3D(inod, GsPt_d, TERM_RES_a)
        end if

    end subroutine momentumStabilitySQRT

    !----------------------------------------------------------------------

    subroutine stressStabilitySQRT(inod, GsPt_b, GsPt_d, TERM_RES_a, TERM_RES_b)

        use PhysicalParameters_mod, only: DfN
        use Tools_mod, only: getStressComponent
        use Residuals_mod, only: LSMEStressSQRT, LSCEStressSQRT, DCSStressSQRT

        implicit none

        PetscInt, intent(in) :: inod
        type(GaussPointQuantities), intent(in) :: GsPt_b, GsPt_d
        PetscScalar, dimension(:), intent(inout) :: TERM_RES_a, TERM_RES_b

        PetscInt  :: i, j, k, istart, Neq_s
        PetscReal :: bifn, jbfn
        PetscReal, dimension(3) :: dFdXi
        PetscScalar :: lsme_term, lsce_term, dcs_term

        Neq_s = Problem_Main_Stability%Neq_s

        bifn = Elem%bfn_p(inod)

        dFdXi(:) = Elem%dFdXi(:,inod)

        istart = size(TERM_RES_a) - Neq_s + 1
        do i = istart, size(TERM_RES_a)

            call getStressComponent(i,istart,j,k)

            lsme_term = LSMEStressSQRT(j,k,dFdXi,GsPt_b%tlsme,GsPt_d%Momentum)

            lsce_term = LSCEStressSQRT(bifn,dFdXi,GsPt_b)

            dcs_term = DCSStressSQRT(j,k,dFdXi,GsPt_b)

            TERM_RES_a(i) = (GsPt_d%Con_EQ(j,k))*bifn &
                            + DfN*dot_product(GsPt_d%dSdXi(j,k,:),dFdXi) &
                            !Stabilization terms
                            + lsme_term + lsce_term*(GsPt_d%Con_EQ(j,k)) &
                            + dcs_term*(GsPt_d%mag_con_EQ)
            
            lsme_term = LSMEStressSQRT(j,k,dFdXi,GsPt_b%tlsme,GsPt_d%Momentum_time)

            TERM_RES_b(i) = (GsPt_d%Con_EQ_time(j,k))*bifn &
                            !Stabilization terms
                            + lsme_term + lsce_term*(GsPt_d%Con_EQ_time(j,k)) &
                            + dcs_term*(GsPt_d%mag_con_EQ_time)

        end do

    end subroutine stressStabilitySQRT

    !----------------------------------------------------------------------

    subroutine momentumOpenBCBoundaryStability(inod, GsPt_b, GsPt_d, TERM_RES_a, TERM_RES_b)

        implicit none

        PetscInt, intent(in) :: inod
        type(GaussPointQuantities), intent(in) :: GsPt_b, GsPt_d
        PetscScalar, dimension(:), intent(inout) :: TERM_RES_a, TERM_RES_b

        PetscReal, dimension(3,3), parameter :: &
            I_tensor = reshape([1,0,0,0,1,0,0,0,1], shape = shape(I_tensor))
        PetscReal :: bifn
        PetscInt  :: i

        bifn = Elem%bfn_p(inod)

        ! TERM_RES_a(1) = (GU_d(1,1) - 0.0d0)*bifn
        ! TERM_RES_a(2) = (U_f_d(2) - 0.0d0)*bifn

        !Open BC
        do i = 1, Problem_Main%Ndim
            TERM_RES_a(i) = -dot_product(GsPt_b%normal,GsPt_d%P_tot(:,i) &
                                        + GsPt_d%P*I_tensor(:,i))*bifn
        end do

        TERM_RES_b(:) = 0.0d0

    end subroutine momentumOpenBCBoundaryStability

    !----------------------------------------------------------------------

    subroutine momentumSlipBoundaryStability(inod, GsPt_b, GsPt_d, TERM_RES_a, TERM_RES_b)

        use PhysicalParameters_mod, only: bslip

        implicit none

        PetscInt, intent(in) :: inod
        type(GaussPointQuantities), intent(in) :: GsPt_b, GsPt_d
        PetscScalar, dimension(:), intent(inout) :: TERM_RES_a, TERM_RES_b

        PetscReal :: bifn, kslip
        PetscScalar :: slip_value, res_value, n_dot_u, ls1, ls2
        PetscScalar, dimension(3,3) :: T_tensor
        PetscReal, dimension(3,3), parameter :: &
            I_tensor = reshape([1,0,0,0,1,0,0,0,1], shape = shape(I_tensor))

        bifn = Elem%bfn_p(inod)

        kslip = 1.0d0/bslip

        T_tensor(:,:) = GsPt_d%P_tot(:,:) !+ P_d*I_tensor(:,:)

        slip_value = -kslip*( dot_product( matmul( GsPt_b%tangent,T_tensor ),GsPt_b%normal ) )
        res_value = dot_product( GsPt_b%tangent,GsPt_d%U_f ) - slip_value
        n_dot_u = dot_product( GsPt_b%normal,GsPt_d%U_f )

        ls1 = n_dot_u*(bifn*GsPt_b%normal(1))
        ls1 = ls1 + res_value*(bifn*GsPt_b%tangent(1))! - slip_value)
        TERM_RES_a(1) = bifn*n_dot_u + ls1

        !Slip in y
        ls2 = n_dot_u*(bifn*GsPt_b%normal(2))
        ls2 = ls2 + res_value*(bifn*GsPt_b%tangent(2))
        TERM_RES_a(2) = bifn*res_value + ls2

        TERM_RES_b(:) = 0.0d0

    end subroutine momentumSlipBoundaryStability

end module ResidualsStability_mod
