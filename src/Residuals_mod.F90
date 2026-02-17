
#include <petsc/finclude/petscksp.h>

module Residuals_mod

    use petscksp

    use PhysicalParameters_mod, only: ProblemParameters
    use FEMParameters_mod, only: FEParameters
    use MeshParameters_mod, only: MeshParameters
    use BoundaryParameters_mod, only: BoundaryParameters
    use GaussParameters_mod, only: GaussIntegration
    use SolutionVariables_mod, only: SolutionArraysType
    use LinearSystemVariables_mod, only: LinearSystemType
    use ElementVariables_mod, only: GaussPointQuantities, NodeArrays

    implicit none

    interface
        subroutine momentum(Neq_f, inod, Elem, GsPt, TERM_RES)
            use ElementVariables_mod, only: GaussPointQuantities, NodeArrays
            implicit none
            PetscInt, intent(in) :: Neq_f, inod
            type(NodeArrays), intent(in) :: Elem
            type(GaussPointQuantities), intent(in) :: GsPt
            PetscReal, dimension(:), intent(inout) :: TERM_RES
        end subroutine momentum

        subroutine stress(Neq_s, inod, Elem, GsPt, TERM_RES)
            use ElementVariables_mod, only: GaussPointQuantities, NodeArrays
            implicit none
            PetscInt, intent(in) :: Neq_s, inod
            type(NodeArrays), intent(in) :: Elem
            type(GaussPointQuantities), intent(in) :: GsPt
            PetscReal, dimension(:), intent(inout) :: TERM_RES
        end subroutine stress
    end interface

    procedure(momentum), pointer :: momentumMain_p
    procedure(stress), pointer :: stressMain_p

    contains

    subroutine momentumMain(Neq_f, inod, Elem, GsPt, TERM_RES)

        use PhysicalParameters_mod, only: ReN
        use Tools_mod, only: doubleDotProduct

        implicit none

        PetscInt, intent(in) :: Neq_f, inod
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscReal, dimension(:), intent(inout) :: TERM_RES

        PetscInt :: i
        PetscReal :: bifn, jbfn
        PetscReal, dimension(3) :: dFdXi
        PetscScalar :: lsme_term, lsic_term, lsce_term
        PetscScalar, dimension(3,3) :: lsce

        bifn = Elem%bfn_p(inod)
        jbfn = Elem%bfn_p(inod)!/dt

        dFdXi(:) = Elem%dFdXi(:,inod)

        do i = 1, Neq_f

            lsme_term = LSMEMomentum(dFdXi,GsPt)
            lsme_term = lsme_term*(GsPt%Momentum(i))

            lsic_term = LSICMomentum(dFdXi(i),GsPt%tlsic)
            lsic_term = lsic_term*(GsPt%Continuity)

            lsce(:,:) = LSCEMomentumT(i,dFdXi,GsPt%tlsce)
            lsce_term = doubleDotProduct(lsce,GsPt%Con_EQ)

            TERM_RES(i) = ReN*(GsPt%Convection(i))*bifn &
                        ! + ReN*Jump(i)*jbfn &
                        + dot_product(GsPt%P_tot(i,:),dFdXi) &
                        !Stabilization terms
                        + lsme_term + lsic_term + lsce_term

        end do

    end subroutine momentumMain

    !----------------------------------------------------------------------

    subroutine continuityBulkMain(Neq_f, inod, Elem, GsPt, TERM_RES)

        implicit none

        PetscInt, intent(in) :: Neq_f, inod
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscReal, dimension(:), intent(inout) :: TERM_RES

        PetscReal :: bifn
        PetscReal, dimension(3) :: dFdXi
        PetscScalar :: lsme_term

        bifn = Elem%bfn_p(inod)

        dFdXi(:) = Elem%dFdXi(:,inod)

        lsme_term = LSMEContinuity(dFdXi,GsPt%tlsme,GsPt%Momentum)
        TERM_RES(Neq_f+1) = (GsPt%Continuity)*bifn + lsme_term

    end subroutine continuityBulkMain

    !----------------------------------------------------------------------

    subroutine stressesBulkMain(Neq_s, inod, Elem, GsPt, TERM_RES)

        use PhysicalParameters_mod, only: DfN
        use Tools_mod, only: getStressComponent

        implicit none

        PetscInt, intent(in) :: Neq_s, inod
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscReal, dimension(:), intent(inout) :: TERM_RES

        PetscInt :: i, j, k, istart
        PetscReal :: bifn, jbfn
        PetscReal, dimension(3) :: dFdXi
        PetscScalar :: lsme_term, lsce_term, dcs_term

        bifn = Elem%bfn_p(inod)
        jbfn = Elem%bfn_p(inod)!/dt

        dFdXi(:) = Elem%dFdXi(:,inod)

        istart = size(TERM_RES) - Neq_s + 1
        do i = istart, size(TERM_RES)

            call getStressComponent(i,istart,j,k)

            lsme_term = LSMEStressT(j,k,dFdXi,GsPt%tlsme,GsPt%Momentum)

            lsce_term = LSCEStressT(bifn,dFdXi,GsPt)
            lsce_term = lsce_term*GsPt%Con_EQ(j,k)

            dcs_term = DCSStressT(j,k,dFdXi,GsPt)
            dcs_term = dcs_term*(GsPt%mag_con_EQ)
            
            TERM_RES(i) = (GsPt%Con_EQ(j,k))*bifn &
                        + DfN*dot_product(GsPt%dTdXi(j,k,:),dFdXi) &
                        ! + Jump(i)*jbfn &
                        !Stabilization terms
                        + lsme_term + lsce_term + dcs_term

        end do

    end subroutine stressesBulkMain

    !----------------------------------------------------------------------

    function LSMEMomentum(dFdXi,GsPt) result(lsme_term)

        use PhysicalParameters_mod, only: ReN

        implicit none

        PetscReal, dimension(:), intent(in) :: dFdXi
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscScalar :: lsme_term

        lsme_term = (GsPt%tlsme)*ReN*dot_product(GsPt%U_f,dFdXi)

    end function LSMEMomentum

    !----------------------------------------------------------------------

    function LSICMomentum(dFdXi,tlsic) result(lsic_term)

        implicit none

        PetscReal, intent(in) :: dFdXi
        PetscScalar, intent(in) :: tlsic
        PetscScalar :: lsic_term

        lsic_term = tlsic*dFdXi

    end function LSICMomentum

    !----------------------------------------------------------------------

    function LSCEMomentumT(i,dFdXi,tlsce) result(lsce)

        use PhysicalParameters_mod, only: Model_Name

        implicit none

        PetscInt, intent(in) :: i
        PetscReal, dimension(:), intent(in) :: dFdXi
        PetscScalar, intent(in) :: tlsce
        PetscScalar, dimension(3,3) :: lsce

        PetscReal, dimension(3,3) :: GUT_lsce, GU_lsce, UCD_lsce

        GUT_lsce(:,:) = 0.0d0 ; GUT_lsce(i,:) = dFdXi(:)

        lsce(:,:) = -tlsce*GUT_lsce(:,:)
        
        ! GU_lsce(:,:) = transpose(GUT_lsce(:,:))
        ! UCD_lsce(:,:) = GU_lsce(:,:)+GUT_lsce(:,:)
        ! lsce (:,:) = -tlsce*UCD_lsce(:,:)

        ! select case(Model_Name)
        ! case ('FENE-CR')!, 'm-L-PTT', 'm-e-PTT')
        !     lsce(:,:) = (GsPt%f_T)*lsce(:,:)
        ! case default
        !     continue
        ! end select

    end function LSCEMomentumT

    !----------------------------------------------------------------------

    function LSMEContinuity(dFdXi,tlsme,Momentum) result(lsme_term)

        implicit none

        PetscReal, dimension(:), intent(in) :: dFdXi
        PetscScalar, intent(in) :: tlsme
        PetscScalar, dimension(3), intent(in) :: Momentum
        PetscScalar :: lsme_term

        lsme_term = tlsme*dot_product(dFdXi,Momentum)

    end function LSMEContinuity

    !----------------------------------------------------------------------

    function LSMEStressT(j,k,dFdXi,tlsme,Momentum) result(lsme_term)

        use PhysicalParameters_mod, only: BvN

        implicit none

        PetscInt, intent(in) :: j, k
        PetscReal, dimension(:), intent(in) :: dFdXi
        PetscScalar, intent(in) :: tlsme
        PetscScalar, dimension(:), intent(in) :: Momentum
        PetscScalar :: lsme_term

        lsme_term = -dFdXi(k)*Momentum(j)
        if (j /= k) then
            lsme_term = lsme_term -dFdXi(j)*Momentum(k)
        end if
        lsme_term = tlsme*lsme_term*(1.0d0-BvN)

    end function LSMEStressT

    !----------------------------------------------------------------------

    function LSCEStressT(bifn,dFdXi,GsPt) result(lsce_term)

        use PhysicalParameters_mod, only: WiN

        implicit none

        PetscReal, intent(in) :: bifn
        PetscReal, dimension(:), intent(in) :: dFdXi
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscScalar :: lsce_term

        lsce_term = WiN*dot_product(GsPt%U_f,dFdXi)
        ! ! lsce_term = lsce_term+(GsPt%f(2))*bifn
        ! lsce_term = lsce_term+(GsPt%f_T)*bifn
        lsce_term = (GsPt%tlsce)*lsce_term

    end function LSCEStressT

    !----------------------------------------------------------------------

    function DCSStressT(j,k,dFdXi,GsPt) result(dcs_term)

        implicit none

        PetscInt, intent(in) :: j, k
        PetscReal, dimension(:), intent(in) :: dFdXi
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscScalar :: dcs_term

        dcs_term = (GsPt%tdcs)*dot_product(GsPt%dTdXi(j,k,:),dFdXi)

    end function DCSStressT

    !----------------------------------------------------------------------

    subroutine momentumMainSQRT(Neq_f, inod, Elem, GsPt, TERM_RES)

        use PhysicalParameters_mod, only: ReN
        use Tools_mod, only: doubleDotProduct

        implicit none

        PetscInt, intent(in) :: Neq_f, inod
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscReal, dimension(:), intent(inout) :: TERM_RES

        PetscInt :: i
        PetscReal :: bifn, jbfn
        PetscReal, dimension(3) :: dFdXi
        PetscScalar :: lsme_term, lsic_term, lsce_term
        PetscScalar, dimension(3,3) :: lsce

        bifn = Elem%bfn_p(inod)
        jbfn = Elem%bfn_p(inod)!/dt

        dFdXi(:) = Elem%dFdXi(:,inod)

        do i = 1, Neq_f

            lsme_term = LSMEMomentum(dFdXi,GsPt)
            lsme_term = lsme_term*(GsPt%Momentum(i))

            lsic_term = LSICMomentum(dFdXi(i),GsPt%tlsic)
            lsic_term = lsic_term*(GsPt%Continuity)
    
            lsce(:,:) = LSCEMomentumSQRT(i,dFdXi,GsPt)
            lsce_term = doubleDotProduct(lsce,GsPt%Con_EQ)
            
            TERM_RES(i) = ReN*(GsPt%Convection(i))*bifn &
                        ! + ReN*Jump(i)*jbfn &
                        + dot_product(GsPt%P_tot(i,:),dFdXi) &
                        !Stabilization terms
                        + lsme_term + lsic_term + lsce_term

        end do

    end subroutine momentumMainSQRT

    !----------------------------------------------------------------------

    subroutine stressesBulkMainSQRT(Neq_s, inod, Elem, GsPt, TERM_RES)

        use PhysicalParameters_mod, only: DfN
        use Tools_mod, only: getStressComponent

        implicit none

        PetscInt, intent(in) :: Neq_s, inod
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscReal, dimension(:), intent(inout) :: TERM_RES

        PetscInt :: i, j, k, istart
        PetscReal :: bifn, jbfn
        PetscReal, dimension(3) :: dFdXi
        PetscScalar :: lsme_term, lsce_term, dcs_term

        bifn = Elem%bfn_p(inod)
        jbfn = Elem%bfn_p(inod)!/dt

        dFdXi(:) = Elem%dFdXi(:,inod)

        istart = size(TERM_RES) - Neq_s + 1
        do i = istart, size(TERM_RES)

            call getStressComponent(i,istart,j,k)

            lsme_term = LSMEStressSQRT(j,k,dFdXi,GsPt%tlsme,GsPt%Momentum)

            lsce_term = LSCEStressSQRT(bifn,dFdXi,GsPt)
            lsce_term = lsce_term*(GsPt%Con_EQ(j,k))

            dcs_term = DCSStressSQRT(j,k,dFdXi,GsPt)
            dcs_term = dcs_term*(GsPt%mag_con_EQ)

            TERM_RES(i) = (GsPt%Con_EQ(j,k))*bifn &
                        + DfN*dot_product(GsPt%dSdXi(j,k,:),dFdXi) &
                        ! + Jump(i)*jbfn &
                        !Stabilization terms
                        + lsme_term + lsce_term + dcs_term

        end do

    end subroutine stressesBulkMainSQRT

    !----------------------------------------------------------------------

    function LSCEMomentumSQRT(i,dFdXi,GsPt) result(lsce)

        use PhysicalParameters_mod, only: BvN, WiN, PM => Problem_Main
        use ConstitutiveModelsSQRT_mod, only: WTensor_p!, WTensor2D, WTensor3D

        implicit none

        PetscInt, intent(in) :: i
        PetscReal, dimension(:), intent(in) :: dFdXi
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscScalar, dimension(3,3) :: lsce

        PetscInt :: j, N
        PetscReal, dimension(3,3) :: GUT_lsce, GU_lsce, W_lsce, UCD_lsce
        PetscReal, dimension(PM%Ndim,PM%Ndim) :: TEMP_S, TEMP_GUT, TEMP_W

        N = PM%Ndim

        GUT_lsce(:,:) = 0.0d0 ; GUT_lsce(i,:) = dFdXi(:)
        GU_lsce(:,:) = transpose(GUT_lsce)

        TEMP_S(:,:) = GsPt%S(1:N,1:N)
        TEMP_GUT(:,:) = GUT_lsce(1:N,1:N)

        W_lsce(:,:) = 0.0d0
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
        !     write(*,'(a)') 'Wrong N in LSCEMomentumSQRT!'
        !     stop
        ! end select

        W_lsce(1:N,1:N) = TEMP_W(:,:)

        UCD_lsce(:,:) = -matmul(GsPt%S,GU_lsce)-matmul(W_lsce,GsPt%S)

        lsce = (GsPt%tlsce)*UCD_lsce(:,:)*(1.0d0-BvN)/WiN

    end function LSCEMomentumSQRT

    !----------------------------------------------------------------------

    function LSMEStressSQRT(j,k,dFdXi,tlsme,Momentum) result(lsme_term)

        use PhysicalParameters_mod, only: WiN

        implicit none

        PetscInt, intent(in) :: j, k
        PetscReal, dimension(:), intent(in) :: dFdXi
        PetscScalar, intent(in) :: tlsme
        PetscScalar, dimension(3), intent(in) :: Momentum
        PetscScalar :: lsme_term

        lsme_term = -dFdXi(k)*Momentum(j)
        if (j /= k) then
            lsme_term = lsme_term - dFdXi(j)*Momentum(k)
        end if
        lsme_term = (tlsme)*lsme_term*WiN

    end function LSMEStressSQRT

    !----------------------------------------------------------------------

    function LSCEStressSQRT(bifn,dFdXi,GsPt) result(lsce_term)

        use PhysicalParameters_mod, only: WiN

        implicit none

        PetscReal, intent(in) :: bifn
        PetscReal, dimension(:), intent(in) :: dFdXi
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscScalar :: lsce_term

        lsce_term = WiN*dot_product(GsPt%U_f,dFdXi)
        ! lsce_term = lsce_term + (GsPt%f(2))*bifn/2.0d0
        lsce_term = lsce_term + (GsPt%f_T)*bifn/2.0d0
        lsce_term = (GsPt%tlsce)*lsce_term

    end function LSCEStressSQRT

    !----------------------------------------------------------------------

    function DCSStressSQRT(j,k,dFdXi,GsPt) result(dcs_term)

        implicit none

        PetscInt, intent(in) :: j, k
        PetscReal, dimension(:), intent(in) :: dFdXi
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscScalar :: dcs_term

        dcs_term = (GsPt%tdcs)*dot_product(GsPt%dSdXi(j,k,:),dFdXi)

    end function DCSStressSQRT

    !----------------------------------------------------------------------

    subroutine momentumInflow(Neq_f, inod, Elem, GsPt, TERM_RES)

        use PhysicalParameters_mod, only: ReN
        use Tools_mod, only: extraEQ, doubleDotProduct

        implicit none

        PetscInt, intent(in) :: Neq_f, inod
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscReal, dimension(:), intent(inout) :: TERM_RES

        PetscInt :: i
        PetscReal :: bifn, dPL
        PetscScalar :: lsme_term, lsic_term, lsce_term
        PetscScalar, dimension(3,3) :: lsce
        PetscReal, dimension(3) :: dFdXi

        bifn = Elem%bfn_p(inod)
        dFdXi(:) = Elem%dFdXi(:,inod)

        dPL = Elem%TEMP_EU(extraEQ('dPL'))

        ! lsme_term = LSMEMomentum(dFdXi,GsPt)
        ! lsme_term = lsme_term*(GsPt%Momentum(1))

        ! lsic_term = LSICMomentum(dFdXi(1),GsPt%tlsic)
        ! lsic_term = lsic_term*(GsPt%Continuity)

        ! ! !T
        ! ! lsce(:,:) = LSCEMomentumT(1,dFdXi,GsPt%tlsce)
        ! ! lsce_term = doubleDotProduct(lsce,GsPt%Con_EQ)
        ! !S
        ! lsce(:,:) = LSCEMomentumSQRT(1,dFdXi,GsPt)
        ! lsce_term = doubleDotProduct(lsce,GsPt%Con_EQ)

        TERM_RES(1) = ReN*(GsPt%dUfdt(1))*bifn + dPL*bifn + &
                    + (GsPt%P_tot(1,2))*dFdXi(2) &
                    + (GsPt%P_tot(1,3))*dFdXi(3) &
                    !Stabilization terms
                    + 0.0d0!lsme_term + lsic_term + lsce_term

        do i = 2, Neq_f
            TERM_RES(i) = (GsPt%U_f(i) - 0.0d0)*bifn
        end do

    end subroutine momentumInflow

    !----------------------------------------------------------------------

    subroutine stressesInflow(Neq_s, inod, Elem, GsPt, TERM_RES)

        use Tools_mod, only: getStressComponent
        use PhysicalParameters_mod, only: DfN
        
        implicit none

        PetscInt, intent(in) :: Neq_s, inod
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscReal, dimension(:), intent(inout) :: TERM_RES

        PetscInt  :: i, j, k, istart
        PetscReal :: bifn
        PetscReal, dimension(3) :: dFdXi
        PetscScalar :: lsme_term, lsce_term, dcs_term

        bifn = Elem%bfn_p(inod)

        dFdXi(:) = Elem%dFdXi(:,inod)

        istart = size(TERM_RES) - Neq_s + 1
        do i = istart, size(TERM_RES)

            call getStressComponent(i,istart,j,k)

            ! !T
            ! lsme_term = LSMEStressT(j,k,dFdXi,GsPt%tlsme,GsPt%Momentum)

            ! lsce_term = LSCEStressT(bifn,dFdXi,GsPt)
            ! lsce_term = lsce_term*GsPt%Con_EQ(j,k)

            ! dcs_term = DCSStressT(j,k,dFdXi,GsPt)
            ! dcs_term = dcs_term*(GsPt%mag_con_EQ)

            ! !S
            ! lsme_term = LSMEStressSQRT(j,k,dFdXi,GsPt%tlsme,GsPt%Momentum)

            ! lsce_term = LSCEStressSQRT(bifn,dFdXi,GsPt)
            ! lsce_term = lsce_term*(GsPt%Con_EQ(j,k))

            ! dcs_term = DCSStressSQRT(j,k,dFdXi,GsPt)
            ! dcs_term = dcs_term*(GsPt%mag_con_EQ)

            TERM_RES(i) = GsPt%Con_EQ(j,k)*bifn &
                        + DfN*dot_product(GsPt%dTdXi(j,k,:),dFdXi) &
                        !Stabilization terms
                        + 0.0d0!lsme_term + lsce_term + dcs_term

        end do

    end subroutine stressesInflow

    !----------------------------------------------------------------------

    subroutine openBC(Neq_f, inod, Elem, GsPt, TERM_RES)

        implicit none

        PetscInt, intent(in) :: Neq_f, inod
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscReal, dimension(:), intent(inout) :: TERM_RES

        PetscReal, dimension(3,3), parameter :: &
            I_tensor = reshape( [1,0,0,0,1,0,0,0,1], shape = shape(I_tensor) )
        PetscInt :: i
        PetscReal :: bifn

        bifn = Elem%bfn_p(inod)

        ! if ( any(GsPt%normal < -1.0d-6) ) write(*,'(3es12.4)') GsPt%normal

        ! TERM_RES(:) = 0.0d0
        do i = 1, Neq_f
            TERM_RES(i) = -dot_product( GsPt%normal,GsPt%P_tot(:,i)&
                                        +GsPt%P*I_tensor(:,i) )*bifn
            ! TERM_RES(i) = -dot_product(GsPt%normal,GsPt%P_tot(:,i))*bifn
        end do

    end subroutine openBC

    !----------------------------------------------------------------------

    subroutine slipBC(inod, Elem, GsPt, TERM_RES)

        use PhysicalParameters_mod, only: bslip

        implicit none

        PetscInt, intent(in) :: inod
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscReal, dimension(:), intent(inout) :: TERM_RES

        PetscReal :: bifn, res_value, ls1, ls2, tlsme_loc
        PetscReal :: kslip, theta, Ur, Uth, slip_value, n_dot_u
        PetscReal, dimension(3,3) :: T_tensor
        PetscReal, dimension(3,3), parameter :: &
            I_tensor = reshape( [1,0,0,0,1,0,0,0,1], shape = shape(I_tensor) )


        bifn = Elem%bfn_p(inod)

        kslip = 1.0d0/bslip

        T_tensor(:,:) = GsPt%P_tot !+ P*I_tensor(:,:)

        slip_value = -kslip*(dot_product(matmul(GsPt%tangent,T_tensor),GsPt%normal))
        res_value = dot_product(GsPt%tangent,GsPt%U_f) - slip_value
        n_dot_u = dot_product(GsPt%normal,GsPt%U_f)

        tlsme_loc = 1.0d0!tlsme

        !No penetration in x
        ls1 = n_dot_u*(bifn*GsPt%normal(1))
        ls1 = ls1 + res_value*(bifn*GsPt%tangent(1))! - slip_value)
        TERM_RES(1) = bifn*n_dot_u + tlsme_loc*ls1

        ! !Slip in x
        ! TERM_RES(1) = kslip*dot_product(tangent,U_f) - tangent(2)*dot_product(normal,P_tot(:,2))
        ! TERM_RES(1) = bifn*TERM_RES(1)/tangent(1)

        !Slip in y
        ls2 = n_dot_u*(bifn*GsPt%normal(2))
        ls2 = ls2 + res_value*(bifn*GsPt%tangent(2))
        TERM_RES(2) = bifn*res_value + tlsme_loc*ls2
        ! TERM_RES(2) = kslip*dot_product(tangent,U_f) - tangent(1)*dot_product(normal,P_tot(:,1))
        ! TERM_RES(2) = bifn*TERM_RES(2)/tangent(2)

        ! theta = atan2(GsPt%Xi_loc(inod,2),GsPt%Xi_loc(inod,1))
        ! !Ur = 0
        ! Ur    =  U_f(1)*cos(theta) + U_f(2)*sin(theta)
        ! TERM_RES(1) = 1.0d8*(Ur - 0.0d0)*bifn

        ! !No slip
        ! TERM_RES(2) = bifn*(dot_product(tangent,U_f))
        ! !Uθ = 0
        ! Uth = -U_f(1)*sin(theta) + U_f(2)*cos(theta)
        ! TERM_RES(2) = 1.0d8*(Uth - 0.0d0)*bifn

        ! !Tangent Momentum Equations
        ! TERM_RES(1) = bifn*dot_product(tangent,Momentum)

    end subroutine slipBC

end module Residuals_mod
