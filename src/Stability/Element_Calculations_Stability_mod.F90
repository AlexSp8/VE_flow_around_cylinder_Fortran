
#include <slepc/finclude/slepceps.h>

module ElementCalculationsStability_mod

    use PhysicalParameters_mod, only: ProblemParameters
    use FEMParameters_mod, only: FEParameters
    use ElementVariables_mod, only: GaussPointQuantities, NodeArrays

    implicit none

    interface setFlowUnknownsStability
        module procedure setFlowUnknownsStabilityReal
        module procedure setFlowUnknownsStabilityComplex
    end interface setFlowUnknownsStability

    contains

    subroutine setFlowQuantitiesStability(Problem, Problem_Stability, Nbf, Elem, GsPt_b, GsPt_d)

        use petscksp, only: PETSC_i
        use ElementCalculations_mod, only: strongFormEquations, setFlowKinematics, &
            strongFormContinuity, setStabilizationTerms
        use Tools_mod, only: doubleDotProduct
        use PhysicalParameters_mod, only: Stress_Reform, QbN, PI, Periodicity
        use ConstitutiveModelsStability_mod, only: strongFormConEQ_Stability
        use ConstitutiveModelsSQRTStability_mod, only: strongFormConEQ_Stability_SQRT

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        PetscInt, intent(in) :: Nbf
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(inout) :: GsPt_b, GsPt_d

        PetscReal :: KwN
        PetscInt :: Neq, Neq_S, N, i
        PetscBool :: check

        Neq = Problem%Neq
        Neq_S = Problem_Stability%Neq

        !Set base problem
        call setFlowUnknownsStability(Problem_Stability, Elem%TEMP_TL_b, Elem, GsPt_b)
        call strongFormEquations(Problem_Stability, Elem%TEMP_EU, GsPt_b)
        call setStabilizationTerms(Nbf, Elem%dFdXi, GsPt_b)

        !Set perturbed problem
        call setFlowUnknownsStability(Problem_Stability, Elem%TEMP_TL_d, Elem, GsPt_d)

        ! KwN = 2.0d0*PI*QbN
        ! do i = 1, size(Periodicity)

        !     select case (Periodicity(i))
        !     case ('X')
        !         N = 1
        !     case ('Y')
        !         N = 2
        !     case ('Z')
        !         N = 3
        !     case default
        !         continue
        !     end select

        !     ! GsPt_d%GU(N,:) = 0.0d0
        !     GsPt_d%GU(N,:) = GsPt_d%GU(N,:) + PETSC_i*KwN*(GsPt_d%U_f(:))           !dUi/dz

        !     ! GsPt_d%G2U(N,:) = 0.0d0
        !     GsPt_d%G2U(N,:) = GsPt_d%G2U(N,:) -(KwN**2)*GsPt_d%U_f(:)

        !     ! GsPt_d%dPdXi(N) = 0.0d0
        !     GsPt_d%dPdXi(N) = GsPt_d%dPdXi(N) + PETSC_i*KwN*(GsPt_d%P)              !dP/dz

        !     ! GsPt_d%dSdXi(:,:,N) = 0.0d0
        !     GsPt_d%dSdXi(:,:,N) = GsPt_d%dSdXi(:,:,N) + PETSC_i*KwN*(GsPt_d%S(:,:)) !dTij/dz

        ! end do

        ! if (Neq_S > Neq) then
        !     call setExtraFlowUnknownsStability(Problem_Stability%Ndim, GsPt_d)
        ! end if

        call setFlowKinematics(Problem_Stability%Ndim,GsPt_d)

        select case (Stress_Reform)
        case ('SQRT')
            call strongFormConEQ_Stability_SQRT(GsPt_b, GsPt_d)
        case default
            call strongFormConEQ_Stability(GsPt_b, GsPt_d)
        end select

        call strongFormMomentumStability(GsPt_b, GsPt_d)

        call strongFormContinuity(GsPt_d)

        GsPt_d%mag_con_EQ = doubleDotProduct(GsPt_b%Con_EQ,GsPt_d%Con_EQ)
        GsPt_d%mag_con_EQ = GsPt_d%mag_con_EQ/(2.0d0*GsPt_b%mag_con_EQ + 1.0d-12)

        GsPt_d%mag_con_EQ_time = doubleDotProduct(GsPt_b%Con_EQ_time,GsPt_d%Con_EQ_time)
        GsPt_d%mag_con_EQ_time = GsPt_d%mag_con_EQ_time/(2.0d0*GsPt_b%mag_con_EQ_time + 1.0d-12)

    end subroutine setFlowQuantitiesStability

    !---------------------------------------------------------------------

    subroutine setFlowUnknownsStabilityReal(Problem_Stability, TEMP_TL, Elem, GsPt)

        use Tools_mod, only: getStressComponent, getTensorComponent

        implicit none

        type(ProblemParameters), intent(in) :: Problem_Stability
        PetscReal, dimension(:,:), intent(in) :: TEMP_TL
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt :: i, j, k, idim, istart, Neq_f, Neq_s
        PetscBool :: check
        
        Neq_f = Problem_Stability%Neq_f

        GsPt%U_f(:) = 0.0d0
        GsPt%GU(:,:) = 0.0d0
        GsPt%G2U(:,:) = 0.0d0

        GsPt%P = 0.0d0
        GsPt%dPdXi(:) = 0.0d0
        
        do i = 1, Neq_f

            GsPt%U_f(i) = dot_product( TEMP_TL(:,i),Elem%bfn_p(:) )

            do j = 1, 3
                GsPt%GU(j,i)  = dot_product( TEMP_TL(:,i),Elem%dFdXi(j,:) )
                GsPt%G2U(j,i) = dot_product( TEMP_TL(:,i),Elem%d2FdXi2(j,:) )
            end do

            GsPt%dPdXi(i) = dot_product( TEMP_TL(:,Neq_f+1),Elem%dFdXi(i,:) )

        end do

        GsPt%P = dot_product( TEMP_TL(:,Neq_f+1),Elem%bfn_p(:) )

        GsPt%S(:,:) = 0.0d0
        GsPt%dSdXi(:,:,:) = 0.0d0

        Neq_s = Problem_Stability%Neq_s
        istart = size(TEMP_TL,2) - Neq_s + 1
        do i = istart, size(TEMP_TL,2)

            call getStressComponent(i,istart,j,k)

            GsPt%S(j,k) = dot_product( TEMP_TL(:,i),Elem%bfn_p(:) )
            GsPt%S(k,j) = GsPt%S(j,k)
            do idim = 1, 3
                GsPt%dSdXi(j,k,idim) = dot_product( TEMP_TL(:,i),Elem%dFdXi(idim,:) )
            end do
            GsPt%dSdXi(k,j,:) = GsPt%dSdXi(j,k,:)

        end do

        check = (Problem_Stability%name == 'Inflow_Stability_1D')
        check = check .or. (Problem_Stability%name == 'Inflow_Stability_1D_2D')
        check = check .or. (Problem_Stability%name == 'Inflow_Stability_2D')
        if (check) then
            ! GsPt%U_f(2)  = 0.0d0  !Uy
            ! GsPt%U_f(3)  = 0.0d0  !Uz
            GsPt%GU(1,:) = 0.0d0  !d/dx
            GsPt%GU(:,2) = 0.0d0  !duy
            GsPt%GU(:,3) = 0.0d0  !duz
            GsPt%P = 0.0d0
            GsPt%dPdXi(:) = 0.0d0
            GsPt%dSdXi(:,:,:) = 0.0d0  !dSjk/dxi
        end if

        ! GsPt%GU_proj(:,:) = 0.0d0
        ! do i = 1, size(TEMP_TL_proj,2)
        !     call getTensorComponent(i,j,k)
        !     GsPt%GU_proj(j,k) = dot_product(TEMP_TL_proj(:,i),Elem%bfn_p(:))
        ! end do

    end subroutine setFlowUnknownsStabilityReal

    !---------------------------------------------------------------------

    subroutine setFlowUnknownsStabilityComplex(Problem_Stability, TEMP_TL, Elem, GsPt)

        use Tools_mod, only: getStressComponent, getTensorComponent

        implicit none

        type(ProblemParameters), intent(in) :: Problem_Stability
        PetscScalar, dimension(:,:), intent(in) :: TEMP_TL
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt :: i, j, k, idim, istart, Neq_f, Neq_s
        PetscBool :: check
        
        Neq_f = Problem_Stability%Neq_f

        GsPt%U_f(:) = 0.0d0
        GsPt%GU(:,:) = 0.0d0
        GsPt%G2U(:,:) = 0.0d0

        GsPt%P = 0.0d0
        GsPt%dPdXi(:) = 0.0d0
        
        do i = 1, Neq_f

            GsPt%U_f(i) = dot_product( TEMP_TL(:,i),Elem%bfn_p(:) )

            do j = 1, 3
                GsPt%GU(j,i)  = dot_product( TEMP_TL(:,i),Elem%dFdXi(j,:) )
                GsPt%G2U(j,i) = dot_product( TEMP_TL(:,i),Elem%d2FdXi2(j,:) )
            end do

            GsPt%dPdXi(i) = dot_product( TEMP_TL(:,Neq_f+1),Elem%dFdXi(i,:) )

        end do

        GsPt%P = dot_product( TEMP_TL(:,Neq_f+1),Elem%bfn_p(:) )

        GsPt%S(:,:) = 0.0d0
        GsPt%dSdXi(:,:,:) = 0.0d0

        Neq_s = Problem_Stability%Neq_s
        istart = size(TEMP_TL,2) - Neq_s + 1
        do i = istart, size(TEMP_TL,2)

            call getStressComponent(i,istart,j,k)

            GsPt%S(j,k) = dot_product( TEMP_TL(:,i),Elem%bfn_p(:) )
            GsPt%S(k,j) = GsPt%S(j,k)
            do idim = 1, 3
                GsPt%dSdXi(j,k,idim) = dot_product( TEMP_TL(:,i),Elem%dFdXi(idim,:) )
            end do
            GsPt%dSdXi(k,j,:) = GsPt%dSdXi(j,k,:)

        end do

        check = (Problem_Stability%name == 'Inflow_Stability_1D')
        check = check .or. (Problem_Stability%name == 'Inflow_Stability_1D_2D')
        check = check .or. (Problem_Stability%name == 'Inflow_Stability_2D')
        if (check) then
            ! GsPt%U_f(2)  = 0.0d0  !Uy
            ! GsPt%U_f(3)  = 0.0d0  !Uz
            GsPt%GU(1,:) = 0.0d0  !d/dx
            GsPt%GU(:,2) = 0.0d0  !duy
            GsPt%GU(:,3) = 0.0d0  !duz
            GsPt%P = 0.0d0
            GsPt%dPdXi(:) = 0.0d0
            GsPt%dSdXi(:,:,:) = 0.0d0  !dSjk/dxi
        end if

        ! GsPt%GU_proj(:,:) = 0.0d0
        ! do i = 1, size(TEMP_TL_proj,2)
        !     call getTensorComponent(i,j,k)
        !     GsPt%GU_proj(j,k) = dot_product(TEMP_TL_proj(:,i),Elem%bfn_p(:))
        ! end do

    end subroutine setFlowUnknownsStabilityComplex

    !---------------------------------------------------------------------

    subroutine setExtraFlowUnknownsStability(N, GsPt_d)

        use petscksp, only: PETSC_i
        use PhysicalParameters_mod, only: KwN
        use ContinuationParameters_mod, only: Stability_Im

        implicit none

        PetscInt, intent(in) :: N
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        !z-derivatives
        if (Stability_Im) then
            GsPt_d%GU(N,:) = PETSC_i*KwN*(GsPt_d%U_f(:))         !dUi/dz
            GsPt_d%dPdXi(N) = PETSC_i*KwN*(GsPt_d%P)             !dP/dz
            GsPt_d%dSdXi(:,:,N) = PETSC_i*KwN*(GsPt_d%S(:,:))    !dTij/dz
        else
            GsPt_d%GU(N,1) = -KwN*(GsPt_d%U_f(1))    !dUx/dz
            GsPt_d%GU(N,2) = -KwN*(GsPt_d%U_f(2))    !dUy/dz
            GsPt_d%GU(N,N) = KwN*(GsPt_d%U_f(N))     !dUz/dz

            GsPt_d%dPdXi(N) = -KwN*GsPt_d%P          !dP/dz
            
            GsPt_d%dSdXi(:,:,N) = -KwN*GsPt_d%S(:,:) !dTij/dz: xx, xy, yx, yy, zz
            GsPt_d%dSdXi(1,N,N) = KwN*GsPt_d%S(1,N)  !dTxz/dz
            GsPt_d%dSdXi(N,1,N) = GsPt_d%dSdXi(1,N,N)   !dTzx/dz
            GsPt_d%dSdXi(2,N,N) = KwN*GsPt_d%S(2,N)  !dTyz/dz
            GsPt_d%dSdXi(N,2,N) = GsPt_d%dSdXi(2,N,N)   !dTzy/dz
        end if

        GsPt_d%G2U(N,:) = -(KwN**2)*GsPt_d%U_f(:)

    end subroutine setExtraFlowUnknownsStability

    !----------------------------------------------------------------------

    subroutine strongFormMomentumStability(GsPt_b, GsPt_d)

        use PhysicalParameters_mod, only: ReN, BvN

        implicit none

        type(GaussPointQuantities), intent(in) :: GsPt_b
        type(GaussPointQuantities), intent(inout) :: GsPt_d

        PetscInt :: i, j
        PetscScalar, dimension(3) :: div_Tve_d, div_Ptot_d

        div_Tve_d(:) = 0.0d0
        do i = 1, 3
            do j = 1, 3
                div_Tve_d(i) = div_Tve_d(i) + GsPt_d%dTdXi(j,i,j)
            end do
        end do

        div_Ptot_d(:) = -GsPt_d%dPdXi(:) + div_Tve_d(:) &
                    + BvN*( (GsPt_d%G2U(1,:))+(GsPt_d%G2U(2,:))+(GsPt_d%G2U(3,:)) )

        GsPt_d%Convection(:) = matmul(GsPt_b%U_f,GsPt_d%GU) + matmul(GsPt_d%U_f,GsPt_b%GU)

        GsPt_d%Momentum(:) = ReN*(GsPt_d%Convection(:)) - div_Ptot_d(:)
        GsPt_d%Momentum_time(:) = ReN*(GsPt_d%U_f(:))

    end subroutine strongFormMomentumStability

end module ElementCalculationsStability_mod
