
#include <slepc/finclude/slepceps.h>

module EquationsStability_mod

    use PhysicalParameters_mod, only: ProblemParameters
    use FEMParameters_mod, only: FEParameters
    use MeshParameters_mod, only: MeshParameters
    use BoundaryParameters_mod, only: BoundaryParameters
    use GaussParameters_mod, only: GaussIntegration
    use SolutionVariables_mod, only: SolutionArraysType
    use ElementVariables_mod, only: GaussPointQuantities, NodeArrays
    use SlepcVariables_mod, only: EVPType
    use ElementCalculations_mod, only: basis

    implicit none

    contains

    subroutine setEquationsBulkStability(Problem, Problem_Stability, FE, Mesh, Boundary, &
        GaussInt, Sol, basis_p, Elem, EVP)

        use ElementCalculations_mod, only: copyToLocal

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        type(EVPType), intent(inout) :: EVP

        PetscInt :: iel_rank, Nel_rank, Nbf
        PetscScalar, dimension(FE%Nbf, Problem_Stability%Neq) :: Res_a, Res_b

        Nel_rank = FE%Nel_Rank
        Nbf = FE%Nbf
        loop_elements:do iel_rank = 1, Nel_rank

            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)

            ! Res_a(:,:) = 0.0d0
            ! Res_b(:,:) = 0.0d0
            call residualBulkStability(Problem, Problem_Stability, Nbf, GaussInt, basis_p, &
                Elem, Res_a, Res_b)

            call storeBulkStability(Problem, Problem_Stability, Nbf, &
                Mesh%Connectivity_Rank(iel_rank,:), Boundary, GaussInt, basis_p, &
                Elem, Res_a, Res_b, EVP)

        end do loop_elements

    end subroutine setEquationsBulkStability

    !---------------------------------------------------------------------

    subroutine residualBulkStability(Problem, Problem_Stability, Nbf, GaussInt, &
        basis_p, Elem, Res_a, Res_b)

        use ElementCalculations_mod, only: setElementBasisFunctions
        use ElementCalculationsStability_mod, only: setFlowQuantitiesStability

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        PetscInt, intent(in) :: Nbf
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, dimension(:,:), intent(out) :: Res_a, Res_b

        PetscInt :: ig
        type(GaussPointQuantities) :: GsPt_b, GsPt_d

        Res_a(:,:) = 0.0d0 ; Res_b(:,:) = 0.0d0
        loop_gauss:do ig = 1, GaussInt%NGauss_b

            call setElementBasisFunctions(ig, GaussInt, Elem)

            call basis_p(Problem%idir, Nbf, Elem, GsPt_b)

            call setFlowQuantitiesStability(Problem, Problem_Stability, Nbf, Elem, GsPt_b, GsPt_d)

            GsPt_b%WET = (GaussInt%Weights_b(ig))*abs(GsPt_b%detJ)

            call equationsResidualsBulkStability(Problem_Stability, Nbf, GsPt_b, GsPt_d, Res_a, Res_b)

        end do loop_gauss

    end subroutine residualBulkStability

    !----------------------------------------------------------------------
    
    subroutine equationsResidualsBulkStability(Problem_Stability, Nbf, &
        GsPt_b, GsPt_d, Res_a, Res_b)

        use ResidualsStability_mod, only: momentumStability_p, &
            continuityStability, stressStability_p

        implicit none

        type(ProblemParameters), intent(in) :: Problem_Stability
        PetscInt, intent(in) :: Nbf
        type(GaussPointQuantities), intent(in) :: GsPt_b, GsPt_d
        PetscScalar, dimension(:,:), intent(out) :: Res_a, Res_b

        PetscInt :: inod
        PetscScalar, dimension(Problem_Stability%Neq) :: TERM_RES_a, TERM_RES_b

        loop_nodes:do inod = 1, Nbf

            TERM_RES_a(:) = 0.0d0 ; TERM_RES_b(:) = 0.0d0

            call momentumStability_p(inod, GsPt_b, GsPt_d, TERM_RES_a, TERM_RES_b)

            call continuityStability(inod, GsPt_b, GsPt_d, TERM_RES_a, TERM_RES_b)

            call stressStability_p(inod, GsPt_b, GsPt_d, TERM_RES_a, TERM_RES_b)

            TERM_RES_b(:) = -TERM_RES_b(:)

            Res_a(inod,:) = Res_a(inod,:) + TERM_RES_a(:)*(GsPt_b%WET)
            Res_b(inod,:) = Res_b(inod,:) + TERM_RES_b(:)*(GsPt_b%WET)

        end do loop_nodes

    end subroutine equationsResidualsBulkStability

    !----------------------------------------------------------------------

    subroutine storeBulkStability(Problem, Problem_Stability, Nbf, Connectivity, Boundary, &
        GaussInt, basis_p, Elem, Res_a, Res_b, EVP)

        use Tools_mod, only: perturbVariable, getRows
        use Storage_mod, only: storeJacobian, storeJacobian_Ac
        use MPIParameters_mod, only: Rank

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        PetscInt, intent(in) :: Nbf
        PetscInt, dimension(:), intent(in) :: Connectivity
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, dimension(:,:), intent(in) :: Res_a, Res_b
        type(EVPType), intent(inout) :: EVP

        PetscInt :: inod, jnod, ieq, inex, jcol, Nrows, Mcols, Neq, Nex
        PetscReal :: eps_jac
        ! PetscScalar :: eps_jac
        PetscScalar, dimension(Nbf, Problem_Stability%Neq) :: dRes_a, dRes_b
        PetscScalar, dimension(Nbf, Problem_Stability%Neq, &
                                Nbf, Problem_Stability%Neq) :: TEMP_A_f, TEMP_B_f
        PetscScalar, dimension(Nbf, Problem_Stability%Neq, &
                                Problem_Stability%Nex) :: TEMP_Ac_f, TEMP_Bc_f
        PetscErrorCode :: ierr
        PetscInt, dimension(Nbf) :: grows, gnodes
        PetscInt, dimension(1), parameter :: ieqs = [1]

        Neq = Problem_Stability%Neq
        Nex = Problem_Stability%Nex

        call MatGetSize(EVP%A_f, Nrows, Mcols, ierr)

        gnodes(:) = Connectivity(:)
        grows(:) = getRows(gnodes,ieqs,Neq)

        !A_f: derivatives of Bulk EQs wrt Bulk unknowns
        TEMP_A_f(:,:,:,:) = 0.0d0 ; TEMP_B_f(:,:,:,:) = 0.0d0
        do inod = 1, Nbf
            do ieq = 1, Neq

                eps_jac = perturbVariable(Elem%TEMP_TL_d(inod,ieq))

                Elem%TEMP_TL_d(inod,ieq) = Elem%TEMP_TL_d(inod,ieq) + eps_jac

                call residualBulkStability(Problem, Problem_Stability, Nbf, GaussInt, &
                    basis_p, Elem, dRes_a, dRes_b)

                Elem%TEMP_TL_d(inod,ieq) = Elem%TEMP_TL_d(inod,ieq) - eps_jac

                do jnod = 1, Nbf
                    TEMP_A_f(inod,ieq,jnod,:) = (dRes_a(jnod,:) - Res_a(jnod,:))/PetscRealPart(eps_jac)
                    TEMP_B_f(inod,ieq,jnod,:) = (dRes_b(jnod,:) - Res_b(jnod,:))/PetscRealPart(eps_jac)
                end do

            end do
        end do

        call storeJacobian(TEMP_A_f, grows, EVP%A_f)
        call storeJacobian(TEMP_B_f, grows, EVP%B_f)

        !Ac_f: Derivatives of Bulk EQs wrt to global unknowns
        TEMP_Ac_f(:,:,:) = 0.0d0 ; TEMP_Bc_f(:,:,:) = 0.0d0
        do inex = 1, Nex

            eps_jac = perturbVariable(Elem%TEMP_EU_d(inex))

            Elem%TEMP_EU_d(inex) = Elem%TEMP_EU_d(inex) + eps_jac

            call residualBulkStability(Problem, Problem_Stability, Nbf, GaussInt, &
                basis_p, Elem, dRes_a, dRes_b)

            Elem%TEMP_EU_d(inex) = Elem%TEMP_EU_d(inex) - eps_jac

            do inod = 1, Nbf
                TEMP_Ac_f(inod,:,inex) = (dRes_a(inod,:) - Res_a(inod,:))/PetscRealPart(eps_jac)
                TEMP_Bc_f(inod,:,inex) = (dRes_b(inod,:) - Res_b(inod,:))/PetscRealPart(eps_jac)
            end do

            jcol = Mcols - Nex + inex -1
            call storeJacobian_Ac(TEMP_Ac_f(:,:,inex), grows, EVP%A_f, jcol)
            call storeJacobian_Ac(TEMP_Bc_f(:,:,inex), grows, EVP%B_f, jcol)

        end do

        call storePeriodicStability(Problem_Stability, Nbf, Connectivity, Boundary, &
            TEMP_A_f, TEMP_B_f, TEMP_Ac_f, TEMP_Bc_f, EVP)

    end subroutine storeBulkStability

    !---------------------------------------------------------------------

    subroutine setEquationsEdgeStability(Problem, Problem_Stability, FE, Mesh, bound, &
        Boundary, GaussInt, Sol, basis_p, Elem, EVP)

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        type(EVPType), intent(inout) :: EVP

        PetscInt :: ieq, ieq_ed

        do ieq = 1, size(bound%eq_ed_S)
            ieq_ed = bound%eq_ed_S(ieq)
            call applyEquationEdgeStability(Problem, Problem_Stability, FE, Mesh, bound, &
                Boundary, GaussInt, Sol, basis_p, Elem, ieq_ed, EVP)
        end do
        
    end subroutine setEquationsEdgeStability

    !---------------------------------------------------------------------

    subroutine applyEquationEdgeStability(Problem, Problem_Stability, FE, Mesh, bound, &
        Boundary, GaussInt, Sol, basis_p, Elem, ieq_ed, EVP)

        use ElementCalculations_mod, only: copyToLocal

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(SolutionArraysType), intent(in) :: Sol
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscInt, intent(in) :: ieq_ed
        type(EVPType), intent(inout) :: EVP

        PetscInt :: iel_rank, iel, ied, Nbf
        PetscScalar, dimension(FE%Nbf, Problem_Stability%Neq) :: Res_a, Res_b

        Nbf = FE%Nbf

        loop_elements:do iel = 1, size(bound%edge_elements_rank)

            iel_rank = bound%edge_elements_rank(iel)

            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)

            ied = bound%edges(iel)

            call residualEdgeStability(Problem, Problem_Stability, FE, GaussInt, &
                basis_p, Elem, ieq_ed, ied, Res_a, Res_b)

            call storeEdgeStability(Problem, Problem_Stability, FE, &
                Mesh%Connectivity_Rank(iel_rank,:), Boundary, GaussInt, basis_p, Elem, &
                ieq_ed, ied, Res_a, Res_b, EVP)

        end do loop_elements

    end subroutine applyEquationEdgeStability

    !---------------------------------------------------------------------

    subroutine residualEdgeStability(Problem, Problem_Stability, FE, GaussInt, &
        basis_p, Elem, ieq_ed, ied, Res_a, Res_b)

        use ElementCalculations_mod, only: setElementBasisFunctionsEdge, setEdgeQuantities
        use ElementCalculationsStability_mod, only: setFlowQuantitiesStability

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscInt, intent(in) :: ieq_ed, ied
        PetscScalar, dimension(:,:), intent(out) :: Res_a, Res_b

        type(GaussPointQuantities) :: GsPt_b, GsPt_d
        PetscInt :: ig, Nbf

        Nbf = FE%Nbf
        Res_a(:,:) = 0.0d0 ; Res_b(:,:) = 0.0d0

        loop_gauss:do ig = 1, GaussInt%NGauss_1D

            call setElementBasisFunctionsEdge(ig, ied, GaussInt, Elem)

            call basis_p(Problem%idir, Nbf, Elem, GsPt_b)

            call setFlowQuantitiesStability(Problem, Problem_Stability, Nbf, Elem, GsPt_b, GsPt_d)

            call setEdgeQuantities(ied, Problem%idir, FE%name, GsPt_b)

            GsPt_b%WET = GaussInt%Weights_1D(ig)*(GsPt_b%dL)

            call equationsResidualsEdgeStability(Problem_Stability, Nbf, ieq_ed, &
                GsPt_b, GsPt_d, Res_a, Res_b)

        end do loop_gauss

    end subroutine residualEdgeStability

    !----------------------------------------------------------------------

    subroutine equationsResidualsEdgeStability(Problem_Stability, Nbf, ieq_ed, &
        GsPt_b, GsPt_d, Res_a, Res_b)

        use ResidualsStability_mod, only: momentumOpenBCBoundaryStability, &
            momentumSlipBoundaryStability

        implicit none

        type(ProblemParameters), intent(in) :: Problem_Stability
        PetscInt, intent(in) :: Nbf, ieq_ed
        type(GaussPointQuantities), intent(in) :: GsPt_b, GsPt_d
        PetscScalar, dimension(:,:), intent(out) :: Res_a, Res_b

        PetscInt :: inod
        PetscScalar, dimension(Problem_Stability%Neq) :: TERM_RES_a, TERM_RES_b

        loop_nodes:do inod = 1, Nbf

            TERM_RES_a(:) = 0.0d0 ; TERM_RES_b(:) = 0.0d0

            select case (ieq_ed)
            case (1)
                call momentumOpenBCBoundaryStability(inod, GsPt_b, GsPt_d, TERM_RES_a, TERM_RES_b)
            case (2)
                call momentumSlipBoundaryStability(inod, GsPt_b, GsPt_d, TERM_RES_a, TERM_RES_b)
            case default
                write(*,'(a)') "Wrong EQ choice in equationsResidualsEdgeStability!"
                stop
            end select

            Res_a(inod,:) = Res_a(inod,:) + TERM_RES_a(:)*GsPt_b%WET
            Res_b(inod,:) = Res_b(inod,:) + TERM_RES_b(:)*GsPt_b%WET

        end do loop_nodes

    end subroutine equationsResidualsEdgeStability

    !----------------------------------------------------------------------

    subroutine storeEdgeStability(Problem, Problem_Stability, FE, Connectivity, &
        Boundary, GaussInt, basis_p, Elem, ieq_ed, ied, Res_a, Res_b, EVP)

        use Tools_mod, only: perturbVariable, getRows
        use Storage_mod, only: storeJacobian, storeJacobian_Ac

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:), intent(in) :: Connectivity
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscInt, intent(in) :: ieq_ed, ied
        PetscScalar, dimension(:,:), intent(in) :: Res_a, Res_b
        type(EVPType), intent(inout) :: EVP

        PetscInt :: inod, jnod, ieq, inex, jcol, Nrows, Mcols, Nbf, Neq, Nex
        PetscReal :: eps_jac
        ! PetscScalar :: eps_jac
        PetscScalar, dimension(FE%Nbf, Problem_Stability%Neq) :: dRes_a, dRes_b
        PetscScalar, dimension(FE%Nbf, Problem_Stability%Neq, &
                                FE%Nbf, Problem_Stability%Neq) :: TEMP_A_f, TEMP_B_f
        PetscScalar, dimension(FE%Nbf, Problem_Stability%Neq, &
                                Problem_Stability%Nex) :: TEMP_Ac_f, TEMP_Bc_f
        PetscErrorCode :: ierr
        PetscInt, dimension(FE%Nbf) :: grows, gnodes
        PetscInt, dimension(1), parameter :: ieqs = [1]

        Nbf = FE%Nbf
        Neq = Problem_Stability%Neq
        Nex = Problem_Stability%Nex

        call MatGetSize(EVP%A_f, Nrows, Mcols, ierr)

        gnodes(:) = Connectivity(:)
        grows(:) = getRows(gnodes(:),ieqs(:),Neq)

        !A_f: derivatives of Bulk EQs wrt Bulk unknowns
        TEMP_A_f(:,:,:,:) = 0.0d0 ; TEMP_B_f(:,:,:,:) = 0.0d0
        do inod = 1, Nbf
            do ieq = 1, Neq

                eps_jac = perturbVariable(Elem%TEMP_TL_d(inod,ieq))

                Elem%TEMP_TL(inod,ieq) = Elem%TEMP_TL_d(inod,ieq) + eps_jac

                call residualEdgeStability(Problem, Problem_Stability, FE, GaussInt, &
                    basis_p, Elem, ieq_ed, ied, dRes_a, dRes_b)

                Elem%TEMP_TL(inod,ieq) = Elem%TEMP_TL_d(inod,ieq) - eps_jac

                do jnod = 1, Nbf
                    TEMP_A_f(inod,ieq,jnod,:) = (dRes_a(jnod,:) - Res_a(jnod,:))/PetscRealPart(eps_jac)
                    TEMP_B_f(inod,ieq,jnod,:) = (dRes_b(jnod,:) - Res_b(jnod,:))/PetscRealPart(eps_jac)
                end do

            end do
        end do

        call storeJacobian(TEMP_A_f, grows, EVP%A_f)
        call storeJacobian(TEMP_B_f, grows, EVP%B_f)

        !Ac_f: Derivatives of Bulk EQs wrt to global unknowns
        TEMP_Ac_f(:,:,:) = 0.0d0 ; TEMP_Bc_f(:,:,:) = 0.0d0
        do inex = 1, Nex

            eps_jac = perturbVariable(Elem%TEMP_EU_d(inex))

            Elem%TEMP_EU(inex) = Elem%TEMP_EU_d(inex) + eps_jac

            call residualEdgeStability(Problem, Problem_Stability, FE, GaussInt, &
                basis_p, Elem, ieq_ed, ied, dRes_a, dRes_b)

            Elem%TEMP_EU(inex) = Elem%TEMP_EU_d(inex) - eps_jac

            do inod = 1, Nbf
                TEMP_Ac_f(inod,:,inex) = (dRes_a(inod,:) - Res_a(inod,:))/PetscRealPart(eps_jac)
                TEMP_Bc_f(inod,:,inex) = (dRes_b(inod,:) - Res_b(inod,:))/PetscRealPart(eps_jac)
            end do

            jcol = Mcols - Nex + inex -1
            call storeJacobian_Ac(TEMP_Ac_f(:,:,inex), grows, EVP%A_f, jcol)
            call storeJacobian_Ac(TEMP_Bc_f(:,:,inex), grows, EVP%B_f, jcol)

        end do

        call storePeriodicStability(Problem_Stability, Nbf, Connectivity, Boundary, &
            TEMP_A_f, TEMP_B_f, TEMP_Ac_f, TEMP_Bc_f, EVP)

    end subroutine storeEdgeStability

    !---------------------------------------------------------------------

    subroutine setEquationsFaceStability(Problem, Problem_Stability, FE, Mesh, bound, &
        Boundary, GaussInt, Sol, basis_p, Elem, EVP)

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        type(EVPType), intent(inout) :: EVP

        PetscInt :: ieq, ieq_fc

        do ieq = 1, size(bound%eq_fc_S)
            ieq_fc = bound%eq_fc_S(ieq)
            call applyEquationFaceStability(Problem, Problem_Stability, FE, Mesh, bound, &
                Boundary, GaussInt, Sol, basis_p, Elem, ieq_fc, EVP)
        end do
        
    end subroutine setEquationsFaceStability

    !---------------------------------------------------------------------

    subroutine applyEquationFaceStability(Problem, Problem_Stability, FE, Mesh, bound, &
        Boundary, GaussInt, Sol, basis_p, Elem, ieq_fc, EVP)

        use ElementCalculations_mod, only: copyToLocal

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscInt, intent(in) :: ieq_fc
        type(EVPType), intent(inout) :: EVP

        PetscInt :: iel_rank, iel, ifc, Nbf
        PetscScalar, dimension(FE%Nbf, Problem_Stability%Neq) :: Res_a, Res_b

        Nbf = FE%Nbf

        loop_elements:do iel = 1, size(bound%face_elements_rank)

            iel_rank = bound%face_elements_rank(iel)

            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)

            ifc = bound%faces(iel)

            call residualFaceStability(Problem, Problem_Stability, FE, GaussInt, &
                basis_p, Elem, ieq_fc, ifc, Res_a, Res_b)

            call storeFaceStability(Problem, Problem_Stability, FE, &
                Mesh%Connectivity_Rank(iel_rank,:), Boundary, GaussInt, basis_p, Elem, &
                ieq_fc, ifc, Res_a, Res_b, EVP)

        end do loop_elements

    end subroutine applyEquationFaceStability

    !---------------------------------------------------------------------

    subroutine residualFaceStability(Problem, Problem_Stability, FE, GaussInt, basis_p, &
        Elem, ieq_fc, ifc, Res_a, Res_b)

        use ElementCalculations_mod, only: setElementBasisFunctionsFace, setFaceQuantities
        use ElementCalculationsStability_mod, only: setFlowQuantitiesStability

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscInt, intent(in) :: ieq_fc, ifc
        PetscScalar, dimension(:,:), intent(out) :: Res_a, Res_b

        type(GaussPointQuantities) :: GsPt_b, GsPt_d
        PetscInt :: ig, Nbf

        Nbf = FE%Nbf
        Res_a(:,:) = 0.0d0 ; Res_b(:,:) = 0.0d0

        loop_gauss:do ig = 1, GaussInt%NGauss_2D

            call setElementBasisFunctionsFace(ig, ifc, GaussInt, Elem)

            call basis_p(Problem%idir, Nbf, Elem, GsPt_b)

            call setFlowQuantitiesStability(Problem, Problem_Stability, Nbf, Elem, GsPt_b, GsPt_d)

            call setFaceQuantities(ifc, FE%name, GsPt_b)

            GsPt_b%WET = GaussInt%Weights_2D(ig)*(GsPt_b%dL)

            call equationsResidualsFaceStability(Problem_Stability, Nbf, ieq_fc, &
                GsPt_b, GsPt_d, Res_a, Res_b)

        end do loop_gauss

    end subroutine residualFaceStability

    !----------------------------------------------------------------------

    subroutine equationsResidualsFaceStability(Problem_Stability, Nbf, ieq_fc, &
        GsPt_b, GsPt_d, Res_a, Res_b)

        use ResidualsStability_mod, only: momentumOpenBCBoundaryStability, momentumSlipBoundaryStability

        implicit none

        type(ProblemParameters), intent(in) :: Problem_Stability
        PetscInt, intent(in) :: Nbf, ieq_fc
        type(GaussPointQuantities), intent(in) :: GsPt_b, GsPt_d
        PetscScalar, dimension(:,:), intent(out) :: Res_a, Res_b

        PetscInt :: inod
        PetscScalar, dimension(Problem_Stability%Neq) :: TERM_RES_a, TERM_RES_b

        loop_nodes:do inod = 1, Nbf

            TERM_RES_a(:) = 0.0d0 ; TERM_RES_b(:) = 0.0d0

            select case (ieq_fc)
            case (1)
                call momentumOpenBCBoundaryStability(inod, GsPt_b, GsPt_d, TERM_RES_a, TERM_RES_b)
            case (2)
                call momentumSlipBoundaryStability(inod, GsPt_b, GsPt_d, TERM_RES_a, TERM_RES_b)
            case default
                write(*,'(a)') "Wrong EQ choice in equationsResidualsFaceStability!"
                stop
            end select

            Res_a(inod,:) = Res_a(inod,:) + TERM_RES_a(:)*(GsPt_b%WET)
            Res_b(inod,:) = Res_b(inod,:) + TERM_RES_b(:)*(GsPt_b%WET)

        end do loop_nodes

    end subroutine equationsResidualsFaceStability

    !----------------------------------------------------------------------

    subroutine storeFaceStability(Problem, Problem_Stability, FE, Connectivity, &
        Boundary, GaussInt, basis_p, Elem, ieq_fc, ifc, Res_a, Res_b, EVP)

        use Tools_mod, only: perturbVariable, getRows
        use Storage_mod, only: storeJacobian, storeJacobian_Ac

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:), intent(in) :: Connectivity
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscInt, intent(in) :: ieq_fc, ifc
        PetscScalar, dimension(:,:), intent(in) :: Res_a, Res_b
        type(EVPType), intent(inout) :: EVP

        PetscInt :: inod, jnod, ieq, inex, jcol, Nrows, Mcols, Nbf, Neq, Nex
        PetscReal :: eps_jac
        ! PetscScalar :: eps_jac
        PetscScalar, dimension(FE%Nbf, Problem_Stability%Neq) :: dRes_a, dRes_b
        PetscScalar, dimension(FE%Nbf, Problem_Stability%Neq, &
                                FE%Nbf, Problem_Stability%Neq) :: TEMP_A_f, TEMP_B_f
        PetscScalar, dimension(FE%Nbf, Problem_Stability%Neq, &
                                Problem_Stability%Nex) :: TEMP_Ac_f, TEMP_Bc_f
        PetscErrorCode :: ierr
        PetscInt, dimension(FE%Nbf) :: grows, gnodes
        PetscInt, dimension(1), parameter :: ieqs = [1]

        Nbf = FE%Nbf
        Neq = Problem_Stability%Neq
        Nex = Problem_Stability%Nex

        call MatGetSize(EVP%A_f, Nrows, Mcols, ierr)

        gnodes(:) = Connectivity(:)
        grows(:) = getRows(gnodes(:),ieqs(:),Neq)

        !A_f: derivatives of Bulk EQs wrt Bulk unknowns
        TEMP_A_f(:,:,:,:) = 0.0d0 ; TEMP_B_f(:,:,:,:) = 0.0d0
        do inod = 1, Nbf
            do ieq = 1, Neq

                eps_jac = perturbVariable(Elem%TEMP_TL_d(inod,ieq))

                Elem%TEMP_TL_d(inod,ieq) = Elem%TEMP_TL_d(inod,ieq) + eps_jac

                call residualFaceStability(Problem, Problem_Stability, FE, GaussInt, &
                    basis_p, Elem, ieq_fc, ifc, dRes_a, dRes_b)

                Elem%TEMP_TL_d(inod,ieq) = Elem%TEMP_TL_d(inod,ieq) - eps_jac

                do jnod = 1, Nbf
                    TEMP_A_f(inod,ieq,jnod,:) = (dRes_a(jnod,:) - Res_a(jnod,:))/PetscRealPart(eps_jac)
                    TEMP_B_f(inod,ieq,jnod,:) = (dRes_b(jnod,:) - Res_b(jnod,:))/PetscRealPart(eps_jac)
                end do

            end do
        end do

        call storeJacobian(TEMP_A_f, grows, EVP%A_f)
        call storeJacobian(TEMP_B_f, grows, EVP%B_f)

        !Ac_f: Derivatives of Bulk EQs wrt to global unknowns
        TEMP_Ac_f(:,:,:) = 0.0d0 ; TEMP_Bc_f(:,:,:) = 0.0d0
        do inex = 1, Nex

            eps_jac = perturbVariable(Elem%TEMP_EU_d(inex))

            Elem%TEMP_EU_d(inex) = Elem%TEMP_EU_d(inex) + eps_jac

            call residualFaceStability(Problem, Problem_Stability, FE, GaussInt, &
                basis_p, Elem, ieq_fc, ifc, dRes_a, dRes_b)

            Elem%TEMP_EU_d(inex) = Elem%TEMP_EU_d(inex) - eps_jac

            do inod = 1, Nbf
                TEMP_Ac_f(inod,:,inex) = (dRes_a(inod,:) - Res_a(inod,:))/PetscRealPart(eps_jac)
                TEMP_Bc_f(inod,:,inex) = (dRes_b(inod,:) - Res_b(inod,:))/PetscRealPart(eps_jac)
            end do

            jcol = Mcols - Nex + inex -1
            call storeJacobian_Ac(TEMP_Ac_f(:,:,inex), grows, EVP%A_f, jcol)
            call storeJacobian_Ac(TEMP_Bc_f(:,:,inex), grows, EVP%B_f, jcol)

        end do

        call storePeriodicStability(Problem_Stability, Nbf, Connectivity, Boundary, &
            TEMP_A_f, TEMP_B_f, TEMP_Ac_f, TEMP_Bc_f, EVP)

    end subroutine storeFaceStability

    !---------------------------------------------------------------------

    subroutine storePeriodicStability(Problem_Stability, Nbf, Connectivity, Boundary, &
        TEMP_A_f, TEMP_B_f, TEMP_Ac_f, TEMP_Bc_f, EVP)

        use petscksp, only: PETSC_i
        use PhysicalParameters_mod, only: QbN, PI
        use Tools_mod, only: indexInArray, getRows
        use Storage_mod, only: storeJacobianPeriodic, storeJacobian_Ac

        implicit none

        type(ProblemParameters), intent(in) :: Problem_Stability
        PetscInt, intent(in) :: Nbf
        PetscInt, dimension(:), intent(in) :: Connectivity
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        PetscScalar, dimension(:,:,:,:), intent(in) :: TEMP_A_f, TEMP_B_f
        PetscScalar, dimension(:,:,:), intent(in) :: TEMP_Ac_f, TEMP_Bc_f
        type(EVPType), intent(inout) :: EVP

        PetscInt, dimension(1), parameter :: ieqs = [1]
        PetscInt :: jnod, inex, Nrows, Mcols, jcol, i, ibnd, gnod, gnod_index
        PetscInt, dimension(Nbf) :: grows, gnodes
        PetscInt, dimension(1) :: grows_p, gnod_p
        PetscErrorCode :: ierr
        PetscScalar :: value

        value = exp(2.0d0*PI*QbN*PETSC_i)

        call MatGetSize(EVP%A_f, Nrows, Mcols, ierr)

        gnodes(:) = Connectivity(:)
        grows(:) = getRows(gnodes,ieqs,Problem_Stability%Neq)

        do i = 1, size(Problem_Stability%PeriodicBNDs,1)

            ibnd = Problem_Stability%PeriodicBNDs(i,1)

            do jnod = 1, size(gnodes)

                gnod = gnodes(jnod)
                gnod_index = indexInArray(Boundary(ibnd)%gnodes_total, gnod)

                if (gnod_index /= 0) then

                    gnod_p(1) = Boundary(ibnd)%gnodes_periodic(gnod_index)
                    grows_p(:) = getRows(gnod_p,ieqs,Problem_Stability%Neq)

                    call storeJacobianPeriodic &
                        (TEMP_A_f(:,:,jnod:jnod,:)*value, grows, grows_p, EVP%A_f)
                    call storeJacobianPeriodic &
                        (TEMP_B_f(:,:,jnod:jnod,:)*value, grows, grows_p, EVP%B_f)
                    do inex = 1, Problem_Stability%Nex
                        jcol = Mcols - Problem_Stability%Nex + inex -1
                        call storeJacobian_Ac &
                            (TEMP_Ac_f(jnod:jnod,:,inex)*value, grows_p, EVP%A_f, jcol)
                        call storeJacobian_Ac &
                            (TEMP_Bc_f(jnod:jnod,:,inex)*value, grows_p, EVP%B_f, jcol)
                    end do

                end if

            end do

        end do

    end subroutine storePeriodicStability

end module EquationsStability_mod
