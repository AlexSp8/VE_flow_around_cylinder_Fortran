
#include <slepc/finclude/slepceps.h>

module Extra_EquationsStability_mod

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

    subroutine setExtraEquationsBulkStability(Problem, Problem_Stability, FE, Mesh, GaussInt, &
        Sol, basis_p, Elem, inex, EVP)

        use MPIParameters_mod, only: Rank
        use ElementCalculations_mod, only: copyToLocal

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscInt, intent(in) :: inex
        type(EVPType), intent(inout) :: EVP
        
        PetscInt :: iel_rank, Nunknowns, Nel, Nbf
        PetscScalar :: value, Res_a_extra, Res_b_extra
        PetscErrorCode :: ierr

        Nel = FE%Nel_Rank
        Nbf = FE%Nbf
        loop_elements:do iel_rank = 1, Nel

            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)

            call extraResidualBulkStability(Problem, Problem_Stability, Nbf, GaussInt, &
                basis_p, Elem, inex, Res_a_extra, Res_b_extra)

            call storeExtraBulkStability(Problem, Problem_Stability, Nbf, &
                Mesh%Connectivity_Rank(iel_rank,:), GaussInt, basis_p, Elem, &
                inex, Res_a_extra, Res_b_extra, EVP)

        end do loop_elements

        value = -Problem_Stability%iex_b_values(inex)

        Nunknowns = Problem_Stability%Nunknowns
        ! if (Rank == 0) call VecSetValues(LS%b_f, 1, Nunknowns+inex-1, value, ADD_VALUES, ierr)

    end subroutine setExtraEquationsBulkStability

    !--------------------------------------------------------------------------

    subroutine extraResidualBulkStability(Problem, Problem_Stability, Nbf, GaussInt, &
        basis_p, Elem, inex, Res_a_extra, Res_b_extra)

        use ElementCalculations_mod, only: setElementBasisFunctions
        use ElementCalculationsStability_mod, only: setFlowQuantitiesStability

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        PetscInt, intent(in) :: Nbf
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscInt, intent(in) :: inex
        PetscScalar, intent(out) :: Res_a_extra, Res_b_extra

        type(GaussPointQuantities) :: GsPt_b, GsPt_d
        PetscInt :: ig

        Res_a_extra = 0.0d0 ; Res_b_extra = 0.0d0
        loop_gauss:do ig = 1, GaussInt%NGauss_b

            call setElementBasisFunctions(ig, GaussInt, Elem)

            call basis_p(Problem%idir, Nbf, Elem, GsPt_b)

            call setFlowQuantitiesStability(Problem, Problem_Stability, Nbf, Elem, GsPt_b, GsPt_d)

            GsPt_d%WET = GaussInt%Weights_b(ig)*abs(GsPt_b%detJ)

            call extraEquationsResidualsBulkStability(inex, GsPt_d, Res_a_extra, Res_b_extra)

        end do loop_gauss

    end subroutine extraResidualBulkStability

    !----------------------------------------------------------------------

    subroutine extraEquationsResidualsBulkStability(inex, GsPt_d, Res_a_extra, Res_b_extra)

        implicit none

        PetscInt, intent(in) :: inex
        type(GaussPointQuantities), intent(in) :: GsPt_d
        PetscScalar, intent(out) :: Res_a_extra, Res_b_extra

        select case(inex)
        case(1)
            !Flow rate
            Res_a_extra = Res_a_extra + (GsPt_d%U_f(1))*(GsPt_d%WET)
            Res_b_extra = 0.0d0
        case default
            write(*,'(a)') "Wrong EQ choice in extraEquationsResidualsBulkStability!"
            stop
        end select

    end subroutine extraEquationsResidualsBulkStability

    !---------------------------------------------------------------------

    subroutine storeExtraBulkStability(Problem, Problem_Stability, Nbf, Connectivity, GaussInt, &
        basis_p, Elem, inex, Res_a_extra, Res_b_extra, EVP)

        use Tools_mod, only: perturbVariable, getRows
        use Storage_mod, only: storeJacobian_Ar, storeJacobian_Ah

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        PetscInt, intent(in) :: Nbf, inex
        PetscInt, dimension(:), intent(in) :: Connectivity
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, intent(in) :: Res_a_extra, Res_b_extra
        type(EVPType), intent(inout) :: EVP

        PetscInt, dimension(1), parameter :: ieqs = [1]
        PetscInt :: inod, ieq, Nrows, Mcols, irow, Neq, Nex
        PetscInt, dimension(Nbf) :: grows, gnodes
        PetscReal :: eps_jac
        PetscScalar :: dRes_a_extra, dRes_b_extra
        PetscScalar, dimension(Nbf, Problem_Stability%Neq) :: TEMP_Ar_f, TEMP_Br_f
        PetscScalar, dimension(Problem_Stability%Nex) :: TEMP_Ah_f, TEMP_Bh_f
        PetscErrorCode :: ierr

        call MatGetSize(EVP%A_f, Nrows, Mcols, ierr)

        Neq = Problem_Stability%Neq
        Nex = Problem_Stability%Nex

        gnodes(:) = Connectivity(:)
        grows(:) = getRows(gnodes(:),ieqs(:),Neq)

        irow = Nrows - Nex + inex - 1

        !Ar_f: derivatives of Extra EQs wrt Bulk unknowns
        TEMP_Ar_f(:,:) = 0.0d0 ; TEMP_Br_f(:,:) = 0.0d0
        do inod = 1, Nbf
            do ieq = 1, Neq

                eps_jac = perturbVariable(Elem%TEMP_TL_d(inod,ieq))

                Elem%TEMP_TL_d(inod,ieq) = Elem%TEMP_TL_d(inod,ieq) + eps_jac

                call extraResidualBulkStability(Problem, Problem_Stability, Nbf, &
                    GaussInt, basis_p, Elem, inex, dRes_a_extra, dRes_b_extra)

                Elem%TEMP_TL_d(inod,ieq) = Elem%TEMP_TL_d(inod,ieq) - eps_jac

                TEMP_Ar_f(inod,ieq) = (dRes_a_extra - Res_a_extra)/eps_jac
                TEMP_Br_f(inod,ieq) = (dRes_b_extra - Res_b_extra)/eps_jac

            end do
        end do

        irow = Nrows - Nex + inex - 1
        call storeJacobian_Ar(grows, TEMP_Ar_f, EVP%A_f, irow)
        call storeJacobian_Ar(grows, TEMP_Br_f, EVP%B_f, irow)

        !Ah_f: Derivatives of Extra EQs wrt to Extra unknowns
        TEMP_Ah_f(:) = 0.0d0 ; TEMP_Bh_f(:) = 0.0d0
        do ieq = 1, Nex

            eps_jac = perturbVariable(Elem%TEMP_EU_d(ieq))

            Elem%TEMP_EU_d(ieq) = Elem%TEMP_EU_d(ieq) + eps_jac

            call extraResidualBulkStability(Problem, Problem_Stability, Nbf, &
                GaussInt, basis_p, Elem, inex, dRes_a_extra, dRes_b_extra)

            Elem%TEMP_EU_d(ieq) = Elem%TEMP_EU_d(ieq) - eps_jac

            TEMP_Ah_f(ieq) = (dRes_a_extra - Res_a_extra)/eps_jac
            TEMP_Bh_f(ieq) = (dRes_b_extra - Res_b_extra)/eps_jac

        end do

        irow = Nrows - Nex + inex - 1
        call storeJacobian_Ah(TEMP_Ah_f, EVP%A_f, inex, irow)
        call storeJacobian_Ah(TEMP_Bh_f, EVP%B_f, inex, irow)

    end subroutine storeExtraBulkStability

    !---------------------------------------------------------------------

    subroutine setExtraEquationsEdgeStability(Problem, Problem_Stability, FE, Mesh, bound, &
        GaussInt, Sol, basis_p, Elem, EVP)

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        type(EVPType), intent(inout) :: EVP

        PetscInt :: inex, iex_ed

        do inex = 1, size(bound%ex_ed_S)
            iex_ed = bound%ex_ed_S(inex)
            call applyExtraEquationEdgeStability(Problem, Problem_Stability, FE, Mesh, bound, &
                GaussInt, Sol, basis_p, Elem, iex_ed, EVP)
        end do
        
    end subroutine setExtraEquationsEdgeStability

    !---------------------------------------------------------------------

    subroutine applyExtraEquationEdgeStability(Problem, Problem_Stability, FE, Mesh, bound, &
        GaussInt, Sol, basis_p, Elem, iex_ed, EVP)

        use ElementCalculations_mod, only: copyToLocal

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(SolutionArraysType), intent(in) :: Sol
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscInt, intent(in) :: iex_ed
        type(EVPType), intent(inout) :: EVP

        PetscInt :: iel_rank, iel, ied, Nbf
        PetscScalar :: Res_a_extra, Res_b_extra

        Nbf = FE%Nbf

        loop_elements:do iel = 1, size(bound%edge_elements_rank)

            iel_rank = bound%edge_elements_rank(iel)

            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)

            ied = bound%edges(iel)

            call extraResidualEdgeStability(Problem, Problem_Stability, FE, GaussInt, &
                basis_p, Elem, iex_ed, ied, Res_a_extra, Res_b_extra)

            call storeExtraEdgeStability(Problem, Problem_Stability, FE, &
                Mesh%Connectivity_Rank(iel_rank,:), GaussInt, basis_p, Elem, &
                iex_ed, ied, Res_a_extra, Res_b_extra, EVP)

        end do loop_elements

    end subroutine applyExtraEquationEdgeStability

    !---------------------------------------------------------------------

    subroutine extraResidualEdgeStability(Problem, Problem_Stability, FE, GaussInt, &
        basis_p, Elem, iex_ed, ied, Res_a_extra, Res_b_extra)

        use ElementCalculations_mod, only: setElementBasisFunctionsEdge, setEdgeQuantities
        use ElementCalculationsStability_mod, only: setFlowQuantitiesStability

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscInt, intent(in) :: iex_ed, ied
        PetscScalar, intent(out) :: Res_a_extra, Res_b_extra

        type(GaussPointQuantities) :: GsPt_b, GsPt_d
        PetscInt :: ig, Nbf

        Nbf = FE%Nbf

        Res_a_extra = 0.0d0 ; Res_b_extra = 0.0d0

        loop_gauss:do ig = 1, GaussInt%NGauss_1D

            call setElementBasisFunctionsEdge(ig, ied, GaussInt, Elem)

            call basis_p(Problem%idir, Nbf, Elem, GsPt_b)

            call setFlowQuantitiesStability(Problem, Problem_Stability, Nbf, Elem, GsPt_b, GsPt_d)

            call setEdgeQuantities(ied, Problem%idir, FE%name, GsPt_d)

            GsPt_d%WET = GaussInt%Weights_1D(ig)*(GsPt_b%dL)

            call extraEquationsResidualsEdgeStability(iex_ed, GsPt_d, Res_a_extra, Res_b_extra)

        end do loop_gauss

    end subroutine extraResidualEdgeStability

    !----------------------------------------------------------------------

    subroutine extraEquationsResidualsEdgeStability(iex_ed, GsPt_d, Res_a_extra, Res_b_extra)

        implicit none

        PetscInt, intent(in) :: iex_ed
        type(GaussPointQuantities), intent(in) :: GsPt_d
        PetscScalar, intent(out) :: Res_a_extra, Res_b_extra

        select case(iex_ed)
        case(1)
            !Flow rate
            Res_a_extra = Res_a_extra + (GsPt_d%U_f(1))*(GsPt_d%WET)
            Res_b_extra = 0.0d0
        case default
            write(*,'(a)') "Wrong EQ choice in extraEquationsResidualsBulkStability!"
            stop
        end select

    end subroutine extraEquationsResidualsEdgeStability

    !----------------------------------------------------------------------

    subroutine storeExtraEdgeStability(Problem, Problem_Stability, FE, Connectivity, &
        GaussInt, basis_p, Elem, iex_ed, ied, Res_a_extra, Res_b_extra, EVP)

        use Tools_mod, only: perturbVariable, getRows
        use Storage_mod, only: storeJacobian_Ar, storeJacobian_Ah

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:), intent(in) :: Connectivity
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscInt, intent(in) :: iex_ed, ied
        PetscScalar, intent(in) :: Res_a_extra, Res_b_extra
        type(EVPType), intent(inout) :: EVP

        PetscInt :: inod, jnod, ieq, inex, irow, Nrows, Mcols, Nbf, Neq, Nex
        PetscReal :: eps_jac
        PetscScalar :: dRes_a_extra, dRes_b_extra
        PetscScalar, dimension(FE%Nbf, Problem_Stability%Neq) :: TEMP_Ar_f, TEMP_Br_f
        PetscScalar, dimension(Problem_Stability%Nex) :: TEMP_Ah_f, TEMP_Bh_f
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
        TEMP_Ar_f(:,:) = 0.0d0 ; TEMP_Br_f(:,:) = 0.0d0
        do inod = 1, Nbf
            do ieq = 1, Neq

                eps_jac = perturbVariable(Elem%TEMP_TL_d(inod,ieq))

                Elem%TEMP_TL_d(inod,ieq) = Elem%TEMP_TL_d(inod,ieq) + eps_jac

                call extraResidualEdgeStability(Problem, Problem_Stability, FE, &
                    GaussInt, basis_p, Elem, iex_ed, ied, dRes_a_extra, dRes_b_extra)

                Elem%TEMP_TL_d(inod,ieq) = Elem%TEMP_TL_d(inod,ieq) - eps_jac

                TEMP_Ar_f(inod,ieq)   = (dRes_a_extra - Res_a_extra)/eps_jac
                TEMP_Br_f(inod,ieq)   = (dRes_b_extra - Res_b_extra)/eps_jac

            end do
        end do

        irow = Nrows - Nex + iex_ed - 1
        call storeJacobian_Ar(grows, TEMP_Ar_f, EVP%A_f, irow)
        call storeJacobian_Ar(grows, TEMP_Br_f, EVP%B_f, irow)

        !Ah_f: Derivatives of Extra EQs wrt to Extra unknowns
        TEMP_Ah_f(:) = 0.0d0 ; TEMP_Bh_f(:) = 0.0d0
        do ieq = 1, Nex

            eps_jac = perturbVariable(Elem%TEMP_EU_d(ieq))

            Elem%TEMP_EU_d(ieq) = Elem%TEMP_EU_d(ieq) + eps_jac

            call extraResidualEdgeStability(Problem, Problem_Stability, FE, &
                GaussInt, basis_p, Elem, iex_ed, ied, dRes_a_extra, dRes_b_extra)

            Elem%TEMP_EU_d(ieq) = Elem%TEMP_EU_d(ieq) - eps_jac

            TEMP_Ah_f(ieq) = (dRes_a_extra - Res_a_extra)/eps_jac
            TEMP_Bh_f(ieq) = (dRes_b_extra - Res_b_extra)/eps_jac

        end do

        irow = Nrows - Nex + iex_ed - 1
        call storeJacobian_Ah(TEMP_Ah_f, EVP%A_f, iex_ed, irow)
        call storeJacobian_Ah(TEMP_Bh_f, EVP%B_f, iex_ed, irow)

    end subroutine storeExtraEdgeStability

    !---------------------------------------------------------------------

    subroutine setExtraEquationsFaceStability(Problem, Problem_Stability, FE, Mesh, bound, &
        GaussInt, Sol, basis_p, Elem, EVP)

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        type(EVPType), intent(inout) :: EVP

        PetscInt :: inex, iex_fc

        do inex = 1, size(bound%ex_fc_S)
            iex_fc = bound%ex_fc_S(inex)
            call applyExtraEquationFaceStability(Problem, Problem_Stability, FE, Mesh, bound, &
                GaussInt, Sol, basis_p, Elem, iex_fc, EVP)
        end do
        
    end subroutine setExtraEquationsFaceStability

    !---------------------------------------------------------------------

    subroutine applyExtraEquationFaceStability(Problem, Problem_Stability, FE, Mesh, bound, &
        GaussInt, Sol, basis_p, Elem, iex_fc, EVP)

        use ElementCalculations_mod, only: copyToLocal

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscInt, intent(in) :: iex_fc
        type(EVPType), intent(inout) :: EVP

        PetscInt :: iel_rank, iel, ied, Nbf
        PetscScalar :: Res_a_extra, Res_b_extra

        Nbf = FE%Nbf
        loop_elements:do iel = 1, size(bound%face_elements_rank)

            iel_rank = bound%face_elements_rank(iel)

            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)
            ied = bound%faces(iel)

            call extraResidualFaceStability(Problem, Problem_Stability, FE, GaussInt, &
                basis_p, Elem, iex_fc, ied, Res_a_extra, Res_b_extra)

            call storeExtraFaceStability(Problem, Problem_Stability, FE, &
                Mesh%Connectivity_Rank(iel_rank,:), GaussInt, basis_p, Elem, &
                iex_fc, ied, Res_a_extra, Res_b_extra, EVP)

        end do loop_elements

    end subroutine applyExtraEquationFaceStability

    !---------------------------------------------------------------------

    subroutine extraResidualFaceStability(Problem, Problem_Stability, FE, &
        GaussInt, basis_p, Elem, iex_fc, ifc, Res_a_extra, Res_b_extra)

        use ElementCalculations_mod, only: setElementBasisFunctionsFace, setFaceQuantities
        use ElementCalculationsStability_mod, only: setFlowQuantitiesStability

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscInt, intent(in) :: iex_fc, ifc
        PetscScalar, intent(out) :: Res_a_extra, Res_b_extra

        type(GaussPointQuantities) :: GsPt_b, GsPt_d
        PetscInt :: ig, Nbf

        Nbf = FE%Nbf

        Res_a_extra = 0.0d0 ; Res_b_extra = 0.0d0

        loop_gauss:do ig = 1, GaussInt%NGauss_2D

            call setElementBasisFunctionsFace(ig, ifc, GaussInt, Elem)

            call basis_p(Problem%idir, Nbf, Elem, GsPt_b)

            call setFlowQuantitiesStability(Problem, Problem_Stability, Nbf, Elem, GsPt_b, GsPt_d)

            call setFaceQuantities(ifc, FE%name, GsPt_d)

            GsPt_d%WET = GaussInt%Weights_2D(ig)*(GsPt_b%dL)

            call extraEquationsResidualsFaceStability(iex_fc, GsPt_d, Res_a_extra, Res_b_extra)

        end do loop_gauss

    end subroutine extraResidualFaceStability

    !----------------------------------------------------------------------

    subroutine extraEquationsResidualsFaceStability(iex_fc, GsPt_d, Res_a_extra, Res_b_extra)

        implicit none

        PetscInt, intent(in) :: iex_fc
        type(GaussPointQuantities), intent(in) :: GsPt_d
        PetscScalar, intent(out) :: Res_a_extra, Res_b_extra

        select case(iex_fc)
        case(1)
            !Flow rate
            Res_a_extra = Res_a_extra + (GsPt_d%U_f(1))*(GsPt_d%WET)
            Res_b_extra = 0.0d0
        case default
            write(*,'(a)') "Wrong EQ choice in extraEquationsResidualsFaceStability!"
            stop
        end select

    end subroutine extraEquationsResidualsFaceStability

    !----------------------------------------------------------------------

    subroutine storeExtraFaceStability(Problem, Problem_Stability, FE, Connectivity, &
        GaussInt, basis_p, Elem, iex_fc, ifc, Res_a_extra, Res_b_extra, EVP)

        use Tools_mod, only: perturbVariable, getRows
        use Storage_mod, only: storeJacobian_Ar, storeJacobian_Ah

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:), intent(in) :: Connectivity
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscInt, intent(in) :: iex_fc, ifc
        PetscScalar, intent(in) :: Res_a_extra, Res_b_extra
        type(EVPType), intent(inout) :: EVP

        PetscInt :: inod, jnod, ieq, inex, irow, Nrows, Mcols, Nbf, Neq, Nex
        PetscReal :: eps_jac
        PetscScalar :: dRes_a_extra, dRes_b_extra
        PetscScalar, dimension(FE%Nbf, Problem_Stability%Neq) :: TEMP_Ar_f, TEMP_Br_f
        PetscScalar, dimension(Problem_Stability%Nex) :: TEMP_Ah_f, TEMP_Bh_f
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
        TEMP_Ar_f(:,:) = 0.0d0 ; TEMP_Br_f(:,:) = 0.0d0
        do inod = 1, Nbf
            do ieq = 1, Neq

                eps_jac = perturbVariable(Elem%TEMP_TL_d(inod,ieq))

                Elem%TEMP_TL_d(inod,ieq) = Elem%TEMP_TL_d(inod,ieq) + eps_jac

                call extraResidualFaceStability(Problem, Problem_Stability, FE, &
                    GaussInt, basis_p, Elem, iex_fc, ifc, dRes_a_extra, dRes_b_extra)

                Elem%TEMP_TL_d(inod,ieq) = Elem%TEMP_TL_d(inod,ieq) - eps_jac

                TEMP_Ar_f(inod,ieq)   = (dRes_a_extra - Res_a_extra)/eps_jac
                TEMP_Br_f(inod,ieq)   = (dRes_b_extra - Res_b_extra)/eps_jac

            end do
        end do

        irow = Nrows - Nex + iex_fc - 1
        call storeJacobian_Ar(grows, TEMP_Ar_f, EVP%A_f, irow)
        call storeJacobian_Ar(grows, TEMP_Br_f, EVP%B_f, irow)

        !Ah_f: Derivatives of Extra EQs wrt to Extra unknowns
        TEMP_Ah_f(:) = 0.0d0 ; TEMP_Bh_f(:) = 0.0d0
        do ieq = 1, Nex

            eps_jac = perturbVariable(Elem%TEMP_EU_d(ieq))

            Elem%TEMP_EU_d(ieq) = Elem%TEMP_EU_d(ieq) + eps_jac

            call extraResidualFaceStability(Problem, Problem_Stability, FE, &
                GaussInt, basis_p, Elem, iex_fc, ifc, dRes_a_extra, dRes_b_extra)

            Elem%TEMP_EU_d(ieq) = Elem%TEMP_EU_d(ieq) - eps_jac

            TEMP_Ah_f(ieq) = (dRes_a_extra - Res_a_extra)/eps_jac
            TEMP_Bh_f(ieq) = (dRes_b_extra - Res_b_extra)/eps_jac

        end do

        irow = Nrows - Nex + iex_fc - 1
        call storeJacobian_Ah(TEMP_Ah_f, EVP%A_f, iex_fc, irow)
        call storeJacobian_Ah(TEMP_Bh_f, EVP%B_f, iex_fc, irow)

    end subroutine storeExtraFaceStability

end module Extra_EquationsStability_mod
