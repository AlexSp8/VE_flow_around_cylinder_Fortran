
#include <petsc/finclude/petscksp.h>

module Extra_Equations_mod

    use petscksp
    use PhysicalParameters_mod, only: ProblemParameters
    use FEMParameters_mod, only: FEParameters
    use MeshParameters_mod, only: MeshParameters
    use BoundaryParameters_mod, only: BoundaryParameters
    use GaussParameters_mod, only: GaussIntegration
    use SolutionVariables_mod, only: SolutionArraysType
    use LinearSystemVariables_mod, only: LinearSystemType
    use ElementVariables_mod, only: GaussPointQuantities, NodeArrays
    use ElementCalculations_mod, only: basis

    implicit none
    
    contains

    subroutine setExtraEquationsBulk(Problem, FE, Mesh, GaussInt, Sol, basis_p, &
        Flag_NR, inex, Elem, LS)

        use ContinuationParameters_mod, only: Continuation_Method
        use PhysicalParameters_mod, only: Total_time
        use MPIParameters_mod, only: Rank
        use ElementCalculations_mod, only: copyToLocal

        implicit none 

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        character(*), intent(in) :: Flag_NR
        PetscInt, intent(in) :: inex
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS
        
        PetscInt :: iel_rank, Nbf, Nunknowns, Nel
        PetscScalar :: value, Res_extra
        PetscErrorCode :: ierr

        Nbf = FE%Nbf
        Nel = FE%Nel_Rank
        loop_elements:do iel_rank = 1, Nel

            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)

            call extraResidualBulk(Problem, Nbf, GaussInt, Sol%name, basis_p, inex, Elem, Res_extra)

            call storeExtraBulk(Problem, Nbf, Mesh%Connectivity_Rank(iel_rank,:), &
                    GaussInt, Sol%name, basis_p, Elem, Res_extra, Flag_NR, inex, LS)

        end do loop_elements

        value = -Problem%iex_b_values(inex)
        if (Continuation_Method == 'Transient') then
            value = value*(1.0d0-exp(-Total_time))
        end if

        Nunknowns = Problem%Nunknowns
        if (Rank == 0) call VecSetValues(LS%b_f, 1, Nunknowns+inex-1, value, ADD_VALUES, ierr)

    end subroutine setExtraEquationsBulk

    !--------------------------------------------------------------------------

    subroutine extraResidualBulk(Problem, Nbf, GaussInt, name, basis_p, inex, Elem, Res_extra)

        use ElementCalculations_mod, only: setElementBasisFunctions, setFlowQuantities

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(GaussIntegration), intent(in) :: GaussInt
        character(*), intent(in) :: name
        procedure(basis), pointer, intent(in) :: basis_p
        PetscInt, intent(in) :: Nbf, inex
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, intent(out) :: Res_extra

        type(GaussPointQuantities) :: GsPt
        PetscInt :: ig

        Res_extra = 0.0d0
        loop_gauss:do ig = 1, GaussInt%NGauss_b

            call setElementBasisFunctions(ig, GaussInt, Elem)

            call basis_p(Problem%idir, Nbf, Elem, GsPt)

            call setFlowQuantities(Problem, Nbf, Elem, GsPt)

            GsPt%WET = GaussInt%Weights_b(ig)*abs(GsPt%detJ)

            select case(name)
            case ('Inflow')
                call extraEquationsResidualsBulkInflow(inex, GsPt, Res_extra)
            case ('Main')
                call extraEquationsResidualsBulkMain(inex, GsPt, Res_extra)
            end select

        end do loop_gauss

    end subroutine extraResidualBulk

    !----------------------------------------------------------------------

    subroutine extraEquationsResidualsBulkMain(inex, GsPt, Res_extra)

        implicit none

        PetscInt, intent(in) :: inex
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscScalar, intent(out) :: Res_extra

        select case(inex)
        case(1)
            !Flow rate
            Res_extra = Res_extra + (GsPt%U_f(1))*(GsPt%WET)
        case default
            write(*,'(a)') "Wrong EQ choice in extraEquationsResidualsBulkMain!"
            stop
        end select

    end subroutine extraEquationsResidualsBulkMain

    !----------------------------------------------------------------------

    subroutine extraEquationsResidualsBulkInflow(inex, GsPt, Res_extra)

        implicit none

        PetscInt, intent(in) :: inex
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscScalar, intent(out) :: Res_extra

        select case(inex)
        case(1)
            !Flow rate
            Res_extra = Res_extra + (GsPt%U_f(1))*(GsPt%WET)
        case default
            write(*,'(a)') "Wrong EQ choice in extraEquationsResidualsBulkInflow!"
            stop
        end select

    end subroutine extraEquationsResidualsBulkInflow

    !---------------------------------------------------------------------

    subroutine storeExtraBulk(Problem, Nbf, Connectivity, GaussInt, name, basis_p, &
        Elem, Res_extra, Flag_NR, inex, LS)

        use Tools_mod, only: perturbVariable, indexInArray, getRows
        use Storage_mod, only: storeResidualGlobal, storeJacobian_Ar, storeJacobian_Ah

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nbf, inex
        PetscInt, dimension(:), intent(in) :: Connectivity
        type(GaussIntegration), intent(in) :: GaussInt
        character(*), intent(in) :: name
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, intent(in) :: Res_extra
        character(*), intent(in) :: Flag_NR
        type(LinearSystemType), intent(inout) :: LS

        PetscInt, dimension(1), parameter :: ieqs = [1]
        PetscInt :: inod, ieq, Nrows, Mcols, irow
        PetscInt, dimension(Nbf) :: grows, gnodes
        PetscReal :: eps_jac
        PetscScalar :: dRes_extra
        PetscScalar, dimension(Nbf, Problem%Neq) :: TEMP_Ar_f
        PetscScalar, dimension(Problem%Nex) :: TEMP_Ah_f
        PetscErrorCode :: ierr

        call MatGetSize(LS%A_f, Nrows, Mcols, ierr)

        gnodes(:) = Connectivity(:)
        grows(:) = getRows(gnodes(:),ieqs(:),Problem%Neq)

        irow = Nrows - Problem%Nex + inex - 1
        call storeResidualGlobal(Res_extra, LS%b_f, irow)

        if (Flag_NR == 'NRP') then

            !Ar_f: derivatives of Extra EQs wrt Bulk unknowns
            TEMP_Ar_f(:,:) = 0.0d0
            do inod = 1, Nbf
                do ieq = 1, Problem%Neq

                    eps_jac = perturbVariable(Elem%TEMP_TL(inod,ieq))
                    Elem%TEMP_TL(inod,ieq) = Elem%TEMP_TL(inod,ieq) + eps_jac

                    call extraResidualBulk(Problem, Nbf, GaussInt, name, basis_p, inex, Elem, dRes_extra)

                    Elem%TEMP_TL(inod,ieq) = Elem%TEMP_TL(inod,ieq) - eps_jac

                    TEMP_Ar_f(inod,ieq) = (dRes_extra - Res_extra)/eps_jac

                end do
            end do

            irow = Nrows - Problem%Nex + inex - 1
            call storeJacobian_Ar(grows, TEMP_Ar_f, LS%A_f, irow)

            !Ah_f: Derivatives of Extra EQs wrt to Extra unknowns
            TEMP_Ah_f(:) = 0.0d0
            do ieq = 1, Problem%Nex

                eps_jac = perturbVariable(Elem%TEMP_EU(ieq))
                Elem%TEMP_EU(ieq) = Elem%TEMP_EU(ieq) + eps_jac

                call extraResidualBulk(Problem, Nbf, GaussInt, name, basis_p, inex, Elem, dRes_extra)

                Elem%TEMP_EU(ieq) = Elem%TEMP_EU(ieq) - eps_jac

                TEMP_Ah_f(ieq) = (dRes_extra - Res_extra)/eps_jac

            end do

            irow = Nrows - Problem%Nex + inex - 1
            call storeJacobian_Ah(TEMP_Ah_f, LS%A_f, inex, irow)

        end if

    end subroutine storeExtraBulk

    !---------------------------------------------------------------------

    subroutine setExtraEquationsEdge(Problem, FE, Mesh, bound, GaussInt, Sol, &
        basis_p, Flag_NR, Elem, LS)

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        character(*), intent(in) :: Flag_NR
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS

        PetscInt :: inex, iex_ed

        do inex = 1, size(bound%ex_ed)
            iex_ed = bound%ex_ed(inex)
            call applyExtraEquationEdge(Problem, FE, Mesh, bound, GaussInt, Sol, &
                basis_p, Flag_NR, iex_ed, Elem, LS)
        end do

    end subroutine setExtraEquationsEdge

    !---------------------------------------------------------------------

    subroutine applyExtraEquationEdge(Problem, FE, Mesh, bound, GaussInt, Sol, &
        basis_p, Flag_NR, iex_ed, Elem, LS)

        use ContinuationParameters_mod, only: Continuation_Method
        use PhysicalParameters_mod, only: Total_time
        use MPIParameters_mod, only: Rank
        use ElementCalculations_mod, only: copyToLocal

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(SolutionArraysType), intent(in) :: Sol
        character(*), intent(in) :: Flag_NR
        PetscInt, intent(in) :: iex_ed
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS

        PetscInt :: iel, iel_rank, ied, Nbf, Nunknowns
        PetscScalar :: value, Res_extra
        PetscErrorCode :: ierr

        Nbf = FE%Nbf
        loop_elements:do iel = 1, size(bound%edge_elements_rank)

            iel_rank = bound%edge_elements_rank(iel)

            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)

            ied = bound%edges(iel)
            call extraResidualEdge(Problem, FE, GaussInt, basis_p, iex_ed, ied, Elem, Res_extra)

            call storeExtraEdge(Problem, FE, Mesh%Connectivity_Rank(iel_rank,:), &
                    GaussInt, basis_p, iex_ed, ied, Elem, Res_extra, Flag_NR, LS)

        end do loop_elements

        value = bound%values_ex_ed(iex_ed)
        ! if (Continuation_Method == 'Transient') then
        !     value = value*(1.0d0-exp(-2.0d0*Total_time))
        ! end if

        Nunknowns = Problem%Nunknowns
        if (Rank == 0) call VecSetValues(LS%b_f, 1, Nunknowns+iex_ed-1, value, ADD_VALUES, ierr)

    end subroutine applyExtraEquationEdge

    !--------------------------------------------------------------------------

    subroutine extraResidualEdge(Problem, FE, GaussInt, basis_p, iex_ed, ied, Elem, Res_extra)

        use ElementCalculations_mod, only: setElementBasisFunctionsEdge, &
            setFlowQuantities, setEdgeQuantities

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        PetscInt, intent(in) :: iex_ed, ied
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, intent(out) :: Res_extra

        type(GaussPointQuantities) :: GsPt
        PetscInt :: ig

        Res_extra = 0.0d0
        loop_gauss:do ig = 1, GaussInt%NGauss_1D

            call setElementBasisFunctionsEdge(ig, ied, GaussInt, Elem)

            call basis_p(Problem%idir, FE%Nbf, Elem, GsPt)

            call setFlowQuantities(Problem, FE%Nbf, Elem, GsPt)

            call setEdgeQuantities(ied, Problem%idir, FE%name, GsPt)

            GsPt%WET = (GaussInt%Weights_1D(ig))*(GsPt%dL)

            call extraEquationsResidualsEdgeMain(iex_ed, GsPt, Res_extra)

        end do loop_gauss

    end subroutine extraResidualEdge

    !----------------------------------------------------------------------

    subroutine extraEquationsResidualsEdgeMain(iex_ed, GsPt, Res_extra)

        implicit none

        PetscInt, intent(in) :: iex_ed
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscScalar, intent(out) :: Res_extra

        PetscReal :: TERM_RES_extra

        TERM_RES_extra = 0.0d0
        select case(iex_ed)
        case(1)
            !Flow rate
            TERM_RES_extra = TERM_RES_extra + (GsPt%U_f(1))*(GsPt%WET)
        case default
            write(*,'(a)') "Wrong EQ choice in extraEquationsResidualsEdgeMain"
            stop
        end select

        Res_extra = Res_extra + TERM_RES_extra

    end subroutine extraEquationsResidualsEdgeMain

    !---------------------------------------------------------------------

    subroutine storeExtraEdge(Problem, FE, Connectivity, GaussInt, basis_p, &
        iex_ed, ied, Elem, Res_extra, Flag_NR, LS)

        use Tools_mod, only: perturbVariable, getRows
        use Storage_mod, only: storeResidualGlobal, storeJacobian_Ar, storeJacobian_Ah

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        PetscInt, intent(in) :: iex_ed, ied
        PetscInt, dimension(:), intent(in) :: Connectivity
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, intent(in) :: Res_extra
        character(*), intent(in) :: Flag_NR
        type(LinearSystemType), intent(inout) :: LS

        PetscInt, dimension(1), parameter :: ieqs = [1]
        PetscInt :: inod, ieq, Nbf, Nrows, Mcols, irow
        PetscInt, dimension(FE%Nbf) :: grows, gnodes
        PetscReal :: eps_jac
        PetscScalar :: dRes_extra
        PetscScalar, dimension(FE%Nbf, Problem%Neq) :: TEMP_Ar_f
        PetscScalar, dimension(Problem%Nex) :: TEMP_Ah_f
        PetscErrorCode :: ierr

        call MatGetSize(LS%A_f, Nrows, Mcols, ierr)

        Nbf = FE%Nbf

        gnodes(:) = Connectivity(:)
        grows(:) = getRows(gnodes(:),ieqs(:),Problem%Neq)

        irow = Nrows - Problem%Nex + iex_ed - 1
        call storeResidualGlobal(Res_extra, LS%b_f, irow)

        if (Flag_NR == 'NRP') then

            !Ar_f: derivatives of Extra EQs wrt Bulk unknowns
            TEMP_Ar_f(:,:) = 0.0d0
            do inod = 1, Nbf
                do ieq = 1, Problem%Neq

                    eps_jac = perturbVariable(Elem%TEMP_TL(inod,ieq))
                    Elem%TEMP_TL(inod,ieq) = Elem%TEMP_TL(inod,ieq) + eps_jac

                    call extraResidualEdge(Problem, FE, GaussInt, basis_p, &
                        iex_ed, ied, Elem, dRes_extra)

                    Elem%TEMP_TL(inod,ieq) = Elem%TEMP_TL(inod,ieq) - eps_jac

                    TEMP_Ar_f(inod,ieq)   = (dRes_extra - Res_extra)/eps_jac

                end do
            end do

            irow = Nrows - Problem%Nex + iex_ed - 1
            call storeJacobian_Ar(grows, TEMP_Ar_f, LS%A_f, irow)

            !Ah_f: Derivatives of Extra EQs wrt to Extra unknowns
            TEMP_Ah_f(:) = 0.0d0
            do ieq = 1, Problem%Nex

                eps_jac = perturbVariable(Elem%TEMP_EU(ieq))
                Elem%TEMP_EU(ieq) = Elem%TEMP_EU(ieq) + eps_jac

                call extraResidualEdge(Problem, FE, GaussInt, basis_p, &
                    iex_ed, ied, Elem, dRes_extra)

                Elem%TEMP_EU(ieq) = Elem%TEMP_EU(ieq) - eps_jac

                TEMP_Ah_f(ieq) = (dRes_extra - Res_extra)/eps_jac

            end do

            irow = Nrows - Problem%Nex + iex_ed - 1
            call storeJacobian_Ah(TEMP_Ah_f, LS%A_f, iex_ed, irow)

        end if

    end subroutine storeExtraEdge

    !---------------------------------------------------------------------

    subroutine setExtraEquationsFace(Problem, FE, Mesh, bound, GaussInt, Sol, &
        basis_p, Flag_NR, Elem, LS)

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        character(*), intent(in) :: Flag_NR
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS

        PetscInt :: inex, iex_fc

        do inex = 1, size(bound%ex_ed)
            iex_fc = bound%ex_ed(inex)
            call applyExtraEquationFace(Problem, FE, Mesh, bound, GaussInt, Sol, &
                basis_p, Flag_NR, iex_fc, Elem, LS)
        end do

    end subroutine setExtraEquationsFace

    !---------------------------------------------------------------------

    subroutine applyExtraEquationFace(Problem, FE, Mesh, bound, GaussInt, Sol, &
        basis_p, Flag_NR, iex_fc, Elem, LS)

        use ContinuationParameters_mod, only: Continuation_Method
        use PhysicalParameters_mod, only: Total_time
        use MPIParameters_mod, only: Rank
        use ElementCalculations_mod, only: copyToLocal

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        character(*), intent(in) :: Flag_NR
        PetscInt, intent(in) :: iex_fc
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS

        PetscInt :: iel, iel_rank, ifc, Nbf, Nunknowns
        PetscScalar :: value, Res_extra
        PetscErrorCode :: ierr

        Nbf = FE%Nbf
        loop_elements:do iel = 1, size(bound%face_elements_rank)

            iel_rank = bound%face_elements_rank(iel)

            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)

            ifc = bound%edges(iel)
            call extraResidualFace(Problem, FE, GaussInt, basis_p, iex_fc, ifc, Elem, Res_extra)

            call storeExtraFace(Problem, FE, Mesh%Connectivity_Rank(iel_rank,:), &
                    GaussInt, basis_p, iex_fc, ifc, Elem, Res_extra, Flag_NR, LS)

        end do loop_elements

        value = bound%values_ex_fc(iex_fc)
        ! if (Continuation_Method == 'Transient') then
        !     value = value*(1.0d0-exp(-2.0d0*Total_time))
        ! end if

        Nunknowns = Problem%Nunknowns
        if (Rank == 0) call VecSetValues(LS%b_f, 1, Nunknowns+iex_fc-1, value, ADD_VALUES, ierr)

    end subroutine applyExtraEquationFace

    !--------------------------------------------------------------------------

    subroutine extraResidualFace(Problem, FE, GaussInt, basis_p, iex_fc, ifc, Elem, Res_extra)

        use ElementCalculations_mod, only: setElementBasisFunctionsFace, &
            setFlowQuantities, setFaceQuantities

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        PetscInt, intent(in) :: iex_fc, ifc
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, intent(out) :: Res_extra

        type(GaussPointQuantities) :: GsPt
        PetscInt :: ig

        Res_extra = 0.0d0
        loop_gauss:do ig = 1, GaussInt%NGauss_2D

            call setElementBasisFunctionsFace(ig, ifc, GaussInt, Elem)

            call basis_p(Problem%idir, FE%Nbf, Elem, GsPt)

            call setFlowQuantities(Problem, FE%Nbf, Elem, GsPt)

            call setFaceQuantities(ifc, FE%name, GsPt)

            GsPt%WET = (GaussInt%Weights_2D(ig))*(GsPt%dL)

            call extraEquationsResidualsFaceMain(iex_fc, GsPt, Res_extra)

        end do loop_gauss

    end subroutine extraResidualFace

    !----------------------------------------------------------------------

    subroutine extraEquationsResidualsFaceMain(iex_fc, GsPt, Res_extra)

        implicit none

        PetscInt, intent(in) :: iex_fc
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscScalar, intent(out) :: Res_extra

        PetscReal :: TERM_RES_extra

        TERM_RES_extra = 0.0d0
        select case(iex_fc)
        case(1)
            !Flow rate
            TERM_RES_extra = TERM_RES_extra + (GsPt%U_f(1))*(GsPt%WET)
        case default
            write(*,'(a)') "Wrong EQ choice in extraEquationsResidualsFaceMain"
            stop
        end select

        Res_extra = Res_extra + TERM_RES_extra

    end subroutine extraEquationsResidualsFaceMain

    !---------------------------------------------------------------------

    subroutine storeExtraFace(Problem, FE, Connectivity, GaussInt, basis_p, &
        iex_fc, ifc, Elem, Res_extra, Flag_NR, LS)

        use Tools_mod, only: perturbVariable, getRows
        use Storage_mod, only: storeResidualGlobal, storeJacobian_Ar, storeJacobian_Ah

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        PetscInt, intent(in) :: iex_fc, ifc
        PetscInt, dimension(:), intent(in) :: Connectivity
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, intent(in) :: Res_extra
        character(*), intent(in) :: Flag_NR
        type(LinearSystemType), intent(inout) :: LS

        PetscInt, dimension(1), parameter :: ieqs = [1]
        PetscInt :: inod, ieq, Nbf, Nrows, Mcols, irow
        PetscInt, dimension(FE%Nbf) :: grows, gnodes
        PetscReal :: eps_jac
        PetscScalar :: dRes_extra
        PetscScalar, dimension(FE%Nbf, Problem%Neq) :: TEMP_Ar_f
        PetscScalar, dimension(Problem%Nex) :: TEMP_Ah_f
        PetscErrorCode :: ierr

        call MatGetSize(LS%A_f, Nrows, Mcols, ierr)

        Nbf = FE%Nbf

        gnodes(:) = Connectivity(:)
        grows(:) = getRows(gnodes(:),ieqs(:),Problem%Neq)

        irow = Nrows - Problem%Nex + iex_fc - 1
        call storeResidualGlobal(Res_extra, LS%b_f, irow)

        if (Flag_NR == 'NRP') then

            !Ar_f: derivatives of Extra EQs wrt Bulk unknowns
            TEMP_Ar_f(:,:) = 0.0d0
            do inod = 1, Nbf
                do ieq = 1, Problem%Neq

                    eps_jac = perturbVariable(Elem%TEMP_TL(inod,ieq))
                    Elem%TEMP_TL(inod,ieq) = Elem%TEMP_TL(inod,ieq) + eps_jac

                    call extraResidualFace(Problem, FE, GaussInt, basis_p, iex_fc, ifc, Elem, dRes_extra)

                    Elem%TEMP_TL(inod,ieq) = Elem%TEMP_TL(inod,ieq) - eps_jac

                    TEMP_Ar_f(inod,ieq)   = (dRes_extra - Res_extra)/eps_jac

                end do
            end do

            irow = Nrows - Problem%Nex + iex_fc - 1
            call storeJacobian_Ar(grows, TEMP_Ar_f, LS%A_f, irow)

            !Ah_f: Derivatives of Extra EQs wrt to Extra unknowns
            TEMP_Ah_f(:) = 0.0d0
            do ieq = 1, Problem%Nex

                eps_jac = perturbVariable(Elem%TEMP_EU(ieq))
                Elem%TEMP_EU(ieq) = Elem%TEMP_EU(ieq) + eps_jac

                call extraResidualFace(Problem, FE, GaussInt, basis_p, iex_fc, ifc, Elem, dRes_extra)

                Elem%TEMP_EU(ieq) = Elem%TEMP_EU(ieq) - eps_jac

                TEMP_Ah_f(ieq) = (dRes_extra - Res_extra)/eps_jac

            end do

            irow = Nrows - Problem%Nex + iex_fc - 1
            call storeJacobian_Ah(TEMP_Ah_f, LS%A_f, iex_fc, irow)

        end if

    end subroutine storeExtraFace

    !----------------------------------------------------------------------

    subroutine arclength(Problem, Nnodes_Total, Sol, Flag_NR, LS)

        use ContinuationParameters_mod, only: LN0, dArN
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Increment

        implicit none 

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nnodes_Total
        type(SolutionArraysType), intent(in) :: Sol
        character(*), intent(in) :: Flag_NR
        type(LinearSystemType), intent(inout) :: LS

        PetscScalar, parameter :: one = 1.0d0
        PetscScalar :: value
        PetscInt, dimension(1) :: irow
        PetscInt :: Nunknowns, Nex, i
        PetscErrorCode :: ierr
        ! PetscInt, dimension(Problem%Nunknowns+Problem%Nex) :: petsc_cols
        ! PetscScalar, dimension(Problem%Nunknowns+Problem%Nex) :: values

        if (Increment <= 2) then

            Nunknowns = Problem%Nunknowns
            Nex = Problem%Nex

            !Continuation EQ: λ-(λ0+(k-1)dλ)
            irow(1) = Nunknowns+Nex-1
            value = Sol%EU(Nex) - (LN0+(Increment-1)*dArN)

            call VecSetValues(LS%b_f, 1, irow(1), value, INSERT_VALUES, ierr)

            if (Flag_NR == 'NRP') then

                call MatAssemblyBegin(LS%A_f, MAT_FINAL_ASSEMBLY, ierr)
                call MatAssemblyEnd(LS%A_f, MAT_FINAL_ASSEMBLY, ierr)
                ! Ac_f(:,Nex)   = 0.0d0   !Derivative of Bulk EQs wrt λ
                ! Ar_f(Nex,:)   = 0.0d0   !Derivative of λ EQ wrt Bulk Unknowns
                ! Ah_f(Nex,Nex) = 1.0d0   !Derivative of λ EQ wrt λ
                call MatZeroRowsColumns(LS%A_f, 1, irow, one, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)

                ! petsc_cols = [ (i, i = 0, size(petsc_cols)-1) ]

                ! ! if (Increment == 1 .and. Iter_f == 1) then
                ! !     call MatSetValues(LS%A_f, size(petsc_cols), petsc_cols, 1, &
                ! !         irow, values, ADD_VALUES, ierr)
                ! !     values(size(values)) = 1.0d0
                ! !     call MatSetValues(LS%A_f, 1, irow, size(petsc_cols), &
                ! !         petsc_cols, values, ADD_VALUES, ierr)
                ! ! else
                !     values(size(values)) = 1.0d0
                !     ! Ac_f(:,Nex)   = 0.0d0   !Derivative of Bulk EQs wrt λ
                !     ! Ah_f(Nex,Nex) = 1.0d0   !Derivative of λ EQ wrt λ
                !     call MatSetValues(LS%A_f, size(petsc_cols), petsc_cols, 1, &
                !         irow, values, INSERT_VALUES, ierr)
                !     ! Ar_f(Nex,:)   = 0.0d0   !Derivative of λ EQ wrt Bulk Unknowns
                !     ! Ah_f(Nex,Nex) = 1.0d0   !Derivative of λ EQ wrt λ
                !     call MatSetValues(LS%A_f, 1, irow, size(petsc_cols), &
                !         petsc_cols, values, INSERT_VALUES, ierr)
                ! ! end if

            end if

            return

        end if

        if (Rank == 0) call pseudoArclength(Problem, Nnodes_Total, Sol, Flag_NR, LS)

    end subroutine arclength

    !--------------------------------------------------------------------------

    subroutine pseudoArclength(Problem, Nnodes_Total, Sol, Flag_NR, LS)

        use ContinuationVariables_mod, only: dSo
        use SolutionVariables_mod, only: dTLodS, dEUodS

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nnodes_Total
        type(SolutionArraysType), intent(in) :: Sol
        character(*), intent(in) :: Flag_NR
        type(LinearSystemType), intent(inout) :: LS

        PetscInt :: i, inod, ieq, irow, Nunknowns, Nex, Neq
        PetscInt, dimension(Problem%Nunknowns+Problem%Nex) :: petsc_cols
        PetscScalar, dimension(Problem%Nunknowns+Problem%Nex) :: values
        PetscScalar :: value
        PetscErrorCode :: ierr

        Nunknowns = Problem%Nunknowns
        Nex = Problem%Nex
        Neq = Problem%Neq

        irow = Nunknowns + Nex - 1

        ! Be_f(Nex) = dot_product(dTLodS,dTL) + dot_product(dEUodS,dEU) - dSo
        value = 0.0d0
        do ieq = 1, Neq
            value = value + dot_product( dTLodS(:,ieq), Sol%TL(:,ieq)-Sol%TLo(:,ieq) )
        end do
        value = value + dot_product( dEUodS(:), Sol%EU(:)-Sol%EUo(:) ) - dSo

        call VecSetValues(LS%b_f, 1, irow, value, INSERT_VALUES, ierr)

        if (Flag_NR == 'NRP') then
            
            petsc_cols = [ (i, i = 0, size(petsc_cols)-1) ]

            !Derivatives of λ EQ wrt Bulk Unknowns
            ! Ar_f(Nex,:) = dTLodS
            i = 0
            do inod = 1, Nnodes_Total
                do ieq = 1, Neq
                    i = i + 1
                    values(i) = dTLodS(inod,ieq)
                end do
            end do
            
            !Derivatives of λ EQ wrt Extra Unknowns
            ! Ah_f(Nex,:) = dEUodS
            values(Nunknowns+1:Nunknowns+Nex) = dEUodS(:)

            call MatSetValues(LS%A_f, 1, irow, size(petsc_cols), petsc_cols, values, INSERT_VALUES, ierr)

        end if

    end subroutine pseudoArclength

end module Extra_Equations_mod
