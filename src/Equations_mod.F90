
#include <petsc/finclude/petscksp.h>

module Equations_mod

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

    subroutine setEquationsBulk(Problem, FE, Mesh, Boundary, GaussInt, Sol, &
        basis_p, Flag_NR, Elem, LS)

        use ElementCalculations_mod, only: copyToLocal

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        character(*), intent(in) :: Flag_NR
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS

        PetscInt :: iel_rank, Nbf, Nel_rank
        PetscScalar, dimension(FE%Nbf, Problem%Neq) :: Res

        Nbf = FE%Nbf
        Nel_rank = FE%Nel_Rank
        loop_elements:do iel_rank = 1, Nel_rank

            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)

            call residualBulk(Problem, Nbf, GaussInt, Sol%name, basis_p, Elem, Res)

            call storeBulk(Problem, Nbf, Mesh%Connectivity_Rank(iel_rank,:), Boundary, &
                GaussInt, Sol%name, basis_p, Elem, Res, Flag_NR, LS)

        end do loop_elements

    end subroutine setEquationsBulk

    !---------------------------------------------------------------------

    subroutine residualBulk(Problem, Nbf, GaussInt, name, basis_p, Elem, Res)

        use ElementCalculations_mod, only: setElementBasisFunctions, setFlowQuantities

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nbf
        type(GaussIntegration), intent(in) :: GaussInt
        character(*), intent(in) :: name
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, dimension(:,:), intent(out) :: Res

        PetscInt :: ig
        type(GaussPointQuantities) :: GsPt

        Res(:,:) = 0.0d0
        loop_gauss:do ig = 1, GaussInt%NGauss_b

            call setElementBasisFunctions(ig, GaussInt, Elem)

            call basis_p(Problem%idir, Nbf, Elem, GsPt)

            call setFlowQuantities(Problem, Nbf, Elem, GsPt)

            GsPt%WET = (GaussInt%Weights_b(ig))*abs(GsPt%detJ)

            select case (name)
            case ('Main')
                call equationsResidualsBulkMain(Problem, Nbf, Elem, GsPt, Res)
            case ('Inflow')
                call equationsResidualsBulkInflow(Problem, Nbf, Elem, GsPt, Res)
            case ('Strong')
                call equationsResidualsBulkStrong(Problem, Nbf, Elem, GsPt, Res)
            end select

        end do loop_gauss

    end subroutine residualBulk

    !----------------------------------------------------------------------

    subroutine equationsResidualsBulkMain(Problem, Nbf, Elem, GsPt, Res)

        use Residuals_mod, only: momentumMain_p, continuityBulkMain, stressMain_p

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nbf
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscScalar, dimension(:,:), intent(out) :: Res

        PetscInt :: inod, ibnd
        PetscBool :: bnd_element
        PetscReal, dimension(Problem%Neq) :: TERM_RES

        loop_nodes:do inod = 1, Nbf

            TERM_RES(:) = 0.0d0

            ! bnd_element = .false.
            ! do ibnd = 1, size(Boundary)
            !     if ( allocated(Boundary(ibnd)%eq_clear) ) then
            !         bnd_element = any( Boundary(ibnd)%gnodes == Mesh%Connectivity(iel,inod) )
            !         if (bnd_element) then
            !             call residualContinuity(inod)
            !             Res(inod,:) = Res(inod,:) + TERM_RES(:)*GsPt%WET
            !             cycle loop_nodes
            !         end if
            !     end if
            ! end do

            call momentumMain_p(Problem%Neq_f, inod, Elem, GsPt, TERM_RES)

            call continuityBulkMain(Problem%Neq_f, inod, Elem, GsPt, TERM_RES)

            call stressMain_p(Problem%Neq_s, inod, Elem, GsPt, TERM_RES)

            Res(inod,:) = Res(inod,:) + TERM_RES(:)*(GsPt%WET)

        end do loop_nodes

    end subroutine equationsResidualsBulkMain

    !----------------------------------------------------------------------

    subroutine equationsResidualsBulkInflow(Problem, Nbf, Elem, GsPt, Res)

        use Residuals_mod, only: momentumInflow, stressesInflow

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nbf
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscScalar, dimension(:,:), intent(out) :: Res

        PetscInt :: inod
        PetscReal, dimension(Problem%Neq) :: TERM_RES

        loop_nodes:do inod = 1, Nbf

            TERM_RES(:) = 0.0d0

            call momentumInflow(Problem%Neq_f, inod, Elem, GsPt, TERM_RES)

            call stressesInflow(Problem%Neq_s, inod, Elem, GsPt, TERM_RES)

            Res(inod,:) = Res(inod,:) + TERM_RES(:)*GsPt%WET

        end do loop_nodes

    end subroutine equationsResidualsBulkInflow

    !---------------------------------------------------------------------

    subroutine equationsResidualsBulkStrong(Problem, Nbf, Elem, GsPt, Res)

        use PhysicalParameters_mod, only: ReN
        use Tools_mod, only: getStressComponent

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nbf
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscScalar, dimension(:,:), intent(out) :: Res

        PetscInt :: inod, Neq_f, Neq_s, i, j, k, istart
        PetscReal :: bifn
        PetscReal, dimension(Problem%Neq) :: TERM_RES

        Neq_f = Problem%Neq_f
        Neq_s = Problem%Neq_s
        loop_nodes:do inod = 1, Nbf

            TERM_RES(:) = 0.0d0

            bifn = Elem%bfn_p(inod)

            do i = 1, Neq_f
                TERM_RES(i) = GsPt%Momentum(i)*bifn
                ! TERM_RES(i) = ReN*(GsPt%Convection(i))*bifn &
                !         + dot_product(GsPt%P_tot(:,i),Elem%dFdXi(:,inod))
            end do

            TERM_RES(Neq_f+1) = (GsPt%Continuity)*bifn

            istart = size(TERM_RES) - Neq_s + 1
            do i = istart, size(TERM_RES)
                call getStressComponent(i,istart,j,k)
                TERM_RES(i) = (GsPt%Con_EQ(j,k))*bifn
            end do

            Res(inod,:) = Res(inod,:) + TERM_RES(:)*(GsPt%WET)

        end do loop_nodes

    end subroutine equationsResidualsBulkStrong

    !----------------------------------------------------------------------

    subroutine storeBulk(Problem, Nbf, Connectivity, Boundary, GaussInt, &
        name, basis_p, Elem, Res, Flag_NR, LS)

        use Tools_mod, only: perturbVariable, indexInArray, getRows
        use Storage_mod, only: storeResidual, storeJacobian, storeJacobianPeriodic, &
            storeJacobian_Ac

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nbf
        PetscInt, dimension(:), intent(in) :: Connectivity
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        character(*), intent(in) :: name
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, dimension(:,:), intent(in) :: Res
        character(*), intent(in) :: Flag_NR
        type(LinearSystemType), intent(inout) :: LS

        PetscInt, dimension(1), parameter :: ieqs = [1]
        PetscInt :: inod, ieq, jnod, inex, Nrows, Mcols, jcol
        PetscInt, dimension(Nbf) :: grows, gnodes
        PetscReal :: eps_jac
        PetscScalar, dimension(Nbf, Problem%Neq) :: dRes
        ! PetscScalar, dimension(Nbf, Nbf, Problem%Neq, Problem%Neq) :: TEMP_A_f
        PetscScalar, dimension(Nbf, Problem%Neq, Nbf, Problem%Neq) :: TEMP_A_f
        PetscScalar, dimension(Nbf, Problem%Neq, Problem%Nex) :: TEMP_Ac_f
        PetscErrorCode :: ierr

        call MatGetSize(LS%A_f, Nrows, Mcols, ierr)

        gnodes(:) = Connectivity(:)
        grows(:) = getRows(gnodes,ieqs,Problem%Neq)

        call storeResidual(Res, grows, LS%b_f)

        if (Flag_NR == 'NRP') then
            !A_f: derivatives of Bulk EQs wrt Bulk unknowns
            TEMP_A_f(:,:,:,:) = 0.0d0 ; dRes(:,:) = 0.0d0
            do inod = 1, Nbf
                do ieq = 1, Problem%Neq

                    eps_jac = perturbVariable(Elem%TEMP_TL(inod,ieq))
                    
                    Elem%TEMP_TL(inod,ieq) = Elem%TEMP_TL(inod,ieq) + eps_jac

                    call residualBulk(Problem, Nbf, GaussInt, name, basis_p, Elem, dRes)

                    Elem%TEMP_TL(inod,ieq) = Elem%TEMP_TL(inod,ieq) - eps_jac

                    do jnod = 1, Nbf
                        ! TEMP_A_f(jnod,inod,ieq,:) = (dRes(jnod,:) - Res(jnod,:))/eps_jac
                        TEMP_A_f(inod,ieq,jnod,:) = (dRes(jnod,:) - Res(jnod,:))/eps_jac
                    end do

                end do
            end do

            call storeJacobian(TEMP_A_f, grows, LS%A_f)

            !Ac_f: Derivatives of Bulk EQs wrt to global unknowns
            TEMP_Ac_f(:,:,:) = 0.0d0 ; dRes(:,:) = 0.0d0
            do inex = 1, Problem%Nex

                eps_jac = perturbVariable(Elem%TEMP_EU(inex))
                Elem%TEMP_EU(inex) = Elem%TEMP_EU(inex) + eps_jac

                call residualBulk(Problem, Nbf, GaussInt, name, basis_p, Elem, dRes)

                Elem%TEMP_EU(inex) = Elem%TEMP_EU(inex) - eps_jac

                do inod = 1, Nbf
                    TEMP_Ac_f(inod,:,inex) = (dRes(inod,:) - Res(inod,:))/eps_jac
                end do

                jcol = Mcols - Problem%Nex + inex -1
                call storeJacobian_Ac(TEMP_Ac_f(:,:,inex), grows, LS%A_f, jcol)

            end do

        end if

        call storePeriodic(Problem, Nbf, Connectivity, Boundary, Res, Flag_NR, &
            TEMP_A_f, TEMP_Ac_f, LS)

    end subroutine storeBulk

    !---------------------------------------------------------------------

    subroutine setEquationsEdge(Problem, FE, Mesh, bound, Boundary, GaussInt, &
        Sol, basis_p, Flag_NR, Elem, LS)

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        character(*), intent(in) :: Flag_NR
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS

        PetscInt :: ieq, ieq_ed

        do ieq = 1, size(bound%eq_ed)
            ieq_ed = bound%eq_ed(ieq)
            call applyEquationEdge(Problem, FE, Mesh, bound, Boundary, GaussInt, &
                Sol, basis_p, Flag_NR, ieq_ed, Elem, LS)
        end do

    end subroutine setEquationsEdge

    !---------------------------------------------------------------------

    subroutine applyEquationEdge(Problem, FE, Mesh, bound, Boundary, GaussInt, Sol, &
        basis_p, Flag_NR, ieq_ed, Elem, LS)

        use ElementCalculations_mod, only: copyToLocal

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        character(*), intent(in) :: Flag_NR
        PetscInt, intent(in) :: ieq_ed
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS

        PetscInt :: iel, iel_rank, ied, Nbf
        PetscScalar, dimension(FE%Nbf, Problem%Neq) :: Res

        Nbf = FE%Nbf
        loop_elements:do iel = 1, size(bound%edge_elements_rank)

            iel_rank = bound%edge_elements_rank(iel)
            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)

            ied = bound%edges(iel)
            call residualEdge(Problem, FE, GaussInt, basis_p, ieq_ed, ied, Elem, Res)

            call storeEdge(Problem, FE, Mesh%Connectivity_Rank(iel_rank,:), &
                    Boundary, GaussInt, basis_p, ieq_ed, ied, Elem, Res, Flag_NR, LS)

        end do loop_elements

    end subroutine applyEquationEdge

    !---------------------------------------------------------------------

    subroutine residualEdge(Problem, FE, GaussInt, basis_p, ieq_ed, ied, Elem, Res)

        use ElementCalculations_mod, only: setElementBasisFunctionsEdge, &
            setFlowQuantities, setEdgeQuantities

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        PetscInt, intent(in) :: ieq_ed, ied
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, dimension(:,:), intent(out) :: Res

        PetscInt :: ig, Nbf
        type(GaussPointQuantities) :: GsPt

        Nbf = FE%Nbf

        Res(:,:) = 0.0d0
        loop_gauss:do ig = 1, GaussInt%NGauss_1D

            call setElementBasisFunctionsEdge(ig, ied, GaussInt, Elem)

            call basis_p(Problem%idir, Nbf, Elem, GsPt)

            call setFlowQuantities(Problem, Nbf, Elem, GsPt)

            call setEdgeQuantities(ied, Problem%idir, FE%name, GsPt)

            GsPt%WET = (GaussInt%Weights_1D(ig))*(GsPt%dL)

            call equationsResidualsEdgeMain(Problem, Nbf, ieq_ed, Elem, GsPt, Res)

        end do loop_gauss

    end subroutine residualEdge

    !----------------------------------------------------------------------

    subroutine equationsResidualsEdgeMain(Problem, Nbf, ieq_ed, Elem, GsPt, Res)

        use Residuals_mod, only: openBC, slipBC

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nbf, ieq_ed
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscScalar, dimension(:,:), intent(out) :: Res

        PetscInt :: inod
        PetscReal, dimension(Problem%Neq) :: TERM_RES

        loop_nodes:do inod = 1, Nbf

            TERM_RES(:) = 0.0d0

            select case (ieq_ed)
            case (1)
                call openBC(Problem%Neq_f, inod, Elem, GsPt, TERM_RES)
            case (2)
                call slipBC(inod, Elem, GsPt, TERM_RES)
            case default
                write(*,'(a)') "Wrong EQ choice in equationsResidualsEdgeMain!"
                stop
            end select

            Res(inod,:) = Res(inod,:) + TERM_RES(:)*(GsPt%WET)

        end do loop_nodes

    end subroutine equationsResidualsEdgeMain

    !---------------------------------------------------------------------

    subroutine storeEdge(Problem, FE, Connectivity, Boundary, GaussInt, basis_p, &
        ieq_ed, ied, Elem, Res, Flag_NR, LS)

        use Tools_mod, only: perturbVariable, indexInArray, getRows
        use Storage_mod, only: storeResidual, storeJacobian, storeJacobianPeriodic, &
            storeJacobian_Ac

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:), intent(in) :: Connectivity
        PetscInt, intent(in) :: ieq_ed, ied
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, dimension(:,:), intent(in) :: Res
        character(*), intent(in) :: Flag_NR
        type(LinearSystemType), intent(inout) :: LS

        PetscInt, dimension(1), parameter :: ieqs = [1]
        PetscInt :: inod, ieq, jnod, inex, Nbf, Nrows, Mcols, jcol
        PetscInt, dimension(FE%Nbf) :: grows, gnodes
        PetscReal :: eps_jac
        PetscScalar, dimension(FE%Nbf, Problem%Neq) :: dRes
        PetscScalar, dimension(FE%Nbf, Problem%Neq, FE%Nbf, Problem%Neq) :: TEMP_A_f
        PetscScalar, dimension(FE%Nbf, Problem%Neq, Problem%Nex) :: TEMP_Ac_f
        PetscErrorCode :: ierr

        call MatGetSize(LS%A_f, Nrows, Mcols, ierr)

        Nbf = FE%Nbf

        gnodes(:) = Connectivity(:)
        grows(:) = getRows(gnodes(:),ieqs(:),Problem%Neq)

        call storeResidual(Res, grows, LS%b_f)

        if (Flag_NR == 'NRP') then
            !A_f: derivatives of Bulk EQs wrt Bulk unknowns
            TEMP_A_f(:,:,:,:) = 0.0d0 ; dRes(:,:) = 0.0d0
            do inod = 1, Nbf
                do ieq = 1, Problem%Neq

                    eps_jac = perturbVariable(Elem%TEMP_TL(inod,ieq))

                    Elem%TEMP_TL(inod,ieq) = Elem%TEMP_TL(inod,ieq) + eps_jac

                    call residualEdge(Problem, FE, GaussInt, basis_p, ieq_ed, ied, Elem, dRes)

                    Elem%TEMP_TL(inod,ieq) = Elem%TEMP_TL(inod,ieq) - eps_jac

                    do jnod = 1, Nbf
                        ! TEMP_A_f(jnod,inod,ieq,:) = (dRes(jnod,:) - Res(jnod,:))/eps_jac
                        TEMP_A_f(inod,ieq,jnod,:) = (dRes(jnod,:) - Res(jnod,:))/eps_jac
                    end do

                end do
            end do

            call storeJacobian(TEMP_A_f, grows, LS%A_f)

            !Ac_f: Derivatives of Bulk EQs wrt to global unknowns
            TEMP_Ac_f(:,:,:) = 0.0d0 ; dRes(:,:) = 0.0d0
            do inex = 1, Problem%Nex

                eps_jac = perturbVariable(Elem%TEMP_EU(inex))
                
                Elem%TEMP_EU(inex) = Elem%TEMP_EU(inex) + eps_jac

                call residualEdge(Problem, FE, GaussInt, basis_p, ieq_ed, ied, Elem, dRes)

                Elem%TEMP_EU(inex) = Elem%TEMP_EU(inex) - eps_jac

                do inod = 1, Nbf
                    TEMP_Ac_f(inod,:,inex) = (dRes(inod,:) - Res(inod,:))/eps_jac
                end do

                jcol = Mcols - Problem%Nex + inex -1
                call storeJacobian_Ac(TEMP_Ac_f(:,:,inex), grows, LS%A_f, jcol)

            end do

        end if

        call storePeriodic(Problem, Nbf, Connectivity, Boundary, Res, Flag_NR, &
            TEMP_A_f, TEMP_Ac_f, LS)

    end subroutine storeEdge

    !---------------------------------------------------------------------

    subroutine setEquationsFace(Problem, FE, Mesh, bound, Boundary, GaussInt, Sol, &
        basis_p, Flag_NR, Elem, LS)

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        character(*), intent(in) :: Flag_NR
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS

        PetscInt :: ieq, ieq_fc

        do ieq = 1, size(bound%eq_fc)
            ieq_fc = bound%eq_fc(ieq)
            call applyEquationFace(Problem, FE, Mesh, bound, Boundary, GaussInt, Sol, &
                basis_p, Flag_NR, ieq_fc, Elem, LS)
        end do

    end subroutine setEquationsFace

    !---------------------------------------------------------------------

    subroutine applyEquationFace(Problem, FE, Mesh, bound, Boundary, GaussInt, Sol, &
        basis_p, Flag_NR, ieq_fc, Elem, LS)

        use ElementCalculations_mod, only: copyToLocal

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), intent(in) :: bound
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        character(*), intent(in) :: Flag_NR
        PetscInt, intent(in) :: ieq_fc
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS

        PetscInt :: iel, iel_rank, ifc, Nbf
        PetscScalar, dimension(FE%Nbf, Problem%Neq) :: Res

        Nbf = FE%Nbf
        loop_elements:do iel = 1, size(bound%face_elements_rank)

            iel_rank = bound%face_elements_rank(iel)
            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)

            ifc = bound%faces(iel)
            call residualFace(Problem, FE, GaussInt, basis_p, ieq_fc, ifc, Elem, Res)

            call storeFace(Problem, FE, Mesh%Connectivity_Rank(iel_rank,:), &
                    Boundary, GaussInt, basis_p, ieq_fc, ifc, Elem, Res, Flag_NR, LS)

        end do loop_elements

    end subroutine applyEquationFace

    !---------------------------------------------------------------------

    subroutine residualFace(Problem, FE, GaussInt, basis_p, ieq_fc, ifc, Elem, Res)

        use ElementCalculations_mod, only: setElementBasisFunctionsFace, &
            setFlowQuantities, setFaceQuantities

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        PetscInt, intent(in) :: ieq_fc, ifc
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, dimension(:,:), intent(out) :: Res

        type(GaussPointQuantities) :: GsPt
        PetscInt :: ig, Nbf

        Nbf = FE%Nbf

        Res(:,:) = 0.0d0
        loop_gauss:do ig = 1, GaussInt%NGauss_2D

            call setElementBasisFunctionsFace(ig, ifc, GaussInt, Elem)

            call basis_p(Problem%idir, Nbf, Elem, GsPt)

            call setFlowQuantities(Problem, Nbf, Elem, GsPt)

            call setFaceQuantities(ifc, FE%name, GsPt)

            GsPt%WET = GaussInt%Weights_2D(ig)*(GsPt%dL)

            ! if ( any(GsPt%normal < -1.0d-6) ) write(*,'(i3,3es12.4)') ifc, GsPt%normal
            
            call equationsResidualsFaceMain(Problem, Nbf, ieq_fc, Elem, GsPt, Res)
            
        end do loop_gauss

    end subroutine residualFace

    !----------------------------------------------------------------------

    subroutine equationsResidualsFaceMain(Problem, Nbf, ieq_fc, Elem, GsPt, Res)

        use Residuals_mod, only: openBC, slipBC

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nbf, ieq_fc
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscScalar, dimension(:,:), intent(out) :: Res

        PetscInt :: inod
        PetscReal, dimension(Problem%Neq) :: TERM_RES

        loop_nodes:do inod = 1, Nbf

            TERM_RES(:) = 0.0d0

            select case (ieq_fc)
            case (1)
                call openBC(Problem%Neq_f, inod, Elem, GsPt, TERM_RES)
            case (2)
                call slipBC(inod, Elem, GsPt, TERM_RES)
            case default
                write(*,'(a)') "Wrong EQ choice in equationsResidualsFaceMain!"
                stop
            end select

            Res(inod,:) = Res(inod,:) + TERM_RES(:)*GsPt%WET

        end do loop_nodes

    end subroutine equationsResidualsFaceMain

    !---------------------------------------------------------------------

    subroutine storeFace(Problem, FE, Connectivity, Boundary, GaussInt, basis_p, &
        ieq_fc, ifc, Elem, Res, Flag_NR, LS)

        use Tools_mod, only: perturbVariable, indexInArray, getRows
        use Storage_mod, only: storeResidual, storeJacobian, storeJacobianPeriodic, &
            storeJacobian_Ac

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:), intent(in) :: Connectivity
        PetscInt, intent(in) :: ieq_fc, ifc
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, dimension(:,:), intent(in) :: Res
        character(*), intent(in) :: Flag_NR
        type(LinearSystemType), intent(inout) :: LS

        PetscInt, dimension(1), parameter :: ieqs = [1]
        PetscInt :: inod, ieq, jnod, inex, Nbf, Nrows, Mcols, jcol
        PetscInt, dimension(FE%Nbf) :: grows, gnodes
        PetscReal :: eps_jac
        PetscScalar, dimension(FE%Nbf, Problem%Neq) :: dRes
        PetscScalar, dimension(FE%Nbf, Problem%Neq, FE%Nbf, Problem%Neq) :: TEMP_A_f
        PetscScalar, dimension(FE%Nbf, Problem%Neq, Problem%Nex) :: TEMP_Ac_f
        PetscErrorCode :: ierr

        call MatGetSize(LS%A_f, Nrows, Mcols, ierr)

        Nbf = FE%Nbf

        gnodes(:) = Connectivity(:)
        grows(:) = getRows(gnodes(:),ieqs(:),Problem%Neq)

        call storeResidual(Res, grows, LS%b_f)

        if (Flag_NR == 'NRP') then
            !A_f: derivatives of Bulk EQs wrt Bulk unknowns
            TEMP_A_f(:,:,:,:) = 0.0d0 ; dRes(:,:) = 0.0d0
            do inod = 1, Nbf
                do ieq = 1, Problem%Neq

                    eps_jac = perturbVariable(Elem%TEMP_TL(inod,ieq))
                    Elem%TEMP_TL(inod,ieq) = Elem%TEMP_TL(inod,ieq) + eps_jac

                    call residualFace(Problem, FE, GaussInt, basis_p, ieq_fc, ifc, Elem, dRes)

                    Elem%TEMP_TL(inod,ieq) = Elem%TEMP_TL(inod,ieq) - eps_jac

                    do jnod = 1, Nbf
                        ! TEMP_A_f(jnod,inod,ieq,:) = (dRes(jnod,:) - Res(jnod,:))/eps_jac
                        TEMP_A_f(inod,ieq,jnod,:) = (dRes(jnod,:) - Res(jnod,:))/eps_jac
                    end do

                end do
            end do

            call storeJacobian(TEMP_A_f, grows, LS%A_f)

            !Ac_f: Derivatives of Bulk EQs wrt to global unknowns
            TEMP_Ac_f(:,:,:) = 0.0d0 ; dRes(:,:) = 0.0d0
            do inex = 1, Problem%Nex

                eps_jac = perturbVariable(Elem%TEMP_EU(inex))
                Elem%TEMP_EU(inex) = Elem%TEMP_EU(inex) + eps_jac

                call residualFace(Problem, FE, GaussInt, basis_p, ieq_fc, ifc, Elem, dRes)

                Elem%TEMP_EU(inex) = Elem%TEMP_EU(inex) - eps_jac

                do inod = 1, Nbf
                    TEMP_Ac_f(inod,:,inex) = (dRes(inod,:) - Res(inod,:))/eps_jac
                end do

                jcol = Mcols - Problem%Nex + inex -1
                call storeJacobian_Ac(TEMP_Ac_f(:,:,inex), grows, LS%A_f, jcol)

            end do

        end if

        call storePeriodic(Problem, Nbf, Connectivity, Boundary, Res, Flag_NR, &
            TEMP_A_f, TEMP_Ac_f, LS)

    end subroutine storeFace

    !---------------------------------------------------------------------

    subroutine storePeriodic(Problem, Nbf, Connectivity, Boundary, Res, Flag_NR, &
        TEMP_A_f, TEMP_Ac_f, LS)

        use Tools_mod, only: indexInArray, getRows
        use Storage_mod, only: storeResidual, storeJacobianPeriodic, storeJacobian_Ac

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nbf
        PetscInt, dimension(:), intent(in) :: Connectivity
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary
        PetscScalar, dimension(:,:), intent(in) :: Res
        character(*), intent(in) :: Flag_NR
        PetscScalar, dimension(:,:,:,:), intent(in) :: TEMP_A_f
        PetscScalar, dimension(:,:,:), intent(in) :: TEMP_Ac_f
        type(LinearSystemType), intent(inout) :: LS

        PetscInt, dimension(1), parameter :: ieqs = [1]
        PetscInt :: jnod, inex, Nrows, Mcols, jcol, i, ibnd, gnod, gnod_index
        PetscInt, dimension(Nbf) :: grows, gnodes
        PetscInt, dimension(1) :: grows_p, gnod_p
        PetscErrorCode :: ierr

        call MatGetSize(LS%A_f, Nrows, Mcols, ierr)

        gnodes(:) = Connectivity(:)
        grows(:) = getRows(gnodes,ieqs,Problem%Neq)

        do i = 1, size(Problem%PeriodicBNDs,1)

            ibnd = Problem%PeriodicBNDs(i,1)

            do jnod = 1, size(gnodes)

                gnod = gnodes(jnod)
                gnod_index = indexInArray(Boundary(ibnd)%gnodes_total, gnod)

                if (gnod_index /= 0) then

                    gnod_p(1) = Boundary(ibnd)%gnodes_periodic(gnod_index)
                    grows_p(:) = getRows(gnod_p,ieqs,Problem%Neq)

                    call storeResidual(Res(jnod:jnod,:), grows_p, LS%b_f)

                    if (Flag_NR == 'NRP') then
                        call storeJacobianPeriodic(TEMP_A_f(:,:,jnod:jnod,:), &
                            grows, grows_p, LS%A_f)
                        do inex = 1, Problem%Nex
                            jcol = Mcols - Problem%Nex + inex -1
                            call storeJacobian_Ac(TEMP_Ac_f(jnod:jnod,:,inex), &
                                grows_p, LS%A_f, jcol)
                        end do
                    end if
                    
                end if

            end do

        end do

    end subroutine storePeriodic

end module Equations_mod
