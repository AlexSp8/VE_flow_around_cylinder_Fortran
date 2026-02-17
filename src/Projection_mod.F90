
#include <petsc/finclude/petscksp.h>

module Projection_mod

    use petscksp
    use PhysicalParameters_mod, only: ProblemParameters
    use FEMParameters_mod, only: FEParameters
    use MeshParameters_mod, only: MeshParameters
    use GaussParameters_mod, only: GaussIntegration
    use SolutionVariables_mod, only: SolutionArraysType
    use LinearSystemVariables_mod, only: LinearSystemType
    use NewtonRaphsonVariables_mod, only: NewtonType
    use ElementVariables_mod, only: GaussPointQuantities, NodeArrays

    implicit none

    contains

    subroutine loopNewtonRaphsonProjection(Problem, FE, Mesh, GaussInt, Sol, Elem, LS)

        use omp_lib, only: omp_get_wtime
        use Tools_mod, only: writeElapsedTime
        use MPIParameters_mod, only: Rank
        use NewtonRaphsonParameters_mod, only: Cor_Norm_eps
        use Petsc_mod, only: initializeLSObjects, solveLinearSystem
        use NR_mod, only: initializeNewtonVariables, checkNewtonConvergence, &
            emergencyNewtonRaphson
        ! use PostProcess_mod, only: writeTecplotFile

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(inout) :: Sol
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS

        type(NewtonType) :: NR
        PetscReal :: tStart, tEnd, tStart2, tStart3
        PetscBool :: check, output, emergency
        PetscErrorCode :: ierr

        tStart = omp_get_wtime()

        NR%name = 'Projection'

        output = (Rank == 0) !.and. (Sol%name == 'Main')

        if (output) then
            write(*, '(a)') '======================================================================'
            write(*, '(25x,a)') 'Projection '//Sol%name
            write(*, '(a)') '======================================================================'
        end if

        NR%xF = 1.0d0
        call initializeNewtonVariables(NR)

        LOOP_NewtonRaphsonProjection:do

            check = (NR%Cor_Norm_New_f < 1.0d-12)!Cor_Norm_eps)
            check = check .and. (NR%Res_Norm < 1.0d-12)!Cor_Norm_eps/100.0d0)
            if (check) then
                exit
            end if

            tStart2 = omp_get_wtime()

            NR%Iter_f = NR%Iter_f + 1

            ! NR%Flag_NR = 'NRP'      !Always Full NR

            call initializeLSObjects(NR%Flag_NR, LS)

            call assemblePetscObjectsProjection(output, Problem, FE, Mesh, &
                GaussInt, Sol, NR%Flag_NR, Elem, LS)

            call VecNormBegin(LS%b_f, NORM_2, NR%Res_Norm, ierr)
            call VecNormEnd(LS%b_f, NORM_2, NR%Res_Norm, ierr)

            tStart3 = omp_get_wtime()
            call solveLinearSystem(NR, LS)
            tEnd = omp_get_wtime()
            if (output) write(*,'(a,f8.3,a)') 'LS Solve time: ', tEnd - tStart3, ' s'

            call VecNormBegin(LS%s_f, NORM_2, NR%Cor_Norm_New_f, ierr)
            call VecNormEnd(LS%s_f, NORM_2, NR%Cor_Norm_New_f, ierr)

            call checkNewtonConvergence(output, NR, emergency)

            if (emergency) then
                call emergencyNewtonRaphson(output, NR, Sol)
                cycle LOOP_NewtonRaphsonProjection
            end if

            call updateNewtonProjection(Problem%Neq_proj, FE, Mesh%Connectivity_Rank, &
                NR%xF, LS%s_f, Sol)

            if (output) then
                ! call writeTecplotFile(Problem, FE, Mesh, Sol%TL, Sol%TL_proj, &
                !     fn, real(NR%Iter_f,8))
                ! pause
            end if

            tEnd = omp_get_wtime()
            if (output) call writeElapsedTime(tStart2, tEnd, 'Newton time')

        end do LOOP_NewtonRaphsonProjection

        call MPI_Barrier(MPI_COMM_WORLD,ierr)

        tEnd = omp_get_wtime()
        if (output) call writeElapsedTime(tStart, tEnd, 'Solution time')

    end subroutine loopNewtonRaphsonProjection

    !--------------------------------------------------------------------------

    subroutine assemblePetscObjectsProjection(output, Problem, FE, Mesh, &
        GaussInt, Sol, Flag_NR, Elem, LS)

        use omp_lib, only: omp_get_wtime
        use ElementCalculations_mod, only: basis, basisInflow_p, basisMain_p
        use Petsc_mod, only: finalAssembly

        implicit none
    
        PetscBool, intent(in) :: output
        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        character(*), intent(in) :: Flag_NR
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS

        PetscReal :: tStart, tEnd
        procedure(basis), pointer :: basis_p

        tStart = omp_get_wtime()

        if (Sol%name == 'Main') then
            basis_p => basisMain_p
        end if

        if (Sol%name == 'Inflow') then
            basis_p => basisInflow_p
        end if

        call setEquationsProjection(Problem, FE, Mesh, GaussInt, Sol, &
            basis_p, Flag_NR, Elem, LS)

        tEnd = omp_get_wtime()
        if (output) write(*,'(a,f8.3,a)', advance='no') 'Projection EQs time: ', tEnd - tStart, ' s'

        call finalAssembly(Flag_NR, LS)

        tEnd = omp_get_wtime()
        if (output) write(*, '(a,f8.3,a)') ' / Assembly time: ', tEnd - tStart, ' s'

    end subroutine assemblePetscObjectsProjection

    !---------------------------------------------------------------------

    subroutine setEquationsProjection(Problem, FE, Mesh, GaussInt, Sol, &
        basis_p, Flag_NR, Elem, LS)

        use ElementCalculations_mod, only: copyToLocal, basis

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(GaussIntegration), intent(in) :: GaussInt
        type(SolutionArraysType), intent(in) :: Sol
        procedure(basis), pointer, intent(in) :: basis_p
        character(*), intent(in) :: Flag_NR
        type(NodeArrays), intent(inout) :: Elem
        type(LinearSystemType), intent(inout) :: LS

        PetscInt :: iel_rank, Nbf, Nel_rank
        PetscScalar, dimension(FE%Nbf, Problem%Neq_proj) :: Res

        Nbf = FE%Nbf
        Nel_rank = FE%Nel_Rank
        loop_elements:do iel_rank = 1, Nel_rank

            call copyToLocal(Nbf, Mesh%Xi_Mesh_Rank, Sol, iel_rank, Elem)

            call residualProjection(Problem, Nbf, GaussInt, basis_p, Elem, Res)

            call storeProjection(Problem, Nbf, Mesh%Connectivity_Rank, &
                GaussInt, basis_p, iel_rank, Elem, Res, Flag_NR, LS)

        end do loop_elements

    end subroutine setEquationsProjection

    !---------------------------------------------------------------------

    subroutine residualProjection(Problem, Nbf, GaussInt, basis_p, Elem, Res)

        use ElementCalculations_mod, only: setElementBasisFunctions, &
            setFlowUnknowns, basis

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nbf
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, dimension(:,:), intent(out) :: Res

        PetscInt :: ig
        type(GaussPointQuantities) :: GsPt

        Res(:,:) = 0.0d0
        loop_gauss:do ig = 1, GaussInt%NGauss_b

            call setElementBasisFunctions(ig, GaussInt, Elem)

            call basis_p(Problem%idir, Nbf, Elem, GsPt)

            call setFlowUnknowns(Problem, Elem, GsPt)

            GsPt%WET = (GaussInt%Weights_b(ig))*abs(GsPt%detJ)

            call equationsResidualsProjection(Elem%bfn_p, GsPt, Res)

        end do loop_gauss

    end subroutine residualProjection

    !----------------------------------------------------------------------

    subroutine equationsResidualsProjection(bfn, GsPt, Res)

        implicit none

        PetscReal, dimension(:), intent(in) :: bfn
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscScalar, dimension(:,:), intent(out) :: Res

        PetscInt :: inod, Nbf
        PetscReal, dimension(size(Res,2)) :: TERM_RES

        Nbf = size(Res,1)
        loop_nodes:do inod = 1, Nbf

            TERM_RES(:) = 0.0d0

            call gradUProjection(bfn(inod), GsPt, TERM_RES)

            Res(inod,:) = Res(inod,:) + TERM_RES(:)*(GsPt%WET)

        end do loop_nodes

    end subroutine equationsResidualsProjection

    !----------------------------------------------------------------------

    subroutine gradUProjection(bifn, GsPt, TERM_RES)

        use Tools_mod, only: getTensorComponent
        
        implicit none

        PetscReal, intent(in) :: bifn
        type(GaussPointQuantities), intent(in) :: GsPt
        PetscReal, dimension(:), intent(inout) :: TERM_RES

        PetscInt :: i, j, k

        do i = 1, size(TERM_RES)
            call getTensorComponent(i,j,k)
            TERM_RES(i) = (GsPt%GU_proj(j,k) - GsPt%GU(j,k))*bifn
        end do

    end subroutine gradUProjection

    !----------------------------------------------------------------------

    subroutine storeProjection(Problem, Nbf, Connectivity_Rank, GaussInt, &
        basis_p, iel_rank, Elem, Res, Flag_NR, LS)

        use ElementCalculations_mod, only: basis
        use Tools_mod, only: perturbVariable, getRows
        use Storage_mod, only: storeResidual, storeJacobian, storeJacobian_Ac

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nbf, iel_rank
        PetscInt, dimension(:,:), intent(in) :: Connectivity_Rank
        type(GaussIntegration), intent(in) :: GaussInt
        procedure(basis), pointer, intent(in) :: basis_p
        type(NodeArrays), intent(inout) :: Elem
        PetscScalar, dimension(:,:), intent(in) :: Res
        character(*), intent(in) :: Flag_NR
        type(LinearSystemType), intent(inout) :: LS

        PetscInt, dimension(1), parameter :: ieqs = [1]
        PetscInt :: inod, ieq, jnod
        PetscInt, dimension(Nbf) :: grows, gnodes
        PetscReal :: eps_jac
        PetscScalar, dimension(Nbf, Problem%Neq_proj) :: dRes
        PetscScalar, dimension(Nbf, Problem%Neq_proj, Nbf, Problem%Neq_proj) :: TEMP_A_f

        gnodes(:) = Connectivity_Rank(iel_rank,:)
        grows(:) = getRows(gnodes(:),ieqs(:),Problem%Neq_proj)

        call storeResidual(Res, grows, LS%b_f)

        if (Flag_NR == 'NRP') then
            !A_f: derivatives of Bulk EQs wrt Bulk unknowns
            TEMP_A_f(:,:,:,:) = 0.0d0 ; dRes(:,:) = 0.0d0
            do inod = 1, Nbf
                do ieq = 1, Problem%Neq_proj

                    eps_jac = perturbVariable(Elem%TEMP_TL_proj(inod,ieq))
                    Elem%TEMP_TL_proj(inod,ieq) = Elem%TEMP_TL_proj(inod,ieq) + eps_jac

                    call residualProjection(Problem, Nbf, GaussInt, basis_p, Elem, dRes)

                    Elem%TEMP_TL_proj(inod,ieq) = Elem%TEMP_TL_proj(inod,ieq) - eps_jac

                    do jnod = 1, Nbf
                        ! TEMP_A_f(jnod,inod,ieq,:) = (dRes(jnod,:) - Res(jnod,:))/eps_jac
                        TEMP_A_f(inod,ieq,jnod,:) = (dRes(jnod,:) - Res(jnod,:))/eps_jac
                    end do

                end do
            end do

            call storeJacobian(TEMP_A_f, grows, LS%A_f)

        end if

    end subroutine storeProjection

    !--------------------------------------------------------------------------

    subroutine updateNewtonProjection(Neq_proj, FE, Connectivity_Rank, xF, s_f_proj, Sol)

        use Tools_mod, only: getRows, getRankNode
        use MPIParameters_mod, only: Rank

        implicit none

        PetscInt, intent(in) :: Neq_proj
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:,:), intent(in) :: Connectivity_Rank
        PetscReal, intent(in) :: xF
        Vec, intent(in) :: s_f_proj
        type(SolutionArraysType), intent(inout) :: Sol

        PetscInt :: inod, ieq, irow, iel_rank, rnod
        PetscInt, dimension(Neq_proj) :: ieqs
        PetscInt, dimension(Neq_proj) :: irows_eqs
        PetscInt, dimension(1) :: gnod
        Vec :: s_f
        PetscScalar, pointer :: s_f_pointer(:)
        VecScatter :: ctx
        PetscErrorCode :: ierr

        call VecScatterCreateToAll(s_f_proj, ctx, s_f, ierr)
        !INSERT_VALUES: any location non-scattered retains its old value
        !SCATTER_fORWARD: scatter from s_f_proj to s_f
        call VecScatterBegin(ctx, s_f_proj, s_f, INSERT_VALUES, SCATTER_FORWARD, ierr)
        call VecScatterEnd(ctx, s_f_proj, s_f, INSERT_VALUES, SCATTER_FORWARD, ierr)

        !Here, we created the pointer to access the values of the solution vector s_f
        call VecGetArrayReadF90(s_f, s_f_pointer, ierr)

        if (Rank == 0) then
            irow = 0
            do inod = 1, FE%Nnodes
                do ieq = 1, Neq_proj
                    irow = irow+1
                    Sol%TL_proj(inod,ieq) = Sol%TL_proj(inod,ieq) - xF*s_f_pointer(irow)
                end do
            end do
        end if

        ieqs = [ (ieq, ieq = 1, Neq_proj) ]
        do iel_rank = 1, FE%Nel_Rank
            do inod = 1, FE%Nbf
                gnod(1) = Connectivity_Rank(iel_rank,inod)
                rnod = getRankNode(iel_rank,inod,FE%Nbf)
                irows_eqs(:) = getRows(gnod(:),ieqs(:),Neq_proj)
                Sol%TL_proj_Rank(rnod,:) = Sol%TL_proj_Rank(rnod,:) - xF*s_f_pointer(irows_eqs(:))
            end do
        end do

        call VecRestoreArrayReadF90(s_f, s_f_pointer, ierr)

        call VecScatterDestroy(ctx, ierr)
        call VecDestroy(s_f, ierr)

    end subroutine updateNewtonProjection

end module Projection_mod
