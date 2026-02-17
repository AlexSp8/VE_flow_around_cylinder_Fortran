
#include <petsc/finclude/petscksp.h>

module Dirichlet_mod

    use petscksp
    use PhysicalParameters_mod, only: ProblemParameters
    use FEMParameters_mod, only: FEParameters
    use MeshParameters_mod, only: MeshParameters
    use BoundaryParameters_mod, only: BoundaryParameters
    use SolutionVariables_mod, only: SolutionArraysType
    use LinearSystemVariables_mod, only: LinearSystemType
    use SlepcVariables_mod, only: EVPType

    implicit none
    
    PetscInt, allocatable, dimension(:) :: Inflow_Nodes_Main

    contains

    subroutine connectInflowToMainArrays(Nnodes, Mesh_Inflow, Mesh_Main, Boundary)

        use Tools_mod, only: allocateArray
        use MeshParameters_mod, only: eps => Mesh_Tolerance

        implicit none

        PetscInt, intent(in) :: Nnodes
        type(MeshParameters), intent(in) :: Mesh_Inflow
        type(MeshParameters), intent(in) :: Mesh_Main
        type(BoundaryParameters), dimension(:), intent(in) :: Boundary

        PetscInt :: inod, jnod, gnod, ibnd, jbnd
        PetscReal, dimension(3) :: dx

        call allocateArray(Inflow_Nodes_Main, Nnodes)

        do ibnd = 1, size(Boundary)
            if (Boundary(ibnd)%name == 'Inlet') then
                jbnd = ibnd
                exit
            end if
        end do

        outer: do inod = 1, Nnodes

            do jnod = 1, size(Boundary(jbnd)%gnodes_total)

                gnod = Boundary(jbnd)%gnodes_total(jnod)
                dx(:) = abs(Mesh_Main%Xi_Mesh(gnod,:) - Mesh_Inflow%Xi_Mesh(inod,:))

                if ( all(dx(:) < eps) ) then
                    Inflow_Nodes_Main(inod) = gnod
                    cycle outer
                end if

            end do

        end do outer

        Inflow_Nodes_Main = pack(Inflow_Nodes_Main, Inflow_Nodes_Main /= 0)

    end subroutine connectInflowToMainArrays

    ! ----------------------------------------------------------------------

    subroutine clearEntries(Neq, bound, Flag_NR, LS)

        use Tools_mod, only: getRows

        implicit none

        PetscInt, intent(in) :: Neq
        type(BoundaryParameters), intent(in) :: bound
        character(*), intent(in) :: Flag_NR
        type(LinearSystemType), intent(inout) :: LS

        PetscScalar, parameter :: zero = 0.0d0
        PetscScalar, dimension(size(bound%gnodes_rank)) :: zeros
        PetscInt, dimension(size(bound%gnodes_rank)) :: irows
        PetscInt, dimension(1) :: jeq
        PetscInt :: ieq
        PetscErrorCode :: ierr

        zeros(:) = 0.0d0
        do ieq = 1, size(bound%eq_clear)

            jeq(1) = bound%eq_clear(ieq)

            irows(:) = getRows(bound%gnodes_rank(:), jeq(:), Neq) - 1

            call VecSetValues(LS%b_f, size(irows), irows, zeros, INSERT_VALUES, ierr)

            if (Flag_NR == 'NRP') then
                call MatZeroRows(LS%A_f, size(irows), irows, zero, Petsc_NULL_VEC, Petsc_NULL_VEC, ierr)
            end if

        end do

    end subroutine clearEntries

    ! ----------------------------------------------------------------------

    subroutine applyDirichletBC(Neq, bound, TL_Rank, Flag_NR, LS)

        use Tools_mod, only: getRows

        implicit none

        PetscInt, intent(in) :: Neq
        type(BoundaryParameters), intent(in) :: bound
        PetscReal, dimension(:,:), intent(in) :: TL_Rank
        character(*), intent(in) :: Flag_NR
        type(LinearSystemType), intent(inout) :: LS

        PetscScalar, parameter :: one = 1.0d0
        PetscScalar, dimension(size(bound%gnodes_rank)) :: res_values
        PetscInt, dimension(size(bound%gnodes_rank)) :: irows
        PetscInt, dimension(1) :: jeq
        PetscInt :: ieq, i
        PetscErrorCode :: ierr

        do ieq = 1, size(bound%values_Dirichlet)

            jeq(1) = bound%eq_Dirichlet(ieq)
            
            irows(:) = getRows(bound%gnodes_rank(:),jeq(:),Neq) - 1

            res_values(:) = TL_Rank(bound%rnodes(:),jeq(1)) - bound%values_Dirichlet(ieq)

            call VecSetValues(LS%b_f, size(irows), irows, res_values, INSERT_VALUES, ierr)

            if (Flag_NR == 'NRP') then
                call MatZeroRows(LS%A_f, size(irows), irows, one, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
                ! !Only insert 1 in the diagonal
                ! call MatSetValues(LS%A_f, size(irows), irows, size(irows), irows, one, INSERT_VALUES, ierr)
                ! !same
                ! do i = 1, size(irows)
                !     call MatSetValues(LS%A_f, 1, irows(i), 1, irows(i), one, INSERT_VALUES, ierr)
                ! end do
            end if

        end do

    end subroutine applyDirichletBC

    ! ----------------------------------------------------------------------

    subroutine applyDirichletBCFullyDeveloped(Flag_NR, LS)

        use Tools_mod, only: getRows
        use PhysicalParameters_mod, only: Problem_Main, Problem_Inflow
        use MPIParameters_mod, only: Rank
        use SolutionVariables_mod, only: Sol_Main, Sol_Inflow

        implicit none

        character(*), intent(in) :: Flag_NR
        type(LinearSystemType), intent(inout) :: LS

        PetscScalar, parameter :: one = 1.0d0
        PetscScalar, dimension(size(Inflow_Nodes_Main)) :: res_values
        PetscInt, dimension(size(Inflow_Nodes_Main)) :: irows
        PetscInt, dimension(Problem_Inflow%Neq) :: jeq
        PetscInt :: ieq
        PetscErrorCode :: ierr

        do ieq = 1, Problem_Inflow%Neq
            jeq(ieq) = ieq
            !Skip P
            if (ieq > Problem_Inflow%Neq_f) jeq(ieq) = ieq+1
        end do

        do ieq = 1, Problem_Inflow%Neq

            irows(:) = getRows(Inflow_Nodes_Main(:),jeq(ieq:ieq),Problem_Main%Neq) - 1

            if (Rank == 0) then
                res_values(:) = Sol_Main%TL(Inflow_Nodes_Main(:),jeq(ieq)) &
                                - Sol_Inflow%TL(:,ieq)
                call VecSetValues(LS%b_f, size(irows), irows, res_values, INSERT_VALUES, ierr)
            end if

            if (Flag_NR == 'NRP') then
                call MatZeroRows(LS%A_f, size(irows), irows, one, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
            end if

        end do

    end subroutine applyDirichletBCFullyDeveloped

    ! ----------------------------------------------------------------------

    subroutine applyPeriodicBC(Neq, bound, Sol, Flag_NR, LS)

        use Tools_mod, only: getRows
        use MPIParameters_mod, only: Rank
        use Petsc_mod, only: flushAssembly, finalAssembly

        implicit none

        PetscInt, intent(in) :: Neq
        type(BoundaryParameters), intent(in) :: bound
        type(SolutionArraysType), intent(in) :: Sol
        character(*), intent(in) :: Flag_NR
        type(LinearSystemType), intent(inout) :: LS
        
        PetscScalar, parameter :: one = 1.0d0
        PetscScalar, dimension(Neq) :: res_values
        PetscInt, dimension(Neq) :: irows, jeqs
        PetscInt, dimension(1) :: jcol_p
        PetscInt :: ieq, inod, gnod, gnod_periodic, i
        PetscErrorCode :: ierr

        jeqs = [ (ieq, ieq = 1, size(jeqs)) ]
        do inod = 1, size(bound%gnodes_total)

            irows(:) = getRows(bound%gnodes_total(inod:inod),jeqs,Neq) - 1

            !Insert TL(gnod,:) - TL(gnod_p,:) = 0
            if (Rank == 0) then
                gnod = bound%gnodes_total(inod)
                gnod_periodic = bound%gnodes_periodic(inod)
                res_values(:) = Sol%TL(gnod,:) - Sol%TL(gnod_periodic,:)
                call VecSetValues(LS%b_f, size(irows), irows, res_values, INSERT_VALUES, ierr)
            end if

            if (Flag_NR == 'NRP') then
                !Insert +1
                call MatZeroRows(LS%A_f, size(irows), irows, one, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
            end if

        end do

        if (Flag_NR == 'NRP') then

            ! call flushAssembly(Flag_NR, LS)

            do inod = 1, size(bound%gnodes_total)

                irows(:) = getRows(bound%gnodes_total(inod:inod),jeqs,Neq) - 1

                do i = 1, size(irows)
                    jcol_p(:) = getRows(bound%gnodes_periodic(inod:inod),jeqs(i:i),Neq) - 1
                    !Insert -1
                    call MatSetValues(LS%A_f, 1, irows(i), 1, jcol_p(1), -one, INSERT_VALUES, ierr)
                end do

            end do

        end if

    end subroutine applyPeriodicBC

    ! ----------------------------------------------------------------------

    subroutine clearEntriesStability(Neq, bound, A_f, B_f)

        use Tools_mod, only: getRows

        implicit none

        PetscInt, intent(in) :: Neq
        type(BoundaryParameters), intent(in) :: bound
        Mat, intent(inout) :: A_f, B_f

        PetscScalar, parameter :: zero = 0.0d0
        PetscInt, dimension(size(bound%gnodes_rank)) :: irows
        PetscInt, dimension(1) :: jeq
        PetscInt :: ieq
        PetscErrorCode :: ierr

        do ieq = 1, size(bound%eq_clear_S)

            jeq(1) = bound%eq_clear_S(ieq)
            irows(:) = getRows(bound%gnodes_rank, jeq, Neq) - 1

            call MatZeroRows(A_f, size(irows), irows, zero, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
            call MatZeroRows(B_f, size(irows), irows, zero, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)

        end do

    end subroutine clearEntriesStability

    ! ----------------------------------------------------------------------

    subroutine applyDirichletBCStability(Neq, bound, A_f, B_f)

        use Tools_mod, only: getRows

        implicit none

        PetscInt, intent(in) :: Neq
        type(BoundaryParameters), intent(in) :: bound
        Mat, intent(inout) :: A_f, B_f

        PetscScalar, parameter :: one = 1.0d0, zero = 0.0d0
        PetscInt, dimension(size(bound%gnodes_rank)) :: irows
        PetscInt, dimension(1) :: jeq
        PetscInt :: ieq
        PetscErrorCode :: ierr

        do ieq = 1, size(bound%eq_Dirichlet_S)

            jeq(1) = bound%eq_Dirichlet_S(ieq)

            irows(:) = getRows(bound%gnodes_rank, jeq, Neq) - 1

            call MatZeroRows(A_f, size(irows), irows,  one, &
                PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
            call MatZeroRows(B_f, size(irows), irows, zero, &
                PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)

        end do

    end subroutine applyDirichletBCStability

    ! ----------------------------------------------------------------------

    subroutine applyPeriodicBCStability(Neq, bound, EVP)

        use petscksp, only: PETSC_i
        use PhysicalParameters_mod, only: PI, QbN
        use Tools_mod, only: getRows
        use Petsc_mod, only: flushAssemblyEVPObjects, finalAssemblyEVPObjects
        use MPIParameters_mod, only: Rank

        implicit none

        PetscInt, intent(in) :: Neq
        type(BoundaryParameters), intent(in) :: bound
        type(EVPType), intent(inout) :: EVP
        
        PetscScalar, parameter :: one = 1.0d0, zero = 0.0d0
        PetscScalar :: value
        PetscInt, dimension(Neq) :: irows, jeqs
        PetscInt, dimension(1) :: jcol_p
        PetscInt :: ieq, inod, i
        PetscErrorCode :: ierr

        value = exp(2.0d0*PI*QbN*PETSC_i)
        if (Rank == 0) then
            write(*,'(a)', advance='no') 'PBC value = '
            write(*,*) value
        end if
        
        jeqs = [ (ieq, ieq = 1, size(jeqs)) ]
        do inod = 1, size(bound%gnodes_total)

            irows(:) = getRows(bound%gnodes_total(inod:inod),jeqs,Neq) - 1

            !Insert +1
            ! call MatZeroRows(EVP%A_f,size(irows), irows,  one, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
            call MatZeroRows(EVP%A_f,size(irows), irows, value, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)
            call MatZeroRows(EVP%B_f,size(irows), irows, zero, PETSC_NULL_VEC, PETSC_NULL_VEC, ierr)

        end do

        call flushAssemblyEVPObjects(EVP%A_f, EVP%B_f)
        ! call finalAssemblyEVPObjects(EVP%A_f, EVP%B_f)

        do inod = 1, size(bound%gnodes_total)

            irows(:) = getRows(bound%gnodes_total(inod:inod),jeqs,Neq) - 1

            do i = 1, size(irows)
                jcol_p(:) = getRows(bound%gnodes_periodic(inod:inod),jeqs(i:i),Neq) - 1
                !Insert -1
                ! call MatSetValues(EVP%A_f, 1, irows(i), 1, jcol_p(1), value, INSERT_VALUES, ierr)
                call MatSetValues(EVP%A_f, 1, irows(i), 1, jcol_p(1), -one, INSERT_VALUES, ierr)
                ! call MatSetValues(EVP%B_f, 1, irows(i), 1, jcol_p(1),   zero, INSERT_VALUES, ierr)
            end do

        end do

    end subroutine applyPeriodicBCStability

end module Dirichlet_mod
