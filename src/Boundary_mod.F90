
#include <petsc/finclude/petscksp.h>

module Boundary_mod

    use petscksp
    use PhysicalParameters_mod, only: ProblemParameters
    use FEMParameters_mod, only: FEParameters
    use MeshParameters_mod, only: MeshParameters
    use BoundaryParameters_mod, only: BoundaryParameters
    use MPIParameters_mod, only: Rank

    implicit none

    interface
        function distanceToBND(ibnd, x) result(dx)
            implicit none
            PetscInt, intent(in) :: ibnd
            PetscReal, dimension(:,:), intent(in) :: x
            PetscReal, dimension(size(x,1)) :: dx
        end function distanceToBND
    end interface

    contains
    
    subroutine setProblemBoundaries(Problem, Problem_Stability, FE, Mesh, Boundary)

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), dimension(:), allocatable, intent(out) :: Boundary

        allocate( Boundary(Problem%Nbd) )

        call setBoundaryNames(Problem%name, Boundary(:)%name)
        call setEquationsOnBoundaries(Problem, Problem_Stability, Boundary)
        call setBoundaryArrays(Problem, FE, Mesh, Boundary)

    end subroutine setProblemBoundaries
    
    ! ----------------------------------------------------------------------
    
    subroutine setBoundaryNames(problem_name, name)

        use PhysicalParameters_mod, only: Periodicity, Symmetry

        implicit none

        character(*), intent(in) :: problem_name
        character(*), dimension(:), intent(out) :: name

        PetscInt :: i

        select case (problem_name)

        case ('Inflow_1D')

            name(:) = 'Wall'

            do i = 1, size(Symmetry)
                select case (Symmetry(i))
                case ('Y')
                    name(:) = ['Y-Symmetry', 'Y-Symmetry']
                case default
                    continue
                end select
            end do

            do i = 1, size(Periodicity)
                select case (Periodicity(i))
                case ('Y')
                    name(:) = ['Periodic', 'Periodic']
                case default
                    continue
                end select
            end do

        case ('Inflow_2D')

            name(:) = 'Wall'

            do i = 1, size(Symmetry)
                select case (Symmetry(i))
                case ('Y')
                    name(:) = ['Y-Symmetry', 'Y-Symmetry', 'Wall', 'Wall']
                case ('Z')
                    name(:) = ['Wall', 'Wall', 'Z-Symmetry', 'Z-Symmetry']
                case default
                    continue
                end select
            end do

            do i = 1, size(Periodicity)
                select case (Periodicity(i))
                case ('Y')
                    name(:) = ['Periodic', 'Periodic', 'Wall', 'Wall']
                case ('Z')
                    name(:) = ['Wall', 'Wall', 'Periodic', 'Periodic']
                case default
                    continue
                end select
            end do

        case ('Cylinder_2D', 'Cylinder_Stability_2D', 'Cylinder_Stability2D_3D')

            name(:) = ['Inlet', 'Outlet', 'Wall', 'Wall', 'Cylinder']

            do i = 1, size(Symmetry)
                select case (Symmetry(i))
                case ('X')
                    name(:) = ['X-Symmetry', 'X-Symmetry', 'Wall', 'Wall', 'Cylinder']
                case ('Y')
                    name(:) = ['Inlet', 'Outlet', 'Y-Symmetry', 'Y-Symmetry', 'Cylinder']
                case default
                    continue
                end select
            end do

            do i = 1, size(Periodicity)
                select case (Periodicity(i))
                case ('X')
                    name(:) = ['Periodic', 'Periodic', 'Wall', 'Wall', 'Cylinder']
                case ('Y')
                    name(:) = ['Inlet', 'Outlet', 'Periodic', 'Periodic', 'Cylinder']
                case default
                    continue
                end select
            end do
        
        case ('Cylinder_3D', 'Cylinder_Stability_3D')

            name(:) = ['Inlet', 'Outlet', 'Wall', 'Wall', 'Cylinder', 'Wall', 'Wall']

            do i = 1, size(Symmetry)
                select case (Symmetry(i))
                case ('X')
                    name(:) = ['X-Symmetry', 'X-Symmetry', 'Wall', 'Wall', &
                                        'Cylinder', 'Wall', 'Wall']
                case ('Y')
                    name(:) = ['Inlet', 'Outlet', 'Y-Symmetry', 'Y-Symmetry', &
                                        'Cylinder', 'Wall', 'Wall']
                case ('Z')
                    name(:) = ['Inlet', 'Outlet', 'Wall', 'Wall', 'Cylinder', &
                                        'Z-Symmetry', 'Z-Symmetry']
                case default
                    continue
                end select
            end do

            do i = 1, size(Periodicity)
                select case (Periodicity(i))
                case ('X')
                    name(:) = ['Periodic', 'Periodic', 'Wall', 'Wall', &
                                        'Cylinder', 'Wall', 'Wall']
                case ('Y')
                    name(:) = ['Inlet', 'Outlet', 'Periodic', 'Periodic', &
                                        'Cylinder', 'Wall', 'Wall']
                case ('Z')
                    name(:) = ['Inlet', 'Outlet', 'Wall', 'Wall', &
                                        'Cylinder', 'Periodic', 'Periodic']
                case default
                    continue
                end select
            end do

        case default
            
        end select

    end subroutine setBoundaryNames
    
    ! ----------------------------------------------------------------------

    subroutine setEquationsOnBoundaries(Problem, Problem_Stability, Boundary)

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(BoundaryParameters), dimension(:), intent(inout) :: Boundary

        PetscInt :: ibnd, Ndim_M
        PetscBool :: check

        do ibnd = 1, size(Boundary)

            call setBoundaryEQs(Problem, Boundary(ibnd))

            check = (Problem_Stability%Ndim > 0)
            if (check) then
                Ndim_M = Problem%Ndim
                call setBoundaryEQsStability(Ndim_M, Problem_Stability, Boundary(ibnd))
            end if

        end do

    end subroutine setEquationsOnBoundaries

    ! ----------------------------------------------------------------------

    subroutine setBoundaryEQs(Problem, bound)

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(BoundaryParameters), intent(inout):: bound

        PetscInt :: i, Neq_f, Neq_s, Neq, Ndim

        Neq_f = Problem%Neq_f
        Neq_s = Problem%Neq_s
        Neq = Problem%Neq
        Ndim = Problem%Ndim

        select case (bound%name)

        case ('Inlet')

        case ('Outlet')

            if (Ndim == 2) then
                allocate( bound%eq_ed(1) )
                bound%eq_ed(:) = [1]      !Open BC
            end if
            if (Ndim == 3) then
                allocate( bound%eq_fc(1) )
                bound%eq_fc(:) = [1]      !Open BC
            end if

            allocate( bound%eq_Dirichlet(1) )
            bound%eq_Dirichlet(:) = [ Neq_f+1 ]    !P

            allocate( bound%values_Dirichlet(1) )
            bound%values_Dirichlet(:) = 0.0d0

        case ('Wall', 'Cylinder')

            allocate( bound%eq_Dirichlet(Neq_f) )
            bound%eq_Dirichlet(:) = [ (i, i = 1, Neq_f) ]    !Ux, Uy, Uz

            allocate( bound%values_Dirichlet(Neq_f) )
            bound%values_Dirichlet(:) = 0.0d0

            if (bound%name == 'Cylinder') then
                if (Ndim == 2) then
                    allocate( bound%epp_ed(1) )
                    bound%epp_ed(:) = [1]   !Fd
                end if
                if (Ndim == 3) then
                    allocate( bound%epp_fc(1) )
                    bound%epp_fc(:) = [1]   !Fd
                end if
            end if

        case ('Y-Symmetry')

            select case (Neq_s)
            case (3)
                allocate( bound%eq_Dirichlet(2) )
                bound%eq_Dirichlet(:) = [ 2, Neq-1 ]    !Uy, Txy
                allocate( bound%values_Dirichlet(2) )
            case (6)
                allocate( bound%eq_Dirichlet(3) )
                bound%eq_Dirichlet(:) = [ 2, Neq-Neq_s+2, Neq-1 ]    !Uy, Txy, Tyz
                allocate( bound%values_Dirichlet(3) )
            case default
                write(*,'(a)') 'Wrong Neq_s in Y-Symmetry in setBoundaryEQs!'
                stop
            end select

            bound%values_Dirichlet(:) = 0.0d0

        case ('Z-Symmetry')

            allocate( bound%eq_Dirichlet(3) )
            bound%eq_Dirichlet(:) = [ Neq_f, Neq-2, Neq-1 ]    !Uz, Txz, Tyz

            allocate( bound%values_Dirichlet(3) )
            bound%values_Dirichlet(:) = 0.0d0

        case ('Slip_Wall', 'Slip_Cylinder')

            allocate( bound%eq_clear(Neq_f) )
            bound%eq_clear(:) = [ (i, i = 1, Neq_f) ]  !Ux, Uy, Uz

            if (Ndim == 2) then
                allocate( bound%eq_ed(1) )
                bound%eq_ed(:) = [2]      !Slip BC
            end if
            if (Ndim == 3) then
                allocate( bound%eq_fc(1) )
                bound%eq_fc(:) = [2]      !Slip BC
            end if

            if (bound%name == 'Slip_Cylinder') then
                if (Ndim == 2) then
                    allocate( bound%epp_ed(1) )
                    bound%epp_ed(:) = [1]   !Fd
                end if
                if (Ndim == 3) then
                    allocate( bound%epp_fc(1) )
                    bound%epp_fc(:) = [1]   !Fd
                end if
            end if

        case ('Periodic')


        case default
            write(*,'(a)') 'Wrong bound%name in setBoundaryEQs!'
            stop
        end select

    end subroutine setBoundaryEQs

    ! ----------------------------------------------------------------------

    subroutine setBoundaryEQsStability(Ndim_M, Problem_Stability, bound)

        implicit none

        PetscInt, intent(in) :: Ndim_M
        type(ProblemParameters), intent(in) :: Problem_Stability
        type(BoundaryParameters), intent(inout) :: bound

        PetscInt :: i, Neq_f, Neq_s, Neq

        Neq_f = Problem_Stability%Neq_f
        Neq_s = Problem_Stability%Neq_s
        Neq = Problem_Stability%Neq

        select case (bound%name)

        case ('Inlet')

            allocate( bound%eq_Dirichlet_S(Neq) )
            bound%eq_Dirichlet_S(:) = [ (i, i = 1, Neq) ]

            allocate( bound%values_Dirichlet_S(Neq) )
            bound%values_Dirichlet_S(:) = 0.0d0

        case ('Outlet')

            ! if (Ndim_M == 2) then
            !     allocate( bound%eq_ed_S(1) )
            !     bound%eq_ed_S(:) = [1]      !Open BC
            !     allocate( bound%eq_fc_S(0) )
            ! end if
            ! if (Ndim_M == 3) then
            !     allocate( bound%eq_ed_S(0) )
            !     allocate( bound%eq_fc_S(1) )
            !     bound%eq_fc_S(:) = [1]      !Open BC
            ! end if

            allocate( bound%eq_Dirichlet_S(Neq) )
            bound%eq_Dirichlet_S(:) = [ (i, i = 1, Neq) ]

            allocate( bound%values_Dirichlet_S(Neq) )
            bound%values_Dirichlet_S(:) = 0.0d0

        case ('Wall', 'Cylinder')

            allocate( bound%eq_Dirichlet_S(Neq_f) )
            bound%eq_Dirichlet_S(:) = [ (i, i = 1, Neq_f) ]      !Ux, Uy, Uz

            allocate( bound%values_Dirichlet_S(Neq_f) )
            bound%values_Dirichlet_S(:) = 0.0d0

            ! if (bound%name == 'Cylinder') then
            !     if (Ndim_M == 2) then
            !         allocate( bound%epp_ed_S(1) )
            !         bound%epp_ed_S(:) = [1]   !Fd
            !     end if
            !     if (Ndim_M == 3) then
            !         allocate( bound%epp_fc_S(1) )
            !         bound%epp_fc_S(:) = [1]   !Fd
            !     end if
            ! end if

        case ('Y-Symmetry')

            select case (Neq_s)
            case (3)
                allocate( bound%eq_Dirichlet_S(2) )
                bound%eq_Dirichlet_S(:) = [ 2, Neq-1 ]    !Uy, Txy
                allocate( bound%values_Dirichlet_S(2) )
            case (6)
                allocate( bound%eq_Dirichlet_S(3) )
                bound%eq_Dirichlet_S(:) = [ 2, Neq-Neq_s+2, Neq-1 ]    !Uy, Txy, Tyz
                allocate( bound%values_Dirichlet_S(3) )
            case default
                write(*,'(a)') 'Wrong Neq_s in Y-Symmetry in setBoundaryEQs!'
                stop
            end select

            bound%values_Dirichlet_S(:) = 0.0d0

        case ('Z-Symmetry')

            allocate( bound%eq_Dirichlet_S(3) )
            bound%eq_Dirichlet_S(:) = [ Neq_f, Neq-2, Neq-1 ]    !Uz, Txz, Tyz

            allocate( bound%values_Dirichlet_S(3) )
            bound%values_Dirichlet_S(:) = 0.0d0

        case ('Slip_Wall', 'Slip_Cylinder')

            allocate( bound%eq_clear_S(Neq_f) ) 
            bound%eq_clear_S(:) = [ (i, i = 1, Neq_f) ]  !Ux, Uy, Uz

            if (Ndim_M == 2) then
                allocate( bound%eq_ed_S(1) )
                bound%eq_ed_S(:) = [2]      !Slip BC
            end if
            if (Ndim_M == 3) then
                allocate( bound%eq_fc_S(1) )
                bound%eq_fc_S(:) = [2]      !Slip BC
            end if

            ! if (bound%name == 'Slip_Cylinder') then
            !     if (Ndim_M == 2) then
            !         allocate( bound%epp_ed_S(1) )
            !         bound%epp_ed_S(:) = [1]   !Fd
            !     end if
            !     if (Ndim_M == 3) then
            !         allocate( bound%epp_fc_S(1) )
            !         bound%epp_fc_S(:) = [1]   !Fd
            !     end if
            ! end if

        case ('Periodic')

        case default
            write(*,'(a)') 'Wrong bound%name in setBoundaryEQsStability!'
            stop
        end select
        
    end subroutine setBoundaryEQsStability

    ! ----------------------------------------------------------------------

    subroutine setBoundaryArrays(Problem, FE, Mesh, Boundary)

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        type(BoundaryParameters), dimension(:), intent(inout) :: Boundary

        procedure(distanceToBND), pointer :: distanceToBND_p => NULL()

        select case (Problem%name)
        case ('Inflow_1D', 'Inflow_2D')
            distanceToBND_p => distanceToBNDInflow
        case ('Cylinder_2D', 'Cylinder_3D')
            distanceToBND_p => distanceToBNDCylinder
        case default
            distanceToBND_p => NULL()
        end select

        call setBoundaryEdges(FE, Mesh, distanceToBND_p, Boundary)

        call setBoundaryFaces(FE, Mesh, distanceToBND_p, Boundary)

        call setBoundaryNodesRank(FE, Mesh, distanceToBND_p, Boundary)

        call setBoundaryCoordinates(FE%Nnodes, Mesh%Xi_Mesh, distanceToBND_p, Boundary)

        call setPeriodicBoundaryNodes(Problem, Boundary)

    end subroutine setBoundaryArrays

    ! ----------------------------------------------------------------------

    subroutine setBoundaryEdges(FE, Mesh, distanceToBND_p, Boundary)

        use Tools_mod, only: allocateArray
        use MeshParameters_mod, only: eps => Mesh_Tolerance

        implicit none

        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        procedure(distanceToBND), pointer, intent(in) :: distanceToBND_p
        type(BoundaryParameters), dimension(:), intent(inout) :: Boundary

        PetscInt, dimension(FE%Nbf_1D) :: nodes
        PetscReal, dimension(FE%Nbf_1D) :: dxi
        PetscReal, dimension(FE%Nbf_1D,3) :: xi
        PetscInt :: ied, ibnd, iel_rank, Nedges_rank

        do ibnd = 1, size(Boundary)
            call allocateArray(Boundary(ibnd)%edges, FE%Nel_Rank*FE%Ned)
            call allocateArray(Boundary(ibnd)%edge_elements_rank, FE%Nel_Rank*FE%Ned)
        end do

        Nedges_rank = 0
        do iel_rank = 1, FE%Nel_Rank
            do ied = 1, FE%Ned

                nodes(:) = setEdgeNodes(iel_rank, ied, FE)
                xi(:,:) = Mesh%Xi_Mesh_Rank(nodes(:),:)

                do ibnd = 1, size(Boundary)

                    dxi(:) = distanceToBND_p(ibnd, xi)

                    if ( all(dxi < eps) ) then

                        Nedges_rank = Nedges_rank + 1

                        Boundary(ibnd)%edges(Nedges_rank) = ied
                        Boundary(ibnd)%edge_elements_rank(Nedges_rank) = iel_rank

                    end if

                end do

            end do
        end do

        do ibnd = 1, size(Boundary)
            associate(B => Boundary(ibnd))
                B%edges = pack(B%edges, B%edges /= 0)
                B%edge_elements_rank = pack(B%edge_elements_rank, B%edge_elements_rank /= 0)
            end associate
        end do

    end subroutine setBoundaryEdges

    ! ----------------------------------------------------------------------

    function setEdgeNodes(iel_rank, ied, FE) result(nodes)

        implicit none

        PetscInt, intent(in) :: iel_rank, ied
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(FE%Nbf_1D) :: nodes

        select case (FE%name)
        case ('Triangle')
            nodes(:) = setEdgeNodesTriangle(iel_rank, ied, FE)
        case ('Quadrangle')
            nodes(:) = setEdgeNodesQuadrangle(iel_rank, ied, FE)
        case ('Tetrahedron')
            nodes(:) = setEdgeNodesTetrahedron(iel_rank, ied, FE)
        case ('Hexahedron')
            nodes(:) = setEdgeNodesHexahedron(iel_rank, ied, FE)
        case default
            write(*,'(a)') 'Wrong FE%name in setEdgeNodes!'
            stop
        end select

    end function setEdgeNodes

    ! ----------------------------------------------------------------------

    function setEdgeNodesTriangle(iel_rank, ied, FE) result(nodes)

        use Tools_mod, only: getRankNode

        implicit none

        PetscInt, intent(in) :: iel_rank, ied
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(FE%Nbf_1D) :: nodes

        PetscInt, dimension(FE%Nbf_1D) :: numbering
        PetscInt :: i, inod

        select case (FE%order)
        case ('Linear')

            select case (ied)
            case (1)
                !η = 0
                numbering(:) = [1, 2]
            case (2)
                !ξ = 0
                numbering(:) = [1, 3]
            case (3)
                !η = 1-ξ
                numbering(:) = [2, 3]
            case default
                write(*,'(A)') 'Wrong ied in Linear setEdgeNodesTriangle!'
                stop
            end select

        case ('Quadratic', 'Serendipity')

            select case (ied)
            case (1)
                !η = 0
                numbering(:) = [1, 2, 4]
            case (2)
                !ξ = 0
                numbering(:) = [1, 3, 6]
            case (3)
                !η = 1-ξ
                numbering(:) = [2, 3, 5]
            case default
                write(*,'(A)') 'Wrong ied in Quadratic/Serendipity setEdgeNodesTriangle!'
                stop
            end select

        case default
            write(*,'(a)') 'Wrong FE%order in setEdgeNodesTriangle!'
            stop
        end select

        do i = 1, size(nodes)
            inod = numbering(i)
            nodes(i) = getRankNode(iel_rank,inod,FE%Nbf)
        end do

    end function setEdgeNodesTriangle

    ! ----------------------------------------------------------------------

    function setEdgeNodesQuadrangle(iel_rank, ied, FE) result(nodes)

        use Tools_mod, only: getRankNode

        implicit none

        PetscInt, intent(in) :: iel_rank, ied
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(FE%Nbf_1D) :: nodes

        PetscInt, dimension(FE%Nbf_1D) :: numbering
        PetscInt :: i, inod

        select case (FE%order)
        case ('Linear')

            select case (ied)
            case (1)
                !η = -1
                numbering(:) = [1, 2]
            case (2)
                !ξ = +1
                numbering(:) = [2, 3]
            case (3)
                !η = +1
                numbering(:) = [3, 4]
            case (4)
                !ξ = -1
                numbering(:) = [4, 1]
            case default
                write(*,'(A)') 'Wrong ied in Linear setEdgeNodesQuadrangle!'
                stop
            end select

        case ('Quadratic', 'Serendipity')

            select case (ied)
            case (1)
                !η = -1
                numbering(:) = [1, 2, 5]
            case (2)
                !ξ = +1
                numbering(:) = [2, 3, 6]
            case (3)
                !η = +1
                numbering(:) = [3, 4, 7]
            case (4)
                !ξ = -1
                numbering(:) = [4, 1, 8]
            case default
                write(*,'(A)') 'Wrong ied in Quadratic/Serendipity setEdgeNodesQuadrangle!'
                stop
            end select

        case default
            write(*,'(a)') 'Wrong FE%order in setEdgeNodesQuadrangle!'
            stop
        end select

        do i = 1, size(nodes)
            inod = numbering(i)
            nodes(i) = getRankNode(iel_rank,inod,FE%Nbf)
        end do

    end function setEdgeNodesQuadrangle

    ! ----------------------------------------------------------------------

    function setEdgeNodesTetrahedron(iel_rank, ied, FE) result(nodes)

        use Tools_mod, only: getRankNode

        implicit none

        PetscInt, intent(in) :: iel_rank, ied
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(FE%Nbf_1D) :: nodes

        PetscInt, dimension(FE%Nbf_1D) :: numbering
        PetscInt :: i, inod

        select case (FE%order)
        case ('Linear')

            select case (ied)
            case (1)
                !η = 0, ζ = 0
                numbering(:) = [1, 2]
            case (2)
                !ξ = 0, ζ = 0
                numbering(:) = [1, 3]
            case (3)
                !ξ = 0, η = 0
                numbering(:) = [1, 4]
            case (4)
                !η = 1-ξ, ζ = 0
                numbering(:) = [2, 3]
            case (5)
                !ξ = 0, ζ = 1-η
                numbering(:) = [3, 4]
            case (6)
                !η = 0, ζ = 1-ξ
                numbering(:) = [2, 4]
            case default
                write(*,'(A)') 'Wrong ied in Linear setEdgeNodesTetrahedron!'
                stop
            end select

        case ('Quadratic', 'Serendipity')

            select case (ied)
            case (1)
                !η = 0, ζ = 0
                numbering(:) = [1, 2, 5]
            case (2)
                !ξ = 0, ζ = 0
                numbering(:) = [1, 3, 7]
            case (3)
                !ξ = 0, η = 0
                numbering(:) = [1, 4, 8]
            case (4)
                !η = 1-ξ, ζ = 0
                numbering(:) = [2, 3, 6]
            case (5)
                !ξ = 0, ζ = 1-η
                numbering(:) = [3, 4, 10]
            case (6)
                !η = 0, ζ = 1-ξ
                numbering(:) = [2, 4, 9]
            case default
                write(*,'(A)') 'Wrong ied in Quadratic/Serendipity setEdgeNodesTetrahedron!'
                stop
            end select

        case default
            write(*,'(a)') 'Wrong FE%order in setEdgeNodesTetrahedron!'
            stop
        end select

        do i = 1, size(nodes)
            inod = numbering(i)
            nodes(i) = getRankNode(iel_rank,inod,FE%Nbf)
        end do

    end function setEdgeNodesTetrahedron

    ! ----------------------------------------------------------------------

    function setEdgeNodesHexahedron(iel_rank, ied, FE) result(nodes)

        use Tools_mod, only: getRankNode

        implicit none

        PetscInt, intent(in) :: iel_rank, ied
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(FE%Nbf_1D) :: nodes

        PetscInt, dimension(FE%Nbf_1D) :: numbering
        PetscInt :: i, inod

        select case (FE%order)
        case ('Linear')

            select case (ied)
            case (1)
                !η = -1, ζ = -1
                numbering(:) = [1, 2]
            case (2)
                !ξ = +1, ζ = -1
                numbering(:) = [2, 3]
            case (3)
                !η = +1, ζ = -1
                numbering(:) = [3, 4]
            case (4)
                !ξ = -1, ζ = -1
                numbering(:) = [4, 1]

            case (5)
                !η = -1, ζ = +1
                numbering(:) = [5, 6]
            case (6)
                !ξ = +1, ζ = +1
                numbering(:) = [6, 7]
            case (7)
                !η = +1, ζ = +1
                numbering(:) = [7, 8]
            case (8)
                !ξ = -1, ζ = +1
                numbering(:) = [8, 5]

            case (9)
                !ξ = -1, η = -1
                numbering(:) = [1, 5]
            case (10)
                !ξ = +1, η = -1
                numbering(:) = [2, 6]
            case (11)
                !ξ = +1, η = +1
                numbering(:) = [3, 7]
            case (12)
                !ξ = -1, η = +1
                numbering(:) = [4, 8]
            case default
                write(*,'(A)') 'Wrong ied in Linear setEdgeNodesHexahedron!'
                stop
            end select

        case ('Quadratic', 'Serendipity')

            select case (ied)
            case (1)
                !η = -1, ζ = -1
                numbering(:) = [1, 2, 9]
            case (2)
                !ξ = +1, ζ = -1
                numbering(:) = [2, 3, 10]
            case (3)
                !η = +1, ζ = -1
                numbering(:) = [3, 4, 11]
            case (4)
                !ξ = -1, ζ = -1
                numbering(:) = [4, 1, 12]

            case (5)
                !η = -1, ζ = +1
                numbering(:) = [5, 6, 13]
            case (6)
                !ξ = +1, ζ = +1
                numbering(:) = [6, 7, 14]
            case (7)
                !η = +1, ζ = +1
                numbering(:) = [7, 8, 15]
            case (8)
                !ξ = -1, ζ = +1
                numbering(:) = [8, 5, 16]

            case (9)
                !ξ = -1, η = -1
                numbering(:) = [1, 5, 17]
            case (10)
                !ξ = +1, η = -1
                numbering(:) = [2, 6, 18]
            case (11)
                !ξ = +1, η = +1
                numbering(:) = [3, 7, 19]
            case (12)
                !ξ = -1, η = +1
                numbering(:) = [4, 8, 20]
            case default
                write(*,'(A)') 'Wrong ied in Quadratic/Serendipity setEdgeNodesHexahedron!'
                stop
            end select

        case default
            write(*,'(a)') 'Wrong FE%order in setEdgeNodesTetrahedron!'
            stop
        end select
        
        do i = 1, size(nodes)
            inod = numbering(i)
            nodes(i) = getRankNode(iel_rank,inod,FE%Nbf)
        end do

    end function setEdgeNodesHexahedron

    ! ----------------------------------------------------------------------

    subroutine setBoundaryFaces(FE, Mesh, distanceToBND_p, Boundary)

        use Tools_mod, only: allocateArray
        use MeshParameters_mod, only: eps => Mesh_Tolerance

        implicit none

        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        procedure(distanceToBND), pointer, intent(in) :: distanceToBND_p
        type(BoundaryParameters), dimension(:), intent(inout) :: Boundary

        PetscInt, dimension(FE%Nbf_2D) :: nodes
        PetscReal, dimension(FE%Nbf_2D) :: dxi
        PetscReal, dimension(FE%Nbf_2D,3) :: xi
        PetscInt :: ifc, ibnd, iel_rank, Nfaces_rank

        do ibnd = 1, size(Boundary)
            call allocateArray(Boundary(ibnd)%faces, FE%Nel_Rank*FE%Nfc)
            call allocateArray(Boundary(ibnd)%face_elements_rank, FE%Nel_Rank*FE%Nfc)
        end do

        Nfaces_rank = 0
        do iel_rank = 1, FE%Nel_Rank
            do ifc = 1, FE%Nfc

                nodes(:) = setFaceNodes(iel_rank, ifc, FE)
                xi(:,:) = Mesh%Xi_Mesh_Rank(nodes(:),:)

                do ibnd = 1, size(Boundary)

                    dxi(:) = distanceToBND_p(ibnd, xi)

                    if ( all(dxi < eps) ) then

                        Nfaces_rank = Nfaces_rank + 1

                        Boundary(ibnd)%faces(Nfaces_rank) = ifc
                        Boundary(ibnd)%face_elements_rank(Nfaces_rank) = iel_rank

                    end if

                end do

            end do
        end do

        do ibnd = 1, size(Boundary)
            associate(B => Boundary(ibnd))
                B%faces = pack(B%faces, B%faces /= 0)
                B%face_elements_rank = pack(B%face_elements_rank, B%face_elements_rank /= 0)
            end associate
        end do

    end subroutine setBoundaryFaces

    ! ----------------------------------------------------------------------

    function setFaceNodes(iel_rank, ifc, FE) result(nodes)

        implicit none

        PetscInt, intent(in) :: iel_rank, ifc
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(FE%Nbf_2D) :: nodes

        select case (FE%name)
        case ('Tetrahedron')
            nodes(:) = setFaceNodesTetrahedron(iel_rank, ifc, FE)
        case ('Hexahedron')
            nodes(:) = setFaceNodesHexahedron(iel_rank, ifc, FE)
        case default
            write(*,'(A)') 'Wrong FE%name in setFaceNodes!'
            stop
        end select

    end function setFaceNodes

    ! ----------------------------------------------------------------------

    function setFaceNodesTetrahedron(iel_rank, ifc, FE) result(nodes)

        use Tools_mod, only: getRankNode

        implicit none

        PetscInt, intent(in) :: iel_rank, ifc
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(FE%Nbf_2D) :: nodes

        PetscInt, dimension(FE%Nbf_2D) :: numbering
        PetscInt :: i, inod

        select case (FE%order)
        case ('Linear')

            select case (ifc)
            case (1)
                numbering(:) = [1, 3, 4]
            case (2)
                numbering(:) = [1, 4, 2]
            case (3)
                numbering(:) = [1, 2, 3]
            case (4)
                numbering(:) = [2, 3, 4]
            case default
                write(*,'(A)') 'Wrong ifc in Linear setFaceNodesTetrahedron!'
                stop
            end select

        case ('Quadratic', 'Serendipity')

            select case (ifc)
            case (1)
                numbering(:) = [1, 3, 4, 7, 10, 8]
            case (2)
                numbering(:) = [1, 2, 4, 5, 9, 8]
            case (3)
                numbering(:) = [1, 2, 3, 5, 6, 7]
            case (4)
                numbering(:) = [2, 4, 3, 6, 9, 10]
            case default
                write(*,'(a)') 'Wrong ifc in Quadratic/Serendipity setFaceNodesTetrahedron!'
                stop
            end select

        case default
            write(*,'(a)') 'Wrong FE%order in setFaceNodesTetrahedron!'
            stop
        end select

        do i = 1, size(nodes)
            inod = numbering(i)
            nodes(i) = getRankNode(iel_rank,inod,FE%Nbf)
        end do

    end function setFaceNodesTetrahedron

    ! ----------------------------------------------------------------------

    function setFaceNodesHexahedron(iel_rank, ifc, FE) result(nodes)

        use Tools_mod, only: getRankNode

        implicit none

        PetscInt, intent(in) :: iel_rank, ifc
        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(FE%Nbf_2D) :: nodes

        PetscInt, dimension(FE%Nbf_2D) :: numbering
        PetscInt :: i, inod

        select case (FE%order)
        case ('Linear')

            select case (ifc)
            case (1)
                numbering(:) = [1, 2, 3, 4]
            case (2)
                numbering(:) = [5, 6, 7, 8]
            case (3)
                numbering(:) = [1, 2, 6, 5]
            case (4)
                numbering(:) = [4, 3, 7, 8]
            case (5)
                numbering(:) = [1, 4, 8, 5]
            case (6)
                numbering(:) = [2, 3, 7, 6]
            case default
                write(*,'(A)') 'Wrong ifc in Linear setFaceNodes!'
                stop
            end select

        case ('Quadratic')

            select case (ifc)
            case (1)
                !ζ = -1
                ! numbering(:) = [1, 2, 3, 4, 9, 10, 11, 12, 25]
                numbering(:) = [1, 2, 3, 4, 9, 10, 11, 12, 21]  !Salome
            case (2)
                !ζ = +1
                numbering(:) = [5, 6, 7, 8, 13, 14, 15, 16, 26]
            case (3)
                !η = -1
                ! numbering(:) = [1, 2, 6, 5, 9, 18, 13, 17, 21]
                numbering(:) = [1, 2, 6, 5, 9, 18, 13, 17, 22]  !Salome
            case (4)
                !η = +1
                ! numbering(:) = [4, 3, 7, 8, 11, 19, 15, 20, 23]
                numbering(:) = [4, 3, 7, 8, 11, 19, 15, 20, 24]  !Salome
            case (5)
                !ξ = -1
                ! numbering(:) = [1, 4, 8, 5, 12, 20, 16, 17, 24]
                numbering(:) = [1, 4, 8, 5, 12, 20, 16, 17, 25]  !Salome
            case (6)
                !ξ = +1
                ! numbering(:) = [2, 3, 7, 6, 10, 19, 14, 18, 22]
                numbering(:) = [2, 3, 7, 6, 10, 19, 14, 18, 23]  !Salome
            case default
                write(*,'(A)') 'Wrong ifc in Quadratic setFaceNodes!'
                stop
            end select

        case ('Serendipity')

            select case (ifc)
            case (1)
                !ζ = -1
                numbering(:) = [1, 2, 3, 4, 9, 10, 11, 12]
            case (2)
                !ζ = +1
                numbering(:) = [5, 6, 7, 8, 13, 14, 15, 16]
            case (3)
                !η = -1
                numbering(:) = [1, 2, 6, 5, 9, 18, 13, 17]
            case (4)
                !η = +1
                numbering(:) = [4, 3, 7, 8, 11, 19, 15, 20]
            case (5)
                !ξ = -1
                numbering(:) = [1, 4, 8, 5, 12, 20, 16, 17]
            case (6)
                !ξ = +1
                numbering(:) = [2, 3, 7, 6, 10, 19, 14, 18]
            case default
                write(*,'(A)') 'Wrong ifc in Serendipity setFaceNodes!'
                stop
            end select

        case default
            write(*,'(a)') 'Wrong FE%order in setFaceNodesHexahedron!'
            stop
        end select

        do i = 1, size(nodes)
            inod = numbering(i)
            nodes(i) = getRankNode(iel_rank,inod,FE%Nbf)
        end do

    end function setFaceNodesHexahedron

    ! ----------------------------------------------------------------------

    subroutine setBoundaryNodesRank(FE, Mesh, distanceToBND_p, Boundary)

        use Tools_mod, only: getRankNode, allocateArray
        use MeshParameters_mod, only: eps => Mesh_Tolerance
        use MPIParameters_mod, only: Rank

        implicit none

        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        procedure(distanceToBND), pointer, intent(in) :: distanceToBND_p
        type(BoundaryParameters), dimension(:), intent(inout) :: Boundary

        PetscInt :: inod, iel_rank, rnod, gnod, ibnd, Nnodes_P_rank, Nnodes_rank
        PetscInt, dimension(1) :: nodes
        PetscReal, dimension(1) :: dxi
        PetscReal, dimension(1,3) :: xi

        do ibnd = 1, size(Boundary)
            call allocateArray(Boundary(ibnd)%gnodes_rank, FE%Nel_Rank*FE%Nbf)
            call allocateArray(Boundary(ibnd)%rnodes, FE%Nel_Rank*FE%Nbf)
            ! call allocateArray(Boundary(ibnd)%node_elements_rank, FE%Nel_Rank*FE%Nbf)
        end do

        Nnodes_rank = 0
        Nnodes_P_rank = 0
        do iel_rank = 1, FE%Nel_Rank
            do inod = 1, FE%Nbf

                rnod = getRankNode(iel_rank,inod,FE%Nbf)
                gnod = Mesh%Connectivity_Rank(iel_rank,inod)
                nodes(1) = rnod
                xi(:,:) = Mesh%Xi_Mesh_Rank(nodes(:),:)

                do ibnd = 1, size(Boundary)

                    dxi(:) = distanceToBND_p(ibnd, xi)

                    if ( all(dxi < eps) ) then

                        Nnodes_rank = Nnodes_rank + 1

                        Boundary(ibnd)%gnodes_rank(Nnodes_rank) = gnod
                        Boundary(ibnd)%rnodes(Nnodes_rank) = rnod
                        ! Boundary(ibnd)%node_elements_rank = iel_rank

                    end if
                end do

            end do
        end do

        do ibnd = 1, size(Boundary)
            associate(B => Boundary(ibnd))
                B%gnodes_rank = pack(B%gnodes_rank, B%gnodes_rank /= 0)
                B%rnodes = pack(B%rnodes, B%rnodes /= 0)
                ! B%node_elements_rank = pack(B%node_elements_rank, B%node_elements_rank /= 0)
            end associate
        end do

    end subroutine setBoundaryNodesRank

    ! ----------------------------------------------------------------------

    subroutine setBoundaryCoordinates(Nnodes_Total, Xi_Mesh, distanceToBND_p, Boundary)

        use Tools_mod, only: allocateArray
        use MeshParameters_mod, only: eps => Mesh_Tolerance

        implicit none

        PetscInt, intent(in) :: Nnodes_Total
        PetscReal, dimension(:,:), intent(in) :: Xi_Mesh
        procedure(distanceToBND), pointer, intent(in) :: distanceToBND_p
        type(BoundaryParameters), dimension(:), intent(inout) :: Boundary

        PetscInt :: gnod, ibnd, Nnodes
        PetscReal, dimension(1) :: dxi
        PetscReal, dimension(1,3) :: xi
        type coordinates
            PetscReal, allocatable, dimension(:) :: x, y, z
        end type coordinates
        type(coordinates), dimension(size(Boundary)) :: xi_dum

        do ibnd = 1, size(Boundary)

            call allocateArray(Boundary(ibnd)%gnodes_total, Nnodes_Total)

            call allocateArray(xi_dum(ibnd)%x, Nnodes_Total)
            xi_dum(ibnd)%x(:) = -huge(1.0d0)
            call allocateArray(xi_dum(ibnd)%y, Nnodes_Total)
            xi_dum(ibnd)%y(:) = -huge(1.0d0)
            call allocateArray(xi_dum(ibnd)%z, Nnodes_Total)
            xi_dum(ibnd)%z(:) = -huge(1.0d0)

        end do

        Nnodes = 0
        do gnod = 1, Nnodes_Total

            xi(1,:) = Xi_Mesh(gnod,:)

            bnd_loop:do ibnd = 1, size(xi_dum)

                dxi(:) = distanceToBND_p(ibnd, xi)

                if ( all(dxi < eps) ) then

                    Nnodes = Nnodes + 1

                    xi_dum(ibnd)%x(Nnodes) = Xi_Mesh(gnod,1)
                    xi_dum(ibnd)%y(Nnodes) = Xi_Mesh(gnod,2)
                    xi_dum(ibnd)%z(Nnodes) = Xi_Mesh(gnod,3)

                    Boundary(ibnd)%gnodes_total(Nnodes) = gnod

                end if
            end do bnd_loop

        end do

        do ibnd = 1, size(xi_dum)

            associate(B => Boundary(ibnd))
                B%gnodes_total = pack(B%gnodes_total, B%gnodes_total /= 0)
            end associate

            associate(dum => xi_dum(ibnd))
                dum%x = pack(dum%x, dum%x > -huge(1.0d0)+1.0d0)
                dum%y = pack(dum%y, dum%y > -huge(1.0d0)+1.0d0)
                dum%z = pack(dum%z, dum%z > -huge(1.0d0)+1.0d0)
            end associate

            Nnodes = size(xi_dum(ibnd)%x)
            call allocateArray(Boundary(ibnd)%xi, Nnodes, 3)
            Boundary(ibnd)%xi(:,1) = xi_dum(ibnd)%x(:)
            Boundary(ibnd)%xi(:,2) = xi_dum(ibnd)%y(:)
            Boundary(ibnd)%xi(:,3) = xi_dum(ibnd)%z(:)

            deallocate( xi_dum(ibnd)%x )
            deallocate( xi_dum(ibnd)%y )
            deallocate( xi_dum(ibnd)%z )

        end do

    end subroutine setBoundaryCoordinates

    ! ----------------------------------------------------------------------

    subroutine setPeriodicBoundaryNodes(Problem, Boundary)

        use Tools_mod, only: allocateArray, getRows
        use PhysicalParameters_mod, only: Periodicity
        use MeshParameters_mod, only: eps => Mesh_Tolerance
        use MPIParameters_mod, only: Rank

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(BoundaryParameters), dimension(:), intent(inout) :: Boundary

        PetscInt :: i, inod, jnod, ibnd, jbnd, lb, ub, Nnodes, Neq, iextra
        PetscInt, dimension(Problem%Neq) :: ieqs
        PetscInt, dimension(1) :: jnodes
        PetscReal :: x1, y1, z1, x2, y2, z2
        PetscBool :: check

        check = .false.
        ieqs = [ (i, i = 1, size(ieqs)) ]
        Neq = Problem%Neq

        do i = 1, size(Problem%PeriodicBNDs,1)

            iextra = 0

            ibnd = Problem%PeriodicBNDs(i,1)
            jbnd = Problem%PeriodicBNDs(i,2)

            Nnodes = size(Boundary(ibnd)%gnodes_total)

            call allocateArray(Boundary(ibnd)%gnodes_periodic, Nnodes)
            call allocateArray(Boundary(jbnd)%gnodes_periodic, Nnodes)

            BND1_loop:do inod = 1, size(Boundary(ibnd)%xi,1)

                x1 = Boundary(ibnd)%xi(inod,1)
                y1 = Boundary(ibnd)%xi(inod,2)
                z1 = Boundary(ibnd)%xi(inod,3)

                BND2_loop:do jnod = 1, size(Boundary(jbnd)%xi,1)

                    x2 = Boundary(jbnd)%xi(jnod,1)
                    y2 = Boundary(jbnd)%xi(jnod,2)
                    z2 = Boundary(jbnd)%xi(jnod,3)

                    if (Periodicity(i) == 'X') then
                        check = ( abs(y1-y2) < eps ) .and. ( abs(z1-z2) < eps )
                    end if

                    if (Periodicity(i) == 'Y') then
                        check = ( abs(x1-x2) < eps ) .and. ( abs(z1-z2) < eps )
                    end if

                    if (Periodicity(i) == 'Z') then
                        check = ( abs(x1-x2) < eps ) .and. ( abs(y1-y2) < eps )
                    end if

                    if (check) then
                        Boundary(ibnd)%gnodes_periodic(inod) = Boundary(jbnd)%gnodes_total(jnod)
                        Boundary(jbnd)%gnodes_periodic(jnod) = Boundary(ibnd)%gnodes_total(inod)
                        exit BND2_loop
                    end if

                end do BND2_loop

                ! if (.not. check) then
                !     if (Rank == 0) then
                !         open(7,file='periodic.dat',position='append')
                !         iextra = iextra+1
                !         write(7,'(a12,3i5,3es22.14)') Problem%name, ibnd, iextra, &
                !             Boundary(ibnd)%gnodes_total(inod), Boundary(ibnd)%xi(inod,:)
                !         close(7)
                !     end if
                ! end if

            end do BND1_loop

            associate(Bi => Boundary(ibnd), Bj => Boundary(jbnd))
                Bi%gnodes_periodic = pack(Bi%gnodes_periodic, Bi%gnodes_periodic /= 0)
                Bj%gnodes_periodic = pack(Bj%gnodes_periodic, Bj%gnodes_periodic /= 0)
            end associate

        end do

        ! if (Rank == 0) then
        !     write(*,*)
        !     do ibnd = 1, size(Boundary)
        !         write(*,'(a,5i5)') Boundary(ibnd)%name, ibnd, size(Boundary(ibnd)%gnodes_total), &
        !             size(Boundary(ibnd)%gnodes_periodic), iextra
        !     end do
        !     write(*,*)
        ! end if

    end subroutine setPeriodicBoundaryNodes

    ! ----------------------------------------------------------------------

    function distanceToBNDCylinder(ibnd, x) result(dx)

        use PhysicalParameters_mod, only: Length, Height, Radius, Width

        implicit none

        PetscInt, intent(in) :: ibnd
        PetscReal, dimension(:,:), intent(in) :: x
        PetscReal, dimension(size(x,1)) :: dx

        select case (ibnd)
        case (1)
            dx(:) = abs(x(:,1) + Length)
        case (2)
            dx(:) = abs(x(:,1) - Length)
        case (3)
            dx(:) = abs(x(:,2) + Height)
        case (4)
            dx(:) = abs(x(:,2) - Height)
        case (5)
            dx(:) = abs(x(:,1)**2 + x(:,2)**2 - Radius**2)
        case (6)
            dx(:) = abs(x(:,3) + Width)
            ! dx(:) = abs(x(:,3) + 0.0d0)
        case (7)
            dx(:) = abs(x(:,3) - Width)
        case default
            write(*,'(a)') 'Wrong ibnd in distanceToBNDCylinder!'
            stop
        end select

    end function distanceToBNDCylinder

    ! ----------------------------------------------------------------------

    function distanceToBNDInflow(ibnd, x) result(dx)

        use PhysicalParameters_mod, only: Height, Width

        implicit none

        PetscInt, intent(in) :: ibnd
        PetscReal, dimension(:,:), intent(in) :: x
        PetscReal, dimension(size(x,1)) :: dx

        select case (ibnd)
        case (1)
            dx(:) = abs(x(:,2) + Height)
        case (2)
            dx(:) = abs(x(:,2) - Height)
        case (3)
            dx(:) = abs(x(:,3) + Width)
        case (4)
            dx(:) = abs(x(:,3) - Width)
        case default
            write(*,'(a)') 'Wrong ibnd in distanceToBNDInflow!'
            stop
        end select

    end function distanceToBNDInflow

end module Boundary_mod
