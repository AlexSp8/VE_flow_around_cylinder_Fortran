
#include <petsc/finclude/petscksp.h>

module FEM_mod

    use PhysicalParameters_mod, only: ProblemParameters
    use FEMParameters_mod, only: FEParameters
    use MeshParameters_mod, only: MeshParameters

    contains

    subroutine setFEParameters(FEType, Problem, FE, Mesh)

        use FEMParameters_mod, only: FEOrder

        implicit none

        character(*), intent(in) :: FEType
        type(ProblemParameters), intent(inout) :: Problem
        type(FEParameters), intent(out) :: FE
        type(MeshParameters), intent(out) :: Mesh
        
        FE%name = FEType
        FE%order = FEOrder
        
        select case (FE%name)

        case ('Line')

            FE%Ndim = 1

            select case (FE%order)
            case ('Linear')
                FE%Nbf_1D = 2
            case ('Quadratic', 'Serendipity')
                FE%Nbf_1D = 3
            case default
                write(*,'(a)') 'Wrong FE%order in Line setFEParameters!'
                stop
            end select
            FE%Nbf_2D = 0
            FE%Nbf_3D = 0
            FE%Nbf = FE%Nbf_1D

            FE%Ned = 0

            FE%Nfc = 0

        case ('Triangle')

            FE%Ndim = 2

            select case (FE%order)
            case ('Linear')
                FE%Nbf_1D = 2
                FE%Nbf_2D = 3
            case ('Quadratic', 'Serendipity')
                FE%Nbf_1D = 3
                FE%Nbf_2D = 6
            case default
                write(*,'(a)') 'Wrong FE%order in Triangle setFEParameters!'
                stop
            end select
            FE%Nbf_3D = 0
            FE%Nbf = FE%Nbf_2D

            FE%Ned = 3

            FE%Nfc = 0

        case ('Quadrangle')

            FE%Ndim = 2

            select case (FE%order)
            case ('Linear')
                FE%Nbf_1D = 2
                FE%Nbf_2D = 4
            case ('Quadratic')
                FE%Nbf_1D = 3
                FE%Nbf_2D = 9
            case ('Serendipity')
                FE%Nbf_1D = 3
                FE%Nbf_2D = 8
            case default
                write(*,'(a)') 'Wrong FE%order in Quadrangle setFEParameters!'
                stop
            end select
            FE%Nbf_3D = 0
            FE%Nbf = FE%Nbf_2D

            FE%Ned = 4

            FE%Nfc = 0

        case ('Tetrahedron')

            FE%Ndim = 3

            select case (FE%order)
            case ('Linear')
                FE%Nbf_1D = 2
                FE%Nbf_2D = 3
                FE%Nbf_3D = 4
            case ('Quadratic', 'Serendipity')
                FE%Nbf_1D = 3
                FE%Nbf_2D = 6
                FE%Nbf_3D = 10
            case default
                write(*,'(a)') 'Wrong FE%order in Tetrahedron setFEParameters!'
                stop
            end select
            FE%Nbf = FE%Nbf_3D

            FE%Ned = 6

            FE%Nfc = 4

        case ('Hexahedron')

            FE%Ndim = 3
            
            select case (FE%order)
            case ('Linear')
                FE%Nbf_1D = 2
                FE%Nbf_2D = 4
                FE%Nbf_3D = 8
            case ('Quadratic')
                FE%Nbf_1D = 3
                FE%Nbf_2D = 9
                FE%Nbf_3D = 27
            case ('Serendipity')
                FE%Nbf_1D = 3
                FE%Nbf_2D = 8
                FE%Nbf_3D = 20
            case default
                write(*,'(a)') 'Wrong FE%order in Hexahedron setFEParameters!'
                stop
            end select
            FE%Nbf = FE%Nbf_3D

            FE%Ned = 12

            FE%Nfc = 6

        case default

            FE%Ndim = 0
            FE%Nbf_1D = 0 ; FE%Nbf_2D = 0 ; FE%Nbf_3D = 0 ; FE%Nbf = 0
            FE%Ned = 0 ; FE%Nfc = 0
            FE%Nel = 0 ; FE%Nel_Rank = 0
            FE%Nnodes = 0 ; FE%Nnodes_Rank = 0

        end select

        call readMeshFile(Problem%name, FE, Mesh)

        call setRankFEParameters(FE)

        ! if (FE%Ndim > 2) then
        !     call reformConnectivity(FE, Mesh)
        ! end if

        call setRankMeshParameters(FE, Mesh)

        Problem%Nunknowns = setFEUnknowns(FE%Nnodes, Problem%Neq)

    end subroutine setFEParameters

    ! ----------------------------------------------------------------------

    subroutine readMeshFile(mesh_name, FE, Mesh)

        implicit none

        character(*), intent(in) :: mesh_name
        type(FEParameters), intent(inout) :: FE
        type(MeshParameters), intent(out) :: Mesh

        PetscInt :: Nel, name_length, ierr

        name_length = len(trim(adjustl(mesh_name)))

        if (name_length <= 0) return
        
        open(8, file=trim(adjustl(mesh_name))//'_MESH.DTA', action='read', iostat=ierr)

        read(8,*) FE%Nnodes, Nel

        call setMeshArrays(FE%Nnodes, Mesh%Xi_Mesh)

        call setConnectivityArrays(Nel, FE, Mesh%Connectivity)

        close(8)

    end subroutine readMeshFile

    !---------------------------------------------------------------------

    subroutine setMeshArrays(Nnodes_Total, xi)

        use PhysicalParameters_mod, only: xDomain
        use Tools_mod, only: allocateArray

        implicit none

        PetscInt, intent(in) :: Nnodes_Total
        PetscReal, dimension(:,:), allocatable, intent(out) :: xi

        PetscInt :: inod, dum, j

        call allocateArray(xi, Nnodes_Total, 3)

        do inod = 1, Nnodes_Total
            read(8,*) dum, ( xi(inod,j), j = 1, size(xi,2) )
        end do

        xi(:,3) = xi(:,3)*xDomain

    end subroutine setMeshArrays

    !---------------------------------------------------------------------

    subroutine setConnectivityArrays(Nel_total, FE, Connectivity)

        use Tools_mod, only: allocateArray

        implicit none

        PetscInt, intent(in) :: Nel_total
        type(FEParameters), intent(inout) :: FE
        PetscInt, dimension(:,:), allocatable, intent(out) :: Connectivity

        PetscInt :: i, dum, iel, inod, istart, iend, idim
        PetscInt, dimension(3) :: Nel
        PetscInt, allocatable, dimension(:,:) :: dum_connectivity

        call allocateArray(dum_connectivity, Nel_total, FE%Nbf)

        Nel(:) = 0
        do i = 1, Nel_total
            read(8,*) dum, iel, ( dum_connectivity(i,inod),inod = 1, mod(iel,100) )
            idim = iel/100
            Nel(idim) = Nel(idim) + 1
        end do

        FE%Nel = Nel( FE%Ndim )

        call allocateArray(Connectivity, FE%Nel, FE%Nbf)

        istart = Nel_total-(FE%Nel)+1 ; iend = Nel_total
        Connectivity(:,:) = dum_connectivity(istart:iend,:)

        deallocate(dum_connectivity)

    end subroutine setConnectivityArrays

    !---------------------------------------------------------------------

    subroutine setRankFEParameters(FE)

        use MPIParameters_mod, only: Rank, NRanks

        implicit none

        type(FEParameters), intent(inout) :: FE

        PetscInt :: nel_missing

        FE%Nel_Rank = (FE%Nel)/NRanks
        if (Rank == NRanks-1) then
            nel_missing = (FE%Nel)-((FE%Nel_Rank)*NRanks)
            FE%Nel_Rank = (FE%Nel_Rank)+nel_missing
        end if

        FE%Nnodes_Rank = (FE%Nel_Rank)*(FE%Nbf)

    end subroutine setRankFEParameters

    !---------------------------------------------------------------------

    subroutine setRankMeshParameters(FE, Mesh)

        use Tools_mod, only: getRankNode, allocateArray
        use MPIParameters_mod, only: Rank, NRanks

        implicit none

        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(inout) :: Mesh

        PetscInt :: iel_rank, iel, Nel_Rank0, inod
        PetscInt, dimension(FE%Nbf) :: gnodes, rnodes

        call allocateArray(Mesh%Connectivity_Rank, FE%Nel_Rank, FE%Nbf)
        call allocateArray(Mesh%Xi_Mesh_Rank, FE%Nnodes_Rank, 3)

        Nel_Rank0 = (FE%Nel)/NRanks

        do iel_rank = 1, FE%Nel_Rank

            iel = iel_rank+Nel_Rank0*Rank
            Mesh%Connectivity_Rank(iel_rank,:) = Mesh%Connectivity(iel,:)

            gnodes(:) = Mesh%Connectivity_Rank(iel_rank,:)
            do inod = 1, FE%Nbf
                rnodes(inod) = getRankNode(iel_rank, inod, FE%Nbf)
            end do
            Mesh%Xi_Mesh_Rank(rnodes(:),:) = Mesh%Xi_Mesh(gnodes(:),:)

        end do

    end subroutine setRankMeshParameters

    !---------------------------------------------------------------------

    function setFEUnknowns(Nnodes, Neq) result(Nunknowns)

        implicit none

        PetscInt, intent(in) :: Nnodes, Neq
        PetscInt :: Nunknowns

        Nunknowns = Nnodes*Neq

    end function setFEUnknowns
    
    !---------------------------------------------------------------------

    subroutine reformConnectivity(FE, Mesh)

        use omp_lib, only: omp_get_wtime
        use mpi, only: MPI_Barrier, MPI_COMM_WORLD, MPI_INTEGER
        use MPIParameters_mod, only: Rank
        use MeshParameters_mod, only: RCM_Reorder

        implicit none

        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(inout) :: Mesh

        PetscReal :: tStart, tEnd
        PetscErrorCode :: ierr

        call calculateCommunications(FE, Mesh%Connectivity)

        tStart = omp_get_wtime()

        call setSurroundingArrays(FE, Mesh)
        ! call setSurroundingConnectivityArrays(FE, Mesh)

        tEnd = omp_get_wtime()
        if (Rank == 0) then
            write(*, '(a,f12.3,a)') 'Surrounding Arrays time: ', (tEnd - tStart)/60.0d0, ' min'
        end if

        tStart = omp_get_wtime()

        Mesh%Connectivity(:,:) = setConnectivityGreedy(FE, Mesh)

        tEnd = omp_get_wtime()

        if (Rank == 0) then
            write(*, '(a,f6.3,a)') 'Greedy Algorithm time: ', tEnd - tStart, ' s'
        end if

        call MPI_Bcast(Mesh%Connectivity, FE%Nel*FE%Nbf, &
            MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        call calculateCommunications(FE, Mesh%Connectivity)

        if ( RCM_Reorder .and. FE%name == 'Tetrahedron' ) then
            if (Rank == 0) call applyRCMReordering(FE, Mesh)
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
            call setNodesCoordinatesRCM(FE%Nnodes, Mesh%Xi_Mesh)
            call setConnectivityArraysRCM(FE, Mesh%Connectivity)
        end if

    end subroutine reformConnectivity

    !-----------------------------------------------------------------------

    subroutine setSurroundingArrays(FE, Mesh)

        use MPIParameters_mod, only: Rank
        use Tools_mod, only: allocateArray

        implicit none

        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(inout) :: Mesh

        PetscInt :: Nnodes, Nel, Nbf, max_Nel_s
        PetscInt :: inod, iel, gnod, gnod_el, elem_size, nodes_size, jel
        PetscInt, dimension(FE%Nbf) :: gnodes
        PetscBool :: check

        ! open(26,file='sur_data.dat')

        Nnodes = FE%Nnodes
        Nel = FE%Nel
        Nbf = FE%Nbf

        !SurN
        allocate(Mesh%SurN(Nnodes))

        max_Nel_s = 1000!Nel
        do gnod = 1, Nnodes

            call allocateArray( Mesh%SurN(gnod)%elem, Nel )
            ! call allocateArray( Mesh%SurN(gnod)%nodes, Nnodes )

            elem_size = 0 ; nodes_size = 0

            do iel = 1, Nel

                gnodes(:) = Mesh%Connectivity(iel,:)

                check = any( gnodes == gnod )
                if (check) then

                    elem_size = elem_size+1
                    Mesh%SurN(gnod)%elem(elem_size) = iel

                    ! do inod = 1, Nbf

                    !     gnod_el = Mesh%Connectivity(iel,inod)

                    !     check = any( Mesh%SurN(gnod)%nodes == gnod_el )
                    !     check = check .or. (gnod_el == gnod)
                    !     if (check) cycle

                    !     nodes_size = nodes_size+1
                    !     Mesh%SurN(gnod)%nodes(nodes_size) = gnod_el

                    ! end do

                end if

            end do

            associate(Sur => Mesh%SurN(gnod))
                Sur%elem = pack(Sur%elem, Sur%elem /= 0)
                ! Sur%nodes = pack(Sur%nodes, Sur%nodes /= 0)
            end associate

            ! if (Rank == 0) then
            !     write(26,'(*(i6))') gnod, size(Mesh%SurN(gnod)%elem), size(Mesh%SurN(gnod)%nodes)
            !     write(26,'(a,*(i6))') 'Elems: ', Mesh%SurN(gnod)%elem
            !     write(26,'(a,*(i6))') 'Nodes: ', Mesh%SurN(gnod)%nodes
            !     write(26,*) '---------------'
            ! end if

        end do

        ! close(26)

        ! !SurE
        ! allocate(Mesh%SurE(Nel))

        ! do iel = 1, Nel

        !     call allocateArray( Mesh%SurE(iel)%elem, Nel )

        !     gnodes(:) = Mesh%Connectivity(iel,:)

        !     elem_size = 0

        !     do jel = 1, Nel
        !         do inod = 1, Nbf
        !             gnod_el = Mesh%Connectivity(jel,inod)
        !             check = any( gnodes == gnod_el )
        !             if (check) then
        !                 elem_size = elem_size+1
        !                 Mesh%SurE(iel)%elem(elem_size) = jel
        !                 exit
        !             end if
        !         end do
        !     end do

        !     associate(Sur => Mesh%SurE(iel))
        !         Sur%elem = pack(Sur%elem, Sur%elem /= 0)
        !         ! Sur%nodes = pack(Sur%nodes, Sur%nodes /= 0)
        !     end associate

        ! end do

    end subroutine setSurroundingArrays

    !-----------------------------------------------------------------------

    function setConnectivityGreedy(FE, Mesh) result(Connectivity_greedy)

        use mpi
        use omp_lib, only: omp_get_wtime
        use MPIParameters_mod, only: NRanks, Rank
        use Tools_mod, only: allocateArray

        implicit none

        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        PetscInt, dimension(FE%Nel,FE%Nbf) :: Connectivity_greedy

        PetscInt :: Nel, Nnodes, Nbf, Nodes_Rank, Nel_Rank, Nel_f, Nel_f0, Nel_checked
        PetscInt :: gnod_min, gnod_max, Nel_assigned
        PetscInt :: iRank, iscore, iel_f, inod_f, iel_s, inod_s, gel, iel, iloop, inod, gnod
        PetscInt :: gel_f, gnod_f, gel_s, gnod_s, score
        PetscInt, dimension(FE%Nbf) :: gnodes_f, gnodes_s, gnodes
        PetscInt, dimension(FE%Nel) :: partition, order
        PetscInt, dimension(0:NRanks-1) :: seeds
        PetscInt, dimension(:), allocatable :: elem_f, elem_f_new
        PetscReal :: tStart, tEnd
        PetscBool :: check
        PetscBool, dimension(FE%Nel) :: assignedElem, checkedElem_s
        PetscErrorCode :: ierr

        tStart = omp_get_wtime()

        iel_f = setSeedElement(FE, Mesh)

        call MPI_Gather(iel_f, 1, MPI_INTEGER, seeds, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        if (Rank /= 0) return

        tEnd = omp_get_wtime()
        write(*, '(a,f6.3,a)') 'Seed Elements time: ', tEnd - tStart, ' s'

        Nel = FE%Nel
        Nnodes = FE%Nnodes
        Nbf = FE%Nbf

        Nodes_Rank = (Nnodes)/NRanks

        assignedElem(:) = .false.
        partition(:) = -1

        ranks_loop:do iRank = 0, NRanks-1

            tStart = omp_get_wtime()

            Nel_assigned = 1

            Nel_f = Nel_assigned

            if (allocated(elem_f)) then
                deallocate(elem_f)
            end if
            call allocateArray(elem_f, Nel_f)
            elem_f(Nel_f) = seeds(iRank)

            assignedElem(elem_f(:)) = .true.

            partition(elem_f(:)) = iRank

            Nel_Rank = FE%Nel_Rank

            gnod_min = iRank*Nodes_Rank+1
            gnod_max = (iRank+1)*Nodes_Rank
            if (iRank == NRanks-1) then
                gnod_max = Nnodes
                Nel_Rank = Nel_Rank + ( Nel-Nel_Rank*NRanks )
            end if

            iscore = Nbf

            iloop = 0
            score_loop:do

                iloop = iloop+1

                if (allocated(elem_f_new)) then
                    deallocate(elem_f_new)
                end if
                call allocateArray(elem_f_new, Nel_Rank)

                Nel_f0 = size(elem_f)
                elem_f_new(1:Nel_f0) = elem_f(:)

                Nel_f = Nel_f0

                Nel_checked = Nel_f

                checkedElem_s(:) = .false.
                ! checkedElem_s(elem_f(:)) = .true.
                iel_f_loop:do iel_f = 1, size(elem_f)

                    gel_f = elem_f(iel_f)

                    gnodes_f(:) = Mesh%Connectivity(gel_f,:)

                    inod_f_loop:do inod_f = 1, size(gnodes_f)

                        gnod_f = gnodes_f(inod_f)

                        iel_s_loop:do iel_s = 1, size(Mesh%SurN(gnod_f)%elem)

                            gel_s = Mesh%SurN(gnod_f)%elem(iel_s)

                            check = checkedElem_s(gel_s)
                            check = check .or. assignedElem(gel_s)
                            if (check) cycle iel_s_loop

                            checkedElem_s(gel_s) = .true.

                            Nel_checked = Nel_checked+1

                            gnodes_s(:) = Mesh%Connectivity(gel_s,:)

                            score = 0
                            inod_s_loop:do inod_s = 1, size(gnodes_s)

                                gnod_s = gnodes_s(inod_s)

                                check = (gnod_s < gnod_min)
                                check = check .or. (gnod_s > gnod_max)
                                if (check) cycle inod_s_loop

                                score = score+1

                            end do inod_s_loop

                            if (score >= iscore) then

                                !Assign element to Rank
                                assignedElem(gel_s) = .true.
                                partition(gel_s) = iRank
                                Nel_assigned = Nel_assigned + 1
                                if (Nel_assigned == Nel_Rank) then
                                    exit score_loop
                                end if

                                !Add element to new frontier
                                Nel_f = Nel_f+1
                                elem_f_new(Nel_f) = gel_s

                                ! if (iRank == 1 .and. iscore == 3) then
                                !     write(*,'(*(i8))') iRank, Nel_assigned, iscore, gel_f, Nel_f, Nel_checked, gel_s, score
                                !     pause
                                ! end if

                            end if

                        end do iel_s_loop

                    end do inod_f_loop

                end do iel_f_loop

                if (Nel_f == Nel_f0) then

                    iscore = iscore-1

                    !Dead end
                    check = (Nel_checked == Nel_f)
                    check = check .or. (iscore < 0)
                    if (check) then
                        iel_loop:do iel = 1, Nel

                            check = assignedElem(iel)
                            if (check) cycle iel_loop

                            gnodes(:) = Mesh%Connectivity(iel,:)

                            score = 0
                            inod_loop:do inod = 1, size(gnodes)

                                gnod = gnodes(inod)

                                check = (gnod < gnod_min)
                                check = check .or. (gnod > gnod_max)
                                if (check) cycle inod_loop

                                score = score+1

                            end do inod_loop

                            if (score >= iscore) then

                                Nel_assigned = Nel_assigned+1

                                Nel_f = 1

                                if (allocated(elem_f)) then
                                    deallocate(elem_f)
                                end if
                                call allocateArray(elem_f, Nel_f)
                                elem_f(Nel_f) = iel

                                assignedElem(elem_f(:)) = .true.

                                partition(elem_f(:)) = iRank

                                exit iel_loop

                            end if

                        end do iel_loop
                    end if

                    ! write(*,'(*(i8))') iloop, Nel_assigned, Nel_f, iscore
                    cycle score_loop

                end if

                elem_f_new = pack(elem_f_new, elem_f_new /= 0)

                if (allocated(elem_f)) then
                    deallocate(elem_f)
                end if
                call allocateArray(elem_f, size(elem_f_new))
                elem_f(:) = elem_f_new(:)

            end do score_loop

            tEnd = omp_get_wtime()
            ! write(*, '(2i8,a,f6.3,a)') iRank, Nel_assigned, ' Rank time: ', tEnd - tStart, ' s'

        end do ranks_loop

        order(:) = -1
        iel = 0
        do iRank = 0, NRanks-1
            do gel = 1, Nel
                if (partition(gel) == -1) partition(gel) = NRanks-1
                if (partition(gel) == iRank) then
                    iel = iel + 1
                    order(iel) = gel
                end if
            end do
        end do

        Connectivity_greedy(:,:) = Mesh%Connectivity(:,:)
        do iel = 1, Nel
            Connectivity_greedy(iel,:) = Mesh%Connectivity(order(iel),:)
        end do

        ! open(13, file='connectivity_new_0.dat', position='append')
        ! do iel = 1, Nel
        !     write(13,*) order(iel), '  -'!Connectivity_greedy(iel,:)
        !     ! write(13,*) Connectivity_greedy(iel,:)
        ! end do
        ! close(13)

        if (allocated(elem_f)) then
            deallocate(elem_f)
        end if
        if (allocated(elem_f_new)) then
            deallocate(elem_f_new)
        end if

    end function setConnectivityGreedy

    !-----------------------------------------------------------------------

    function setSeedElement(FE, Mesh) result(seedElement)

        use MPIParameters_mod, only: NRanks, Rank
        use Tools_mod, only: sortArray

        implicit none

        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh
        PetscInt :: seedElement

        PetscInt, parameter :: Nref = 1
        PetscInt :: Nel, Nnodes, Nbf
        PetscInt :: iel, inod, gnod, jnod, sur_gnod
        PetscInt :: nodes_Rank0, iSeed, gnod_max, gnod_min, score
        PetscInt, dimension(FE%Nel) :: score_el, order_el
        PetscBool, dimension(FE%Nel) :: assigned_el
        PetscBool :: output, check

        Nel = FE%Nel
        Nnodes = FE%Nnodes
        Nbf = FE%Nbf

        assigned_el(:) = .false.

        nodes_Rank0 = (Nnodes)/NRanks

        ! output = (FE%Ndim > 1) !.and. (Rank == 0)

        gnod_min = Rank*nodes_Rank0+1
        gnod_max = (Rank+1)*nodes_Rank0

        if (Rank == NRanks-1) then
            gnod_max = Nnodes
        end if

        score_el(:) = 0
        order_el(:) = [ (iel, iel = 1, Nel) ]
        do iel = 1, Nel
            do inod = 1, Nbf

                gnod = Mesh%Connectivity(iel,inod)

                check = (gnod >= gnod_min)
                check = check .and. (gnod <= gnod_max)
                if (check) then
                    score = 10**(Nref+2)
                    score_el(iel) = score_el(iel) + score
                end if

                score = 10**(Nref-1)

                if (score > 0) then

                    score_el(iel) = score_el(iel)
                    do jnod = 1, size(Mesh%SurN(gnod)%nodes)

                        sur_gnod = Mesh%SurN(gnod)%nodes(jnod)

                        call surroundingScore(Mesh, sur_gnod, gnod_min, gnod_max, &
                            score, score_el(iel))

                    end do
                end if

            end do
        end do

        call sortArray(score_el, order_el)

        seedElement = 1
        if ( score_el(Nel) > score_el(1) ) then

            do iel = Nel, 1, -1
                iSeed = order_el(iel)
                if ( .not. assigned_el(iSeed) ) then

                    seedElement = iSeed
                    assigned_el(iSeed) = .true.

                    ! if (output) then
                    !     iSeed = seedElement
                    !     write(*,'(*(i12))') Rank, gnod_max, &
                    !         iSeed, Mesh%Connectivity(iSeed,:), score_el(iel)
                    ! end if

                    exit

                end if
            end do

        end if

    end function setSeedElement

    !-----------------------------------------------------------------------

    recursive subroutine surroundingScore(Mesh, gnod, gnod_min, gnod_max, &
        weight, score)

        use MPIParameters_mod, only: Rank

        implicit none

        type(MeshParameters), intent(in) :: Mesh
        PetscInt, intent(in) :: gnod, gnod_min, gnod_max, weight
        PetscInt, intent(inout) :: score

        PetscInt, parameter :: decay = 10
        PetscInt :: jnod, sur_gnod
        PetscBool :: check
        
        if (weight == 0) return

        check = (gnod >= gnod_min)
        check = check .and. (gnod <= gnod_max)
        if (check) then
            score = score + weight
        end if

        ! if (Rank == 0) then
        !     print*, gnod, weight, score
        !     pause
        ! end if

        do jnod = 1, size(Mesh%SurN(gnod)%nodes)

            sur_gnod = Mesh%SurN(gnod)%nodes(jnod)

            call surroundingScore(Mesh, sur_gnod, gnod_min, gnod_max, &
                weight/decay, score)

        end do

    end subroutine surroundingScore

    !-----------------------------------------------------------------------

    subroutine calculateCommunications(FE, Connectivity)

        use mpi
        use MPIParameters_mod, only: NRanks, Rank
        use Tools_mod, only: sortFile

        implicit none

        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:,:), intent(in) :: Connectivity

        PetscInt :: iel_rank, iel, inod, gnod, Nnodes_Rank, Nnodes_off_rank
        PetscInt :: nel_Rank0, nnodes_Rank0, gnod_min, gnod_max
        PetscBool :: check

        PetscBool, save :: first_call = .true.
        PetscInt :: ierr, fh, amode
        PetscInt :: status(MPI_STATUS_SIZE)
        character(100) :: fn, line

        nel_Rank0 = (FE%Nel)/NRanks
        nnodes_Rank0 = (FE%Nnodes)/NRanks

        gnod_min = Rank*nnodes_Rank0+1
        gnod_max = (Rank+1)*nnodes_Rank0

        if (Rank == NRanks-1) then
            gnod_max = FE%Nnodes
        end if

        Nnodes_Rank = 0
        Nnodes_off_rank = 0
        do iel_rank = 1, FE%Nel_Rank
            iel = iel_rank+nel_Rank0*Rank
            do inod = 1, FE%Nbf
                gnod = Connectivity(iel,inod)
                check = (gnod >= gnod_min)
                check = check .and. (gnod <= gnod_max)
                if (check) then
                    Nnodes_Rank = Nnodes_Rank + 1
                else
                    Nnodes_off_rank = Nnodes_off_rank + 1
                end if
            end do
        end do

        write(line,'(4i8,f16.4)') Rank, FE%Nnodes_Rank, Nnodes_Rank, Nnodes_off_rank, &
            (100.0d0*Nnodes_Rank/FE%Nnodes_Rank)!, (100.0d0*Nnodes_off_rank/FE%Nnodes_Rank)
        line = trim(line) // new_line('a')

        fn = 'comm.dat'

        if (first_call) then
            amode = MPI_MODE_WRONLY + MPI_MODE_CREATE
        else
            amode = MPI_MODE_WRONLY + MPI_MODE_APPEND
        end if

        call MPI_File_open(MPI_COMM_WORLD, fn, amode, MPI_INFO_NULL, fh, ierr)

        call MPI_File_write_ordered(fh, line, len_trim(line), MPI_CHARACTER, status, ierr)

        call MPI_File_close(fh, ierr)

        first_call = .false.

    end subroutine calculateCommunications

    !-----------------------------------------------------------------------

    subroutine setSurroundingConnectivityArrays(FE, Mesh)

        use MPIParameters_mod, only: Rank

        implicit none

        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(inout) :: Mesh

        PetscInt :: Nnodes, Nel, Nbf
        PetscInt :: iel, jel, inod, jnod, knod, gnod, iNsuN0, iNsuNf, Nsun_tot, iNsuN
        PetscInt :: i, J, K, L, II, JJ, KK, LL
        PetscInt :: IDEsuN, IDNsuN, maxIDNsuN
        PetscInt, allocatable, dimension(:) :: LPN, ITEMP

        Nnodes = FE%Nnodes
        Nel = FE%Nel
        Nbf = FE%Nbf
        !-------------------------------------------------------------
        !    ELEMENTS SURROUNDING NODES
        !-------------------------------------------------------------
        !Mesh%EsuN1 --> ARRAY HOLDING ELEMENTS THAT SURROUND EACH NODE
        !Mesh%EsuN2 --> INDEXING ARRAY

        !allocate ARRAY Mesh%EsuN2 & LPN
        allocate( Mesh%EsuN2(Nnodes+1) )
        allocate( LPN(Nnodes) )
        allocate( ITEMP(Nnodes) )

        !COUNT NUMBER OF ELEMENTS CONNECTED TO EACH NODE
        Mesh%EsuN2(:) = 0

        do iel = 1, Nel
            do inod = 1, Nbf
                gnod = Mesh%Connectivity(iel,inod)
                Mesh%EsuN2(gnod) = Mesh%EsuN2(gnod) + 1
            end do
        end do

        !RESHUFFLING
        do i = 2, Nnodes+1
            Mesh%EsuN2(i) = Mesh%EsuN2(i) + Mesh%EsuN2(i-1)
        end do

        !DEFINE IDEsuN
        IDEsuN = Mesh%EsuN2(Nnodes+1)

        !RESHUFFLING
        do i = Nnodes+1, 2, -1
            Mesh%EsuN2(i) = Mesh%EsuN2(i-1)
        end do
        Mesh%EsuN2(1) = 0

        !allocate ARRAY Mesh%EsuN1
        allocate(Mesh%EsuN1(IDEsuN))

        Mesh%EsuN1(:) = 0

        !STORE THE ELEMENTS IN ARRAY Mesh%EsuN1
        do iel = 1, Nel
            do inod = 1, Nbf
                gnod = Mesh%Connectivity(iel,inod)
                jel = Mesh%EsuN2(gnod)
                LOOP_IN: do 
                    jel = jel + 1
                    if (Mesh%EsuN1(jel) == 0) THEN
                        Mesh%EsuN1(jel) = iel
                    EXIT LOOP_IN
                    end if 
                end do LOOP_IN
            end do
        end do

        !-------------------------------------------------------------
        !    NODES SURROUNDING NODES
        !-------------------------------------------------------------
        !Mesh%NsuN1 --> ARRAY HOLDING NODES THAT SURROUND EACH NODE
        !Mesh%NsuN2 --> INDEXING ARRAY

        !DEFINE maxIDNsuN
        maxIDNsuN = (Nbf-1)*IDEsuN

        !allocate ARRAY Mesh%NsuN2
        allocate(Mesh%NsuN2(Nnodes+1))
        allocate(Mesh%NsuN1(maxIDNsuN))

        if (Rank == 0) then
            print*, IDEsuN, maxIDNsuN, Nnodes, Nel
        end if

        !INITIALIZATION
        Mesh%NsuN1(:) = 0
        Mesh%NsuN2(:) = 0

        !ASSEMBLE NsuN
        JJ = 0
        LOOP_1: do i = 1, Nnodes

            LPN = 0
            !For all elements that surround node i
            LOOP_2: do jel = Mesh%EsuN2(i)+1, Mesh%EsuN2(i+1)

                iel = Mesh%EsuN1(jel)

                ! LOOP_3: do inod = Nbf, 1, -1
                    LOOP_3: do inod = 1, NBF

                    gnod = Mesh%Connectivity(iel,inod)
                    knod = gnod
                    if (LPN(gnod) == 0) then
                        JJ = JJ + 1
                        Mesh%NsuN1(JJ) = knod
                        ! if (knod == i) then
                        !         if (Rank == 0) then
                        !             print*, i, JJ, iel
                        !             print *, Mesh%NsuN1(Mesh%NsuN2(i)+1)
                        !         endif
                        ! endif
                        LPN(gnod) = inod
                    end if

                end do LOOP_3

            end do LOOP_2

            Mesh%NsuN2(i+1) = JJ
            ! if (i /= Mesh%NsuN1(Mesh%NsuN2(i)+1)) then
            !     print*, 'MY MESH IS FUCKED UP! POURYA WAS RIGHT!'
            !     stop
            ! end if

            ! if (i == 3) then
            !     !NSUN1 returns global number, NSUN2 returns indexing
            !     if (Rank == 0) then
            !         iNsuN0 = Mesh%NsuN2(i)+1
            !         iNsuNf = Mesh%NsuN2(i+1)
            !         Nsun_tot = iNsuNf - iNsuN0
            !         print*, "total surrounding nodes", Nsun_tot!, iNsuNf, iNsuN0
            !         print*, Mesh%NsuN1(iNsuN0), i ! ALEX NODE NUMBER has to be equal to i
            !         do iNsuN = 1, Nsun_tot
            !             print *, iNsuN, Mesh%NsuN1(iNsuN0 + iNsuN)
            !         end do
            !     end if
            ! end if

        end do LOOP_1

        !DEFINE IDNsuN
        IDNsuN = Mesh%NsuN2(Nnodes+1)

        !RESHUFFLING
        LOOPNODES: do I = 1, Nnodes
            L = 0
            do J = Mesh%NsuN2(I)+1, Mesh%NsuN2(I+1)
                L = L + 1
                ITEMP(L) = Mesh%NsuN1(J)
            end do
            LL = 0
            do J = Mesh%NsuN2(I)+1, Mesh%NsuN2(I+1)
                LL = LL + 1
                jnod = 0
                jnod = MINVAL(ITEMP(1:L), MASK = ITEMP(1:L) > 0)
                do KK = 1, L
                    if (jnod == ITEMP(KK)) then
                        ITEMP(KK) = 0
                    end if
                end do
                Mesh%NsuN1(J) = jnod
            end do
        end do LOOPNODES

    end subroutine setSurroundingConnectivityArrays

    !-----------------------------------------------------------------------

    subroutine applyRCMReordering(FE, Mesh)

        use OMP_LIB, only: omp_get_wtime
        use Tools_mod, only: writeElapsedTime
        use Reverse_Cuthill_McKee_Node_Reordering, only: Reverce_CM

        implicit none

        type(FEParameters), intent(in) :: FE
        type(MeshParameters), intent(in) :: Mesh

        PetscInt  :: inod, idim, iel
        PetscReal :: tStart, tEnd

        tStart = omp_get_wtime()

        open(31, file='Salome_nodes.txt')
        open(32, file='Salome_elements.txt')

        do inod = 1, FE%Nnodes
            write(31,*) ( Mesh%Xi_Mesh(inod,idim), idim = 1, size(Mesh%Xi_Mesh,2) )
        end do

        do iel = 1, FE%Nel
            write(32,'(*(i8))') ( Mesh%Connectivity(iel,inod),inod = 1, FE%Nbf )
        end do

        close(31)
        close(32)

        !RCM (Reverse Cuthill-McKee): turn a sparse matrix into banded with reversed index numbers
        call Reverce_CM()

        tEnd = omp_get_wtime()

        call writeElapsedTime(tStart, tEnd, 'RCM time')

    end subroutine applyRCMReordering

    !----------------------------------------------------------------------- 

    subroutine setNodesCoordinatesRCM(Nnodes, xi)

        implicit none

        PetscInt, intent(in) :: Nnodes
        PetscReal, dimension(:,:), intent(out) :: xi

        PetscInt :: inod, j

        open(31, file='Salome_rcm_nodes.txt')
        do inod = 1, Nnodes
            read(31,*) ( xi(inod,j), j = 1, size(xi,2) )
        end do
        close(31)

    end subroutine setNodesCoordinatesRCM

    ! ----------------------------------------------------------------------

    subroutine setConnectivityArraysRCM(FE, Connectivity)

        implicit none

        type(FEParameters), intent(in) :: FE
        PetscInt, dimension(:,:), intent(out) :: Connectivity

        PetscInt :: inod, iel

        open(34, file='Salome_rcm_elements.txt')

        do iel = 1, FE%Nel
            read(34,*) ( Connectivity(iel,inod), inod = 1, FE%Nbf )
        end do

        close(34)

    end subroutine setConnectivityArraysRCM

end module FEM_mod
