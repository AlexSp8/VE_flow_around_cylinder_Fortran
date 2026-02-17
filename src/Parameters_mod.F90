
#include <petsc/finclude/petscksp.h>

module ContinuationParameters_mod

    implicit none

    !Arclength, Stability, Transient
    ! character(*), parameter :: Continuation_Method = 'Transient'
    character(*), parameter :: Continuation_Method = 'Arclength'

    !WiN, ReN, Time, BvN, VeN, DfN
    ! character(*), parameter :: Cvar1_Name = 'Time'
    character(*), parameter :: Cvar1_Name = 'WiN'

    !BvN, VeN, ReN, DfN
    ! character(*), dimension(3), parameter :: Cvar2_Names = ['BvN', 'VeN', 'ReN']
    ! character(*), dimension(1), parameter :: Cvar2_Names = ['BvN']
    character(*), dimension(0), parameter :: Cvar2_Names = ''

    PetscBool, parameter :: adjust_dt = .false.

    ! PetscReal, parameter :: Dx_write = 1.0d0
    PetscReal, parameter :: Dx_write = -0.2d0

    PetscReal, parameter :: LN0 = 0.668713921947d0, Eps_Pr = 1.0d-3, dArN = -1.0d-3, &
        dS0 = 1.0d0, dS_max = 50.0d0

    PetscReal, parameter :: Cont_var0_Stability = +huge(1.0d0)
    PetscReal, parameter :: Stability_step = 0.2d0!*Cont_var0_Stability
    PetscBool, parameter :: Stability_Im = .false.

end module ContinuationParameters_mod

!----------------------------------------------------------------------------------------

module PhysicalParameters_mod

    implicit none

    ! !2D test
    ! PetscReal, parameter :: PI = 4.0d0*atan(1.0d0)

    ! PetscReal, parameter :: Time0 = 0.0d0, dt0 = 0.1d0, Time_f = 1.0d6
    ! PetscReal, target :: ReN, WiN, Total_time, VeN, BvN, DfN

    ! !L-PTT, e-PTT, m-L-PTT, m-e-PTT, FENE-CR, FENE-P, GIESEKUS
    ! character(*), parameter :: Model_Name = 'L-PTT'

    ! PetscReal, parameter :: ReN0 = 0.0d0, ReN_f = 0.02d0, dReN = 0.01d0
    ! PetscReal, parameter :: WiN0 = 0.05d0, WiN_f = 35.0d0, dWiN = 0.5d0
    ! ! PetscReal, parameter :: WiN0 = 23.55d0, WiN_f = 25.0d0, dWiN = 0.5d0
    ! PetscReal, parameter :: BvN0 = 0.05d0, BvN_f = 0.052d0, dBvN = 0.001d0
    ! PetscReal, parameter :: VeN0 = 0.05d0, VeN_f = 0.052d0, dVeN = 0.001d0
    ! PetscReal, parameter :: DfN0 = 0.0d0, DfN_f = 0.0d0, dDfN = 1.0d-4
    ! PetscReal, parameter :: QbN0 = 0.0d0, QbN_f = 0.0d0, dQbN = 0.05d0
    ! PetscReal :: QbN

    ! !Domain
    ! PetscInt, parameter :: xDomain = 1
    ! PetscReal, parameter :: BR = 0.5d0, Radius = 1.0d0
    ! PetscReal, parameter :: Height = Radius/BR, Length = 12.5d0*Height, Width = 0.5d0!1.0d0
    ! !SQRT
    ! character(*), parameter :: Stress_Reform = 'SQRT'

    ! PetscReal, parameter :: bslip = 1.0d0

    ! PetscReal, parameter :: LwN = 2.38d0/3.175d0, KwN = (2.0d0*PI)/LwN

    ! !Inflow_1D, Inflow_2D
    ! character(*), parameter :: Inflow_Name = 'Inflow_1D'
    ! !1D: 4, 2D: 9
    ! PetscInt, parameter :: Inflow_Neq_proj = 0
    ! !Inflow_Stability_1D, Inflow_Stability_1D_2D, Inflow_Stability_2D
    ! character(*), parameter :: Inflow_Name_Stability = ''
    ! PetscInt, parameter :: Inflow_Neq_proj_S = 0

    ! !Cylinder_2D, Cylinder_3D
    ! character(*), parameter :: Main_Name = 'Cylinder_2D'
    ! !2D: 4, 3D: 9
    ! PetscInt, parameter :: Main_Neq_proj = 0
    ! !Cylinder_Stability_2D, Cylinder_Stability_2D_3D, Cylinder_Stability_3D
    ! character(*), parameter :: Main_Name_Stability = 'Cylinder_Stability_2D'
    ! PetscInt, parameter :: Main_Neq_proj_S = 0

    ! PetscInt, parameter :: Nbd_symmetric = 0
    ! character(*), parameter, dimension(Nbd_symmetric):: Symmetry = ''

    ! PetscInt, parameter :: Nbd_periodic = 0
    ! character(*), parameter, dimension(Nbd_periodic) :: Periodicity  = ''
    ! PetscInt, parameter, dimension(Nbd_periodic,2) :: BNDs_p_Inflow = 0![3,4]
    ! PetscInt, parameter, dimension(Nbd_periodic,2) :: BNDs_p_Main  = 0![6,7]

    !3D stability
    PetscReal, parameter :: PI = 4.0d0*atan(1.0d0)

    PetscReal, parameter :: Time0 = 0.0d0, dt0 = 0.1d0, Time_f = 1.0d6
    PetscReal, target :: ReN, WiN, Total_time, VeN, BvN, DfN, dum

    !L-PTT, e-PTT, m-L-PTT, m-e-PTT, FENE-CR, FENE-P, GIESEKUS
    character(*), parameter :: Model_Name = 'L-PTT'

    PetscReal, parameter :: ReN0 = 0.0d0, ReN_f = 4.0d0, dReN = 0.2d0
    ! PetscReal, parameter :: WiN0 = 0.1d0, WiN_f = 4.0d0, dWiN = 0.1d0
    PetscReal, parameter :: WiN0 = 3.0d0, WiN_f = +huge(1.0d0), dWiN = -0.02d0
    PetscReal, parameter :: BvN0 = 0.1d0, BvN_f = 0.8d0, dBvN = -0.01d0
    ! PetscReal, parameter :: BvN0 = 9.22d0/13.75d0, BvN_f = 0.8d0, dBvN = 0.05d0
    ! PetscReal, parameter :: BvN0 = 0.59d0, BvN_f = 0.8d0, dBvN = 0.05d0
    ! PetscReal, parameter :: VeN0 = 0.01d0, VeN_f = 0.05d0, dVeN = -0.0002d0
    PetscReal, parameter :: VeN0 = 1.0d0/150.0d0, VeN_f = 0.03d0, dVeN = 0.005d0
    PetscReal, parameter :: DfN0 = 1.0d-3, DfN_f = 0.0d0, dDfN = 1.0d-4
    PetscReal, parameter :: QbN0 = 0.0d0, QbN_f = 0.0d0, dQbN = 0.05d0
    PetscReal :: QbN

    !Domain
    PetscInt, parameter :: xDomain = 1
    PetscReal, parameter :: BR = 0.5d0, Radius = 1.0d0
    PetscReal, parameter :: Height = Radius/BR, Length = 12.5d0*Height, &
        Width = 0.5d0*xDomain/2

    !SQRT
    character(*), parameter :: Stress_Reform = ''

    PetscReal, parameter :: bslip = 1.0d0

    PetscReal, parameter :: LwN = 2.38d0/3.175d0, KwN = (2.0d0*PI)/LwN

    !Inflow_1D, Inflow_2D
    character(*), parameter :: Inflow_Name = 'Inflow_2D'
    !1D: 4, 2D: 9
    PetscInt, parameter :: Inflow_Neq_proj = 0
    !Inflow_Stability_1D, Inflow_Stability_1D_2D, Inflow_Stability_2D
    character(*), parameter :: Inflow_Name_Stability = ''
    PetscInt, parameter :: Inflow_Neq_proj_S = 0

    !Cylinder_2D, Cylinder_3D
    character(*), parameter :: Main_Name = 'Cylinder_3D'
    !2D: 4, 3D: 9
    PetscInt, parameter :: Main_Neq_proj = 0
    !Cylinder_Stability_2D, Cylinder_Stability_2D_3D, Cylinder_Stability_3D
    character(*), parameter :: Main_Name_Stability = ''
    PetscInt, parameter :: Main_Neq_proj_S = 0

    PetscInt, parameter :: Nbd_symmetric = 1
    character(*), parameter, dimension(Nbd_symmetric):: Symmetry = 'Z'

    ! PetscInt, parameter :: Nbd_periodic = 1
    ! character(*), parameter, dimension(Nbd_periodic) :: Periodicity  = 'Z'
    ! PetscInt, parameter, dimension(Nbd_periodic,2) :: BNDs_p_Inflow = [3,4]
    ! PetscInt, parameter, dimension(Nbd_periodic,2) :: BNDs_p_Main  = [6,7]
    PetscInt, parameter :: Nbd_periodic = 0
    character(*), parameter, dimension(Nbd_periodic) :: Periodicity  = ''
    PetscInt, parameter, dimension(Nbd_periodic,2) :: BNDs_p_Inflow = 0![3,4]
    PetscInt, parameter, dimension(Nbd_periodic,2) :: BNDs_p_Main  = 0![6,7]

    type ProblemParameters

        character(30) :: name
        PetscInt :: Ndim, Neq_f, Neq_s, Neq, Nex_b, Nex_ed, Nex_fc, Nex, Neq_proj
        PetscInt :: Nbd, Nunknowns, Nunknowns_proj
        PetscInt, allocatable, dimension(:) :: iex_b_values
        PetscInt, allocatable, dimension(:,:) :: PeriodicBNDs
        PetscInt, allocatable, dimension(:) :: idir

    end type ProblemParameters

    type(ProblemParameters) :: Problem_Inflow, Problem_Inflow_Stability, &
        Problem_Main, Problem_Main_Stability

end module PhysicalParameters_mod

!----------------------------------------------------------------------------------------

module FEMParameters_mod

    implicit none

    !Linear, Quadratic (tetra Salome), Serendipity (visualization)
    character(*), parameter :: FEOrder = 'Linear'
    !Line, Triangle, Quadrangle
    ! character(*), parameter :: FEType_Inflow = 'Line'
    character(*), parameter :: FEType_Inflow = 'Quadrangle'
    !Triangle, Quadrangle, Tetrahedron, Hexahedron
    ! character(*), parameter :: FEType_Main = 'Quadrangle'
    character(*), parameter :: FEType_Main = 'Hexahedron'

    type FEParameters

        character(20) :: name, order
        PetscInt :: Ndim
        PetscInt :: Nbf_1D, Nbf_2D, Nbf_3D, Nbf
        PetscInt :: Ned, Nfc
        PetscInt :: Nel, Nel_Rank
        PetscInt :: Nnodes, Nnodes_Rank

    end type FEParameters

    type(FEParameters) :: FE_Inflow, FE_Main

end module FEMParameters_mod

!----------------------------------------------------------------------------------------

module NewtonRaphsonParameters_mod

    implicit none

    PetscBool, parameter :: Iterative_Solver = .false., Monitor_Convergence = .true.
    !KSPGMRES, KSPFGMRES, KSPPREONLY
    character(*), parameter :: Solver_Name = 'KSPFGMRES'
    !PCFIELDSPLIT, PCASM, PCLU
    character(*), parameter :: PC_Name = 'PCFIELDSPLIT'
    PetscInt, parameter :: Niter = 50, LS_maxit = 1000
    PetscReal, parameter :: Cor_Norm_eps = 1.0d-8, MNR_eps = 1.0d-1, &
        Cor_Norm_div = 1.0d8, LS_eps = 1.0d-12

end module NewtonRaphsonParameters_mod

!----------------------------------------------------------------------------------------

module MeshParameters_mod

    implicit none

    PetscBool, parameter :: RCM_Reorder = .false.
    PetscReal, parameter :: Mesh_Tolerance = 1.0d-6

    type SurroundingArrays
        PetscInt, allocatable, dimension(:) :: elem, nodes
    end type SurroundingArrays

    type MeshParameters
        PetscInt,  allocatable, dimension(:,:) :: Connectivity, Connectivity_Rank
        PetscReal, allocatable, dimension(:,:) :: Xi_Mesh, Xi_Mesh_Rank
        PetscInt,  allocatable, dimension(:)   :: EsuN1, NsuN1, EsuN2, NsuN2
        type(SurroundingArrays), allocatable, dimension(:) :: SurN, SurE
    end type MeshParameters

    type(MeshParameters) :: Mesh_Inflow, Mesh_Main

end module MeshParameters_mod

!----------------------------------------------------------------------------------------

module BoundaryParameters_mod

    implicit none

    type BoundaryParameters

        character(len=20) :: name

        PetscInt, allocatable, dimension(:) :: eq_ed, ex_ed
        PetscReal, allocatable, dimension(:) :: values_ex_ed

        PetscInt, allocatable, dimension(:) :: eq_fc, ex_fc
        PetscReal, allocatable, dimension(:) :: values_ex_fc

        PetscInt, allocatable, dimension(:) :: eq_clear

        PetscInt, allocatable, dimension(:) :: eq_Dirichlet
        PetscReal, allocatable, dimension(:) :: values_Dirichlet

        PetscInt, allocatable, dimension(:) :: epp_ed, epp_fc

        PetscInt, allocatable, dimension(:) :: eq_ed_S, ex_ed_S
        PetscReal, allocatable, dimension(:) :: values_ex_ed_S

        PetscInt, allocatable, dimension(:) :: eq_fc_S, ex_fc_S
        PetscReal, allocatable, dimension(:) :: values_ex_fc_S

        PetscInt, allocatable, dimension(:) :: eq_clear_S

        PetscInt, allocatable, dimension(:) :: eq_Dirichlet_S
        PetscReal, allocatable, dimension(:) :: values_Dirichlet_S

        PetscInt, allocatable, dimension(:) :: epp_ed_S, epp_fc_S

        PetscInt, allocatable, dimension(:) :: gnodes_rank, rnodes, &
                                               edges, edge_elements_rank, &
                                               faces, face_elements_rank!, node_elements_rank

        PetscInt, allocatable, dimension(:) :: gnodes_total, gnodes_periodic
        PetscReal, allocatable, dimension(:,:) :: xi

    end type BoundaryParameters

    type(BoundaryParameters), dimension(:), allocatable :: Boundary_Main, Boundary_Inflow

end module BoundaryParameters_mod

!----------------------------------------------------------------------------------------

module GaussParameters_mod

    implicit none

    type GaussIntegration

        character(20) :: order

        PetscInt :: NGauss_1D, NGauss_2D, NGauss_3D, NGauss_LC, NGauss_b
        PetscInt :: Nbf, Ned, Nfc

        PetscReal, dimension(:), allocatable :: Weights_1D
        PetscReal, dimension(:), allocatable :: Weights_2D
        PetscReal, dimension(:), allocatable :: Weights_3D
        PetscReal, dimension(:), allocatable :: Weights_LC
        PetscReal, dimension(:), allocatable :: Weights_b

        PetscReal, dimension(:),   allocatable :: Points_1D
        PetscReal, dimension(:,:), allocatable :: Points_2D
        PetscReal, dimension(:,:), allocatable :: Points_3D
        PetscReal, dimension(:),   allocatable :: Points_LC

        PetscReal, dimension(:,:,:), allocatable :: Bfn_E, Dfdc_E, Dfde_E, Dfds_E, &
                                                    D2fdc2_E, D2fde2_E, D2fds2_E, &
                                                    D2fdce_E, D2fdcs_E, D2fdes_E

        PetscReal, dimension(:,:,:), allocatable :: Bfn_F, Dfdc_F, Dfde_F, Dfds_F, &
                                                    D2fdc2_F, D2fde2_F, D2fds2_F, &
                                                    D2fdce_F, D2fdcs_F, D2fdes_F

        PetscReal, dimension(:,:), allocatable :: Bfn, Dfdc, Dfde, Dfds, &
                                                    D2fdc2, D2fde2, D2fds2, &
                                                    D2fdce, D2fdcs, D2fdes

        PetscReal, dimension(:,:,:), allocatable :: dfdci
        PetscReal, dimension(:,:,:,:), allocatable :: dfdci_E, dfdci_F

    end type GaussIntegration

    type(GaussIntegration) :: GaussInt_Inflow, GaussInt_Main

end module GaussParameters_mod

!----------------------------------------------------------------------------------------

module MPIParameters_mod

    implicit none

    PetscMPIInt :: Rank, NRanks

end module MPIParameters_mod

!----------------------------------------------------------------------------------------

module SlepcParameters_mod

    implicit none

    !nev: number of requested eigenvalues to compute: default nev = 1
    !maxit: max iterations (default 100)
    PetscInt, parameter :: Nshifts = 1, nev = 50, maxit = 100
    !tolerance (default 1.0d-8)
    PetscReal, parameter :: EVP_tol = 1.0d-10

end module SlepcParameters_mod
