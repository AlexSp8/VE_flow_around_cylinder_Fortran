
#include <petsc/finclude/petscksp.h>
#include <slepc/finclude/slepceps.h>

module ContinuationVariables_mod

    implicit none

    PetscInt  :: Increment
    PetscReal :: dt
    PetscReal :: dSo

    type ContVarType
        PetscReal, pointer :: p
        PetscReal :: dvar, var_f, var_o, var0
        character(20) :: name
    end type ContVarType
    type(ContVarType) :: Cvar1, Cvar2

end module ContinuationVariables_mod

!----------------------------------------------------------------------------------------

module SolutionVariables_mod

    implicit none

    type SolutionArraysType
        character(20) :: name
        PetscReal, allocatable, dimension(:,:) :: TL, TLo, TLb, TLp
        PetscReal, allocatable, dimension(:)   :: EU, EUo, EUb, EUp
        PetscReal, allocatable, dimension(:,:) :: TL_Rank, TLo_Rank, TLb_Rank, TLp_Rank
        PetscReal, allocatable, dimension(:,:) :: TL_proj, TL_proj_Rank
        PetscReal, allocatable, dimension(:,:) :: TL_b
        PetscScalar, allocatable, dimension(:,:) :: TL_d
        PetscReal, allocatable, dimension(:)   :: EU_b
        PetscScalar, allocatable, dimension(:) :: EU_d
        PetscReal, allocatable, dimension(:,:) :: TL_b_Rank
        PetscScalar, allocatable, dimension(:,:) :: TL_d_Rank
    end type SolutionArraysType

    type(SolutionArraysType) :: Sol_Inflow, Sol_Main

    PetscReal, allocatable, dimension(:)   :: dEUodS
    PetscReal, allocatable, dimension(:,:) :: dTLodS
    PetscReal, allocatable, dimension(:,:) :: dTLodS_Rank

end module SolutionVariables_mod

!----------------------------------------------------------------------------------------

module ElementVariables_mod

    implicit none

    type GaussPointQuantities
        
        PetscScalar, dimension(3,3) :: GU, GUT, G_dot, GU_tr, GUT_tr, G_dot_tr, GU_proj
        PetscScalar, dimension(3,3) :: P_tot, G2U
        PetscScalar, dimension(3) :: U_f, U_fo, U_fb, dPdXi, dUfdt
        PetscScalar, dimension(3) :: Momentum, Convection, Momentum_time
        PetscScalar :: Continuity, P

        PetscScalar, dimension(3,3,3) :: dTdXi, dCdXi, dSdXi, dGU_proj
        PetscScalar, dimension(3,3) :: S, So, Sb, dSdt
        PetscScalar, dimension(3,3) :: Con_EQ, T_ve, Con_EQ_time, C_Tensor
        PetscScalar, dimension(3) :: f
        PetscScalar, dimension(2) :: g
        PetscScalar, dimension(2,3) :: dgdxi
        PetscScalar :: f_T

        PetscReal, dimension(3) :: normal, tangent, tangent2
        PetscReal :: dL

        PetscReal, dimension(3,3) :: dXidci
        PetscReal, dimension(3) :: Xi

        PetscReal :: detJ, WET
        ! PetscReal :: X0, dX0dC, dX0dE, Y0, dY0dC, dY0dE, CJAC0, AJAC0

        PetscScalar :: mag_con_EQ, mag_con_EQ_time
        PetscScalar :: tlsic, tlsme, tlsce, tdcs

    end type GaussPointQuantities

    type NodeArrays

        PetscReal, allocatable, dimension(:,:) :: TEMP_TL, TEMP_TLo, TEMP_TLb
        PetscReal, allocatable, dimension(:)   :: TEMP_EU

        PetscReal, allocatable, dimension(:,:) :: TEMP_TL_proj

        PetscReal, allocatable, dimension(:,:) :: TEMP_TL_b
        PetscReal, allocatable, dimension(:)   :: TEMP_EU_b

        PetscScalar, allocatable, dimension(:,:) :: TEMP_TL_d
        PetscScalar, allocatable, dimension(:) :: TEMP_EU_d

        PetscReal, allocatable, dimension(:,:) :: dFdXi, d2FdXi2
        PetscReal, allocatable, dimension(:,:) :: Xi_loc

        PetscReal, dimension(:), pointer :: bfn_p, dfdc_p, dfde_p, dfds_p
        ! PetscReal, dimension(:,:), pointer :: dfdci_p, d2fdci2_p
        PetscReal, dimension(:), pointer :: d2fdc2_p, d2fde2_p, d2fds2_p, &
            d2fdce_p, d2fdcs_p, d2fdes_p

        ! PetscReal, allocatable, dimension(Neq) :: Jump
        ! PetscReal, allocatable, dimension(:,:) :: Xi0_loc, dFdXi0

    end type NodeArrays

    type(NodeArrays) :: NodeArrays_Inflow, NodeArrays_Main

end module ElementVariables_mod

!----------------------------------------------------------------------------------------

module NewtonRaphsonVariables_mod

    type NewtonType

        character(20) :: name
        PetscInt :: Iter_f, xNiter
        PetscReal :: Res_Norm, Cor_Norm_Old_f, Cor_Norm_New_f, xF, Cor_Norm0
        character(len=3) :: Flag_NR
        PetscBool :: Emergency

    end type NewtonType

end module NewtonRaphsonVariables_mod

!----------------------------------------------------------------------------------------

module LinearSystemVariables_mod

    use petscksp

    type LinearSystemType

        character(20) :: name
        Mat :: A_f, A_f_factored
        Vec :: b_f, s_f
        KSP :: ksp
        PC  :: pc
        PetscInt :: maxit, iter
        PetscBool :: reCreate

    end type LinearSystemType

    type(LinearSystemType) :: LS_Inflow, LS_Main, LS_Inflow_proj, LS_Main_proj

end module LinearSystemVariables_mod

!----------------------------------------------------------------------------------------

module SlepcVariables_mod

    use slepceps
    
    implicit none

    type EVPType

        KSP :: ksp
        EPS :: eps
        Mat :: A_f, B_f
        Vec :: Xr
        PetscScalar :: K_Lead
        PetscReal :: target_real, Kr_Leado, Ki_max, xCr, x0, LwN_Lead
        PetscBool :: critical

    end type EVPType

    type(EVPType) :: EVP_Inflow, EVP_Main

end module SlepcVariables_mod
