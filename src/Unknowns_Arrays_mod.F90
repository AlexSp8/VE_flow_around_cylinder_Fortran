
#include <petsc/finclude/petscksp.h>

module UnknownsArrays_mod

    use PhysicalParameters_mod, only: ProblemParameters
    use FEMParameters_mod, only: FEParameters
    use SolutionVariables_mod, only: SolutionArraysType

    implicit none
    
    contains

    subroutine allocateUnknownsArrays(Problem, Problem_Stability, FE, Sol)

        use Tools_mod, only: allocateArray
        use ContinuationParameters_mod, only: Continuation_Method

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(SolutionArraysType), intent(inout) :: Sol

        PetscInt :: Neq, Neq_S, Nex, Nex_s, Neq_proj, Neq_proj_S
        PetscBool :: check

        Neq = Problem%Neq
        if (Neq == 0) return
        Nex = Problem%Nex
        Neq_proj = Problem%Neq_proj
        Neq_proj_S = Problem_Stability%Neq_proj
        Neq_S = Problem_Stability%Neq
        Nex_S = Problem_Stability%Nex

        select case (Problem%name)
        case('Inflow_1D', 'Inflow_2D')
            Sol%name = 'Inflow'
        case('Cylinder_2D', 'Cylinder_3D')
            Sol%name = 'Main'
        case default
            Sol%name = ''
        end select

        call allocateArray(Sol%TL, FE%Nnodes, Neq)
        call allocateArray(Sol%TLo, FE%Nnodes, Neq)
        call allocateArray(Sol%TLb, FE%Nnodes, Neq)
        call allocateArray(Sol%TLp, FE%Nnodes, Neq)

        call allocateArray(Sol%TL_Rank, FE%Nnodes_Rank, Neq)
        call allocateArray(Sol%TLo_Rank, FE%Nnodes_Rank, Neq)
        call allocateArray(Sol%TLb_Rank, FE%Nnodes_Rank, Neq)
        call allocateArray(Sol%TLp_Rank, FE%Nnodes_Rank, Neq)

        call allocateArray(Sol%EU, Nex)
        call allocateArray(Sol%EUo, Nex)
        call allocateArray(Sol%EUb, Nex)
        call allocateArray(Sol%EUp, Nex)

        call allocateArray(Sol%TL_proj, FE%Nnodes, Neq_proj)
        call allocateArray(Sol%TL_proj_Rank, FE%Nnodes_Rank, Neq_proj)

        call allocateArray(Sol%TL_b, FE%Nnodes, Neq_S)
        call allocateArray(Sol%TL_d, FE%Nnodes, Neq_S)
        call allocateArray(Sol%TL_b_Rank, FE%Nnodes_Rank, Neq_S)
        call allocateArray(Sol%TL_d_Rank, FE%Nnodes_Rank, Neq_S)
        call allocateArray(Sol%EU_b, Nex_S)
        call allocateArray(Sol%EU_d, Nex_S)

        check = (Sol%name == 'Main') .and. (Continuation_Method == 'Arclength')
        if (check) then
            call allocateUnknownsArraysArclength(Problem, FE)
        end if

    end subroutine allocateUnknownsArrays

    ! ----------------------------------------------------------------------

    subroutine allocateUnknownsArraysArclength(Problem, FE)

        use Tools_mod, only: allocateArray
        use MPIParameters_mod, only: Rank
        use SolutionVariables_mod, only: dTLodS, dEUodS, dTLodS_Rank

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(FEParameters), intent(in) :: FE

        if (Rank == 0) then
            call allocateArray(dTLodS, FE%Nnodes, Problem%Neq)
        end if
        call allocateArray(dEUodS, Problem%Nex)
        call allocateArray(dTLodS_Rank, FE%Nnodes_Rank, Problem%Neq)

    end subroutine allocateUnknownsArraysArclength

end module UnknownsArrays_mod
