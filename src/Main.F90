! ----------------------------------------------------------------------
!   PROGRAM FOR SOLVING SYSTEMS OF PDEs IN A CARTESIAN DOMAIN X,Y,Z
! ----------------------------------------------------------------------
!   PROGRAM WRITTEN BY J.A.T. & A.N.B. (1985, MIT    ) 
!   REVISED BY Y.D.                    (2013, UPATRAS)
!   REVISED BY S.E.V.                  (2018, UPATRAS)
!   REMIXED & REMASTERED BY A.S.       (2024, UPATRAS)
! ----------------------------------------------------------------------

#include <petsc/finclude/petscksp.h>

program FEM2D
    
    use petscksp
    use omp_lib, only: omp_get_wtime

    use ContinuationParameters_mod, only: Continuation_Method
    use PhysicalParameters_mod, only: &
        Inflow_Name, Inflow_Neq_proj, BNDs_p_Inflow, Problem_Inflow, &
        Inflow_Name_Stability, Inflow_Neq_proj_S, Problem_Inflow_Stability, &
        Main_Name, Main_Neq_proj, BNDs_p_Main, Problem_Main, &
        Main_Name_Stability, Main_Neq_proj_S, Problem_Main_Stability
    use FEMParameters_mod, only: FEType_Inflow, FE_Inflow, FEType_Main, FE_Main
    use MeshParameters_mod, only: Mesh_Inflow, Mesh_Main
    use BoundaryParameters_mod, only: Boundary_Inflow, Boundary_Main
    use GaussParameters_mod, only: GaussInt_Inflow, GaussInt_Main

    use Physical_mod, only: setProblemParameters, setDimensionlessNumbers
    use FEM_mod, only: setFEParameters, setFEUnknowns
    use Boundary_mod, only: setProblemBoundaries
    use Dirichlet_mod, only: connectInflowToMainArrays
    use Gauss_mod, only: setGaussParameters

    use LinearSystemVariables_mod, only: LS_Inflow, LS_Main, LS_Inflow_proj, LS_Main_proj
    use Petsc_mod, only: initializePetsc, finalizePetsc, writeMemoryUsage, &
        createPetscObjects, createPetscObjectsProjection, createSlepcObjects

    use UnknownsArrays_mod, only: allocateUnknownsArrays
    use SolutionVariables_mod, only: Sol_Inflow, Sol_Main

    use SlepcVariables_mod, only: EVP_Inflow, EVP_Main

    use ElementVariables_mod, only: NodeArrays_Inflow, NodeArrays_Main
    use ElementCalculations_mod, only: allocateLocalArrays

    use Continuation_mod, only: loopContinuation

    implicit none

    PetscReal :: tStart
    PetscErrorCode :: ierr

    tStart = omp_get_wtime()

    call initializePetsc(Continuation_Method)

    call setProblemParameters(Inflow_Name, Inflow_Neq_proj, BNDs_p_Inflow, Problem_Inflow)
    call setProblemParameters(Inflow_Name_Stability, Inflow_Neq_proj_S, BNDs_p_Inflow, &
        Problem_Inflow_Stability)

    call setProblemParameters(Main_Name, Main_Neq_proj, BNDs_p_Main, Problem_Main)
    call setProblemParameters(Main_Name_Stability, Main_Neq_proj_S, BNDs_p_Main, Problem_Main_Stability)

    call setFEParameters(FEType_Inflow, Problem_Inflow, FE_Inflow, Mesh_Inflow)
    call setFEParameters(FEType_Main, Problem_Main, FE_Main, Mesh_Main)

    Problem_Inflow%Nunknowns_proj = setFEUnknowns(FE_Inflow%Nnodes, Problem_Inflow%Neq_proj)
    Problem_Inflow_Stability%Nunknowns = setFEUnknowns(FE_Inflow%Nnodes, Problem_Inflow_Stability%Neq)
    Problem_Inflow_Stability%Nunknowns_proj = setFEUnknowns(FE_Inflow%Nnodes, &
        Problem_Inflow_Stability%Neq_proj)

    Problem_Main%Nunknowns_proj = setFEUnknowns(FE_Main%Nnodes, Problem_Main%Neq_proj)
    Problem_Main_Stability%Nunknowns = setFEUnknowns(FE_Main%Nnodes, Problem_Main_Stability%Neq)
    Problem_Main_Stability%Nunknowns_proj = setFEUnknowns(FE_Main%Nnodes, Problem_Main_Stability%Neq_proj)

    call setProblemBoundaries(Problem_Inflow, Problem_Inflow_Stability, &
        FE_inflow, Mesh_Inflow, Boundary_Inflow)
    call setProblemBoundaries(Problem_Main, Problem_Main_Stability, FE_Main, &
        Mesh_Main, Boundary_Main)

    call connectInflowToMainArrays(FE_Inflow%Nnodes, Mesh_Inflow, Mesh_Main, Boundary_Main)
    
    call setGaussParameters(FE_Inflow, GaussInt_Inflow)
    call setGaussParameters(FE_Main, GaussInt_Main)

    call createPetscObjects(Problem_Inflow, FE_Inflow, Mesh_Inflow, &
        Boundary_Inflow, LS_Inflow)
    call createPetscObjectsProjection(Problem_Inflow, FE_Inflow%Ndim, LS_Inflow_proj)

    call createPetscObjects(Problem_Main, FE_Main, Mesh_Main, &
        Boundary_Main, LS_Main)
    call createPetscObjectsProjection(Problem_Main, FE_Main%Ndim, LS_Main_proj)

    call createSlepcObjects(Problem_Inflow_Stability, FE_Inflow, Mesh_Inflow, &
        Boundary_Inflow, EVP_Inflow)
    call createSlepcObjects(Problem_Main_Stability, FE_Main, Mesh_Main, &
        Boundary_Main, EVP_Main)

    call allocateUnknownsArrays(Problem_Inflow, Problem_Inflow_Stability, FE_Inflow, Sol_Inflow)
    call allocateUnknownsArrays(Problem_Main, Problem_Main_Stability, FE_Main, Sol_Main)

    call allocateLocalArrays(Problem_Inflow, Problem_Inflow_Stability, FE_Inflow, NodeArrays_Inflow)
    call allocateLocalArrays(Problem_Main, Problem_Main_Stability, FE_Main, NodeArrays_Main)
    
    call setDimensionlessNumbers()
    
    call writeMemoryUsage()

    call loopContinuation(tStart)

    call finalizePetsc(Continuation_Method)

end program FEM2D
