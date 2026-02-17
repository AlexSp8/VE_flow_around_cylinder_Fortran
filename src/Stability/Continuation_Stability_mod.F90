
#include <slepc/finclude/slepceps.h>

module ContinuationStability_mod

    use slepceps
    use PhysicalParameters_mod, only: ProblemParameters
    use SolutionVariables_mod, only: SolutionArraysType
    use SlepcVariables_mod, only: EVPType

    implicit none

    contains

    subroutine loopStability()

        use omp_lib, only: omp_get_wtime
        use ContinuationParameters_mod, only: Cvar2_Names, adjust_dt
        use PhysicalParameters_mod, only: dt0
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Increment, Cvar2, dt
        use Physical_mod, only: setCVar, initializeDimensionlessNumbers
        use Stability_mod, only: startingSolution
        use SlepcVariables_mod, only: EVP_Inflow, EVP_Main
        use Tools_mod, only: writeElapsedTime

        implicit none

        PetscInt :: ivar2, N
        PetscReal :: tStart, tEnd, tStart2, Kr_Leado

        tStart = omp_get_wtime()

        ! ! if (adjust_dt) then
        !     dt = dt0
        ! ! end if

        Increment = 0
        
        call initializeVariablesStability(EVP_Inflow)
        call initializeVariablesStability(EVP_Main)

        N = max(1, size(Cvar2_Names))
        do ivar2 = 1, N

            tStart2 = omp_get_wtime()

            if (size(Cvar2_Names) > 0) then
                call setCVar(Cvar2_Names(ivar2),Cvar2)
            else
                call setCVar('',Cvar2)                
            end if

            if (ivar2 > 1) then
                Cvar2%p = Cvar2%p + Cvar2%dvar
            end if

            if (Rank == 0) then
                write(*,*)
                write(*,'(a)') '**********************************************************************'
                write(*,'(27x,2a)') trim(adjustl(Cvar2%name)), ' Continuation'
                write(*,'(a)') '**********************************************************************'
            end if

            call var2Continuation(tStart)

            tEnd = omp_get_wtime()
            if (Rank == 0) then
                call writeElapsedTime(tStart2, tEnd, trim(adjustl(Cvar2%name))//' Continuation time')
                call writeElapsedTime(tStart, tEnd, 'Total elapsed time')
            end if

            call initializeDimensionlessNumbers()
            call startingSolution(EVP_Main%K_Lead, 'read')

        end do

    end subroutine loopStability

    !---------------------------------------------------------------------

    subroutine initializeVariablesStability(EVP)

        implicit none

        type(EVPType), intent(inout) :: EVP
        
        EVP%Kr_Leado = -1.0d8
        EVP%K_Lead = -1.0d8 + PETSC_I*0.0d0
        EVP%Ki_max = 0.0d0

        EVP%target_real = 0.0d0

        EVP%xCr = -1.0d0
        ! EVP%x0 = 0.0d0
        EVP%LwN_Lead = 0.0d0

        EVP%critical = .false.

    end subroutine initializeVariablesStability

    !---------------------------------------------------------------------

    subroutine var2Continuation(tStart)

        use omp_lib, only: omp_get_wtime
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Cvar2
        use Stability_mod, only: previousSolution
        use SlepcVariables_mod, only: EVP_Main
        use Tools_mod, only: writeElapsedTime

        implicit none

        PetscReal, intent(in) :: tStart

        PetscReal, parameter :: eps = 1.0d-12
        PetscReal :: Kr_Leado, tStart2, tEnd
        PetscInt :: i, N
        character(50) :: str

        N = nint(1+(Cvar2%var_f - Cvar2%p)/(Cvar2%dvar+eps))

        do i = 1, N

            tStart2 = omp_get_wtime()

            if (Cvar2%p < 0.0d0) Cvar2%p = 0.0d0

            if (Rank == 0) then
                call makeVar2Directories(Cvar2%name, Cvar2%p)
                write(*,*)
                write(str,'(3a,f8.4)') 'QbN Continuation on ', trim(adjustl(Cvar2%name)), ' = ', &
                    Cvar2%p
                write(*,'(a)')'//////////////////////////////////////////////////////////////////////'
                write(*,'(18x,a)') str
                write(*,'(a)')'//////////////////////////////////////////////////////////////////////'
            end if

            call loopQContinuation(tStart)

            if (Rank == 0) then
                call var2CriticalConditions(Cvar2%name, Cvar2%p)
                tEnd = omp_get_wtime()
                write(str,'(3a,f8.4)') 'QbN Continuation on ', trim(adjustl(Cvar2%name)), ' = ', Cvar2%p
                call writeElapsedTime(tStart2, tEnd, trim(adjustl(str))//' time')
                call writeElapsedTime(tStart, tEnd, 'Total elapsed time')
            end if

            call previousSolution(EVP_Main%K_Lead, 'read')
            Cvar2%var_o = Cvar2%p
            Cvar2%p = Cvar2%p + Cvar2%dvar

        end do

    end subroutine var2Continuation

    !---------------------------------------------------------------------

    subroutine makeVar2Directories(name, value)

        use ContinuationParameters_mod, only: Continuation_Method

        implicit none

        character(*), intent(in) :: name
        PetscReal, intent(in) :: value

        character(200) :: folder, str
        character(len=:), allocatable :: val

        write(str,'(f12.4)') value
        val = trim(adjustl(str))

        folder = 'Results/Base/'//trim(adjustl(name))//'/'
        call execute_command_line('mkdir -p '//folder)
        folder = trim(adjustl(folder))//trim(adjustl(name))//'_'//trim(adjustl(val))
        call execute_command_line('mkdir -p '//folder)

        call execute_command_line('mkdir -p '//trim(adjustl(folder))//'/Sol/')
        call execute_command_line('mkdir -p '//trim(adjustl(folder))//'/Tec/')
        call execute_command_line('mkdir -p '//trim(adjustl(folder))//'/Flow_data/')
        call execute_command_line('mkdir -p '//trim(adjustl(folder))//'/Line/')
        call execute_command_line('mkdir -p '//trim(adjustl(folder))//'/Txx/')
        call execute_command_line('mkdir -p '//trim(adjustl(folder))//'/Convergence/')

        if (Continuation_Method == 'Stability') then
            folder = 'Results/Stability/'//trim(adjustl(name))//'/'
            call execute_command_line('mkdir -p '//folder)
            folder = trim(adjustl(folder))//trim(adjustl(name))//'_'//trim(adjustl(val))
            call execute_command_line('mkdir -p '//folder)
        end if

    end subroutine makeVar2Directories

    !---------------------------------------------------------------------

    subroutine var2CriticalConditions(name, value)

        implicit none

        character(*), intent(in) :: name
        PetscReal, intent(in) :: value

        PetscInt :: ioerr
        PetscReal :: QbN, xCr, LwN, QbNc, xCr_min, LwNc
        character(200) :: fn
        character(len=:), allocatable :: val

        write(fn,'(f12.4)') value
        val = trim(adjustl(fn))

        fn = 'Results/Stability/'//trim(adjustl(name))//'/'
        fn = trim(adjustl(fn))//trim(adjustl(name))//'_'//trim(adjustl(val))
        fn = trim(adjustl(fn))//'/Critical_conditions.dat'

        QbNc = 0.0d0
        xCr_min = huge(1.0d0)
        LwNc = 0.0d0

        open(4, file=fn, iostat=ioerr, position='rewind')
        if (ioerr /= 0) return
        do
            read(4,*,iostat=ioerr) QbN, xCr, LwN
            if (ioerr /= 0) exit
            if (xCr < xCr_min) then
                xCr_min = xCr
                QbNc = QbN
                LwNc = LwN
            end if
        end do
        close(4)

        fn = 'Results/Stability/'//trim(adjustl(name))//'/Critical_conditions.dat'
        open(4, file=fn, iostat=ioerr, position='append')
        write(4,'(2f12.4,2f16.8)') value, QbNc, xCr_min, LwNc
        close(4)

    end subroutine var2CriticalConditions

    !---------------------------------------------------------------------

    subroutine loopQContinuation(tStart)

        use omp_lib, only: omp_get_wtime
        use ContinuationParameters_mod, only: Cvar2_Names
        use PhysicalParameters_mod, only: QbN_f, QbN0, dQbN, Nbd_periodic, QbN
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Cvar2, Cvar1, Increment
        use Stability_mod, only: previousSolution
        use SlepcVariables_mod, only: EVP_Main
        use Tools_mod, only: writeElapsedTime

        implicit none

        PetscReal, intent(in) :: tStart

        PetscReal, parameter :: eps = 1.0d-12
        PetscBool, parameter :: check = (size(Cvar2_Names) > 0)
        PetscReal :: tStart2, tEnd, QbNo
        PetscInt :: iQ, NQ
        character(100) :: str

        QbNo = QbN0
        ! if (Increment == 0) QbNo = 0.45d0

        NQ = nint(1+(QbN_f-QbNo)/(dQbN+eps))
        ! if ( Nbd_periodic < 1 ) NQ = 1

        QbN = QbNo
        do iQ = 1, NQ

            tStart2 = omp_get_wtime()

            if (QbN < 0.0d0) QbN = 0.0d0

            if (Rank == 0) then

                call makeQbNDirectories(Cvar2%name, Cvar2%p, QbN)

                write(*,*)
                write(str,'(10x,2a,f8.4,3a,f8.4)') trim(adjustl(Cvar1%name)), &
                    ' Continuation on QbN = ', QbN, ', ', trim(adjustl(Cvar2%name)), ' = ', Cvar2%p
                write(*,'(a)') '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
                write(*,'(a)') str
                write(*,'(a)') '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'

            end if

            call var1Continuation(tStart)

            if (check .and. iQ == 1) then
                call previousSolution(EVP_Main%K_Lead, 'write')
            end if

            tEnd = omp_get_wtime()
            if (Rank == 0) then
                write(str,'(2a,f8.4,3a,f8.4)') trim(adjustl(Cvar1%name)), &
                    ' Continuation on QbN = ', QbN, ', ', trim(adjustl(Cvar2%name)), ' = ', Cvar2%p
                call writeElapsedTime(tStart2, tEnd, trim(adjustl(str))//' time')
                call writeElapsedTime(tStart, tEnd, 'Total elapsed time')
            end if

            QbN = QbN + dQbN

        end do

    end subroutine loopQContinuation

    !---------------------------------------------------------------------

    subroutine makeQbNDirectories(name, value, QbN)

        implicit none

        character(*), intent(in) :: name
        PetscReal, intent(in) :: value, QbN

        character(200) :: folder, str
        character(len=:), allocatable :: val, Q_val

        write(str,'(f12.4)') value
        val = trim(adjustl(str))

        write(str,'(f12.4)') QbN
        Q_val = trim(adjustl(str))

        folder = 'Results/Stability/'//trim(adjustl(name))//'/'
        folder = trim(adjustl(folder))//trim(adjustl(name))//'_'//trim(adjustl(val))
        folder = trim(adjustl(folder))//'/QbN_'//trim(adjustl(Q_val))
        call execute_command_line('mkdir -p '//folder)

        call execute_command_line('mkdir -p '//trim(adjustl(folder))//'/Eigenvectors/')
        call execute_command_line('mkdir -p '//trim(adjustl(folder))//'/Eigenvalues/')
        call execute_command_line('mkdir -p '//trim(adjustl(folder))//'/Line/')
        call execute_command_line('mkdir -p '//trim(adjustl(folder))//'/Final_solution/')
        call execute_command_line('mkdir -p '//trim(adjustl(folder))//'/Stability_data/')

        deallocate(val) ; deallocate(Q_val)

    end subroutine makeQbNDirectories

    !---------------------------------------------------------------------

    subroutine var1Continuation(tStart)

        use omp_lib, only: omp_get_wtime
        use ContinuationParameters_mod, only: adjust_dt, Cvar2_Names
        use PhysicalParameters_mod, only: QbN, ReN, WiN, BvN, VeN, DfN
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Increment, Cvar2, Cvar1
        use Tools_mod, only: writeElapsedTime
        use Petsc_mod, only: writeMemoryUsage
        use SlepcVariables_mod, only: EVP_Main
        use Physical_mod, only: initializeDimensionlessNumbers

        implicit none

        PetscReal, intent(in) :: tStart

        PetscReal, parameter :: eps = 1.0d-8
        PetscReal :: tStart2, tEnd, Kr_Lead
        PetscBool :: check

        loop_continuation:do

            tStart2 = omp_get_wtime()

            Increment = Increment + 1

            if (Rank == 0) then

                call printStepInfoStability()

                open(10, file='Results/Base/Info/Variables.dat',position='append')

                Kr_Lead = PetscRealPart(EVP_Main%K_Lead)

                write(10,'(3f8.3,2es16.6,(*(f8.3)))') Cvar2%p, QbN, Cvar1%p, &
                    Kr_Lead, EVP_Main%Kr_Leado, ReN, WiN, BvN, VeN, DfN

                close(10)

            end if
            
            call readVar1Solution()

            call continuationStepStability()

            check = (EVP_Main%critical) .and. (size(Cvar2_Names) > 0)
            if (check) then
                EVP_Main%target_real = 0.0d0
                EVP_Main%xCr = -1.0d0
                EVP_Main%critical = .false.
                exit
            end if

            check = (Cvar1%p+eps > Cvar1%var_f)
            if (check) then
                ! call initializeDimensionlessNumbers()
                Cvar1%p = EVP_Main%x0
                call initializeVariablesStability(EVP_Main)
                exit
            end if

            if (adjust_dt) then
                call updateTimeStep()
            else
                call updateContinuationParameter()
            end if

            tEnd = omp_get_wtime()
            if (Rank == 0) then
                call writeElapsedTime(tStart2, tEnd, 'Continuation step time')
                call writeElapsedTime(tStart, tEnd, 'Total elapsed time')
            end if

            call writeMemoryUsage()

        end do loop_continuation

    end subroutine var1Continuation

    !---------------------------------------------------------------------

    subroutine readVar1Solution()

        use ContinuationVariables_mod, only: Cvar2, Cvar1
        use PhysicalParameters_mod, only: Problem_Main, Problem_Inflow
        use FEMParameters_mod, only: FE_Inflow, FE_Main
        use MeshParameters_mod, only: Mesh_Inflow, Mesh_Main
        use SolutionVariables_mod, only: Sol_Inflow, Sol_Main
        use PostProcess_mod, only: readSolution

        implicit none

        PetscReal :: iter_dump
        character(100) :: fn, str, folder
        character(len=:), allocatable :: val1, val2

        write(str,'(f12.4)') Cvar1%p
        ! write(str,'(f12.4)') nint(Cvar1%p * 10.0**4) / 10.0**4
        val1 = trim(adjustl(str))

        write(str,'(f12.4)') Cvar2%p
        val2 = trim(adjustl(str))

        folder = 'Results/Base/'//trim(adjustl(Cvar2%name))//'/' &
        ! folder = 'Sol/'//trim(adjustl(Cvar2%name))//'/' &
            //trim(adjustl(Cvar2%name))//'_'//trim(adjustl(val2))//'/Sol/'
        fn = trim(adjustl(folder))//trim(adjustl(Problem_Inflow%name))//'_'// &
                trim(adjustl(Cvar1%name))//'_'//val1
        call readSolution(FE_Inflow, Mesh_Inflow%Connectivity, iter_dump, &
                            fn, Sol_Inflow)

        folder = 'Results/Base/'//trim(adjustl(Cvar2%name))//'/' &
        ! folder = 'Sol/'//trim(adjustl(Cvar2%name))//'/' &
            //trim(adjustl(Cvar2%name))//'_'//trim(adjustl(val2))//'/Sol/'
        fn = trim(adjustl(folder))//trim(adjustl(Problem_Main%name))//'_'// &
                trim(adjustl(Cvar1%name))//'_'//val1
        call readSolution(FE_Main, Mesh_Main%Connectivity, iter_dump, fn, Sol_Main)

    end subroutine readVar1Solution

    !---------------------------------------------------------------------

    subroutine continuationStepStability()

        use PhysicalParameters_mod, only: Problem_Inflow, Problem_Main, &
            Problem_Inflow_Stability, Problem_Main_Stability
        use FEMParameters_mod, only: FE_Inflow, FE_Main
        use MeshParameters_mod, only: Mesh_Inflow, Mesh_Main
        use BoundaryParameters_mod, only: Boundary_Inflow, Boundary_Main
        use GaussParameters_mod, only: GaussInt_Inflow, GaussInt_Main
        use ContinuationVariables_mod, only: Cvar1
        use ElementVariables_mod, only: NodeArrays_Inflow, NodeArrays_Main
        use SolutionVariables_mod, only: Sol_Inflow, Sol_Main
        use LinearSystemVariables_mod, only: LS_Inflow, LS_Main, &
            LS_Main_proj, LS_Inflow_proj
        use SlepcVariables_mod, only: EVP_Inflow, EVP_Main
        use NR_mod, only: loopNewtonRaphson
        use Stability_mod, only: stabilityProblem
        use PostProcess_mod, only: postProcess

        implicit none

        PetscReal, save :: twrite0_Main = 0.0d0, twrite0_Inflow = 0.0d0
        PetscInt, save :: iStab = 0
        PetscBool :: check

        !Inflow
        if (Problem_Inflow%Ndim > 0) then

            call loopNewtonRaphson(Problem_Inflow, FE_Inflow, Mesh_Inflow, &
                Boundary_Inflow, GaussInt_Inflow, Sol_Inflow, NodeArrays_Inflow, LS_Inflow)

            call postProcess(Problem_Inflow, FE_Inflow, Mesh_Inflow, Boundary_Inflow, &
                GaussInt_Inflow, Sol_Inflow, NodeArrays_Inflow, LS_Inflow_proj, twrite0_Inflow)

            call updateBaseSolution(Sol_Inflow)
            if (Problem_Inflow_Stability%Neq > 0) then
                call updateStabilityArrays(Problem_Inflow, Problem_Inflow_Stability, Sol_Inflow)
            end if

            check = (Problem_Inflow_Stability%Neq > 0)
            check = check .and. checkStabilityIncrement()
            ! check = check .or. (iStab > 0)
            if (check) then
                call stabilityProblem(Problem_Inflow, Problem_Inflow_Stability, &
                    FE_Inflow, Mesh_Inflow, Boundary_Inflow, GaussInt_Inflow, &
                    Sol_Inflow, NodeArrays_Inflow, EVP_Inflow)
            end if

        end if

        call loopNewtonRaphson(Problem_Main, FE_Main, Mesh_Main, &
            Boundary_Main, GaussInt_Main, Sol_Main, NodeArrays_Main, LS_Main)

        call postProcess(Problem_Main, FE_Main, Mesh_Main, &
            Boundary_Main, GaussInt_Main, Sol_Main, NodeArrays_Main, LS_Main_proj, twrite0_Main)

        call updateBaseSolution(Sol_Main)

        if (Problem_Main_Stability%Neq > 0) then
            call updateStabilityArrays(Problem_Main, Problem_Main_Stability, Sol_Main)
        end if

        check = (Problem_Main_Stability%Neq > 0)
        check = check .and. checkStabilityIncrement()
        check = check .or. (iStab > 0)
        if (.not. check) return

        iStab = iStab+1
        if (iStab == 1) EVP_Main%x0 = Cvar1%p

        call stabilityProblem(Problem_Main, Problem_Main_Stability, FE_Main, Mesh_Main, &
            Boundary_Main, GaussInt_Main, Sol_Main, NodeArrays_Main, EVP_Main)

    end subroutine continuationStepStability

    !---------------------------------------------------------------------

    subroutine printStepInfoStability()

        use ContinuationVariables_mod, only: Increment, Cvar1
        use SolutionVariables_mod, only: Sol_Main

        implicit none

        write(*,*)''
        write(*,'(a)') '######################################################################'
        write(*,60) Increment, trim(adjustl(Cvar1%name)), Cvar1%p, Sol_Main%EU(:)
        write(*,'(a)') '######################################################################'
        write(*,*)''

        60 format('STEP = ',i5,5x,a,' = ',f12.4,5x,' EU = ',*(f10.3))

    end subroutine printStepInfoStability

    !---------------------------------------------------------------------

    subroutine updateBaseSolution(Sol)

        use MPIParameters_mod, only: Rank

        implicit none

        type(SolutionArraysType), intent(inout) :: Sol

        PetscInt :: Neq_S

        if (Rank == 0) then
            Sol%TLp(:,:) = Sol%TL(:,:)
            Sol%TLb(:,:) = Sol%TLo(:,:)
            Sol%TLo(:,:) = Sol%TL(:,:)
            Sol%TL(:,:) = Sol%TLp(:,:)
        end if

        Sol%TLp_rank(:,:) = Sol%TL_rank(:,:)
        Sol%TLb_rank(:,:) = Sol%TLo_rank(:,:)
        Sol%TLo_rank(:,:) = Sol%TL_rank(:,:)
        Sol%TL_rank(:,:)  = Sol%TLp_rank(:,:)

        Sol%EUp(:) = Sol%EU(:)
        Sol%EUb(:) = Sol%EUo(:)
        Sol%EUo(:) = Sol%EU(:)
        Sol%EU(:) = Sol%EUp(:)

    end subroutine updateBaseSolution

    !--------------------------------------------------------------------------

    subroutine updateStabilityArrays(Problem, Problem_Stability, Sol)

        use MPIParameters_mod, only: Rank

        implicit none
        
        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(SolutionArraysType), intent(inout) :: Sol

        PetscInt :: Ndim, Neq, Ndim_s, Nex_s

        Ndim = Problem%Ndim
        Neq = Problem%Neq
        Ndim_s = Problem_Stability%Ndim
        Nex_s = Problem_Stability%Nex

        if (Rank == 0) then

            Sol%TL_b(:,:) = 0.0d0
            if (Ndim_s == Ndim+1) then
               Sol%TL_b(:,1:Ndim) = Sol%TL(:,1:Ndim)
               !Uz stored at Ndim+1
               Sol%TL_b(:,Ndim+2:Neq+1) = Sol%TL(:,Ndim+1:Neq)
            else
                Sol%TL_b(:,1:Neq) = Sol%TL(:,:)
            end if

        end if

        Sol%TL_b_rank(:,:) = 0.0d0
        if (Ndim_s == Ndim+1) then
           Sol%TL_b_rank(:,1:Ndim) = Sol%TL_rank(:,1:Ndim)
            !Uz stored at Ndim+1
           Sol%TL_b_rank(:,Ndim+2:Neq+1) = Sol%TL_rank(:,Ndim+1:Neq)
        else
            Sol%TL_b_rank(:,1:Neq) = Sol%TL_rank(:,:)
        end if

        Sol%EU_b(:) = 0.0d0
        Sol%EU_b(:) = Sol%EU(1:Nex_s)

    end subroutine updateStabilityArrays

    !--------------------------------------------------------------------------

    function checkStabilityIncrement() result(check)

        use ContinuationParameters_mod, only: Cont_var0_Stability
        use ContinuationVariables_mod, only: Cvar1

        implicit none

        PetscBool :: check

        PetscReal, parameter :: eps = 1.0d-8

        check = ( Cvar1%p > (Cont_var0_Stability - eps) )

    end function checkStabilityIncrement

    !--------------------------------------------------------------------------

    subroutine updateContinuationParameter()

        use ContinuationVariables_mod, only: Cvar1
        
        implicit none

        PetscReal :: step

        step = Cvar1%dvar
        ! if (Cvar1%p < 0.5d0) step = 2*step

        Cvar1%p = Cvar1%p + step

    end subroutine updateContinuationParameter

    !--------------------------------------------------------------------------

    subroutine updateTimeStep()

        use PhysicalParameters_mod, only: dt0
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: dt, Increment
        
        implicit none

        PetscInt, save :: k = 0
        PetscReal :: factor
        PetscBool :: check
        PetscErrorCode :: ierr

        check = (mod(Increment-1,9) == 0)
        if (check) k = k+1

        factor = 10.0d0**k
        dt = dt + dt0*factor

        ! factor = (1.0d0/10.0d0**k)
        ! dt = dt - dt0*factor

        if (Rank == 0) then
            open(12,file='Results/Base/Info/dt.dat',position='append')
            write(12,'(i8,2es20.8)') k, dt0*factor, dt
            close(12,status='keep')

            write(*,'(a,es20.12)') 'dt = ', dt
            write(*, '(a)') '--------------------------------------------------------------'
        end if

        if (dt > 1.0d8 .or. dt < 1.0d-12) then
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
            stop
        end if

    end subroutine updateTimeStep

end module ContinuationStability_mod
