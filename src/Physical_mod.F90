
#include <petsc/finclude/petscksp.h>

module Physical_mod

    use PhysicalParameters_mod, only: ProblemParameters

    contains

    subroutine setProblemParameters(Problem_name, Neq_proj, BNDs_p, Problem)

        use Tools_mod, only: allocateArray
        use ContinuationParameters_mod, only: Continuation_Method
        use PhysicalParameters_mod, only: Height, Width
        
        implicit none

        character(*), intent(in) :: Problem_name
        PetscInt, intent(in) :: Neq_proj
        PetscInt, dimension(:,:), intent(in) :: BNDs_p
        type(ProblemParameters), intent(out) :: Problem

        PetscReal :: flow_rate

        Problem%name = Problem_name

        Problem%Neq_proj = Neq_proj

        select case (Problem%name)

        case ('Inflow_1D')

            Problem%Ndim = 1

            call allocateArray(Problem%idir, Problem%Ndim)
            Problem%idir(1) = 2

            Problem%Neq_f = 2
            Problem%Neq_s = 3
            Problem%Neq = 5

            Problem%Nex_b = 1
            Problem%Nex_ed = 0
            Problem%Nex_fc = 0
            Problem%Nex = Problem%Nex_b + Problem%Nex_ed + Problem%Nex_fc

            call allocateArray(Problem%iex_b_values, Problem%Nex_b)
            flow_rate = 2.0d0*Height
            Problem%iex_b_values(1) = flow_rate

            Problem%Nbd = 2

        case ('Inflow_Stability_1D')

            Problem%Ndim = 1

            call allocateArray(Problem%idir, Problem%Ndim)
            Problem%idir(1) = 2

            Problem%Neq_f = 2
            Problem%Neq_s = 3
            Problem%Neq = 5

            Problem%Nex_b = 0
            Problem%Nex_ed = 0
            Problem%Nex_fc = 0
            Problem%Nex = Problem%Nex_b + Problem%Nex_ed + Problem%Nex_fc

            call allocateArray(Problem%iex_b_values, Problem%Nex_b)

            Problem%Nbd = 2

        case ('Inflow_Stability_1D_2D')

            Problem%Ndim = 2

            call allocateArray(Problem%idir, Problem%Ndim)
            Problem%idir(:) = [2, 3]

            Problem%Neq_f = 3
            Problem%Neq_s = 6
            Problem%Neq = 9

            Problem%Nex_b = 0
            Problem%Nex_ed = 0
            Problem%Nex_fc = 0
            Problem%Nex = Problem%Nex_b + Problem%Nex_ed + Problem%Nex_fc

            call allocateArray(Problem%iex_b_values, Problem%Nex_b)

            Problem%Nbd = 2

        case ('Cylinder_2D')

            Problem%Ndim = 2

            call allocateArray(Problem%idir, Problem%Ndim)
            Problem%idir(:) = [1, 2]

            Problem%Neq_f = 2
            Problem%Neq_s = 3
            Problem%Neq = 6

            Problem%Nex_b = 0
            Problem%Nex_ed = 0
            Problem%Nex_fc = 0
            Problem%Nex = Problem%Nex_b + Problem%Nex_ed + Problem%Nex_fc
            if (Continuation_Method == 'Arclength') Problem%Nex = Problem%Nex + 1

            call allocateArray(Problem%iex_b_values, Problem%Nex_b)

            Problem%Nbd = 5

        case ('Cylinder_Stability_2D')

            Problem%Ndim = 2

            call allocateArray(Problem%idir, Problem%Ndim)
            Problem%idir(:) = [1, 2]

            Problem%Neq_f = 2
            Problem%Neq_s = 3
            Problem%Neq = 6

            Problem%Nex_b = 0
            Problem%Nex_ed = 0
            Problem%Nex_fc = 0
            Problem%Nex = Problem%Nex_b + Problem%Nex_ed + Problem%Nex_fc

            call allocateArray(Problem%iex_b_values, Problem%Nex_b)

            Problem%Nbd = 5

        case ('Cylinder_Stability_2D_3D')

            Problem%Ndim = 3

            call allocateArray(Problem%idir, 2)
            Problem%idir(:) = [1, 2]

            Problem%Neq_f = 3
            Problem%Neq_s = 6
            Problem%Neq = 10

            Problem%Nex_b = 0
            Problem%Nex_ed = 0
            Problem%Nex_fc = 0
            Problem%Nex = Problem%Nex_b + Problem%Nex_ed + Problem%Nex_fc

            call allocateArray(Problem%iex_b_values, Problem%Nex_b)

            Problem%Nbd = 5

        case ('Inflow_2D')

            Problem%Ndim = 2

            call allocateArray(Problem%idir, Problem%Ndim)
            Problem%idir(:) = [2, 3]

            Problem%Neq_f = 3
            Problem%Neq_s = 6
            Problem%Neq = 9

            Problem%Nex_b = 1
            Problem%Nex_ed = 0
            Problem%Nex_fc = 0
            Problem%Nex = Problem%Nex_b + Problem%Nex_ed + Problem%Nex_fc

            call allocateArray(Problem%iex_b_values, Problem%Nex_b)
            flow_rate = 2.0d0*Height*2.0d0*Width
            Problem%iex_b_values(1) = flow_rate

            Problem%Nbd = 4

        case ('Inflow_Stability_2D')

            Problem%Ndim = 2

            call allocateArray(Problem%idir, Problem%Ndim)
            Problem%idir(:) = [2, 3]

            Problem%Neq_f = 3
            Problem%Neq_s = 6
            Problem%Neq = 9

            Problem%Nex_b = 0
            Problem%Nex_ed = 0
            Problem%Nex_fc = 0
            Problem%Nex = Problem%Nex_b + Problem%Nex_ed + Problem%Nex_fc

            call allocateArray(Problem%iex_b_values, Problem%Nex_b)

            Problem%Nbd = 4

        case ('Cylinder_3D')

            Problem%Ndim = 3

            call allocateArray(Problem%idir, Problem%Ndim)
            Problem%idir(:) = [1, 2, 3]

            Problem%Neq_f = 3
            Problem%Neq_s = 6
            Problem%Neq = 10

            Problem%Nex_b = 0
            Problem%Nex_ed = 0
            Problem%Nex_fc = 0
            Problem%Nex = Problem%Nex_b + Problem%Nex_ed + Problem%Nex_fc
            if (Continuation_Method == 'Arclength') Problem%Nex = Problem%Nex + 1

            call allocateArray(Problem%iex_b_values, Problem%Nex_b)

            Problem%Nbd = 7

        case ('Cylinder_Stability_3D')

            Problem%Ndim = 3

            call allocateArray(Problem%idir, Problem%Ndim)
            Problem%idir(:) = [1, 2, 3]

            Problem%Neq_f = 3
            Problem%Neq_s = 6
            Problem%Neq = 10

            Problem%Nex_b = 0
            Problem%Nex_ed = 0
            Problem%Nex_fc = 0
            Problem%Nex = Problem%Nex_b + Problem%Nex_ed + Problem%Nex_fc

            call allocateArray(Problem%iex_b_values, Problem%Nex_b)
            
            Problem%Nbd = 7

        case default

            Problem%Ndim = 0

            call allocateArray(Problem%idir, Problem%Ndim)

            Problem%Neq_f = 0
            Problem%Neq_s = 0
            Problem%Neq = 0

            Problem%Nex_b = 0
            Problem%Nex_ed = 0
            Problem%Nex_fc = 0
            Problem%Nex = 0

            call allocateArray(Problem%iex_b_values, Problem%Nex_b)

            Problem%Nbd = 0

            call allocateArray(Problem%PeriodicBNDs, 0, 0)

            Problem%Nunknowns = 0

            return

        end select
        
        if (size(BNDs_p,1) > 0) then
            call setProblemPeriodicBoundaries(BNDs_p, Problem)
        end if

    end subroutine setProblemParameters

    ! ----------------------------------------------------------------------
    
    subroutine setProblemPeriodicBoundaries(BNDs_p, Problem)

        use Tools_mod, only: allocateArray

        implicit none

        PetscInt, dimension(:,:), intent(in) :: BNDs_p
        type(ProblemParameters), intent(inout) :: Problem

        PetscInt :: Nbd_periodic

        Nbd_periodic = size(BNDs_p,1)

        call allocateArray(Problem%PeriodicBNDs, Nbd_periodic, 2)

        Problem%PeriodicBNDs(:,:) = BNDs_p(:,:)

    end subroutine setProblemPeriodicBoundaries

    !---------------------------------------------------------------------

    subroutine setDimensionlessNumbers()

        use ContinuationParameters_mod, only: Continuation_Method, &
            LN0, Cvar2_Names, Cvar1_Name
        use MPIParameters_mod, only: Rank
        use ContinuationVariables_mod, only: Cvar1, Cvar2

        implicit none

        PetscInt :: ivar2
        character(20) :: name

        call initializeDimensionlessNumbers()
        
        call setCVar(Cvar1_Name,Cvar1)

        if (Continuation_Method == 'Arclength') then
            Cvar1%p = (Cvar1%p)*LN0
        end if

        ! Cvar2 = Cvar1
        ! if (size(Cvar2_Names) > 0) then
        !     name = Cvar2_Names(1)
        !     call setCVar(name,Cvar2)
        ! end if

        if (Rank == 0) then
            call writeInitialValues()
        end if

    end subroutine setDimensionlessNumbers

    !---------------------------------------------------------------------

    subroutine initializeDimensionlessNumbers()

        use ContinuationParameters_mod, only: Continuation_Method
        use PhysicalParameters_mod, only: ReN, WiN, BvN, VeN, DfN, Total_time, &
            ReN0, WiN0, BvN0, VeN0, Time0, DfN0, dt0
        use ContinuationVariables_mod, only: dt

        implicit none

        PetscBool :: check

        ReN = ReN0
        WiN = WiN0
        Total_time = Time0
        BvN = BvN0
        VeN = VeN0
        DfN = DfN0

        dt = +huge(1.0d0)
        check = (Continuation_Method == 'Transient')
        if (check) then
            dt = dt0
        end if

    end subroutine initializeDimensionlessNumbers

    !---------------------------------------------------------------------

    subroutine setCVar(name, Cvar)

        use PhysicalParameters_mod, only: ReN, ReN0, ReN_f, dReN, &
            WiN, WiN0, WiN_f, dWiN, BvN, BvN0, BvN_f, dBvN, &
            VeN, VeN0, VeN_f, dVeN, DfN, DfN0, DfN_f, dDfN, &
            Total_time, dt0, Time0, Time_f, dum
        use ContinuationVariables_mod, only: ContVarType

        implicit none

        character(*), intent(in) :: name
        type(ContVarType), intent(inout) :: Cvar

        Cvar%name = name

        select case (name)
        case ('BvN')
            Cvar%p => BvN
            Cvar%dvar = dBvN
            Cvar%var_f = BvN_f
            Cvar%var0 = BvN0
        case ('VeN')
            Cvar%p => VeN
            Cvar%dvar = dVeN
            Cvar%var_f = VeN_f
            Cvar%var0 = VeN0
        case ('ReN')
            Cvar%p => ReN
            Cvar%dvar = dReN
            Cvar%var_f = ReN_f
            Cvar%var0 = ReN0
        case ('DfN')
            Cvar%p => DfN
            Cvar%dvar = dDfN
            Cvar%var_f = DfN_f
            Cvar%var0 = DfN0
        case ('WiN')
            Cvar%p => WiN
            Cvar%dvar = dWiN
            Cvar%var_f = WiN_f
            Cvar%var0 = WiN0
        case ('Time')
            Cvar%p => Total_time
            Cvar%dvar = dt0
            Cvar%var_f = Time_f
            Cvar%var0 = Time0
        case default
            dum = 0.1d0
            Cvar%p => dum!NULL()
            Cvar%dvar = 0.0d0
            Cvar%var_f = Cvar%p
            Cvar%var0 = Cvar%p
            Cvar%name = 'Dum'
        end select

        Cvar%var_o = Cvar%p

    end subroutine setCVar

    !---------------------------------------------------------------------

    subroutine writeInitialValues()

        use PhysicalParameters_mod, only: Length, Height, Radius, Width, BR, Model_Name, &
            BvN, VeN, KwN, LwN, Main_Name_Stability, ReN, WiN, Stress_Reform, DfN
        use ContinuationParameters_mod, only: Continuation_Method, Eps_Pr, Stability_Im
        use FEMParameters_mod, only: FEOrder, FEType_Main
        use ContinuationVariables_mod, only: dt

        implicit none

        character(50) :: date, host

        write(*,*)
        write(*,'(a)') '----------------------------------------------------------------------'

        call fdate(date)
        write(*,'(a)') date
        write(*,*)

        call execute_command_line('hostname >> host.txt')
        open(20, file='host.txt')
        read(20,*) host
        close(20, status='delete')
        write(*,'(2a)') 'Running in: ', host
        write(*,*)
        
        select case (Stress_Reform)
        case ('SQRT')
            write(*,'(a)') "Solving for: S"
        case default
            write(*,'(a)') "Solving for: T"
        end select
        write(*,*)

        write(*,'(2a,1x,a)') "Elements: ", FEOrder, FEType_Main
        write(*,*)

        write(*,'(4(a,f8.4))') "Domain: L = ", Length, ", H = ", Height, &
                                ", W = ", Width, ", BR = ", BR
        write(*,*)

        write(*,'(2(a,f8.4),a,es12.4)') "ReN = ", ReN, ", WiN = ", WiN, ", DfN = ", DfN
        write(*,*)

        write(*,'(2a,2(a,f8.4))') "Model: ", Model_Name, ": BvN = ", BvN, ", VeN = ", VeN
        write(*,*)

        write(*,'(2a)') "Continuation Method: ", Continuation_Method

        write(*,'(a,es12.4)') "dt = ", dt
        write(*,*)

        if (Continuation_Method == 'Arclength') then
            write(*,'(a,es12.4)') "dS eps = ", Eps_Pr
        end if

        if (Main_Name_Stability == 'Cylinder_Stability_2D_3D') then
            write(*,*)
            write(*,'(2(a,f8.4),a,L)') "KwN = ", KwN, ", LwN = ", LwN, &
                ", Imaginary = ", Stability_Im
        end if

        write(*,'(a)') '----------------------------------------------------------------------'
        write(*,*)

    end subroutine writeInitialValues

end module Physical_mod
