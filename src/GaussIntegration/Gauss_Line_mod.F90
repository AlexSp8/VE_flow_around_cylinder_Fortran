
#include <petsc/finclude/petscksp.h>

module GaussLine_mod

    use GaussParameters_mod, only: GaussIntegration

    implicit none

    PetscInt, parameter :: NGauss_1D = 2

    contains

    subroutine setGaussParametersLine(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt
        
        GaussInt%NGauss_1D = NGauss_1D
        if (GaussInt%order == 'Quadratic' .or. GaussInt%order == 'Serendipity') then
            GaussInt%NGauss_1D = 3
        end if
        GaussInt%NGauss_2D = 0
        GaussInt%NGauss_3D = 0
        GaussInt%NGauss_LC = 0

        GaussInt%NGauss_b = GaussInt%NGauss_1D

        call setGaussPointsLine(GaussInt)
        call setGaussWeightsLine(GaussInt)
        call setGaussBasisFunctionsLine(GaussInt)

    end subroutine setGaussParametersLine

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsLine(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        call setGaussPointsLine1D(GaussInt)

    end subroutine setGaussPointsLine

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsLine1D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        allocate( GaussInt%Points_1D(GaussInt%NGauss_1D) )

        ![-1,1]
        select case (GaussInt%NGauss_1D)
        case (1)
            GaussInt%Points_1D(1) = 0.0d0+1.0d-8
        case (2)
            GaussInt%Points_1D(1) = -1.0d0/sqrt(3.0d0)
            GaussInt%Points_1D(2) = +1.0d0/sqrt(3.0d0)
        case (3)
            GaussInt%Points_1D(1) = 0.0d0
            GaussInt%Points_1D(2) = +sqrt(3.0d0/5.0d0)
            GaussInt%Points_1D(3) = -sqrt(3.0d0/5.0d0)
        case (4)
            GaussInt%Points_1D(1) = +sqrt( (3.0d0/7.0d0)-(2.0d0/7.0d0)*sqrt(6.0d0/5.0d0) )
            GaussInt%Points_1D(2) = -sqrt( (3.0d0/7.0d0)-(2.0d0/7.0d0)*sqrt(6.0d0/5.0d0) )
            GaussInt%Points_1D(3) = +sqrt( (3.0d0/7.0d0)+(2.0d0/7.0d0)*sqrt(6.0d0/5.0d0) )
        case (5)
            GaussInt%Points_1D(1) = 0.0d0
            GaussInt%Points_1D(2) = +(1.0d0/3.0d0)*sqrt( 5.0d0-2.0d0*sqrt(10.0d0/7.0d0) )
            GaussInt%Points_1D(3) = -(1.0d0/3.0d0)*sqrt( 5.0d0-2.0d0*sqrt(10.0d0/7.0d0) )
            GaussInt%Points_1D(4) = +(1.0d0/3.0d0)*sqrt( 5.0d0+2.0d0*sqrt(10.0d0/7.0d0) )
            GaussInt%Points_1D(5) = -(1.0d0/3.0d0)*sqrt( 5.0d0+2.0d0*sqrt(10.0d0/7.0d0) )
        case default
            write(*,'(a)') 'Wrong NGauss_1D in setGaussPointsLine1D!'
            stop
        end select

        GaussInt%Points_1D(:) = GaussInt%Points_1D(:)

    end subroutine setGaussPointsLine1D

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsLine(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        call setGaussWeightsLine1D(GaussInt)

        allocate( GaussInt%Weights_b(GaussInt%NGauss_1D) )
        GaussInt%Weights_b(:) = GaussInt%Weights_1D(:)

    end subroutine setGaussWeightsLine

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsLine1D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        allocate( GaussInt%Weights_1D(GaussInt%NGauss_1D) )

        select case (GaussInt%NGauss_1D)
        case (1)
            GaussInt%Weights_1D(1) = 2.0d0
        case (2)
            GaussInt%Weights_1D(1) = 1.0d0
            GaussInt%Weights_1D(2) = 1.0d0
        case (3)
            GaussInt%Weights_1D(1) = 8.0d0/9.0d0
            GaussInt%Weights_1D(2) = 5.0d0/9.0d0
            GaussInt%Weights_1D(3) = 5.0d0/9.0d0
        case (4)
            GaussInt%Weights_1D(1) = (18.0d0+sqrt(30.0d0))/36.0d0
            GaussInt%Weights_1D(2) = (18.0d0+sqrt(30.0d0))/36.0d0
            GaussInt%Weights_1D(3) = (18.0d0-sqrt(30.0d0))/36.0d0
            GaussInt%Weights_1D(4) = (18.0d0-sqrt(30.0d0))/36.0d0
        case (5)
            GaussInt%Weights_1D(1) = 128.0d0/225.0d0
            GaussInt%Weights_1D(1) = (332.0d0+13.0d0*sqrt(70.0d0))/900.0d0
            GaussInt%Weights_1D(2) = (332.0d0+13.0d0*sqrt(70.0d0))/900.0d0
            GaussInt%Weights_1D(3) = (332.0d0-13.0d0*sqrt(70.0d0))/900.0d0
            GaussInt%Weights_1D(4) = (332.0d0-13.0d0*sqrt(70.0d0))/900.0d0
        case default
            write(*,'(a)') 'Wrong NGauss_1D in setGaussWeightsLine1D!'
            stop
        end select

    end subroutine setGaussWeightsLine1D

    !----------------------------------------------------------------------- 

    subroutine setGaussBasisFunctionsLine(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        call setGaussBasisFunctionsLine1D(GaussInt)

    end subroutine setGaussBasisFunctionsLine

    !----------------------------------------------------------------------- 

    subroutine setGaussBasisFunctionsLine1D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        PetscInt :: ig
        PetscReal :: c
        PetscReal, dimension(GaussInt%Nbf) :: bfn, dfdc, d2fdc2    !C-> ξ

        allocate( GaussInt%Bfn(GaussInt%Nbf, GaussInt%NGAUSS_1D) )
        allocate( GaussInt%dfdc(GaussInt%Nbf, GaussInt%NGAUSS_1D) )
        allocate( GaussInt%d2fdc2(GaussInt%Nbf, GaussInt%NGAUSS_1D) )

        do ig = 1, GaussInt%NGauss_1D

            c = GaussInt%Points_1D(ig)

            select case (GaussInt%order)
            case ('Linear')
                call setLinearBasisFunctions1D(bfn, dfdc, c)
                d2fdc2(:) = 0.0d0
            case ('Quadratic', 'Serendipity')
                call setQuadraticBasisFunctions1D(bfn, dfdc, d2fdc2, c)
            case default
                write(*,'(a)') 'Wrong GaussInt%order in setGaussBasisFunctionsLine1D!'
                stop
            end select
            GaussInt%Bfn(:,ig) = bfn(:)
            GaussInt%dfdc(:,ig) = dfdc(:)
            GaussInt%d2fdc2(:,ig) = d2fdc2(:)

        end do

    end subroutine setGaussBasisFunctionsLine1D

    !----------------------------------------------------------------------- 

    subroutine setLinearBasisFunctions1D(bfn, dfdc, c)

        implicit none

        PetscReal, intent(in)  :: c
        PetscReal, dimension(:), intent(out) :: bfn, dfdc

        !Basis function, φ
        bfn(1) = (1.0d0-c)/2.0d0
        bfn(2) = (1.0d0+c)/2.0d0

        !dφ/dξ
        dfdc(1) = -1.0d0/2.0d0
        dfdc(2) =  1.0d0/2.0d0

    end subroutine setLinearBasisFunctions1D

    !----------------------------------------------------------------------- 

    subroutine setQuadraticBasisFunctions1D(bfn, dfdc, d2fdc2, c)

        implicit none

        PetscReal, intent(in)  :: c
        PetscReal, dimension(:), intent(out) :: bfn, dfdc, d2fdc2

        !Basis function, φ
        bfn(1)  = -(1.0d0-c)*c/2.0d0
        bfn(2)  = +(1.0d0+c)*c/2.0d0
        bfn(3)  = +(1.0d0-c)*(1.0d0+c)

        !dφ/dξ
        dfdc(1) = -(1.0d0-2.0d0*c)/2.0d0
        dfdc(2) = +(1.0d0+2.0d0*c)/2.0d0
        dfdc(3) = -2.0d0*c
        
        !d2φ/dξ2
        d2fdc2(1) = +1.0d0
        d2fdc2(2) = +1.0d0
        d2fdc2(3) = -2.0d0

    end subroutine setQuadraticBasisFunctions1D

end module GaussLine_mod
