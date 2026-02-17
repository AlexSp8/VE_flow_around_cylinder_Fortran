
#include <petsc/finclude/petscksp.h>

module GaussTriangle_mod

    use GaussParameters_mod, only: GaussIntegration

    implicit none

    PetscInt, parameter :: NGauss_1D = 2, NGauss_2D = 4, NGauss_LC = 6

    contains

    subroutine setGaussParametersTriangle(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt
        
        GaussInt%NGauss_1D = NGauss_1D
        GaussInt%NGauss_2D = NGauss_2D
        if (GaussInt%order == 'Quadratic' .or. GaussInt%order == 'Serendipity') then
            GaussInt%NGauss_1D = 3
            GaussInt%NGauss_2D = 7
        end if
        GaussInt%NGauss_3D = 0
        GaussInt%NGauss_LC = NGauss_LC

        GaussInt%NGauss_b = GaussInt%NGauss_2D

        call setGaussPointsTriangle(GaussInt)
        call setGaussWeightsTriangle(GaussInt)
        call setGaussBasisFunctionsTriangle(GaussInt)

    end subroutine setGaussParametersTriangle

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsTriangle(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        call setGaussPointsTriangle1D(GaussInt)

        call setGaussPointsTriangle2D(GaussInt)

        call setGaussPointsTriangleLC(GaussInt)

    end subroutine setGaussPointsTriangle

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsTriangle1D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        allocate( GaussInt%Points_1D(GaussInt%NGauss_1D) )

        !Not in [-1,1], but in [0,1]: ξ -> (ξ+1)/2, η -> (η+1)/2, ζ -> (ζ+1)/2
        select case (GaussInt%NGauss_1D)
        case (1)
            GaussInt%Points_1D(1) = 0.0d0+1.0d-8
        case (2)
            GaussInt%Points_1D(1) = 1.0d0-(1.0d0/sqrt(3.0d0))
            GaussInt%Points_1D(2) = 1.0d0+(1.0d0/sqrt(3.0d0))
        case (3)
            GaussInt%Points_1D(1) = 1.0d0-0.0d0
            GaussInt%Points_1D(2) = 1.0d0+sqrt(3.0d0/5.0d0)
            GaussInt%Points_1D(3) = 1.0d0-sqrt(3.0d0/5.0d0)
        case default
            write(*,'(a)') 'Wrong NGauss_1D in setGaussPointsTriangle1D!'
            stop
        end select

        GaussInt%Points_1D(:) = GaussInt%Points_1D(:)/2.0d0

    end subroutine setGaussPointsTriangle1D

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsTriangle2D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        allocate( GaussInt%Points_2D(GaussInt%NGauss_2D,2) )

        select case (GaussInt%NGauss_2D)
        case (1)
            !Precision 1
            GaussInt%Points_2D(1,1) = 1.0d0/3.0d0

            GaussInt%Points_2D(1,2) = 1.0d0/3.0d0

        case (3)
            !Precision 2
            GaussInt%Points_2D(1,1) = 0.5d0
            GaussInt%Points_2D(2,1) = 0.5d0
            GaussInt%Points_2D(3,1) = 0.0d0 + 1.0d-8

            GaussInt%Points_2D(1,2) = 0.0d0 + 1.0d-8
            GaussInt%Points_2D(2,2) = 0.5d0
            GaussInt%Points_2D(3,2) = 0.5d0

        case (4)
            !Precision 3
            GaussInt%Points_2D(1,1) = 1.0d0/3.0d0
            GaussInt%Points_2D(2,1) = 0.2d0
            GaussInt%Points_2D(3,1) = 0.2d0
            GaussInt%Points_2D(4,1) = 0.6d0

            GaussInt%Points_2D(1,2) = 1.0d0/3.0d0
            GaussInt%Points_2D(2,2) = 0.2d0
            GaussInt%Points_2D(3,2) = 0.6d0
            GaussInt%Points_2D(4,2) = 0.2d0

        case (6)
            !Precision 4
            GaussInt%Points_2D(1,1) = 0.816847572980459d0
            GaussInt%Points_2D(2,1) = 0.091576213509771d0
            GaussInt%Points_2D(3,1) = GaussInt%Points_2D(2,1)
            GaussInt%Points_2D(4,1) = 0.108103018168070d0
            GaussInt%Points_2D(5,1) = 0.445948490915965d0
            GaussInt%Points_2D(6,1) = GaussInt%Points_2D(5,1)

            GaussInt%Points_2D(1,2) = 0.091576213509771d0
            GaussInt%Points_2D(2,2) = 0.816847572980459d0
            GaussInt%Points_2D(3,2) = GaussInt%Points_2D(1,2)
            GaussInt%Points_2D(4,2) = 0.445948490915965d0
            GaussInt%Points_2D(5,2) = 0.108103018168070d0
            GaussInt%Points_2D(6,2) = GaussInt%Points_2D(4,2)

        case (7)
            !Precision 5
            GaussInt%Points_2D(1,1) = 1.0d0/3.0d0
            GaussInt%Points_2D(2,1) = 0.10128650732345633d0
            GaussInt%Points_2D(3,1) = 0.79742698535308720d0
            GaussInt%Points_2D(4,1) = GaussInt%Points_2D(2,1)
            GaussInt%Points_2D(5,1) = 0.47014206410511505d0
            GaussInt%Points_2D(6,1) = 0.05971587178976981d0
            GaussInt%Points_2D(7,1) = GaussInt%Points_2D(5,1)

            GaussInt%Points_2D(1,2) = 1.0d0/3.0d0
            GaussInt%Points_2D(2,2) = GaussInt%Points_2D(2,1)
            GaussInt%Points_2D(3,2) = GaussInt%Points_2D(2,1)
            GaussInt%Points_2D(4,2) = GaussInt%Points_2D(3,1)
            GaussInt%Points_2D(5,2) = GaussInt%Points_2D(5,1)
            GaussInt%Points_2D(6,2) = GaussInt%Points_2D(5,1)
            GaussInt%Points_2D(7,2) = GaussInt%Points_2D(6,1)

        case default
            write(*,'(a)') 'Wrong NGauss_2D in setGaussPointsTriangle2D!'
            stop
        end select

    end subroutine setGaussPointsTriangle2D

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsTriangleLC(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        allocate( GaussInt%Points_LC(GaussInt%NGauss_LC) )

        select case (GaussInt%NGauss_LC)

        case (6)

            GaussInt%Points_LC(1) = 1.0d0-0.9324695142031520278123d0
            GaussInt%Points_LC(2) = 1.0d0-0.6612093864662645136610d0
            GaussInt%Points_LC(3) = 1.0d0-0.2386191860831969086305d0
            GaussInt%Points_LC(4) = 1.0d0+0.2386191860831969086305d0
            GaussInt%Points_LC(5) = 1.0d0+0.6612093864662645136610d0
            GaussInt%Points_LC(6) = 1.0d0+0.9324695142031520278123d0

        case default
            write(*,'(a)') 'Wrong NGauss_LC in setGaussPointsTriangleLC!'
            stop
        end select

        GaussInt%Points_LC(:) = GaussInt%Points_LC(:)/2.0d0

    end subroutine setGaussPointsTriangleLC

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsTriangle(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        call setGaussWeightsTriangle1D(GaussInt)

        call setGaussWeightsTriangle2D(GaussInt)

        call setGaussWeightsTriangleLC(GaussInt)

        allocate( GaussInt%Weights_b(GaussInt%NGauss_2D) )
        GaussInt%Weights_b(:) = GaussInt%Weights_2D(:)

    end subroutine setGaussWeightsTriangle

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsTriangle1D(GaussInt)

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
        case default
            write(*,'(a)') 'Wrong NGauss_1D in setGaussWeightsTriangle1D!'
            stop
        end select

        GaussInt%Weights_1D(:) = GaussInt%Weights_1D(:)/2.0d0

    end subroutine setGaussWeightsTriangle1D

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsTriangle2D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        allocate( GaussInt%Weights_2D(GaussInt%NGauss_2D) )

        select case (GaussInt%NGauss_2D)

        case (1)
            !Precision 1
            GaussInt%Weights_2D(1) = 1.0d0

        case (3)
            !Precision 2
            GaussInt%Weights_2D(1) = 1.0d0/3.0d0
            GaussInt%Weights_2D(2) = 1.0d0/3.0d0
            GaussInt%Weights_2D(3) = 1.0d0/3.0d0

        case (4)
            !Precision 3
            GaussInt%Weights_2D(1) = -27.0d0/48.0d0
            GaussInt%Weights_2D(2) =  25.0d0/48.0d0
            GaussInt%Weights_2D(3) =  25.0d0/48.0d0
            GaussInt%Weights_2D(4) =  25.0d0/48.0d0

        case (6)
            !Precision 4
            GaussInt%Weights_2D(1:3) = 0.109951743655322d0
            GaussInt%Weights_2D(4:6) = 0.223381589678011d0

        case (7)
            !Precision 5
            GaussInt%Weights_2D(1) = 0.225d0
            GaussInt%Weights_2D(2:4) = 0.12593918054482717d0
            GaussInt%Weights_2D(5:7) = 0.13239415278850616d0

        case default
            write(*,'(a)') 'Wrong NGauss_2D in setGaussWeightsTriangle2D!'
            stop
        end select
        
        GaussInt%Weights_2D(:) = GaussInt%Weights_2D(:)/2.0d0

    end subroutine setGaussWeightsTriangle2D

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsTriangleLC(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        allocate( GaussInt%Weights_LC(GaussInt%NGauss_LC) )

        select case (GaussInt%NGauss_LC)

        case (6)

            GaussInt%Weights_LC(1) = 0.1713244923791703450403d0
            GaussInt%Weights_LC(2) = 0.3607615730481386075698d0
            GaussInt%Weights_LC(3) = 0.4679139345726910473899d0
            GaussInt%Weights_LC(4) = 0.4679139345726910473899d0
            GaussInt%Weights_LC(5) = 0.3607615730481386075698d0
            GaussInt%Weights_LC(6) = 0.1713244923791703450403d0

        case default
            write(*,'(a)') 'Wrong NGauss_LC in setGaussWeightsTriangleLC!'
            stop
        end select

        GaussInt%Weights_LC(:) = GaussInt%Weights_LC(:)/2.0d0

    end subroutine setGaussWeightsTriangleLC

    !----------------------------------------------------------------------- 

    subroutine setGaussBasisFunctionsTriangle(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        call setGaussBasisFunctionsTriangleEdge(GaussInt)
        call setGaussBasisFunctionsTriangle2D(GaussInt)

    end subroutine setGaussBasisFunctionsTriangle

    !----------------------------------------------------------------------- 

    subroutine setGaussBasisFunctionsTriangleEdge(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        PetscInt :: ig, ied
        PetscReal :: c, e
        PetscReal, dimension(2) :: xi
        PetscReal, dimension(GaussInt%Nbf) :: bfn, dfdc, dfde, d2fdc2, d2fde2, d2fdce   !C-> ξ, E-> η

        allocate( GaussInt%Bfn_E(GaussInt%Nbf, GaussInt%NGAUSS_1D, GaussInt%Ned) )
        allocate( GaussInt%dfdc_E(GaussInt%Nbf, GaussInt%NGAUSS_1D, GaussInt%Ned) )
        allocate( GaussInt%dfde_E(GaussInt%Nbf, GaussInt%NGAUSS_1D, GaussInt%Ned) )
        allocate( GaussInt%d2fdc2_E(GaussInt%Nbf, GaussInt%NGAUSS_1D, GaussInt%Ned) )
        allocate( GaussInt%d2fde2_E(GaussInt%Nbf, GaussInt%NGAUSS_1D, GaussInt%Ned) )
        allocate( GaussInt%d2fdce_E(GaussInt%Nbf, GaussInt%NGAUSS_1D, GaussInt%Ned) )

        do ig = 1, GaussInt%NGauss_1D
            do ied = 1, GaussInt%Ned

                select case (ied)
                case (1)
                    c = GaussInt%Points_1D(ig)
                    e = 0.0d0
                case (2)
                    c = 0.0d0
                    e = GaussInt%Points_1D(ig)
                case (3)
                    c = GaussInt%Points_1D(ig)
                    e = 1.0d0-c
                case default
                    write(*,'(a)') 'Wrong ied in setGaussBasisFunctionsTriangleEdge!'
                    stop
                end select

                xi = [c, e]

                select case (GaussInt%order)
                case ('Linear')
                    call setLinearTriangle(bfn, dfdc, dfde, c, e)
                    ! call setLinearTriangle_Bfn(xi, bfn)
                    ! call setTriangle_dBfn(setLinearTriangle_Bfn, xi, dfdc, dfde)
                    d2fdc2(:) = 0.0d0 ; d2fdce(:) = 0.0d0 ; d2fde2(:) = 0.0d0
                case ('Quadratic', 'Serendipity')
                    call setQuadraticTriangle(bfn, dfdc, dfde, d2fdc2, d2fde2, d2fdce, c, e)
                    ! call setQuadraticTriangle_Bfn(xi, bfn)
                    ! call setTriangle_dBfn(setQuadraticTriangle_Bfn, xi, dfdc, dfde)
                    ! call setTriangle_d2Bfn(setQuadraticTriangle_Bfn, xi, d2fdc2, d2fde2, d2fdce)
                case default
                    write(*,'(a)') 'Wrong GaussInt%order in setGaussBasisFunctionsTriangleEdge!'
                    stop
                end select
                GaussInt%Bfn_E(:,ig,ied)  = bfn(:)
                GaussInt%dfdc_E(:,ig,ied) = dfdc(:)
                GaussInt%dfde_E(:,ig,ied) = dfde(:)
                GaussInt%d2fdc2_E(:,ig,ied) = d2fdc2(:)
                GaussInt%d2fde2_E(:,ig,ied) = d2fde2(:)
                GaussInt%d2fdce_E(:,ig,ied) = d2fdce(:)

            end do
        end do

    end subroutine setGaussBasisFunctionsTriangleEdge

    !----------------------------------------------------------------------- 

    subroutine setGaussBasisFunctionsTriangle2D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        PetscInt :: ig
        PetscReal :: c, e
        PetscReal, dimension(2) :: xi
        PetscReal, dimension(GaussInt%Nbf) :: bfn, dfdc, dfde, d2fdc2, d2fde2, d2fdce   !C-> ξ, E-> η

        allocate( GaussInt%Bfn(GaussInt%Nbf, GaussInt%NGAUSS_2D) )
        allocate( GaussInt%dfdc(GaussInt%Nbf, GaussInt%NGAUSS_2D) )
        allocate( GaussInt%dfde(GaussInt%Nbf, GaussInt%NGAUSS_2D) )
        allocate( GaussInt%d2fdc2(GaussInt%Nbf, GaussInt%NGAUSS_2D) )
        allocate( GaussInt%d2fde2(GaussInt%Nbf, GaussInt%NGAUSS_2D) )
        allocate( GaussInt%d2fdce(GaussInt%Nbf, GaussInt%NGAUSS_2D) )

        do ig = 1, GaussInt%NGauss_2D

            c = GaussInt%Points_2D(ig,1)
            e = GaussInt%Points_2D(ig,2)

            xi = [c, e]

            select case (GaussInt%order)
            case ('Linear')
                call setLinearTriangle(bfn, dfdc, dfde, c, e)
                ! call setLinearTriangle_Bfn(xi, bfn)
                ! call setTriangle_dBfn(setLinearTriangle_Bfn, xi, dfdc, dfde)
                d2fdc2(:) = 0.0d0 ; d2fdce(:) = 0.0d0 ; d2fde2(:) = 0.0d0
            case ('Quadratic', 'Serendipity')
                call setQuadraticTriangle(bfn, dfdc, dfde, d2fdc2, d2fde2, d2fdce, c, e)
                ! call setQuadraticTriangle_Bfn(xi, bfn)
                ! call setTriangle_dBfn(setQuadraticTriangle_Bfn, xi, dfdc, dfde)
                ! call setTriangle_d2Bfn(setQuadraticTriangle_Bfn, xi, d2fdc2, d2fde2, d2fdce)
            case default
                write(*,'(a)') 'Wrong GaussInt%order in setGaussBasisFunctionsTriangle2D!'
                stop
            end select
            GaussInt%Bfn(:,ig)  = bfn(:)
            GaussInt%dfdc(:,ig) = dfdc(:)
            GaussInt%dfde(:,ig) = dfde(:)
            GaussInt%d2fdc2(:,ig) = d2fdc2(:)
            GaussInt%d2fde2(:,ig) = d2fde2(:)
            GaussInt%d2fdce(:,ig) = d2fdce(:)

        end do

    end subroutine setGaussBasisFunctionsTriangle2D

    !----------------------------------------------------------------------- 

    subroutine setLinearTriangle(bfn, dfdc, dfde, c, e)

        implicit none

        PetscReal, intent(in)  :: c, e
        PetscReal, dimension(:), intent(out) :: bfn, dfdc, dfde

        !Basis function, φ
        bfn(1)  = 1.0d0-c-e
        bfn(2)  = c
        bfn(3)  = e

        !dφ/dξ
        dfdc(1) = -1.0d0
        dfdc(2) = +1.0d0
        dfdc(3) =  0.0d0

        !dφ/dη
        dfde(1) = -1.0d0
        dfde(2) =  0.0d0
        dfde(3) = +1.0d0

    end subroutine setLinearTriangle

    !----------------------------------------------------------------------- 

    subroutine setQuadraticTriangle(bfn, dfdc, dfde, d2fdc2, d2fde2, d2fdce, c, e)

        implicit none

        PetscReal, intent(in)  :: c, e
        PetscReal, dimension(:), intent(out) :: bfn, dfdc, dfde, d2fdc2, d2fde2, d2fdce

        !Basis function, φ
        bfn(1)  = (1.0d0-c-e)*(1.0d0-2.0d0*c-2.0d0*e)
        bfn(2)  = c*(2.0d0*c-1.0d0)
        bfn(3)  = e*(2.0d0*e-1.0d0)
        bfn(4)  = 4.0d0*c*(1.0d0-c-e)
        bfn(5)  = 4.0d0*c*e
        bfn(6)  = 4.0d0*e*(1.0d0-c-e)

        !dφ/dξ
        dfdc(1) = -3.0d0+4.0d0*c+4.0d0*e
        dfdc(2) = +4.0d0*c-1.0d0
        dfdc(3) = +0.0d0
        dfdc(4) = +4.0d0-8.0d0*c-4.0d0*e
        dfdc(5) = +4.0d0*e
        dfdc(6) = -4.0d0*e

        !d2φ/dξ2
        d2fdc2(1) = +4.0d0
        d2fdc2(2) = +4.0d0
        d2fdc2(3) = +0.0d0
        d2fdc2(4) = -8.0d0
        d2fdc2(5) = +0.0d0
        d2fdc2(6) = +0.0d0

        !d2φ/dξdη
        d2fdce(1) = +4.0d0
        d2fdce(2) = +0.0d0
        d2fdce(3) = +0.0d0
        d2fdce(4) = -4.0d0
        d2fdce(5) = +4.0d0
        d2fdce(6) = -4.0d0

        !dφ/dη
        dfde(1) = -3.0d0+4.0d0*c+4.0d0*e
        dfde(2) = +0.0d0
        dfde(3) = +4.0d0*e-1.0d0
        dfde(4) = -4.0d0*c
        dfde(5) = +4.0d0*c
        dfde(6) = +4.0d0-4.0d0*c-8.0d0*e

        !d2φ/dη2
        d2fde2(1) = +4.0d0
        d2fde2(2) = +0.0d0
        d2fde2(3) = +4.0d0
        d2fde2(4) = +0.0d0
        d2fde2(5) = +0.0d0
        d2fde2(6) = -8.0d0

    end subroutine setQuadraticTriangle

    !----------------------------------------------------------------------- 

    subroutine setLinearTriangle_Bfn(xi, bfn)

        implicit none

        PetscReal, dimension(:), intent(in)  :: xi
        PetscReal, dimension(:), intent(out) :: bfn

        PetscReal :: c, e

        c = xi(1) ; e = xi(2)

        !Basis function, φ
        bfn(1)  = 1.0d0-c-e
        bfn(2)  = c
        bfn(3)  = e

    end subroutine setLinearTriangle_Bfn

    !----------------------------------------------------------------------- 

    subroutine setQuadraticTriangle_Bfn(xi, bfn)

        implicit none

        PetscReal, dimension(:), intent(in)  :: xi
        PetscReal, dimension(:), intent(out) :: bfn

        PetscReal :: c, e

        c = xi(1) ; e = xi(2)

        !Basis function, φ
        bfn(1)  = (1.0d0-c-e)*(1.0d0-2.0d0*c-2.0d0*e)
        bfn(2)  = c*(2.0d0*c-1.0d0)
        bfn(3)  = e*(2.0d0*e-1.0d0)
        bfn(4)  = 4.0d0*c*(1.0d0-c-e)
        bfn(5)  = 4.0d0*c*e
        bfn(6)  = 4.0d0*e*(1.0d0-c-e)

    end subroutine setQuadraticTriangle_Bfn

    !-----------------------------------------------------------------------

    subroutine setTriangle_dBfn(f_sub, xi, dfdc, dfde)

        implicit none

        interface
            subroutine f_sub(xi, f)
                PetscReal, dimension(:), intent(in)  :: xi
                PetscReal, dimension(:), intent(out) :: f
            end subroutine f_sub
        end interface
        PetscReal, dimension(:), intent(in)  :: xi
        PetscReal, dimension(:), intent(out) :: dfdc, dfde

        call derivative(f_sub, xi, 1, dfdc)
        call derivative(f_sub, xi, 2, dfde)

    end subroutine setTriangle_dBfn

    !-----------------------------------------------------------------------

    subroutine setTriangle_d2Bfn(f_sub, xi, d2fdc2, d2fde2, d2fdce)

        implicit none

        interface
            subroutine f_sub(xi, f)
                PetscReal, dimension(:), intent(in)  :: xi
                PetscReal, dimension(:), intent(out) :: f
            end subroutine f_sub
        end interface
        PetscReal, dimension(:), intent(in)  :: xi
        PetscReal, dimension(:), intent(out) :: d2fdc2, d2fde2, d2fdce

        call derivative2(f_sub, xi, 1, d2fdc2)
        call derivative2(f_sub, xi, 2, d2fde2)

        call derivative2_mixed(f_sub, xi, 1, 2, d2fdce)

    end subroutine setTriangle_d2Bfn

    !----------------------------------------------------------------------- 

    subroutine derivative(f_sub, xi, k, df)

        implicit none

        interface
            subroutine f_sub(xi, f)
                PetscReal, dimension(:), intent(in)  :: xi
                PetscReal, dimension(:), intent(out) :: f
            end subroutine f_sub
        end interface
        PetscReal, dimension(:), intent(in) :: xi
        PetscInt, intent(in) :: k
        PetscReal, dimension(:), intent(out) :: df

        PetscReal, parameter :: eps = 1.0d-8
        PetscReal, dimension(size(df)) :: f1, f2
        PetscReal, dimension(size(xi)) :: dx

        call f_sub(xi,f1)

        dx(:) = xi(:)
        dx(k) = dx(k) + eps

        call f_sub(dx, f2)

        df(:) = (f2(:)-f1(:))/eps

    end subroutine derivative

    !----------------------------------------------------------------------- 

    subroutine derivative2(f_sub, xi, k, d2f)

        implicit none

        interface
            subroutine f_sub(xi, f)
                PetscReal, dimension(:), intent(in)  :: xi
                PetscReal, dimension(:), intent(out) :: f
            end subroutine f_sub
        end interface
        PetscReal, dimension(:), intent(in) :: xi
        PetscInt, intent(in) :: k
        PetscReal, dimension(:), intent(out) :: d2f

        PetscReal, parameter :: eps = 1.0d-5
        PetscReal, dimension(size(d2f)) :: f1, f2, f3, f4, f5
        PetscReal, dimension(size(xi)) :: dx

        !Forward
        call f_sub(xi,f1)

        dx(:) = xi(:)
        dx(k) = dx(k) + eps
        call f_sub(dx, f2)

        dx(k) = dx(k) + eps
        call f_sub(dx, f3)

        dx(k) = dx(k) + eps
        call f_sub(dx, f4)

        ! d2f(:) = (f3(:)-2.0d0*f2(:)+f1(:))/(eps**2)
        d2f(:) = (-f4(:)+4.0d0*f3(:)-5.0d0*f2(:)+2.0d0*f1(:))/(eps**2)

    end subroutine derivative2
    
    !----------------------------------------------------------------------- 

    subroutine derivative2_mixed(f_sub, xi, j, k, d2f)

        implicit none

        interface
            subroutine f_sub(xi, f)
                PetscReal, dimension(:), intent(in)  :: xi
                PetscReal, dimension(:), intent(out) :: f
            end subroutine f_sub
        end interface
        PetscReal, dimension(:), intent(in) :: xi
        PetscInt, intent(in) :: j, k
        PetscReal, dimension(:), intent(out) :: d2f

        PetscReal, parameter :: eps = 1.0d-5
        PetscReal, dimension(size(d2f)) :: f1, f2, f3, f4
        PetscReal, dimension(size(xi)) :: dx

        !Forward
        call f_sub(xi,f1)

        dx(:) = xi(:)

        dx(j) = dx(j) + eps
        call f_sub(dx, f2)

        dx(:) = xi(:)
        dx(k) = dx(k) + eps
        call f_sub(dx, f3)

        dx(:) = xi(:)
        dx(j) = dx(j) + eps
        dx(k) = dx(k) + eps
        call f_sub(dx, f4)

        d2f(:) = (f4(:)-f3(:)-f2(:)+f1(:))/(eps**2)

    end subroutine derivative2_mixed

end module GaussTriangle_mod
