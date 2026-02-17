
#include <petsc/finclude/petscksp.h>

module GaussTetrahedron_mod

    use GaussParameters_mod, only: GaussIntegration

    implicit none

    PetscInt, parameter :: NGauss_1D = 2, NGauss_2D = 4, NGauss_3D = 5, NGauss_LC = 6

    contains

    subroutine setGaussParametersTetrahedron(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt
        
        GaussInt%NGauss_1D = NGauss_1D
        GaussInt%NGauss_2D = NGauss_2D
        GaussInt%NGauss_3D = NGauss_3D
        if (GaussInt%order == 'Quadratic' .or. GaussInt%order == 'Serendipity') then
            GaussInt%NGauss_1D = 3
            GaussInt%NGauss_2D = 7
            GaussInt%NGauss_3D = 11
        end if
        GaussInt%NGauss_LC = NGauss_LC

        GaussInt%NGauss_b = GaussInt%NGauss_3D

        call setGaussPointsTetrahedron(GaussInt)
        call setGaussWeightsTetrahedron(GaussInt)
        call setGaussBasisFunctionsTetrahedron(GaussInt)

    end subroutine setGaussParametersTetrahedron

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsTetrahedron(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        call setGaussPointsTetrahedron1D(GaussInt)

        call setGaussPointsTetrahedron2D(GaussInt)

        call setGaussPointsTetrahedron3D(GaussInt)

        call setGaussPointsTetrahedronLC(GaussInt)

    end subroutine setGaussPointsTetrahedron

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsTetrahedron1D(GaussInt)

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
            write(*,'(a)') 'Wrong NGauss_1D in setGaussPointsTetrahedron1D!'
            stop
        end select

        GaussInt%Points_1D(:) = GaussInt%Points_1D(:)/2.0d0

    end subroutine setGaussPointsTetrahedron1D

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsTetrahedron2D(GaussInt)

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
            GaussInt%Points_2D(2,1) = 0.79742698535308720d0
            GaussInt%Points_2D(3,1) = 0.10128650732345633d0
            GaussInt%Points_2D(4,1) = GaussInt%Points_2D(3,1)
            GaussInt%Points_2D(5,1) = 0.05971587178976981d0
            GaussInt%Points_2D(6,1) = 0.47014206410511505d0
            GaussInt%Points_2D(7,1) = GaussInt%Points_2D(6,1)

            GaussInt%Points_2D(1,2) = 1.0d0/3.0d0
            GaussInt%Points_2D(2,2) = 0.10128650732345633d0
            GaussInt%Points_2D(3,2) = 0.79742698535308720d0
            GaussInt%Points_2D(4,2) = GaussInt%Points_2D(2,2)
            GaussInt%Points_2D(5,2) = 0.47014206410511505d0
            GaussInt%Points_2D(6,2) = 0.05971587178976981d0
            GaussInt%Points_2D(7,2) = GaussInt%Points_2D(5,2)

        case default
            write(*,'(a)') 'Wrong NGauss_2D in setGaussPointsTetrahedron2D!'
            stop
        end select

    end subroutine setGaussPointsTetrahedron2D

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsTetrahedron3D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        allocate( GaussInt%Points_3D(GaussInt%NGauss_3D,3) )

        select case (GaussInt%NGauss_3D)

        case (4)
            !Precision 2
            GaussInt%Points_3D(:,1) = [0.5854101966249685d0, 0.1381966011250105d0, &
                0.1381966011250105d0, 0.1381966011250105d0]

            GaussInt%Points_3D(:,2) = [0.1381966011250105d0, 0.5854101966249685d0, &
                0.1381966011250105d0, 0.1381966011250105d0]
   
            GaussInt%Points_3D(:,3) = [0.1381966011250105d0, 0.1381966011250105d0, &
                0.5854101966249685d0, 0.1381966011250105d0]

        case (5)
            !Precision 3
            GaussInt%Points_3D(:,1) = [0.25d0, 1.0d0/6.0d0, 1.0d0/6.0d0, 1.0d0/6.0d0, 0.50d0]

            GaussInt%Points_3D(:,2) = [0.25d0, 1.0d0/6.0d0, 1.0d0/6.0d0, 0.5d0, 1.0d0/6.0d0]
   
            GaussInt%Points_3D(:,3) = [0.25d0, 1.0d0/6.0d0, 0.50d0, 1.0d0/6.0d0, 1.0d0/6.0d0]

        case (10)
            !Precision 3
            GaussInt%Points_3D(:,1) = [0.5684305841968444d0, 0.1438564719343852d0, &
                0.1438564719343852d0, 0.1438564719343852d0, 0.0d0, 0.5d0, 0.5d0, &
                0.5d0, 0.0d0, 0.0d0]

            GaussInt%Points_3D(:,2) = [0.1438564719343852d0, 0.1438564719343852d0, &
                0.1438564719343852d0, 0.5684305841968444d0, 0.5d0, 0.0d0, 0.5d0, &
                0.0d0, 0.5d0, 0.0d0]

            GaussInt%Points_3D(:,3) = [0.1438564719343852d0, 0.1438564719343852d0, &
                0.5684305841968444d0, 0.1438564719343852d0, 0.5d0, 0.5d0, 0.0d0, &
                0.0d0, 0.0d0, 0.5d0]

        case (11)
            !Precision 4
            GaussInt%Points_3D(:,1) = [0.25d0, 0.7857142857142857d0, 0.0714285714285714d0, &
                0.0714285714285714d0, 0.0714285714285714d0, 0.1005964238332008d0, &
                0.3994035761667992d0, 0.3994035761667992d0, 0.3994035761667992d0, &
                0.1005964238332008d0, 0.1005964238332008d0]

            GaussInt%Points_3D(:,2) = [0.25d0, 0.0714285714285714d0, 0.0714285714285714d0, &
                0.0714285714285714d0, 0.7857142857142857d0, 0.3994035761667992d0, &
                0.1005964238332008d0, 0.3994035761667992d0, 0.1005964238332008d0, &
                0.3994035761667992d0, 0.1005964238332008d0]

            GaussInt%Points_3D(:,3) = [0.25d0, 0.0714285714285714d0, 0.0714285714285714d0, &
                0.7857142857142857d0, 0.0714285714285714d0, 0.3994035761667992d0, &
                0.3994035761667992d0, 0.1005964238332008d0, 0.1005964238332008d0, &
                0.1005964238332008d0, 0.3994035761667992d0]

        case (15)
            !Precision 5
            GaussInt%Points_3D(:,1) = [0.25d0, 0.00d0, 1.0d0/3.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0, &
                0.7272727272727273d0, 1.0d0/11.0d0, 1.0d0/11.0d0, 1.0d0/11.0d0, &
                0.4334498464263357d0, 0.0665501535736643d0, 0.0665501535736643d0, &
                0.0665501535736643d0, 0.4334498464263357d0, 0.4334498464263357d0]

            GaussInt%Points_3D(:,2) = [0.25d0, 1.0d0/3.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0, 0.0d0, &
                1.0d0/11.0d0, 1.0d0/11.0d0, 1.0d0/11.0d0, 0.7272727272727273d0, &
                0.0665501535736643d0, 0.4334498464263357d0, 0.0665501535736643d0, &
                0.4334498464263357d0, 0.0665501535736643d0, 0.4334498464263357d0]

            GaussInt%Points_3D(:,3) = [0.25d0, 1.0d0/3.0d0, 1.0d0/3.0d0, 0.0d0, 1.0d0/3.0d0, &
                1.0d0/11.0d0, 1.0d0/11.0d0, 0.7272727272727273d0, 1.0d0/11.0d0, &
                0.0665501535736643d0, 0.0665501535736643d0, 0.4334498464263357d0, &
                0.4334498464263357d0, 0.4334498464263357d0, 0.0665501535736643d0]

        case default
            write(*,'(a)') 'Wrong GaussInt%Points_3D in setGaussPointsTetrahedron3D!'
            stop
        end select
    
    end subroutine setGaussPointsTetrahedron3D

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsTetrahedronLC(GaussInt)

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
            write(*,'(a)') 'Wrong GaussInt%NGauss_LC in setGaussPointsTetrahedronLC!'
            stop
        end select

        GaussInt%Points_LC(:) = GaussInt%Points_LC(:)/2.0d0

    end subroutine setGaussPointsTetrahedronLC

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsTetrahedron(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        call setGaussWeightsTetrahedron1D(GaussInt)

        call setGaussWeightsTetrahedron2D(GaussInt)

        call setGaussWeightsTetrahedron3D(GaussInt)

        call setGaussWeightsTetrahedronLC(GaussInt)

        allocate( GaussInt%Weights_b(GaussInt%NGauss_3D) )
        GaussInt%Weights_b(:) = GaussInt%Weights_3D(:)

    end subroutine setGaussWeightsTetrahedron

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsTetrahedron1D(GaussInt)

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
            write(*,'(a)') 'Wrong GaussInt%NGauss_1D in setGaussWeightsTetrahedron1D!'
            stop
        end select

        GaussInt%Weights_1D(:) = GaussInt%Weights_1D(:)/2.0d0

    end subroutine setGaussWeightsTetrahedron1D

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsTetrahedron2D(GaussInt)

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
            write(*,'(a)') 'Wrong GaussInt%NGauss_2D in setGaussWeightsTetrahedron2D!'
            stop
        end select
        
        GaussInt%Weights_2D(:) = GaussInt%Weights_2D(:)/2.0d0

    end subroutine setGaussWeightsTetrahedron2D

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsTetrahedron3D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        allocate( GaussInt%Weights_3D(GaussInt%NGauss_3D) )

        select case (GaussInt%NGauss_3D)

        case (4)
            !Precision 2
            GaussInt%Weights_3D(:) = 0.25d0

        case (5)
            !Precision 3
            GaussInt%Weights_3D(:) = [-0.80d0, 0.45d0, 0.45d0, 0.45d0, 0.45d0]

        case (10)
            !Precision 3
            GaussInt%Weights_3D(:) = [0.2177650698804054d0, 0.2177650698804054d0, &
                0.2177650698804054d0, 0.2177650698804054d0, 0.0214899534130631d0, &
                0.0214899534130631d0, 0.0214899534130631d0, 0.0214899534130631d0, &
                0.0214899534130631d0, 0.0214899534130631d0]

        case (11)
            !Precision 4
            GaussInt%Weights_3D(:) = [-0.0789333333333333d0, 0.0457333333333333d0, &
                0.0457333333333333d0, 0.0457333333333333d0, 0.0457333333333333d0, &
                0.1493333333333333d0, 0.1493333333333333d0, 0.1493333333333333d0, &
                0.1493333333333333d0, 0.1493333333333333d0, 0.1493333333333333d0]

        case (15)
            !Precision 5
            GaussInt%Weights_3D(:) = [0.1817020685825351d0, 0.0361607142857143d0, &
                0.0361607142857143d0, 0.0361607142857143d0, 0.0361607142857143d0, &
                0.0698714945161738d0, 0.0698714945161738d0, 0.0698714945161738d0, &
                0.0698714945161738d0, 0.0656948493683187d0, 0.0656948493683187d0, &
                0.0656948493683187d0, 0.0656948493683187d0, 0.0656948493683187d0, &
                0.0656948493683187d0]

        case default
            write(*,'(a)') 'Wrong GaussInt%NGauss_3D in setGaussWeightsTetrahedron3D!'
            stop
        end select

        GaussInt%Weights_3D(:) = GaussInt%Weights_3D(:)/6.0d0

    end subroutine setGaussWeightsTetrahedron3D

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsTetrahedronLC(GaussInt)

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
            write(*,'(a)') 'Wrong GaussInt%NGauss_LC in setGaussWeightsTetrahedronLC!'
            stop
        end select

        GaussInt%Weights_LC(:) = GaussInt%Weights_LC(:)/2.0d0

    end subroutine setGaussWeightsTetrahedronLC

    !----------------------------------------------------------------------- 

    subroutine setGaussBasisFunctionsTetrahedron(GaussInt)

        implicit none
        
        type(GaussIntegration), intent(inout) :: GaussInt

        call setGaussBasisFunctionsTetrahedronEdge(GaussInt)
        call setGaussBasisFunctionsTetrahedronFace(GaussInt)
        call setGaussBasisFunctionsTetrahedron3D(GaussInt)

    end subroutine setGaussBasisFunctionsTetrahedron

    !----------------------------------------------------------------------- 

    subroutine setGaussBasisFunctionsTetrahedronEdge(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        PetscInt :: ig, ied
        PetscReal :: c, e, s    !C-> ξ, E-> η, s -> ζ
        PetscReal, dimension(3) :: xi
        PetscReal, dimension(GaussInt%Nbf) :: bfn, dfdc, dfde, dfds, &
            d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes

        allocate( GaussInt%Bfn_E(GaussInt%Nbf, GaussInt%NGAUSS_1D, GaussInt%Ned) )
        allocate( GaussInt%dfdc_E(GaussInt%Nbf, GaussInt%NGAUSS_1D, GaussInt%Ned) )
        allocate( GaussInt%dfde_E(GaussInt%Nbf, GaussInt%NGAUSS_1D, GaussInt%Ned) )
        allocate( GaussInt%dfds_E(GaussInt%Nbf, GaussInt%NGAUSS_1D, GaussInt%Ned) )
        allocate( GaussInt%d2fdc2_E(GaussInt%Nbf, GaussInt%NGAUSS_1D, GaussInt%Ned) )
        allocate( GaussInt%d2fde2_E(GaussInt%Nbf, GaussInt%NGAUSS_1D, GaussInt%Ned) )
        allocate( GaussInt%d2fds2_E(GaussInt%Nbf, GaussInt%NGAUSS_1D, GaussInt%Ned) )
        allocate( GaussInt%d2fdce_E(GaussInt%Nbf, GaussInt%NGAUSS_1D, GaussInt%Ned) )
        allocate( GaussInt%d2fdcs_E(GaussInt%Nbf, GaussInt%NGAUSS_1D, GaussInt%Ned) )
        allocate( GaussInt%d2fdes_E(GaussInt%Nbf, GaussInt%NGAUSS_1D, GaussInt%Ned) )

        do ig = 1, GaussInt%NGauss_1D
            do ied = 1, GaussInt%Ned

                select case (ied)
                case (1)
                    c  = GaussInt%Points_1D(ig)
                    e  = 0.0d0
                    s = 0.0d0
                case (2)
                    c  = 0.0d0
                    e  = GaussInt%Points_1D(ig)
                    s = 0.0d0
                case (3)
                    c  = 0.0d0
                    e  = 0.0d0
                    s = GaussInt%Points_1D(ig)
                case (4)
                    c  = GaussInt%Points_1D(ig)
                    e  = 1.0d0-c
                    s = 0.0d0
                case (5)
                    c  = 0.0d0
                    e  = GaussInt%Points_1D(ig)
                    s = 1.0d0-e
                case (6)
                    c  = GaussInt%Points_1D(ig)
                    e  = 0.0d0
                    s = 1.0d0-c
                case default
                    write(*,'(a)') 'Wrong ied in setGaussBasisFunctionsTetrahedronEdge!'
                    stop
                end select

                xi = [c, e, s]

                select case (GaussInt%order)
                case ('Linear')
                    call setLinearTetrahedron(bfn, dfdc, dfde, dfds, c, e, s)
                    ! call setLinearTetrahedron_Bfn(xi, bfn)
                    ! call setTetrahedron_dBfn(setLinearTetrahedron_Bfn, xi, dfdc, dfde, dfds)
                    d2fdc2(:) = 0.0d0 ; d2fde2(:) = 0.0d0 ; d2fds2(:) = 0.0d0
                    d2fdce(:) = 0.0d0 ; d2fdcs(:) = 0.0d0 ; d2fdes(:) = 0.0d0
                case ('Quadratic', 'Serendipity')
                    call setQuadraticTetrahedron(bfn, dfdc, dfde, dfds, &
                                    d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes, c, e, s)
                    ! call setQuadraticTetrahedron_Bfn(xi, bfn)
                    ! call setTetrahedron_dBfn(setQuadraticTetrahedron_Bfn, xi, dfdc, dfde, dfds)
                    ! call setTetrahedron_d2Bfn(setQuadraticTetrahedron_Bfn, xi, d2fdc2, d2fde2, d2fds2, &
                    !                                 d2fdce, d2fdcs, d2fdes)
                case default
                    write(*,'(a)') 'Wrong GaussInt%order in setGaussBasisFunctionsTetrahedronEdge!'
                    stop
                end select
                GaussInt%Bfn_E(:,ig,ied)  = bfn(:)
                GaussInt%dfdc_E(:,ig,ied) = dfdc(:)
                GaussInt%dfde_E(:,ig,ied) = dfde(:)
                GaussInt%dfds_E(:,ig,ied) = dfds(:)
                GaussInt%d2fdc2_E(:,ig,ied) = d2fdc2(:)
                GaussInt%d2fde2_E(:,ig,ied) = d2fde2(:)
                GaussInt%d2fds2_E(:,ig,ied) = d2fds2(:)
                GaussInt%d2fdce_E(:,ig,ied) = d2fdce(:)
                GaussInt%d2fdcs_E(:,ig,ied) = d2fdcs(:)
                GaussInt%d2fdes_E(:,ig,ied) = d2fdes(:)

            end do
        end do

    end subroutine setGaussBasisFunctionsTetrahedronEdge

    !----------------------------------------------------------------------- 

    subroutine setGaussBasisFunctionsTetrahedronFace(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        PetscInt :: ig, ifc
        PetscReal :: c, e, s    !C-> ξ, E-> η, s -> ζ
        PetscReal, dimension(3) :: xi
        PetscReal, dimension(GaussInt%Nbf) :: bfn, dfdc, dfde, dfds, &
            d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes

        allocate( GaussInt%Bfn_F(GaussInt%Nbf, GaussInt%NGAUSS_2D, GaussInt%Nfc) )
        allocate( GaussInt%dfdc_F(GaussInt%Nbf, GaussInt%NGAUSS_2D, GaussInt%Nfc) )
        allocate( GaussInt%dfde_F(GaussInt%Nbf, GaussInt%NGAUSS_2D, GaussInt%Nfc) )
        allocate( GaussInt%dfds_F(GaussInt%Nbf, GaussInt%NGAUSS_2D, GaussInt%Nfc) )
        allocate( GaussInt%d2fdc2_F(GaussInt%Nbf, GaussInt%NGAUSS_2D, GaussInt%Nfc) )
        allocate( GaussInt%d2fde2_F(GaussInt%Nbf, GaussInt%NGAUSS_2D, GaussInt%Nfc) )
        allocate( GaussInt%d2fds2_F(GaussInt%Nbf, GaussInt%NGAUSS_2D, GaussInt%Nfc) )
        allocate( GaussInt%d2fdce_F(GaussInt%Nbf, GaussInt%NGAUSS_2D, GaussInt%Nfc) )
        allocate( GaussInt%d2fdcs_F(GaussInt%Nbf, GaussInt%NGAUSS_2D, GaussInt%Nfc) )
        allocate( GaussInt%d2fdes_F(GaussInt%Nbf, GaussInt%NGAUSS_2D, GaussInt%Nfc) )

        do ig = 1, GaussInt%NGauss_2D
            do ifc = 1, GaussInt%Nfc

                select case (ifc)
                case (1)
                    c = 0.0d0
                    e = GaussInt%Points_2D(ig,1)
                    s = GaussInt%Points_2D(ig,2)
                case (2)
                    c = GaussInt%Points_2D(ig,1)
                    e = 0.0d0
                    s = GaussInt%Points_2D(ig,2)
                case (3)
                    c = GaussInt%Points_2D(ig,1)
                    e = GaussInt%Points_2D(ig,2)
                    s = 0.0d0
                case (4)
                    c = GaussInt%Points_2D(ig,1)
                    e = GaussInt%Points_2D(ig,2)
                    s = 1.0d0-c-e
                case default
                    write(*,'(a)') 'Wrong ifc in setGaussBasisFunctionsTetrahedronFace!'
                    stop
                end select

                xi = [c, e, s]

                select case (GaussInt%order)
                case ('Linear')
                    call setLinearTetrahedron(bfn, dfdc, dfde, dfds, c, e, s)
                    ! call setLinearTetrahedron_Bfn(xi, bfn)
                    ! call setTetrahedron_dBfn(setLinearTetrahedron_Bfn, xi, dfdc, dfde, dfds)
                    d2fdc2(:) = 0.0d0 ; d2fde2(:) = 0.0d0 ; d2fds2(:) = 0.0d0
                    d2fdce(:) = 0.0d0 ; d2fdcs(:) = 0.0d0 ; d2fdes(:) = 0.0d0
                case ('Quadratic', 'Serendipity')
                    call setQuadraticTetrahedron(bfn, dfdc, dfde, dfds, &
                                    d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes, c, e, s)
                    ! call setQuadraticTetrahedron_Bfn(xi, bfn)
                    ! call setTetrahedron_dBfn(setQuadraticTetrahedron_Bfn, xi, dfdc, dfde, dfds)
                    ! call setTetrahedron_d2Bfn(setQuadraticTetrahedron_Bfn, xi, d2fdc2, d2fde2, d2fds2, &
                    !                                 d2fdce, d2fdcs, d2fdes)
                case default
                    write(*,'(a)') 'Wrong GaussInt%order in setGaussBasisFunctionsTetrahedronFace!'
                    stop
                end select
                GaussInt%Bfn_F(:,ig,ifc)  = bfn(:)
                GaussInt%dfdc_F(:,ig,ifc) = dfdc(:)
                GaussInt%dfde_F(:,ig,ifc) = dfde(:)
                GaussInt%dfds_F(:,ig,ifc) = dfds(:)
                GaussInt%d2fdc2_F(:,ig,ifc) = d2fdc2(:)
                GaussInt%d2fde2_F(:,ig,ifc) = d2fde2(:)
                GaussInt%d2fds2_F(:,ig,ifc) = d2fds2(:)
                GaussInt%d2fdce_F(:,ig,ifc) = d2fdce(:)
                GaussInt%d2fdcs_F(:,ig,ifc) = d2fdcs(:)
                GaussInt%d2fdes_F(:,ig,ifc) = d2fdes(:)

            end do
        end do

    end subroutine setGaussBasisFunctionsTetrahedronFace

    !----------------------------------------------------------------------- 

    subroutine setGaussBasisFunctionsTetrahedron3D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        PetscInt :: ig
        PetscReal :: c, e, s    !C-> ξ, E-> η, s -> ζ
        PetscReal, dimension(3) :: xi
        PetscReal, dimension(GaussInt%Nbf) :: bfn, dfdc, dfde, dfds, &
            d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes

        allocate( GaussInt%Bfn(GaussInt%Nbf, GaussInt%NGAUSS_3D) )
        allocate( GaussInt%dfdc(GaussInt%Nbf, GaussInt%NGAUSS_3D) )
        allocate( GaussInt%dfde(GaussInt%Nbf, GaussInt%NGAUSS_3D) )
        allocate( GaussInt%dfds(GaussInt%Nbf, GaussInt%NGAUSS_3D) )
        allocate( GaussInt%d2fdc2(GaussInt%Nbf, GaussInt%NGAUSS_3D) )
        allocate( GaussInt%d2fde2(GaussInt%Nbf, GaussInt%NGAUSS_3D) )
        allocate( GaussInt%d2fds2(GaussInt%Nbf, GaussInt%NGAUSS_3D) )
        allocate( GaussInt%d2fdce(GaussInt%Nbf, GaussInt%NGAUSS_3D) )
        allocate( GaussInt%d2fdcs(GaussInt%Nbf, GaussInt%NGAUSS_3D) )
        allocate( GaussInt%d2fdes(GaussInt%Nbf, GaussInt%NGAUSS_3D) )

        do ig = 1, GaussInt%NGauss_3D

            c = GaussInt%Points_3D(ig,1)
            e = GaussInt%Points_3D(ig,2)
            s = GaussInt%Points_3D(ig,3)

            xi = [c, e, s]

            select case (GaussInt%order)
            case ('Linear')
                call setLinearTetrahedron(bfn, dfdc, dfde, dfds, c, e, s)
                ! call setLinearTetrahedron_Bfn(xi, bfn)
                ! call setTetrahedron_dBfn(setLinearTetrahedron_Bfn, xi, dfdc, dfde, dfds)
                d2fdc2(:) = 0.0d0 ; d2fde2(:) = 0.0d0 ; d2fds2(:) = 0.0d0
                d2fdce(:) = 0.0d0 ; d2fdcs(:) = 0.0d0 ; d2fdes(:) = 0.0d0
            case ('Quadratic', 'Serendipity')
                call setQuadraticTetrahedron(bfn, dfdc, dfde, dfds, &
                    d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes, c, e, s)
                ! call setQuadraticTetrahedron_Bfn(xi, bfn)
                ! call setTetrahedron_dBfn(setQuadraticTetrahedron_Bfn, xi, dfdc, dfde, dfds)
                ! call setTetrahedron_d2Bfn(setQuadraticTetrahedron_Bfn, xi, d2fdc2, d2fde2, d2fds2, &
                !                                 d2fdce, d2fdcs, d2fdes)
            case default
                write(*,'(a)') 'Wrong GaussInt%order in setGaussBasisFunctionsTetrahedron3D!'
                stop
            end select
            GaussInt%Bfn(:,ig)   = bfn(:)
            GaussInt%dfdc(:,ig)  = dfdc(:)
            GaussInt%dfde(:,ig)  = dfde(:)
            GaussInt%dfds(:,ig) = dfds(:)
            GaussInt%d2fdc2(:,ig) = d2fdc2(:)
            GaussInt%d2fde2(:,ig) = d2fde2(:)
            GaussInt%d2fds2(:,ig) = d2fds2(:)
            GaussInt%d2fdce(:,ig) = d2fdce(:)
            GaussInt%d2fdcs(:,ig) = d2fdcs(:)
            GaussInt%d2fdes(:,ig) = d2fdes(:)

        end do

    end subroutine setGaussBasisFunctionsTetrahedron3D

    !----------------------------------------------------------------------- 

    subroutine setLinearTetrahedron(bfn, dfdc, dfde, dfds, c, e, s)

        implicit none

        PetscReal, intent(in)  :: c, e, s
        PetscReal, dimension(:), intent(out) :: bfn, dfdc, dfde, dfds

        !Basis function, φ
        bfn(1)  = 1.0d0-c-e-s
        bfn(2)  = c
        bfn(3)  = e
        bfn(4)  = s
        !dφ/dξ
        dfdc(1) = -1.0d0
        dfdc(2) = +1.0d0
        dfdc(3) =  0.0d0
        dfdc(4) =  0.0d0
        !dφ/dη
        dfde(1) = -1.0d0
        dfde(2) =  0.0d0
        dfde(3) = +1.0d0
        dfde(4) =  0.0d0
        !dφ/dζ
        dfds(1) = -1.0d0
        dfds(2) =  0.0d0
        dfds(3) =  0.0d0
        dfds(4) = +1.0d0

    end subroutine setLinearTetrahedron

    !----------------------------------------------------------------------- 

    subroutine setQuadraticTetrahedron(bfn, dfdc, dfde, dfds, &
        d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes, c, e, s)

        implicit none

        PetscReal, intent(in)  :: c, e, s
        PetscReal, dimension(:), intent(out) :: bfn, dfdc, dfde, dfds, &
            d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes

        !Basis function, φ
        bfn(1)  = +(1.0d0-c-e-s)*(1.0d0-2.0d0*c-2.0d0*e-2.0d0*s)
        bfn(2)  = +c*(2.0d0*c-1.0d0)
        bfn(3)  = +e*(2.0d0*e-1.0d0)
        bfn(4)  = +s*(2.0d0*s-1.0d0)
        bfn(5)  = +4.0d0*(1.0d0-c-e-s)*c
        bfn(6)  = +4.0d0*c*e
        bfn(7)  = +4.0d0*(1.0d0-c-e-s)*e
        bfn(8)  = +4.0d0*(1.0d0-c-e-s)*s
        bfn(9)  = +4.0d0*c*s
        bfn(10) = +4.0d0*e*s

        !dφ/dξ
        dfdc(1)  = -3.0d0+4.0d0*(c+e+s)
        dfdc(2)  = +4.0d0*c-1.0d0
        dfdc(3)  = +0.0d0
        dfdc(4)  = +0.0d0
        dfdc(5)  = +4.0d0*(1.0d0-2.0d0*c-e-s)
        dfdc(6)  = +4.0d0*e
        dfdc(7)  = -4.0d0*e
        dfdc(8)  = -4.0d0*s
        dfdc(9)  = +4.0d0*s
        dfdc(10) = +0.0d0

        !d2φ/dξ2
        d2fdc2(1)  = +4.0d0
        d2fdc2(2)  = +4.0d0
        d2fdc2(3)  = +0.0d0
        d2fdc2(4)  = +0.0d0
        d2fdc2(5)  = -8.0d0
        d2fdc2(6)  = +0.0d0
        d2fdc2(7)  = +0.0d0
        d2fdc2(8)  = +0.0d0
        d2fdc2(9)  = +0.0d0
        d2fdc2(10) = +0.0d0

        !d2φ/dξdη
        d2fdce(1)  = +4.0d0
        d2fdce(2)  =  0.0d0
        d2fdce(3)  = +0.0d0
        d2fdce(4)  = +0.0d0
        d2fdce(5)  = -4.0d0
        d2fdce(6)  = +4.0d0
        d2fdce(7)  = -4.0d0
        d2fdce(8)  =  0.0d0
        d2fdce(9)  =  0.0d0
        d2fdce(10) = +0.0d0

        !d2φ/dξdζ
        d2fdcs(1)  = +4.0d0
        d2fdcs(2)  =  0.0d0
        d2fdcs(3)  = +0.0d0
        d2fdcs(4)  = +0.0d0
        d2fdcs(5)  = -4.0d0
        d2fdcs(6)  =  0.0d0
        d2fdcs(7)  =  0.0d0
        d2fdcs(8)  = -4.0d0
        d2fdcs(9)  = +4.0d0
        d2fdcs(10) = +0.0d0

        !dφ/dη
        dfde(1)  = -3.0d0+4.0d0*(c+e+s)
        dfde(2)  = +0.0d0
        dfde(3)  = +4.0d0*e-1.0d0
        dfde(4)  = +0.0d0
        dfde(5)  = -4.0d0*c
        dfde(6)  = +4.0d0*c
        dfde(7)  = +4.0d0*(1.0d0-c-2.0d0*e-s)
        dfde(8)  = -4.0d0*s
        dfde(9)  = +0.0d0
        dfde(10) = +4.0d0*s

        !d2φ/dη2
        d2fde2(1)  = +4.0d0
        d2fde2(2)  = +0.0d0
        d2fde2(3)  = +4.0d0
        d2fde2(4)  = +0.0d0
        d2fde2(5)  = +0.0d0
        d2fde2(6)  = +0.0d0
        d2fde2(7)  = -8.0d0
        d2fde2(8)  = +0.0d0
        d2fde2(9)  = +0.0d0
        d2fde2(10) = +0.0d0

        !d2φ/dηdζ
        d2fdes(1)  = +4.0d0
        d2fdes(2)  = +0.0d0
        d2fdes(3)  =  0.0d0
        d2fdes(4)  = +0.0d0
        d2fdes(5)  =  0.0d0
        d2fdes(6)  =  0.0d0
        d2fdes(7)  = -4.0d0
        d2fdes(8)  = -4.0d0
        d2fdes(9)  = +0.0d0
        d2fdes(10) = +4.0d0

        !dφ/dζ
        dfds(1)  = -3.0d0+4.0d0*(c+e+s)
        dfds(2)  = +0.0d0
        dfds(3)  = +0.0d0
        dfds(4)  = +4.0d0*s-1.0d0
        dfds(5)  = -4.0d0*c
        dfds(6)  = +0.0d0
        dfds(7)  = -4.0d0*e
        dfds(8)  = +4.0d0*(1.0d0-c-e-2.0d0*s)
        dfds(9)  = +4.0d0*c
        dfds(10) = +4.0d0*e

        !d2φ/dζ2
        d2fds2(1)  = +4.0d0
        d2fds2(2)  = +0.0d0
        d2fds2(3)  = +0.0d0
        d2fds2(4)  = +4.0d0
        d2fds2(5)  = +0.0d0
        d2fds2(6)  = +0.0d0
        d2fds2(7)  = +0.0d0
        d2fds2(8)  = -8.0d0
        d2fds2(9)  = +0.0d0
        d2fds2(10) = +0.0d0

    end subroutine setQuadraticTetrahedron

    !----------------------------------------------------------------------- 

    subroutine setLinearTetrahedron_Bfn(xi, bfn)

        implicit none

        PetscReal, dimension(:), intent(in)  :: xi
        PetscReal, dimension(:), intent(out) :: bfn

        PetscReal :: c, e, s

        c = xi(1) ; e = xi(2) ; s = xi(3)

        !Basis function, φ
        bfn(1)  = 1.0d0-c-e-s
        bfn(2)  = c
        bfn(3)  = e
        bfn(4)  = s

    end subroutine setLinearTetrahedron_Bfn

    !----------------------------------------------------------------------- 

    subroutine setQuadraticTetrahedron_Bfn(xi, bfn)

        implicit none

        PetscReal, dimension(:), intent(in)  :: xi
        PetscReal, dimension(:), intent(out) :: bfn

        PetscReal :: c, e, s

        c = xi(1) ; e = xi(2) ; s = xi(3)

        !Basis function, φ
        bfn(1)  = +(1.0d0-c-e-s)*(1.0d0-2.0d0*c-2.0d0*e-2.0d0*s)
        bfn(2)  = +c*(2.0d0*c-1.0d0)
        bfn(3)  = +e*(2.0d0*e-1.0d0)
        bfn(4)  = +s*(2.0d0*s-1.0d0)
        bfn(5)  = +4.0d0*(1.0d0-c-e-s)*c
        bfn(6)  = +4.0d0*c*e
        bfn(7)  = +4.0d0*(1.0d0-c-e-s)*e
        bfn(8)  = +4.0d0*(1.0d0-c-e-s)*s
        bfn(9)  = +4.0d0*c*s
        bfn(10) = +4.0d0*e*s

    end subroutine setQuadraticTetrahedron_Bfn
    
    !-----------------------------------------------------------------------

    subroutine setTetrahedron_dBfn(f_sub, xi, dfdc, dfde, dfds)

        implicit none

        interface
            subroutine f_sub(xi, f)
                PetscReal, dimension(:), intent(in)  :: xi
                PetscReal, dimension(:), intent(out) :: f
            end subroutine f_sub
        end interface
        PetscReal, dimension(:), intent(in)  :: xi
        PetscReal, dimension(:), intent(out) :: dfdc, dfde, dfds

        call derivative(f_sub, xi, 1, dfdc)
        call derivative(f_sub, xi, 2, dfde)
        call derivative(f_sub, xi, 3, dfds)

    end subroutine setTetrahedron_dBfn

    !-----------------------------------------------------------------------

    subroutine setTetrahedron_d2Bfn(f_sub, xi, d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes)

        implicit none

        interface
            subroutine f_sub(xi, f)
                PetscReal, dimension(:), intent(in)  :: xi
                PetscReal, dimension(:), intent(out) :: f
            end subroutine f_sub
        end interface
        PetscReal, dimension(:), intent(in)  :: xi
        PetscReal, dimension(:), intent(out) :: d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes

        call derivative2(f_sub, xi, 1, d2fdc2)
        call derivative2(f_sub, xi, 2, d2fde2)
        call derivative2(f_sub, xi, 3, d2fds2)

        call derivative2_mixed(f_sub, xi, 1, 2, d2fdce)
        call derivative2_mixed(f_sub, xi, 1, 3, d2fdcs)
        call derivative2_mixed(f_sub, xi, 2, 3, d2fdes)

    end subroutine setTetrahedron_d2Bfn

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

end module GaussTetrahedron_mod
