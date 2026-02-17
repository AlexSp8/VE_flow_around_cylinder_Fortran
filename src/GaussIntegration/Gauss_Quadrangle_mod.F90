
#include <petsc/finclude/petscksp.h>

module GaussQuadrangle_mod

    use GaussParameters_mod, only: GaussIntegration

    implicit none

    PetscInt, parameter :: NGauss_1D = 2, NGauss_2D = NGauss_1D**2, NGauss_LC = 6

    contains

    subroutine setGaussParametersQuadrangle(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt
        
        GaussInt%NGauss_1D = NGauss_1D
        if (GaussInt%order == 'Quadratic' .or. GaussInt%order == 'Serendipity') then
            GaussInt%NGauss_1D = 3
        end if
        GaussInt%NGauss_2D = GaussInt%NGauss_1D**2
        GaussInt%NGauss_3D = 0
        GaussInt%NGauss_LC = NGauss_LC

        GaussInt%NGauss_b = GaussInt%NGauss_2D

        call setGaussPointsQuadrangle(GaussInt)
        call setGaussWeightsQuadrangle(GaussInt)
        call setGaussBasisFunctionsQuadrangle(GaussInt)

    end subroutine setGaussParametersQuadrangle

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsQuadrangle(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        call setGaussPointsQuadrangle1D(GaussInt)

        call setGaussPointsQuadrangle2D(GaussInt)

        call setGaussPointsQuadrangleLC(GaussInt)

    end subroutine setGaussPointsQuadrangle

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsQuadrangle1D(GaussInt)

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
            GaussInt%Points_1D(1) = 0.0d0+1.0d-8
            GaussInt%Points_1D(2) = +sqrt(3.0d0/5.0d0)
            GaussInt%Points_1D(3) = -sqrt(3.0d0/5.0d0)
        case default
            write(*,'(a)') 'Wrong GaussInt%NGauss_1D in setGaussPointsQuadrangle1D!'
            stop
        end select

    end subroutine setGaussPointsQuadrangle1D

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsQuadrangle2D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        PetscInt :: i, j, Ng_total

        allocate( GaussInt%Points_2D(GaussInt%NGauss_2D,2) )

        Ng_total = 0
        do i = 1, GaussInt%NGauss_1D
            do j = 1, GaussInt%NGauss_1D
                Ng_total = Ng_total+1
                GaussInt%Points_2D(Ng_total,1) = GaussInt%Points_1D(i)
                GaussInt%Points_2D(Ng_total,2) = GaussInt%Points_1D(j)
            end do
        end do

    end subroutine setGaussPointsQuadrangle2D

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsQuadrangleLC(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        allocate( GaussInt%Points_LC(GaussInt%NGauss_LC) )

        select case (GaussInt%NGauss_LC)

        case (6)

            GaussInt%Points_LC(1) = -0.9324695142031520278123d0
            GaussInt%Points_LC(2) = -0.6612093864662645136610d0
            GaussInt%Points_LC(3) = -0.2386191860831969086305d0
            GaussInt%Points_LC(4) = +0.2386191860831969086305d0
            GaussInt%Points_LC(5) = +0.6612093864662645136610d0
            GaussInt%Points_LC(6) = +0.9324695142031520278123d0

        case default
            write(*,'(a)') 'Wrong NGauss_LC in setGaussPointsQuadrangleLC!'
            stop
        end select

    end subroutine setGaussPointsQuadrangleLC

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsQuadrangle(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        call setGaussWeightsQuadrangle1D(GaussInt)

        call setGaussWeightsQuadrangle2D(GaussInt)

        call setGaussWeightsQuadrangleLC(GaussInt)

        allocate( GaussInt%Weights_b(GaussInt%NGauss_2D) )
        GaussInt%Weights_b(:) = GaussInt%Weights_2D(:)

    end subroutine setGaussWeightsQuadrangle

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsQuadrangle1D(GaussInt)

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
            write(*,'(a)') 'Wrong NGauss_1D in setGaussWeightsQuadrangle1D!'
            stop
        end select

    end subroutine setGaussWeightsQuadrangle1D

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsQuadrangle2D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        PetscInt :: i, j, Ng_total

        allocate( GaussInt%Weights_2D(GaussInt%NGauss_2D) )

        Ng_total = 0
        do i = 1, GaussInt%NGauss_1D
            do j = 1, GaussInt%NGauss_1D
                Ng_total = Ng_total+1
                GaussInt%Weights_2D(Ng_total) = GaussInt%Weights_1D(i)&
                                                        *GaussInt%Weights_1D(j)
            end do
        end do

    end subroutine setGaussWeightsQuadrangle2D

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsQuadrangleLC(GaussInt)

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
            write(*,'(a)') 'Wrong NGauss_LC in setGaussWeightsQuadrangleLC!'
            stop
        end select

    end subroutine setGaussWeightsQuadrangleLC

    !----------------------------------------------------------------------- 

    subroutine setGaussBasisFunctionsQuadrangle(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        call setGaussBasisFunctionsQuadrangleEdge(GaussInt)
        call setGaussBasisFunctionsQuadrangle2D(GaussInt)

    end subroutine setGaussBasisFunctionsQuadrangle

    !----------------------------------------------------------------------- 

    subroutine setGaussBasisFunctionsQuadrangleEdge(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        PetscInt :: ig, ied
        PetscReal :: c, e
        PetscReal, dimension(2) :: xi
        PetscReal, dimension(GaussInt%Nbf) :: bfn, dfdc, dfde, d2fdc2, d2fde2, d2fdce    !C-> ξ, E-> η

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
                    e = -1.0d0
                case (2)
                    c = +1.0d0
                    e = GaussInt%Points_1D(ig)
                case (3)
                    c = GaussInt%Points_1D(ig)
                    e = +1.0d0
                case (4)
                    c = -1.0d0
                    e = GaussInt%Points_1D(ig)
                case default
                    write(*,'(A)') 'Wrong ied in setGaussBasisFunctionsQuadrangleEdge!'
                    stop
                end select

                xi = [c, e]

                select case (GaussInt%order)
                case ('Linear')
                    call setLinearQuadrangle(bfn, dfdc, dfde, c, e)
                    ! call setLinearQuadrangle_Bfn(xi, bfn)
                    ! call setQuadrangle_dBfn(setLinearQuadrangle_Bfn, xi, dfdc, dfde)
                    d2fdc2(:) = 0.0d0 ; d2fde2(:) = 0.0d0 ; d2fdce(:) = 0.0d0
                case ('Quadratic')
                    call setQuadraticQuadrangle(bfn, dfdc, dfde, d2fdc2, d2fde2, d2fdce, c, e)
                    ! call setQuadraticQuadrangle_Bfn(xi, bfn)
                    ! call setQuadrangle_dBfn(setQuadraticQuadrangle_Bfn, xi, dfdc, dfde)
                    ! call setQuadrangle_d2Bfn(setQuadraticQuadrangle_Bfn, xi, d2fdc2, d2fde2, d2fdce)
                case ('Serendipity')
                    call setSerendipityQuadrangle(bfn, dfdc, dfde, d2fdc2, d2fde2, d2fdce, c, e)
                    ! call setSerendipityQuadrangle_Bfn(xi, bfn)
                    ! call setQuadrangle_dBfn(setSerendipityQuadrangle_Bfn, xi, dfdc, dfde)
                    ! call setQuadrangle_d2Bfn(setSerendipityQuadrangle_Bfn, xi, d2fdc2, d2fde2, d2fdce)
                case default
                    write(*,'(a)') 'Wrong GaussInt%order in setGaussBasisFunctionsQuadrangleEdge!'
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

    end subroutine setGaussBasisFunctionsQuadrangleEdge

    !----------------------------------------------------------------------- 

    subroutine setGaussBasisFunctionsQuadrangle2D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        PetscInt :: ig, i
        PetscReal :: c, e
        PetscReal, dimension(2) :: xi
        PetscReal, dimension(GaussInt%Nbf) :: bfn, dfdc, dfde, d2fdc2, d2fde2, d2fdce    !C-> ξ, E-> η

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
                call setLinearQuadrangle(bfn, dfdc, dfde, c, e)
                ! call setLinearQuadrangle_Bfn(xi, bfn)
                ! call setQuadrangle_dBfn(setLinearQuadrangle_Bfn, xi, dfdc, dfde)
                d2fdc2(:) = 0.0d0 ; d2fde2(:) = 0.0d0 ; d2fdce(:) = 0.0d0
            case ('Quadratic')
                call setQuadraticQuadrangle(bfn, dfdc, dfde, d2fdc2, d2fde2, d2fdce, c, e)
                ! call setQuadraticQuadrangle_Bfn(xi, bfn)
                ! call setQuadrangle_dBfn(setQuadraticQuadrangle_Bfn, xi, dfdc, dfde)
                ! call setQuadrangle_d2Bfn(setQuadraticQuadrangle_Bfn, xi, d2fdc2, d2fde2, d2fdce)
            case ('Serendipity')
                call setSerendipityQuadrangle(bfn, dfdc, dfde, d2fdc2, d2fde2, d2fdce, c, e)
                ! call setSerendipityQuadrangle_Bfn(xi, bfn)
                ! call setQuadrangle_dBfn(setSerendipityQuadrangle_Bfn, xi, dfdc, dfde)
                ! call setQuadrangle_d2Bfn(setSerendipityQuadrangle_Bfn, xi, d2fdc2, d2fde2, d2fdce)
            case default
                write(*,'(a)') 'Wrong GaussInt%order in setGaussBasisFunctionsQuadrangle2D!'
                stop
            end select
            GaussInt%Bfn(:,ig)  = bfn(:)
            GaussInt%dfdc(:,ig) = dfdc(:)
            GaussInt%dfde(:,ig) = dfde(:)
            GaussInt%d2fdc2(:,ig) = d2fdc2(:)
            GaussInt%d2fde2(:,ig) = d2fde2(:)
            GaussInt%d2fdce(:,ig) = d2fdce(:)

        end do

    end subroutine setGaussBasisFunctionsQuadrangle2D

    !----------------------------------------------------------------------- 

    subroutine setLinearQuadrangle(bfn, dfdc, dfde, c, e)

        implicit none

        PetscReal, intent(in)  :: c, e
        PetscReal, dimension(:), intent(out) :: bfn, dfdc, dfde

        !Basis function, φ
        bfn(1)  = (1.0d0-c)*(1.0d0-e)/4.0d0
        bfn(2)  = (1.0d0+c)*(1.0d0-e)/4.0d0
        bfn(3)  = (1.0d0+c)*(1.0d0+e)/4.0d0
        bfn(4)  = (1.0d0-c)*(1.0d0+e)/4.0d0

        !dφ/dξ
        dfdc(1) = -(1.0d0-e)/4.0d0
        dfdc(2) = +(1.0d0-e)/4.0d0
        dfdc(3) = +(1.0d0+e)/4.0d0
        dfdc(4) = -(1.0d0+e)/4.0d0

        !dφ/dη
        dfde(1) = -(1.0d0-c)/4.0d0
        dfde(2) = -(1.0d0+c)/4.0d0
        dfde(3) = +(1.0d0+c)/4.0d0
        dfde(4) = +(1.0d0-c)/4.0d0

    end subroutine setLinearQuadrangle

    !----------------------------------------------------------------------- 

    subroutine setQuadraticQuadrangle(bfn, dfdc, dfde, d2fdc2, d2fde2, d2fdce, c, e)

        implicit none

        PetscReal, intent(in)  :: c, e
        PetscReal, dimension(:), intent(out) :: bfn, dfdc, dfde, d2fdc2, d2fde2, d2fdce

        !Basis function, φ
        bfn(1)  = +(1.0d0-c)*c*(1.0d0-e)*e/4.0d0
        bfn(2)  = -(1.0d0+c)*c*(1.0d0-e)*e/4.0d0
        bfn(3)  = +(1.0d0+c)*c*(1.0d0+e)*e/4.0d0
        bfn(4)  = -(1.0d0-c)*c*(1.0d0+e)*e/4.0d0
        bfn(5)  = -(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e/2.0d0
        bfn(6)  = +(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)/2.0d0
        bfn(7)  = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e/2.0d0
        bfn(8)  = -(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)/2.0d0
        bfn(9)  = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)

        !dφ/dξ
        dfdc(1) = +(1.0d0-2.0d0*c)*(1.0d0-e)*e/4.0d0
        dfdc(2) = -(1.0d0+2.0d0*c)*(1.0d0-e)*e/4.0d0
        dfdc(3) = +(1.0d0+2.0d0*c)*(1.0d0+e)*e/4.0d0
        dfdc(4) = -(1.0d0-2.0d0*c)*(1.0d0+e)*e/4.0d0
        dfdc(5) = -(-2.0d0*c)*(1.0d0-e)*e/2.0d0
        dfdc(6) = +(1.0d0+2.0d0*c)*(1.0d0+e)*(1.0d0-e)/2.0d0
        dfdc(7) = +(-2.0d0*c)*(1.0d0+e)*e/2.0d0
        dfdc(8) = -(1.0d0-2.0d0*c)*(1.0d0+e)*(1.0d0-e)/2.0d0
        dfdc(9) = +(-2.0d0*c)*(1.0d0+e)*(1.0d0-e)

        !d2φ/dξ2
        d2fdc2(1) = +(-2.0d0)*(1.0d0-e)*e/4.0d0
        d2fdc2(2) = -(+2.0d0)*(1.0d0-e)*e/4.0d0
        d2fdc2(3) = +(+2.0d0)*(1.0d0+e)*e/4.0d0
        d2fdc2(4) = -(-2.0d0)*(1.0d0+e)*e/4.0d0
        d2fdc2(5) = -(-2.0d0)*(1.0d0-e)*e/2.0d0
        d2fdc2(6) = +(+2.0d0)*(1.0d0+e)*(1.0d0-e)/2.0d0
        d2fdc2(7) = +(-2.0d0)*(1.0d0+e)*e/2.0d0
        d2fdc2(8) = -(-2.0d0)*(1.0d0+e)*(1.0d0-e)/2.0d0
        d2fdc2(9) = +(-2.0d0)*(1.0d0+e)*(1.0d0-e)

        !d2φ/dξdη
        d2fdce(1) = +(1.0d0-2.0d0*c)*(1.0d0-2.0d0*e)/4.0d0
        d2fdce(2) = -(1.0d0+2.0d0*c)*(1.0d0-2.0d0*e)/4.0d0
        d2fdce(3) = +(1.0d0+2.0d0*c)*(1.0d0+2.0d0*e)/4.0d0
        d2fdce(4) = -(1.0d0-2.0d0*c)*(1.0d0+2.0d0*e)/4.0d0
        d2fdce(5) = -(-2.0d0*c)*(1.0d0-2.0d0*e)/2.0d0
        d2fdce(6) = +(1.0d0+2.0d0*c)*(-2.0d0*e)/2.0d0
        d2fdce(7) = +(-2.0d0*c)*(1.0d0+2.0d0*e)/2.0d0
        d2fdce(8) = -(1.0d0-2.0d0*c)*(-2.0d0*e)/2.0d0
        d2fdce(9) = +(-2.0d0*c)*(-2.0d0*e)

        !dφ/dη
        dfde(1) = +(1.0d0-c)*c*(1.0d0-2.0d0*e)/4.0d0
        dfde(2) = -(1.0d0+c)*c*(1.0d0-2.0d0*e)/4.0d0
        dfde(3) = +(1.0d0+c)*c*(1.0d0+2.0d0*e)/4.0d0
        dfde(4) = -(1.0d0-c)*c*(1.0d0+2.0d0*e)/4.0d0
        dfde(5) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-2.0d0*e)/2.0d0
        dfde(6) = +(1.0d0+c)*c*(-2.0d0*e)/2.0d0
        dfde(7) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+2.0d0*e)/2.0d0
        dfde(8) = -(1.0d0-c)*c*(-2.0d0*e)/2.0d0
        dfde(9) = +(1.0d0+c)*(1.0d0-c)*(-2.0d0*e)

        !d2φ/dη2
        d2fde2(1) = +(1.0d0-c)*c*(-2.0d0)/4.0d0
        d2fde2(2) = -(1.0d0+c)*c*(-2.0d0)/4.0d0
        d2fde2(3) = +(1.0d0+c)*c*(+2.0d0)/4.0d0
        d2fde2(4) = -(1.0d0-c)*c*(+2.0d0)/4.0d0
        d2fde2(5) = -(1.0d0+c)*(1.0d0-c)*(-2.0d0)/2.0d0
        d2fde2(6) = +(1.0d0+c)*c*(-2.0d0)/2.0d0
        d2fde2(7) = +(1.0d0+c)*(1.0d0-c)*(+2.0d0)/2.0d0
        d2fde2(8) = -(1.0d0-c)*c*(-2.0d0)/2.0d0
        d2fde2(9) = +(1.0d0+c)*(1.0d0-c)*(-2.0d0)

    end subroutine setQuadraticQuadrangle

    !----------------------------------------------------------------------- 

    subroutine setSerendipityQuadrangle(bfn, dfdc, dfde, d2fdc2, d2fde2, d2fdce, c, e)

        implicit none

        PetscReal, intent(in)  :: c, e
        PetscReal, dimension(:), intent(out) :: bfn, dfdc, dfde, d2fdc2, d2fde2, d2fdce

        !Basis function, φ
        bfn(1) = -(1.0d0-c)*(1.0d0-e)*(1.0d0+c+e)/4.0d0
        bfn(2) = -(1.0d0+c)*(1.0d0-e)*(1.0d0-c+e)/4.0d0
        bfn(3) = -(1.0d0+c)*(1.0d0+e)*(1.0d0-c-e)/4.0d0
        bfn(4) = -(1.0d0-c)*(1.0d0+e)*(1.0d0+c-e)/4.0d0
        bfn(5) = +(1.0d0+c)*(1.0d0-c)*(1.0d0-e)/2.0d0
        bfn(6) = +(1.0d0+c)*(1.0d0+e)*(1.0d0-e)/2.0d0
        bfn(7) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)/2.0d0
        bfn(8) = +(1.0d0-c)*(1.0d0+e)*(1.0d0-e)/2.0d0

        !dφ/dξ
        dfdc(1) = -((-1.0d0)*(1.0d0+c+e)+(1.0d0-c)*(+1.0d0))*(1.0d0-e)/4.0d0
        dfdc(2) = -((+1.0d0)*(1.0d0-c+e)+(1.0d0+c)*(-1.0d0))*(1.0d0-e)/4.0d0
        dfdc(3) = -((+1.0d0)*(1.0d0-c-e)+(1.0d0+c)*(-1.0d0))*(1.0d0+e)/4.0d0
        dfdc(4) = -((-1.0d0)*(1.0d0+c-e)+(1.0d0-c)*(+1.0d0))*(1.0d0+e)/4.0d0
        dfdc(5) = +(-2.0d0*c)*(1.0d0-e)/2.0d0
        dfdc(6) = +(+1.0d0)*(1.0d0+e)*(1.0d0-e)/2.0d0
        dfdc(7) = +(-2.0d0*c)*(1.0d0+e)/2.0d0
        dfdc(8) = +(-1.0d0)*(1.0d0+e)*(1.0d0-e)/2.0d0

        !d2φ/dξ2
        d2fdc2(1) = -((-1.0d0)*(+1.0d0)+(-1.0d0)*(+1.0d0))*(1.0d0-e)/4.0d0
        d2fdc2(2) = -((+1.0d0)*(-1.0d0)+(+1.0d0)*(-1.0d0))*(1.0d0-e)/4.0d0
        d2fdc2(3) = -((+1.0d0)*(-1.0d0)+(+1.0d0)*(-1.0d0))*(1.0d0+e)/4.0d0
        d2fdc2(4) = -((-1.0d0)*(+1.0d0)+(-1.0d0)*(+1.0d0))*(1.0d0+e)/4.0d0
        d2fdc2(5) = +(-2.0d0)*(1.0d0-e)/2.0d0
        d2fdc2(6) = +0.0d0
        d2fdc2(7) = +(-2.0d0)*(1.0d0+e)/2.0d0
        d2fdc2(8) = +0.0d0

        !d2φ/dξdη
        d2fdce(1) = -( -(+1.0d0)*(1.0d0-e)+(-(1.0d0+c+e)+(1.0d0-c))*(-1.0d0) )/4.0d0
        d2fdce(2) = -( +(+1.0d0)*(1.0d0-e)+(+(1.0d0-c+e)-(1.0d0+c))*(-1.0d0) )/4.0d0
        d2fdce(3) = -( +(-1.0d0)*(1.0d0+e)+(+(1.0d0-c-e)-(1.0d0+c))*(+1.0d0) )/4.0d0
        d2fdce(4) = -( -(-1.0d0)*(1.0d0+e)+(-(1.0d0+c-e)+(1.0d0-c))*(+1.0d0) )/4.0d0
        d2fdce(5) = +(-2.0d0*c)*(-1.0d0)/2.0d0
        d2fdce(6) = +(+1.0d0)*(-2.0d0*e)/2.0d0
        d2fdce(7) = +(-2.0d0*c)*(+1.0d0)/2.0d0
        d2fdce(8) = +(-1.0d0)*(-2.0d0*e)/2.0d0

        !dφ/dη
        dfde(1) = -(1.0d0-c)*((-1.0d0)*(1.0d0+c+e)+(1.0d0-e)*(+1.0d0))/4.0d0
        dfde(2) = -(1.0d0+c)*((-1.0d0)*(1.0d0-c+e)+(1.0d0-e)*(+1.0d0))/4.0d0
        dfde(3) = -(1.0d0+c)*((+1.0d0)*(1.0d0-c-e)+(1.0d0+e)*(-1.0d0))/4.0d0
        dfde(4) = -(1.0d0-c)*((+1.0d0)*(1.0d0+c-e)+(1.0d0+e)*(-1.0d0))/4.0d0
        dfde(5) = +(1.0d0+c)*(1.0d0-c)*(-1.0d0)/2.0d0
        dfde(6) = +(1.0d0+c)*(-2.0d0*e)/2.0d0
        dfde(7) = +(1.0d0+c)*(1.0d0-c)*(+1.0d0)/2.0d0
        dfde(8) = +(1.0d0-c)*(-2.0d0*e)/2.0d0

        !d2φ/dη2
        d2fde2(1) = -(1.0d0-c)*(-(+1.0d0)+(-1.0d0))/4.0d0
        d2fde2(2) = -(1.0d0+c)*(-(+1.0d0)+(-1.0d0))/4.0d0
        d2fde2(3) = -(1.0d0+c)*(+(-1.0d0)-(+1.0d0))/4.0d0
        d2fde2(4) = -(1.0d0-c)*(+(-1.0d0)-(+1.0d0))/4.0d0
        d2fde2(5) = +0.0d0
        d2fde2(6) = +(1.0d0+c)*(-2.0d0)/2.0d0
        d2fde2(7) = +0.0d0
        d2fde2(8) = +(1.0d0-c)*(-2.0d0)/2.0d0

    end subroutine setSerendipityQuadrangle
    
    !----------------------------------------------------------------------- 

    subroutine setLinearQuadrangle_Bfn(xi, bfn)

        implicit none

        PetscReal, dimension(:), intent(in)  :: xi
        PetscReal, dimension(:), intent(out) :: bfn

        PetscReal :: c, e

        c = xi(1) ; e = xi(2)

        !Basis function, φ
        bfn(1)  = (1.0d0-c)*(1.0d0-e)/4.0d0
        bfn(2)  = (1.0d0+c)*(1.0d0-e)/4.0d0
        bfn(3)  = (1.0d0+c)*(1.0d0+e)/4.0d0
        bfn(4)  = (1.0d0-c)*(1.0d0+e)/4.0d0

    end subroutine setLinearQuadrangle_Bfn

    !----------------------------------------------------------------------- 

    subroutine setQuadraticQuadrangle_Bfn(xi, bfn)

        implicit none

        PetscReal, dimension(:), intent(in)  :: xi
        PetscReal, dimension(:), intent(out) :: bfn

        PetscReal :: c, e

        c = xi(1) ; e = xi(2)

        !Basis function, φ
        bfn(1)  = +(1.0d0-c)*c*(1.0d0-e)*e/4.0d0
        bfn(2)  = -(1.0d0+c)*c*(1.0d0-e)*e/4.0d0
        bfn(3)  = +(1.0d0+c)*c*(1.0d0+e)*e/4.0d0
        bfn(4)  = -(1.0d0-c)*c*(1.0d0+e)*e/4.0d0
        bfn(5)  = -(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e/2.0d0
        bfn(6)  = +(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)/2.0d0
        bfn(7)  = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e/2.0d0
        bfn(8)  = -(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)/2.0d0
        bfn(9)  = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)

    end subroutine setQuadraticQuadrangle_Bfn

    !----------------------------------------------------------------------- 

    subroutine setSerendipityQuadrangle_Bfn(xi, bfn)

        implicit none

        PetscReal, dimension(:), intent(in)  :: xi
        PetscReal, dimension(:), intent(out) :: bfn

        PetscReal :: c, e

        c = xi(1) ; e = xi(2)

        !Basis function, φ
        bfn(1) = -(1.0d0-c)*(1.0d0-e)*(1.0d0+c+e)/4.0d0
        bfn(2) = -(1.0d0+c)*(1.0d0-e)*(1.0d0-c+e)/4.0d0
        bfn(3) = -(1.0d0+c)*(1.0d0+e)*(1.0d0-c-e)/4.0d0
        bfn(4) = -(1.0d0-c)*(1.0d0+e)*(1.0d0+c-e)/4.0d0
        bfn(5) = +(1.0d0+c)*(1.0d0-c)*(1.0d0-e)/2.0d0
        bfn(6) = +(1.0d0+c)*(1.0d0+e)*(1.0d0-e)/2.0d0
        bfn(7) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)/2.0d0
        bfn(8) = +(1.0d0-c)*(1.0d0+e)*(1.0d0-e)/2.0d0

    end subroutine setSerendipityQuadrangle_Bfn

    !-----------------------------------------------------------------------

    subroutine setQuadrangle_dBfn(f_sub, xi, dfdc, dfde)

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

    end subroutine setQuadrangle_dBfn

    !-----------------------------------------------------------------------

    subroutine setQuadrangle_d2Bfn(f_sub, xi, d2fdc2, d2fde2, d2fdce)

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

    end subroutine setQuadrangle_d2Bfn

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

end module GaussQuadrangle_mod
