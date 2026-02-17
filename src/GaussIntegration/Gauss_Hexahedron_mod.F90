
#include <petsc/finclude/petscksp.h>

module GaussHexahedron_mod

    use GaussParameters_mod, only: GaussIntegration

    implicit none

    PetscInt, parameter :: NGauss_1D = 2, NGauss_2D = NGauss_1D**2, NGauss_3D = NGauss_1D**3, &
                            NGauss_LC = 6

    contains

    subroutine setGaussParametersHexahedron(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        GaussInt%NGauss_1D = NGauss_1D
        if (GaussInt%order == 'Quadratic' .or. GaussInt%order == 'Serendipity') then
            GaussInt%NGauss_1D = 3
        end if
        GaussInt%NGauss_2D = GaussInt%NGauss_1D**2
        GaussInt%NGauss_3D = GaussInt%NGauss_1D**3
        GaussInt%NGauss_LC = NGauss_LC

        GaussInt%NGauss_b = GaussInt%NGauss_3D
        
        call setGaussPointsHexahedron(GaussInt)
        call setGaussWeightsHexahedron(GaussInt)
        call setGaussBasisFunctionsHexahedron(GaussInt)

    end subroutine setGaussParametersHexahedron

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsHexahedron(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        call setGaussPointsHexahedron1D(GaussInt)

        call setGaussPointsHexahedron2D(GaussInt)

        call setGaussPointsHexahedron3D(GaussInt)

        call setGaussPointsHexahedronLC(GaussInt)

    end subroutine setGaussPointsHexahedron

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsHexahedron1D(GaussInt)

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
        case (4)
            GaussInt%Points_1D(1) = +sqrt(3.0d0/7.0d0 - (2.0d0/7.0d0)*sqrt(6.0d0/5.0d0))
            GaussInt%Points_1D(2) = -sqrt(3.0d0/7.0d0 - (2.0d0/7.0d0)*sqrt(6.0d0/5.0d0))
            GaussInt%Points_1D(3) = +sqrt(3.0d0/7.0d0 + (2.0d0/7.0d0)*sqrt(6.0d0/5.0d0))
            GaussInt%Points_1D(4) = -sqrt(3.0d0/7.0d0 + (2.0d0/7.0d0)*sqrt(6.0d0/5.0d0))
        case default
            write(*,'(a)') 'Wrong GaussInt%NGauss_1D in setGaussPointsHexahedron1D!'
            stop
        end select

    end subroutine setGaussPointsHexahedron1D

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsHexahedron2D(GaussInt)

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

    end subroutine setGaussPointsHexahedron2D

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsHexahedron3D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        PetscInt :: i, j, k, Ng_total

        allocate( GaussInt%Points_3D(GaussInt%NGauss_3D,3) )

        Ng_total = 0
        do i = 1, GaussInt%NGauss_1D
            do j = 1, GaussInt%NGauss_1D
                do k = 1, GaussInt%NGauss_1D
                    Ng_total = Ng_total+1
                    GaussInt%Points_3D(Ng_total,1) = GaussInt%Points_1D(i)
                    GaussInt%Points_3D(Ng_total,2) = GaussInt%Points_1D(j)
                    GaussInt%Points_3D(Ng_total,3) = GaussInt%Points_1D(k)
                end do
            end do
        end do
    
    end subroutine setGaussPointsHexahedron3D

    !----------------------------------------------------------------------- 

    subroutine setGaussPointsHexahedronLC(GaussInt)

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
            write(*,'(a)') 'Wrong GaussInt%NGauss_LC in setGaussPointsHexahedronLC!'
            stop
        end select

    end subroutine setGaussPointsHexahedronLC

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsHexahedron(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        call setGaussWeightsHexahedron1D(GaussInt)

        call setGaussWeightsHexahedron2D(GaussInt)

        call setGaussWeightsHexahedron3D(GaussInt)

        call setGaussWeightsHexahedronLC(GaussInt)

        allocate( GaussInt%Weights_b(GaussInt%NGauss_3D) )
        GaussInt%Weights_b(:) = GaussInt%Weights_3D(:)

    end subroutine setGaussWeightsHexahedron

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsHexahedron1D(GaussInt)

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
        case default
            write(*,'(a)') 'Wrong GaussInt%NGauss_1D in setGaussWeightsHexahedron1D!'
            stop
        end select

    end subroutine setGaussWeightsHexahedron1D

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsHexahedron2D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        PetscInt :: i, j, Ng_total

        allocate( GaussInt%Weights_2D(GaussInt%NGauss_2D) )

        Ng_total = 0
        do i = 1, GaussInt%NGauss_1D
            do j = 1, GaussInt%NGauss_1D
                Ng_total = Ng_total+1
                GaussInt%Weights_2D(Ng_total) = GaussInt%Weights_1D(i)*&
                                                        GaussInt%Weights_1D(j)
            end do
        end do
        
    end subroutine setGaussWeightsHexahedron2D

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsHexahedron3D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        PetscInt :: i, j, k, Ng_total

        allocate( GaussInt%Weights_3D(GaussInt%NGauss_3D) )

        Ng_total = 0
        do i = 1, GaussInt%NGauss_1D
            do j = 1, GaussInt%NGauss_1D
                do k = 1, GaussInt%NGauss_1D
                    Ng_total = Ng_total+1
                    GaussInt%Weights_3D(Ng_total) = GaussInt%Weights_1D(i)* &
                        GaussInt%Weights_1D(j)*GaussInt%Weights_1D(k)
                end do
            end do
        end do

    end subroutine setGaussWeightsHexahedron3D

    !----------------------------------------------------------------------- 

    subroutine setGaussWeightsHexahedronLC(GaussInt)

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
            write(*,'(a)') 'Wrong GaussInt%NGauss_LC in setGaussWeightsHexahedronLC!'
            stop
        end select

    end subroutine setGaussWeightsHexahedronLC

    !----------------------------------------------------------------------- 

    subroutine setGaussBasisFunctionsHexahedron(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        call setGaussBasisFunctionsHexahedronEdge(GaussInt)
        call setGaussBasisFunctionsHexahedronFace(GaussInt)
        call setGaussBasisFunctionsHexahedron3D(GaussInt)

    end subroutine setGaussBasisFunctionsHexahedron

    !----------------------------------------------------------------------- 

    subroutine setGaussBasisFunctionsHexahedronEdge(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        PetscInt :: ig, ied, Nbf, Ng, Ned
        PetscReal :: c, e, s
        PetscReal, dimension(3) :: xi
        PetscReal, dimension(GaussInt%Nbf) :: bfn, dfdc, dfde, dfds, &
            d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes    !C-> ξ, E-> η, s -> ζ

        Nbf = GaussInt%Nbf
        Ng = GaussInt%NGauss_1D
        Ned = GaussInt%Ned

        allocate( GaussInt%Bfn_E(Nbf, Ng, Ned) )
        allocate( GaussInt%dfdc_E(Nbf, Ng, Ned) )
        allocate( GaussInt%dfde_E(Nbf, Ng, Ned) )
        allocate( GaussInt%dfds_E(Nbf, Ng, Ned) )
        allocate( GaussInt%d2fdc2_E(Nbf, Ng, Ned) )
        allocate( GaussInt%d2fde2_E(Nbf, Ng, Ned) )
        allocate( GaussInt%d2fds2_E(Nbf, Ng, Ned) )
        allocate( GaussInt%d2fdce_E(Nbf, Ng, Ned) )
        allocate( GaussInt%d2fdcs_E(Nbf, Ng, Ned) )
        allocate( GaussInt%d2fdes_E(Nbf, Ng, Ned) )

        do ig = 1, Ng
            do ied = 1, Ned

                select case (ied)

                case (1)
                    c = GaussInt%Points_1D(ig)
                    e = -1.0d0
                    s = -1.0d0
                case (2)
                    c = +1.0d0
                    e = GaussInt%Points_1D(ig)
                    s = -1.0d0
                case (3)
                    c = GaussInt%Points_1D(ig)
                    e = +1.0d0
                    s = -1.0d0
                case (4)
                    c = -1.0d0
                    e = GaussInt%Points_1D(ig)
                    s = -1.0d0

                case (5)
                    c = GaussInt%Points_1D(ig)
                    e = -1.0d0
                    s = +1.0d0
                case (6)
                    c = +1.0d0
                    e = GaussInt%Points_1D(ig)
                    s = +1.0d0
                case (7)
                    c = GaussInt%Points_1D(ig)
                    e = +1.0d0
                    s = +1.0d0
                case (8)
                    c = -1.0d0
                    e = GaussInt%Points_1D(ig)
                    s = +1.0d0

                case (9)
                    c = -1.0d0
                    e = -1.0d0
                    s = GaussInt%Points_1D(ig)
                case (10)
                    c = +1.0d0
                    e = -1.0d0
                    s = GaussInt%Points_1D(ig)
                case (11)
                    c = +1.0d0
                    e = +1.0d0
                    s = GaussInt%Points_1D(ig)
                case (12)
                    c = -1.0d0
                    e = +1.0d0
                    s = GaussInt%Points_1D(ig)

                case default
                    write(*,'(a)') 'Wrong ied in setGaussBasisFunctionsHexahedronEdge!'
                    stop
                end select

                xi(:) = [c, e, s]

                select case (GaussInt%order)
                case ('Linear')
                    call setLinearHexahedron(bfn, dfdc, dfde, dfds, c, e, s)
                    ! call setLinearHexahedron_Bfn(xi, bfn)
                    ! call setHexahedron_dBfn(setLinearHexahedron_Bfn, xi, dfdc, dfde, dfds)
                    d2fdc2(:) = 0.0d0 ; d2fde2(:) = 0.0d0 ; d2fds2(:) = 0.0d0
                    d2fdce(:) = 0.0d0 ; d2fdcs(:) = 0.0d0 ; d2fdes(:) = 0.0d0
                case ('Quadratic')
                    call setQuadraticHexahedron(bfn, dfdc, dfde, dfds, &
                                d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes, c, e, s)
                    ! call setQuadraticHexahedron_Bfn(xi, bfn)
                    ! call setHexahedron_dBfn(setQuadraticHexahedron_Bfn, xi, dfdc, dfde, dfds)
                    ! call setHexahedron_d2Bfn(setQuadraticHexahedron_Bfn, xi, d2fdc2, d2fde2, d2fds2, &
                    !                                 d2fdce, d2fdcs, d2fdes)
                case ('Serendipity')
                    call setSerendipityHexahedron(bfn, dfdc, dfde, dfds, &
                                        d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes, c, e, s)
                    ! call setSerendipityHexahedron_Bfn(xi, bfn)
                    ! call setHexahedron_dBfn(setSerendipityHexahedron_Bfn, xi, dfdc, dfde, dfds)
                    ! call setHexahedron_d2Bfn(setSerendipityHexahedron_Bfn, xi, d2fdc2, d2fde2, d2fds2, &
                    !                                 d2fdce, d2fdcs, d2fdes)
                case default
                    write(*,'(a)') 'Wrong GaussInt%order in setGaussBasisFunctionsHexahedronEdge!'
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

    end subroutine setGaussBasisFunctionsHexahedronEdge

    !----------------------------------------------------------------------- 

    subroutine setGaussBasisFunctionsHexahedronFace(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        PetscInt :: ig, ifc
        PetscReal :: c, e, s
        PetscReal, dimension(3) :: xi
        PetscReal, dimension(GaussInt%Nbf) :: bfn, dfdc, dfde, dfds, &
            d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes    !C-> ξ, E-> η, s -> ζ

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
                    c = GaussInt%Points_2D(ig,1)
                    e = GaussInt%Points_2D(ig,2)
                    s = -1.0d0
                case (2)
                    c = GaussInt%Points_2D(ig,1)
                    e = GaussInt%Points_2D(ig,2)
                    s = +1.0d0
                case (3)
                    c = GaussInt%Points_2D(ig,1)
                    e = -1.0d0
                    s = GaussInt%Points_2D(ig,2)
                case (4)
                    c = GaussInt%Points_2D(ig,1)
                    e = +1.0d0
                    s = GaussInt%Points_2D(ig,2)
                case (5)
                    c = -1.0d0
                    e = GaussInt%Points_2D(ig,1)
                    s = GaussInt%Points_2D(ig,2)
                case (6)
                    c = +1.0d0
                    e = GaussInt%Points_2D(ig,1)
                    s = GaussInt%Points_2D(ig,2)
                case default
                    write(*,'(a)') 'Wrong ifc in setGaussBasisFunctionsHexahedronFace!'
                    stop
                end select

                xi(:) = [c, e, s]

                select case (GaussInt%order)
                case ('Linear')
                    call setLinearHexahedron(bfn, dfdc, dfde, dfds, c, e, s)
                    ! call setLinearHexahedron_Bfn(xi, bfn)
                    ! call setHexahedron_dBfn(setLinearHexahedron_Bfn, xi, dfdc, dfde, dfds)
                    d2fdc2(:) = 0.0d0 ; d2fde2(:) = 0.0d0 ; d2fds2(:) = 0.0d0
                    d2fdce(:) = 0.0d0 ; d2fdcs(:) = 0.0d0 ; d2fdes(:) = 0.0d0
                case ('Quadratic')
                    call setQuadraticHexahedron(bfn, dfdc, dfde, dfds, &
                                        d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes, c, e, s)
                    ! call setQuadraticHexahedron_Bfn(xi, bfn)
                    ! call setHexahedron_dBfn(setQuadraticHexahedron_Bfn, xi, dfdc, dfde, dfds)
                    ! call setHexahedron_d2Bfn(setQuadraticHexahedron_Bfn, xi, d2fdc2, d2fde2, d2fds2, &
                    !                                 d2fdce, d2fdcs, d2fdes)
                case ('Serendipity')
                    call setSerendipityHexahedron(bfn, dfdc, dfde, dfds, &
                        d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes, c, e, s)
                    ! call setSerendipityHexahedron_Bfn(xi, bfn)
                    ! call setHexahedron_dBfn(setSerendipityHexahedron_Bfn, xi, dfdc, dfde, dfds)
                    ! call setHexahedron_d2Bfn(setSerendipityHexahedron_Bfn, xi, d2fdc2, d2fde2, d2fds2, &
                    !                                 d2fdce, d2fdcs, d2fdes)
                case default
                    write(*,'(a)') 'Wrong GaussInt%order in setGaussBasisFunctionsHexahedronFace!'
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

    end subroutine setGaussBasisFunctionsHexahedronFace

    !----------------------------------------------------------------------- 

    subroutine setGaussBasisFunctionsHexahedron3D(GaussInt)

        implicit none

        type(GaussIntegration), intent(inout) :: GaussInt

        PetscInt :: ig, i
        PetscReal :: c, e, s           !C-> ξ, E-> η, s -> ζ
        PetscReal, dimension(3) :: xi
        PetscReal, dimension(GaussInt%Nbf) :: bfn, dfdc, dfde, dfds, &
            d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes

        allocate(    GaussInt%Bfn(GaussInt%Nbf, GaussInt%NGauss_3D) )
        allocate(   GaussInt%dfdc(GaussInt%Nbf, GaussInt%NGauss_3D) )
        allocate(   GaussInt%dfde(GaussInt%Nbf, GaussInt%NGauss_3D) )
        allocate(   GaussInt%dfds(GaussInt%Nbf, GaussInt%NGauss_3D) )
        allocate( GaussInt%d2fdc2(GaussInt%Nbf, GaussInt%NGauss_3D) )
        allocate( GaussInt%d2fde2(GaussInt%Nbf, GaussInt%NGauss_3D) )
        allocate( GaussInt%d2fds2(GaussInt%Nbf, GaussInt%NGauss_3D) )
        allocate( GaussInt%d2fdce(GaussInt%Nbf, GaussInt%NGauss_3D) )
        allocate( GaussInt%d2fdcs(GaussInt%Nbf, GaussInt%NGauss_3D) )
        allocate( GaussInt%d2fdes(GaussInt%Nbf, GaussInt%NGauss_3D) )

        do ig = 1, GaussInt%NGauss_3D

            c = GaussInt%Points_3D(ig,1)
            e = GaussInt%Points_3D(ig,2)
            s = GaussInt%Points_3D(ig,3)

            xi(:) = [c, e, s]

            select case (GaussInt%order)
            case ('Linear')
                call setLinearHexahedron(bfn, dfdc, dfde, dfds, c, e, s)
                ! call setLinearHexahedron_Bfn(xi, bfn)
                ! call setHexahedron_dBfn(setLinearHexahedron_Bfn, xi, dfdc, dfde, dfds)
                d2fdc2(:) = 0.0d0 ; d2fde2(:) = 0.0d0 ; d2fds2(:) = 0.0d0
                d2fdce(:) = 0.0d0 ; d2fdcs(:) = 0.0d0 ; d2fdes(:) = 0.0d0
            case ('Quadratic')
                call setQuadraticHexahedron(bfn, dfdc, dfde, dfds, &
                    d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes, c, e, s)
                ! call setQuadraticHexahedron_Bfn(xi, bfn)
                ! call setHexahedron_dBfn(setQuadraticHexahedron_Bfn, xi, dfdc, dfde, dfds)
                ! call setHexahedron_d2Bfn(setQuadraticHexahedron_Bfn, xi, d2fdc2, d2fde2, d2fds2, &
                !                                 d2fdce, d2fdcs, d2fdes)
            case ('Serendipity')
                call setSerendipityHexahedron(bfn, dfdc, dfde, dfds, &
                    d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes, c, e, s)
                ! call setSerendipityHexahedron_Bfn(xi, bfn)
                ! call setHexahedron_dBfn(setSerendipityHexahedron_Bfn, xi, dfdc, dfde, dfds)
                ! call setHexahedron_d2Bfn(setSerendipityHexahedron_Bfn, xi, d2fdc2, d2fde2, &
                !                         d2fds2, d2fdce, d2fdcs, d2fdes)
            case default
                write(*,'(a)') 'Wrong GaussInt%order in setGaussBasisFunctionsHexahedron3D!'
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

    end subroutine setGaussBasisFunctionsHexahedron3D

    !----------------------------------------------------------------------- 

    subroutine setLinearHexahedron(bfn, dfdc, dfde, dfds, c, e, s)

        implicit none

        PetscReal,               intent(in)  :: c, e, s
        PetscReal, dimension(:), intent(out) :: bfn, dfdc, dfde, dfds

        !Basis function, φ
        bfn(1)  = (1.0d0-c)*(1.0d0-e)*(1.0d0-s)/8.0d0
        bfn(2)  = (1.0d0+c)*(1.0d0-e)*(1.0d0-s)/8.0d0
        bfn(3)  = (1.0d0+c)*(1.0d0+e)*(1.0d0-s)/8.0d0
        bfn(4)  = (1.0d0-c)*(1.0d0+e)*(1.0d0-s)/8.0d0
        bfn(5)  = (1.0d0-c)*(1.0d0-e)*(1.0d0+s)/8.0d0
        bfn(6)  = (1.0d0+c)*(1.0d0-e)*(1.0d0+s)/8.0d0
        bfn(7)  = (1.0d0+c)*(1.0d0+e)*(1.0d0+s)/8.0d0
        bfn(8)  = (1.0d0-c)*(1.0d0+e)*(1.0d0+s)/8.0d0

        !dφ/dξ
        dfdc(1) = -(1.0d0-e)*(1.0d0-s)/8.0d0
        dfdc(2) = +(1.0d0-e)*(1.0d0-s)/8.0d0
        dfdc(3) = +(1.0d0+e)*(1.0d0-s)/8.0d0
        dfdc(4) = -(1.0d0+e)*(1.0d0-s)/8.0d0
        dfdc(5) = -(1.0d0-e)*(1.0d0+s)/8.0d0
        dfdc(6) = +(1.0d0-e)*(1.0d0+s)/8.0d0
        dfdc(7) = +(1.0d0+e)*(1.0d0+s)/8.0d0
        dfdc(8) = -(1.0d0+e)*(1.0d0+s)/8.0d0

        !dφ/dη
        dfde(1) = -(1.0d0-c)*(1.0d0-s)/8.0d0
        dfde(2) = -(1.0d0+c)*(1.0d0-s)/8.0d0
        dfde(3) = +(1.0d0+c)*(1.0d0-s)/8.0d0
        dfde(4) = +(1.0d0-c)*(1.0d0-s)/8.0d0
        dfde(5) = -(1.0d0-c)*(1.0d0+s)/8.0d0
        dfde(6) = -(1.0d0+c)*(1.0d0+s)/8.0d0
        dfde(7) = +(1.0d0+c)*(1.0d0+s)/8.0d0
        dfde(8) = +(1.0d0-c)*(1.0d0+s)/8.0d0

        !dφ/dζ
        dfds(1) = -(1.0d0-c)*(1.0d0-e)/8.0d0
        dfds(2) = -(1.0d0+c)*(1.0d0-e)/8.0d0
        dfds(3) = -(1.0d0+c)*(1.0d0+e)/8.0d0
        dfds(4) = -(1.0d0-c)*(1.0d0+e)/8.0d0
        dfds(5) = +(1.0d0-c)*(1.0d0-e)/8.0d0
        dfds(6) = +(1.0d0+c)*(1.0d0-e)/8.0d0
        dfds(7) = +(1.0d0+c)*(1.0d0+e)/8.0d0
        dfds(8) = +(1.0d0-c)*(1.0d0+e)/8.0d0

    end subroutine setLinearHexahedron

    !----------------------------------------------------------------------- 

    subroutine setQuadraticHexahedron(bfn, dfdc, dfde, dfds, &
        d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes, c, e, s)

        implicit none

        PetscReal, intent(in) :: c, e, s
        PetscReal, dimension(:), intent(out) :: bfn, dfdc, dfde, dfds, &
            d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes

        !Basis function, φ
        !Corners
        bfn(1)  = -(1.0d0-c)*c*(1.0d0-e)*e*(1.0d0-s)*s/8.0d0
        bfn(2)  = +(1.0d0+c)*c*(1.0d0-e)*e*(1.0d0-s)*s/8.0d0
        bfn(3)  = -(1.0d0+c)*c*(1.0d0+e)*e*(1.0d0-s)*s/8.0d0
        bfn(4)  = +(1.0d0-c)*c*(1.0d0+e)*e*(1.0d0-s)*s/8.0d0
        bfn(5)  = +(1.0d0-c)*c*(1.0d0-e)*e*(1.0d0+s)*s/8.0d0
        bfn(6)  = -(1.0d0+c)*c*(1.0d0-e)*e*(1.0d0+s)*s/8.0d0
        bfn(7)  = +(1.0d0+c)*c*(1.0d0+e)*e*(1.0d0+s)*s/8.0d0
        bfn(8)  = -(1.0d0-c)*c*(1.0d0+e)*e*(1.0d0+s)*s/8.0d0
        !Mid-edge ξ = 0
        bfn(9)  = +(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e*(1.0d0-s)*s/4.0d0
        bfn(11) = -(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e*(1.0d0-s)*s/4.0d0
        bfn(13) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e*(1.0d0+s)*s/4.0d0
        bfn(15) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e*(1.0d0+s)*s/4.0d0
        !Mid-edge η = 0
        bfn(10) = -(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)*s/4.0d0
        bfn(12) = +(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)*s/4.0d0
        bfn(14) = +(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*s/4.0d0
        bfn(16) = -(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*s/4.0d0
        !Mid-edge ζ = 0
        bfn(17) = +(1.0d0-c)*c*(1.0d0-e)*e*(1.0d0+s)*(1.0d0-s)/4.0d0
        bfn(18) = -(1.0d0+c)*c*(1.0d0-e)*e*(1.0d0+s)*(1.0d0-s)/4.0d0
        bfn(19) = +(1.0d0+c)*c*(1.0d0+e)*e*(1.0d0+s)*(1.0d0-s)/4.0d0
        bfn(20) = -(1.0d0-c)*c*(1.0d0+e)*e*(1.0d0+s)*(1.0d0-s)/4.0d0
        !Mid-face ξ = +-1
        ! bfn(22) = +(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/2.0d0
        bfn(23) = +(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, ξ = +1
        ! bfn(24) = -(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/2.0d0
        bfn(25) = -(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, ξ = -1
        !Mid-face η = +-1
        ! bfn(21) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e*(1.0d0+s)*(1.0d0-s)/2.0d0
        bfn(22) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, η = -1
        ! bfn(23) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e*(1.0d0+s)*(1.0d0-s)/2.0d0
        bfn(24) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, η = +1
        !Mid-face ζ = +-1
        ! bfn(25) = -(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)*s/2.0d0
        bfn(21) = -(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)*s/2.0d0  !Salome, ζ = -1
        bfn(26) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*s/2.0d0
        !Central
        bfn(27) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)

        !dφ/dξ
        !Corners
        dfdc(1)  = -(1.0d0-2.0d0*c)*(1.0d0-e)*e*(1.0d0-s)*s/8.0d0
        dfdc(2)  = +(1.0d0+2.0d0*c)*(1.0d0-e)*e*(1.0d0-s)*s/8.0d0
        dfdc(3)  = -(1.0d0+2.0d0*c)*(1.0d0+e)*e*(1.0d0-s)*s/8.0d0
        dfdc(4)  = +(1.0d0-2.0d0*c)*(1.0d0+e)*e*(1.0d0-s)*s/8.0d0
        dfdc(5)  = +(1.0d0-2.0d0*c)*(1.0d0-e)*e*(1.0d0+s)*s/8.0d0
        dfdc(6)  = -(1.0d0+2.0d0*c)*(1.0d0-e)*e*(1.0d0+s)*s/8.0d0
        dfdc(7)  = +(1.0d0+2.0d0*c)*(1.0d0+e)*e*(1.0d0+s)*s/8.0d0
        dfdc(8)  = -(1.0d0-2.0d0*c)*(1.0d0+e)*e*(1.0d0+s)*s/8.0d0
        !Mid-edge ξ = 0
        dfdc(9)  = +(-2.0d0*c)*(1.0d0-e)*e*(1.0d0-s)*s/4.0d0
        dfdc(11) = -(-2.0d0*c)*(1.0d0+e)*e*(1.0d0-s)*s/4.0d0
        dfdc(13) = -(-2.0d0*c)*(1.0d0-e)*e*(1.0d0+s)*s/4.0d0
        dfdc(15) = +(-2.0d0*c)*(1.0d0+e)*e*(1.0d0+s)*s/4.0d0
        !Mid-edge η = 0
        dfdc(10) = -(1.0d0+2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)*s/4.0d0
        dfdc(12) = +(1.0d0-2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)*s/4.0d0
        dfdc(14) = +(1.0d0+2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*s/4.0d0
        dfdc(16) = -(1.0d0-2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*s/4.0d0
        !Mid-edge ζ = 0
        dfdc(17) = +(1.0d0-2.0d0*c)*(1.0d0-e)*e*(1.0d0+s)*(1.0d0-s)/4.0d0
        dfdc(18) = -(1.0d0+2.0d0*c)*(1.0d0-e)*e*(1.0d0+s)*(1.0d0-s)/4.0d0
        dfdc(19) = +(1.0d0+2.0d0*c)*(1.0d0+e)*e*(1.0d0+s)*(1.0d0-s)/4.0d0
        dfdc(20) = -(1.0d0-2.0d0*c)*(1.0d0+e)*e*(1.0d0+s)*(1.0d0-s)/4.0d0
        !Mid-face ξ = +-1
        ! dfdc(22) = +(1.0d0+2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/2.0d0
        dfdc(23) = +(1.0d0+2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, ξ = +1
        ! dfdc(24) = -(1.0d0-2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/2.0d0
        dfdc(25) = -(1.0d0-2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, ξ = -1
        !Mid-face η = +-1
        ! dfdc(21) = -(-2.0d0*c)*(1.0d0-e)*e*(1.0d0+s)*(1.0d0-s)/2.0d0
        dfdc(22) = -(-2.0d0*c)*(1.0d0-e)*e*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, η = -1
        ! dfdc(23) = +(-2.0d0*c)*(1.0d0+e)*e*(1.0d0+s)*(1.0d0-s)/2.0d0
        dfdc(24) = +(-2.0d0*c)*(1.0d0+e)*e*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, η = +1
        !Mid-face ζ = +-1
        ! dfdc(25) = -(-2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)*s/2.0d0
        dfdc(21) = -(-2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)*s/2.0d0  !Salome, ζ = -1
        dfdc(26) = +(-2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*s/2.0d0
        !Central
        dfdc(27) = +(-2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)

        !d2φ/dξ2
        !Corners
        d2fdc2(1)  = -(-2.0d0)*(1.0d0-e)*e*(1.0d0-s)*s/8.0d0
        d2fdc2(2)  = +(+2.0d0)*(1.0d0-e)*e*(1.0d0-s)*s/8.0d0
        d2fdc2(3)  = -(+2.0d0)*(1.0d0+e)*e*(1.0d0-s)*s/8.0d0
        d2fdc2(4)  = +(-2.0d0)*(1.0d0+e)*e*(1.0d0-s)*s/8.0d0
        d2fdc2(5)  = +(-2.0d0)*(1.0d0-e)*e*(1.0d0+s)*s/8.0d0
        d2fdc2(6)  = -(+2.0d0)*(1.0d0-e)*e*(1.0d0+s)*s/8.0d0
        d2fdc2(7)  = +(+2.0d0)*(1.0d0+e)*e*(1.0d0+s)*s/8.0d0
        d2fdc2(8)  = -(-2.0d0)*(1.0d0+e)*e*(1.0d0+s)*s/8.0d0
        !Mid-edge ξ = 0
        d2fdc2(9)  = +(-2.0d0)*(1.0d0-e)*e*(1.0d0-s)*s/4.0d0
        d2fdc2(11) = -(-2.0d0)*(1.0d0+e)*e*(1.0d0-s)*s/4.0d0
        d2fdc2(13) = -(-2.0d0)*(1.0d0-e)*e*(1.0d0+s)*s/4.0d0
        d2fdc2(15) = +(-2.0d0)*(1.0d0+e)*e*(1.0d0+s)*s/4.0d0
        !Mid-edge η = 0
        d2fdc2(10) = -(+2.0d0)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)*s/4.0d0
        d2fdc2(12) = +(-2.0d0)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)*s/4.0d0
        d2fdc2(14) = +(+2.0d0)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*s/4.0d0
        d2fdc2(16) = -(-2.0d0)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*s/4.0d0
        !Mid-edge ζ = 0
        d2fdc2(17) = +(-2.0d0)*(1.0d0-e)*e*(1.0d0+s)*(1.0d0-s)/4.0d0
        d2fdc2(18) = -(+2.0d0)*(1.0d0-e)*e*(1.0d0+s)*(1.0d0-s)/4.0d0
        d2fdc2(19) = +(+2.0d0)*(1.0d0+e)*e*(1.0d0+s)*(1.0d0-s)/4.0d0
        d2fdc2(20) = -(-2.0d0)*(1.0d0+e)*e*(1.0d0+s)*(1.0d0-s)/4.0d0
        !Mid-face ξ = +-1
        ! d2fdc2(22) = +(+2.0d0)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/2.0d0
        d2fdc2(23) = +(+2.0d0)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, ξ = +1
        ! d2fdc2(24) = -(-2.0d0)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/2.0d0
        d2fdc2(25) = -(-2.0d0)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, ξ = -1
        !Mid-face η = +-1
        ! d2fdc2(21) = -(-2.0d0)*(1.0d0-e)*e*(1.0d0+s)*(1.0d0-s)/2.0d0
        d2fdc2(22) = -(-2.0d0)*(1.0d0-e)*e*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, η = -1
        ! d2fdc2(23) = +(-2.0d0)*(1.0d0+e)*e*(1.0d0+s)*(1.0d0-s)/2.0d0
        d2fdc2(24) = +(-2.0d0)*(1.0d0+e)*e*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, η = +1
        !Mid-face ζ = +-1
        ! d2fdc2(25) = -(-2.0d0)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)*s/2.0d0
        d2fdc2(21) = -(-2.0d0)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)*s/2.0d0  !Salome, ζ = -1
        d2fdc2(26) = +(-2.0d0)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*s/2.0d0
        !Central
        d2fdc2(27) = +(-2.0d0)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)

        !d2φ/dξdη
        !Corners
        d2fdce(1)  = -(1.0d0-2.0d0*c)*(1.0d0-2.0d0*e)*(1.0d0-s)*s/8.0d0
        d2fdce(2)  = +(1.0d0+2.0d0*c)*(1.0d0-2.0d0*e)*(1.0d0-s)*s/8.0d0
        d2fdce(3)  = -(1.0d0+2.0d0*c)*(1.0d0+2.0d0*e)*(1.0d0-s)*s/8.0d0
        d2fdce(4)  = +(1.0d0-2.0d0*c)*(1.0d0+2.0d0*e)*(1.0d0-s)*s/8.0d0
        d2fdce(5)  = +(1.0d0-2.0d0*c)*(1.0d0-2.0d0*e)*(1.0d0+s)*s/8.0d0
        d2fdce(6)  = -(1.0d0+2.0d0*c)*(1.0d0-2.0d0*e)*(1.0d0+s)*s/8.0d0
        d2fdce(7)  = +(1.0d0+2.0d0*c)*(1.0d0+2.0d0*e)*(1.0d0+s)*s/8.0d0
        d2fdce(8)  = -(1.0d0-2.0d0*c)*(1.0d0+2.0d0*e)*(1.0d0+s)*s/8.0d0
        !Mid-edge ξ = 0
        d2fdce(9)  = +(-2.0d0*c)*(1.0d0-2.0d0*e)*(1.0d0-s)*s/4.0d0
        d2fdce(11) = -(-2.0d0*c)*(1.0d0+2.0d0*e)*(1.0d0-s)*s/4.0d0
        d2fdce(13) = -(-2.0d0*c)*(1.0d0-2.0d0*e)*(1.0d0+s)*s/4.0d0
        d2fdce(15) = +(-2.0d0*c)*(1.0d0+2.0d0*e)*(1.0d0+s)*s/4.0d0
        !Mid-edge η = 0
        d2fdce(10) = -(1.0d0+2.0d0*c)*(-2.0d0*e)*(1.0d0-s)*s/4.0d0
        d2fdce(12) = +(1.0d0-2.0d0*c)*(-2.0d0*e)*(1.0d0-s)*s/4.0d0
        d2fdce(14) = +(1.0d0+2.0d0*c)*(-2.0d0*e)*(1.0d0+s)*s/4.0d0
        d2fdce(16) = -(1.0d0-2.0d0*c)*(-2.0d0*e)*(1.0d0+s)*s/4.0d0
        !Mid-edge ζ = 0
        d2fdce(17) = +(1.0d0-2.0d0*c)*(1.0d0-2.0d0*e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        d2fdce(18) = -(1.0d0+2.0d0*c)*(1.0d0-2.0d0*e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        d2fdce(19) = +(1.0d0+2.0d0*c)*(1.0d0+2.0d0*e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        d2fdce(20) = -(1.0d0-2.0d0*c)*(1.0d0+2.0d0*e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        !Mid-face ξ = +-1
        ! d2fdce(22) = +(1.0d0+2.0d0*c)*(-2.0d0*e)*(1.0d0+s)*(1.0d0-s)/2.0d0
        d2fdce(23) = +(1.0d0+2.0d0*c)*(-2.0d0*e)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, ξ = +1
        ! d2fdce(24) = -(1.0d0-2.0d0*c)*(-2.0d0*e)*(1.0d0+s)*(1.0d0-s)/2.0d0
        d2fdce(25) = -(1.0d0-2.0d0*c)*(-2.0d0*e)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, ξ = -1
        !Mid-face η = +-1
        ! d2fdce(21) = -(-2.0d0*c)*(1.0d0-2.0d0*e)*(1.0d0+s)*(1.0d0-s)/2.0d0
        d2fdce(22) = -(-2.0d0*c)*(1.0d0-2.0d0*e)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, η = -1
        ! dfdc(23) = +(-2.0d0*c)*(1.0d0+2.0d0*e)*(1.0d0+s)*(1.0d0-s)/2.0d0
        d2fdce(24) = +(-2.0d0*c)*(1.0d0+2.0d0*e)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, η = +1
        !Mid-face ζ = +-1
        ! d2fdce(25) = -(-2.0d0*c)*(-2.0d0*e)*(1.0d0-s)*s/2.0d0
        d2fdce(21) = -(-2.0d0*c)*(-2.0d0*e)*(1.0d0-s)*s/2.0d0  !Salome, ζ = -1
        d2fdce(26) = +(-2.0d0*c)*(-2.0d0*e)*(1.0d0+s)*s/2.0d0
        !Central
        d2fdce(27) = +(-2.0d0*c)*(-2.0d0*e)*(1.0d0+s)*(1.0d0-s)

        !d2φ/dξdζ
        !Corners
        d2fdcs(1)  = -(1.0d0-2.0d0*c)*(1.0d0-e)*e*(1.0d0-2.0d0*s)/8.0d0
        d2fdcs(2)  = +(1.0d0+2.0d0*c)*(1.0d0-e)*e*(1.0d0-2.0d0*s)/8.0d0
        d2fdcs(3)  = -(1.0d0+2.0d0*c)*(1.0d0+e)*e*(1.0d0-2.0d0*s)/8.0d0
        d2fdcs(4)  = +(1.0d0-2.0d0*c)*(1.0d0+e)*e*(1.0d0-2.0d0*s)/8.0d0
        d2fdcs(5)  = +(1.0d0-2.0d0*c)*(1.0d0-e)*e*(1.0d0+2.0d0*s)/8.0d0
        d2fdcs(6)  = -(1.0d0+2.0d0*c)*(1.0d0-e)*e*(1.0d0+2.0d0*s)/8.0d0
        d2fdcs(7)  = +(1.0d0+2.0d0*c)*(1.0d0+e)*e*(1.0d0+2.0d0*s)/8.0d0
        d2fdcs(8)  = -(1.0d0-2.0d0*c)*(1.0d0+e)*e*(1.0d0+2.0d0*s)/8.0d0
        !Mid-edge ξ = 0
        d2fdcs(9)  = +(-2.0d0*c)*(1.0d0-e)*e*(1.0d0-2.0d0*s)/4.0d0
        d2fdcs(11) = -(-2.0d0*c)*(1.0d0+e)*e*(1.0d0-2.0d0*s)/4.0d0
        d2fdcs(13) = -(-2.0d0*c)*(1.0d0-e)*e*(1.0d0+2.0d0*s)/4.0d0
        d2fdcs(15) = +(-2.0d0*c)*(1.0d0+e)*e*(1.0d0+2.0d0*s)/4.0d0
        !Mid-edge η = 0
        d2fdcs(10) = -(1.0d0+2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-2.0d0*s)/4.0d0
        d2fdcs(12) = +(1.0d0-2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-2.0d0*s)/4.0d0
        d2fdcs(14) = +(1.0d0+2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+2.0d0*s)/4.0d0
        d2fdcs(16) = -(1.0d0-2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+2.0d0*s)/4.0d0
        !Mid-edge ζ = 0
        d2fdcs(17) = +(1.0d0-2.0d0*c)*(1.0d0-e)*e*(-2.0d0*s)/4.0d0
        d2fdcs(18) = -(1.0d0+2.0d0*c)*(1.0d0-e)*e*(-2.0d0*s)/4.0d0
        d2fdcs(19) = +(1.0d0+2.0d0*c)*(1.0d0+e)*e*(-2.0d0*s)/4.0d0
        d2fdcs(20) = -(1.0d0-2.0d0*c)*(1.0d0+e)*e*(-2.0d0*s)/4.0d0
        !Mid-face ξ = +-1
        ! d2fdcs(22) = +(1.0d0+2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(-2.0d0*s)/2.0d0
        d2fdcs(23) = +(1.0d0+2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(-2.0d0*s)/2.0d0  !Salome, ξ = +1
        ! d2fdcs(24) = -(1.0d0-2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(-2.0d0*s)/2.0d0
        d2fdcs(25) = -(1.0d0-2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(-2.0d0*s)/2.0d0  !Salome, ξ = -1
        !Mid-face η = +-1
        ! d2fdcs(21) = -(-2.0d0*c)*(1.0d0-e)*e*(-2.0d0*s)/2.0d0
        d2fdcs(22) = -(-2.0d0*c)*(1.0d0-e)*e*(-2.0d0*s)/2.0d0  !Salome, η = -1
        ! d2fdcs(23) = +(-2.0d0*c)*(1.0d0+e)*e*(-2.0d0*s)/2.0d0
        d2fdcs(24) = +(-2.0d0*c)*(1.0d0+e)*e*(-2.0d0*s)/2.0d0  !Salome, η = +1
        !Mid-face ζ = +-1
        ! d2fdcs(25) = -(-2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-2.0d0*s)/2.0d0
        d2fdcs(21) = -(-2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-2.0d0*s)/2.0d0  !Salome, ζ = -1
        d2fdcs(26) = +(-2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+2.0d0*s)/2.0d0
        !Central
        d2fdcs(27) = +(-2.0d0*c)*(1.0d0+e)*(1.0d0-e)*(-2.0d0*s)

        !dφ/dη
        !Corners
        dfde(1)  = -(1.0d0-c)*c*(1.0d0-2.0d0*e)*(1.0d0-s)*s/8.0d0
        dfde(2)  = +(1.0d0+c)*c*(1.0d0-2.0d0*e)*(1.0d0-s)*s/8.0d0
        dfde(3)  = -(1.0d0+c)*c*(1.0d0+2.0d0*e)*(1.0d0-s)*s/8.0d0
        dfde(4)  = +(1.0d0-c)*c*(1.0d0+2.0d0*e)*(1.0d0-s)*s/8.0d0
        dfde(5)  = +(1.0d0-c)*c*(1.0d0-2.0d0*e)*(1.0d0+s)*s/8.0d0
        dfde(6)  = -(1.0d0+c)*c*(1.0d0-2.0d0*e)*(1.0d0+s)*s/8.0d0
        dfde(7)  = +(1.0d0+c)*c*(1.0d0+2.0d0*e)*(1.0d0+s)*s/8.0d0
        dfde(8)  = -(1.0d0-c)*c*(1.0d0+2.0d0*e)*(1.0d0+s)*s/8.0d0
        !Mid-edge ξ = 0
        dfde(9)  = +(1.0d0+c)*(1.0d0-c)*(1.0d0-2.0d0*e)*(1.0d0-s)*s/4.0d0
        dfde(11) = -(1.0d0+c)*(1.0d0-c)*(1.0d0+2.0d0*e)*(1.0d0-s)*s/4.0d0
        dfde(13) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-2.0d0*e)*(1.0d0+s)*s/4.0d0
        dfde(15) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+2.0d0*e)*(1.0d0+s)*s/4.0d0
        !Mid-edge η = 0
        dfde(10) = -(1.0d0+c)*c*(-2.0d0*e)*(1.0d0-s)*s/4.0d0
        dfde(12) = +(1.0d0-c)*c*(-2.0d0*e)*(1.0d0-s)*s/4.0d0
        dfde(14) = +(1.0d0+c)*c*(-2.0d0*e)*(1.0d0+s)*s/4.0d0
        dfde(16) = -(1.0d0-c)*c*(-2.0d0*e)*(1.0d0+s)*s/4.0d0
        !Mid-edge ζ = 0
        dfde(17) = +(1.0d0-c)*c*(1.0d0-2.0d0*e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        dfde(18) = -(1.0d0+c)*c*(1.0d0-2.0d0*e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        dfde(19) = +(1.0d0+c)*c*(1.0d0+2.0d0*e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        dfde(20) = -(1.0d0-c)*c*(1.0d0+2.0d0*e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        !Mid-face ξ = +-1
        ! dfde(22) = +(1.0d0+c)*c*(-2.0d0*e)*(1.0d0+s)*(1.0d0-s)/2.0d0
        dfde(23) = +(1.0d0+c)*c*(-2.0d0*e)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, ξ = +1
        ! dfde(24) = -(1.0d0-c)*c*(-2.0d0*e)*(1.0d0+s)*(1.0d0-s)/2.0d0
        dfde(25) = -(1.0d0-c)*c*(-2.0d0*e)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, ξ = -1
        !Mid-face η = +-1
        ! dfde(21) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-2.0d0*e)*(1.0d0+s)*(1.0d0-s)/2.0d0
        dfde(22) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-2.0d0*e)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, η = -1
        ! dfde(23) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+2.0d0*e)*(1.0d0+s)*(1.0d0-s)/2.0d0
        dfde(24) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+2.0d0*e)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, η = +1
        !Mid-face ζ = +-1
        ! dfde(25) = -(1.0d0+c)*(1.0d0-c)*(-2.0d0*e)*(1.0d0-s)*s/2.0d0
        dfde(21) = -(1.0d0+c)*(1.0d0-c)*(-2.0d0*e)*(1.0d0-s)*s/2.0d0  !Salome, ζ = -1
        dfde(26) = +(1.0d0+c)*(1.0d0-c)*(-2.0d0*e)*(1.0d0+s)*s/2.0d0
        !Central
        dfde(27) = +(1.0d0+c)*(1.0d0-c)*(-2.0d0*e)*(1.0d0+s)*(1.0d0-s)

        !d2φ/dη2
        !Corners
        d2fde2(1)  = -(1.0d0-c)*c*(-2.0d0)*(1.0d0-s)*s/8.0d0
        d2fde2(2)  = +(1.0d0+c)*c*(-2.0d0)*(1.0d0-s)*s/8.0d0
        d2fde2(3)  = -(1.0d0+c)*c*(+2.0d0)*(1.0d0-s)*s/8.0d0
        d2fde2(4)  = +(1.0d0-c)*c*(+2.0d0)*(1.0d0-s)*s/8.0d0
        d2fde2(5)  = +(1.0d0-c)*c*(-2.0d0)*(1.0d0+s)*s/8.0d0
        d2fde2(6)  = -(1.0d0+c)*c*(-2.0d0)*(1.0d0+s)*s/8.0d0
        d2fde2(7)  = +(1.0d0+c)*c*(+2.0d0)*(1.0d0+s)*s/8.0d0
        d2fde2(8)  = -(1.0d0-c)*c*(+2.0d0)*(1.0d0+s)*s/8.0d0
        !Mid-edge ξ = 0
        d2fde2(9)  = +(1.0d0+c)*(1.0d0-c)*(-2.0d0)*(1.0d0-s)*s/4.0d0
        d2fde2(11) = -(1.0d0+c)*(1.0d0-c)*(+2.0d0)*(1.0d0-s)*s/4.0d0
        d2fde2(13) = -(1.0d0+c)*(1.0d0-c)*(-2.0d0)*(1.0d0+s)*s/4.0d0
        d2fde2(15) = +(1.0d0+c)*(1.0d0-c)*(+2.0d0)*(1.0d0+s)*s/4.0d0
        !Mid-edge η = 0
        d2fde2(10) = -(1.0d0+c)*c*(-2.0d0)*(1.0d0-s)*s/4.0d0
        d2fde2(12) = +(1.0d0-c)*c*(-2.0d0)*(1.0d0-s)*s/4.0d0
        d2fde2(14) = +(1.0d0+c)*c*(-2.0d0)*(1.0d0+s)*s/4.0d0
        d2fde2(16) = -(1.0d0-c)*c*(-2.0d0)*(1.0d0+s)*s/4.0d0
        !Mid-edge ζ = 0
        d2fde2(17) = +(1.0d0-c)*c*(-2.0d0)*(1.0d0+s)*(1.0d0-s)/4.0d0
        d2fde2(18) = -(1.0d0+c)*c*(-2.0d0)*(1.0d0+s)*(1.0d0-s)/4.0d0
        d2fde2(19) = +(1.0d0+c)*c*(+2.0d0)*(1.0d0+s)*(1.0d0-s)/4.0d0
        d2fde2(20) = -(1.0d0-c)*c*(+2.0d0)*(1.0d0+s)*(1.0d0-s)/4.0d0
        !Mid-face ξ = +-1
        ! d2fde2(22) = +(1.0d0+c)*c*(-2.0d0)*(1.0d0+s)*(1.0d0-s)/2.0d0
        d2fde2(23) = +(1.0d0+c)*c*(-2.0d0)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, ξ = +1
        ! d2fde2(24) = -(1.0d0-c)*c*(-2.0d0)*(1.0d0+s)*(1.0d0-s)/2.0d0
        d2fde2(25) = -(1.0d0-c)*c*(-2.0d0)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, ξ = -1
        !Mid-face η = +-1
        ! d2fde2(21) = -(1.0d0+c)*(1.0d0-c)*(-2.0d0)*(1.0d0+s)*(1.0d0-s)/2.0d0
        d2fde2(22) = -(1.0d0+c)*(1.0d0-c)*(-2.0d0)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, η = -1
        ! d2fde2(23) = +(1.0d0+c)*(1.0d0-c)*(+2.0d0)*(1.0d0+s)*(1.0d0-s)/2.0d0
        d2fde2(24) = +(1.0d0+c)*(1.0d0-c)*(+2.0d0)*(1.0d0+s)*(1.0d0-s)/2.0d0  !Salome, η = +1
        !Mid-face ζ = +-1
        ! d2fde2(25) = -(1.0d0+c)*(1.0d0-c)*(-2.0d0)*(1.0d0-s)*s/2.0d0
        d2fde2(21) = -(1.0d0+c)*(1.0d0-c)*(-2.0d0)*(1.0d0-s)*s/2.0d0  !Salome, ζ = -1
        d2fde2(26) = +(1.0d0+c)*(1.0d0-c)*(-2.0d0)*(1.0d0+s)*s/2.0d0
        !Central
        d2fde2(27) = +(1.0d0+c)*(1.0d0-c)*(-2.0d0)*(1.0d0+s)*(1.0d0-s)

        !d2φ/dηζ
        !Corners
        d2fdes(1)  = -(1.0d0-c)*c*(1.0d0-2.0d0*e)*(1.0d0-2.0d0*s)/8.0d0
        d2fdes(2)  = +(1.0d0+c)*c*(1.0d0-2.0d0*e)*(1.0d0-2.0d0*s)/8.0d0
        d2fdes(3)  = -(1.0d0+c)*c*(1.0d0+2.0d0*e)*(1.0d0-2.0d0*s)/8.0d0
        d2fdes(4)  = +(1.0d0-c)*c*(1.0d0+2.0d0*e)*(1.0d0-2.0d0*s)/8.0d0
        d2fdes(5)  = +(1.0d0-c)*c*(1.0d0-2.0d0*e)*(1.0d0+2.0d0*s)/8.0d0
        d2fdes(6)  = -(1.0d0+c)*c*(1.0d0-2.0d0*e)*(1.0d0+2.0d0*s)/8.0d0
        d2fdes(7)  = +(1.0d0+c)*c*(1.0d0+2.0d0*e)*(1.0d0+2.0d0*s)/8.0d0
        d2fdes(8)  = -(1.0d0-c)*c*(1.0d0+2.0d0*e)*(1.0d0+2.0d0*s)/8.0d0
        !Mid-edge ξ = 0
        d2fdes(9)  = +(1.0d0+c)*(1.0d0-c)*(1.0d0-2.0d0*e)*(1.0d0-2.0d0*s)/4.0d0
        d2fdes(11) = -(1.0d0+c)*(1.0d0-c)*(1.0d0+2.0d0*e)*(1.0d0-2.0d0*s)/4.0d0
        d2fdes(13) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-2.0d0*e)*(1.0d0+2.0d0*s)/4.0d0
        d2fdes(15) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+2.0d0*e)*(1.0d0+2.0d0*s)/4.0d0
        !Mid-edge η = 0
        d2fdes(10) = -(1.0d0+c)*c*(-2.0d0*e)*(1.0d0-2.0d0*s)/4.0d0
        d2fdes(12) = +(1.0d0-c)*c*(-2.0d0*e)*(1.0d0-2.0d0*s)/4.0d0
        d2fdes(14) = +(1.0d0+c)*c*(-2.0d0*e)*(1.0d0+2.0d0*s)/4.0d0
        d2fdes(16) = -(1.0d0-c)*c*(-2.0d0*e)*(1.0d0+2.0d0*s)/4.0d0
        !Mid-edge ζ = 0
        d2fdes(17) = +(1.0d0-c)*c*(1.0d0-2.0d0*e)*(-2.0d0*s)/4.0d0
        d2fdes(18) = -(1.0d0+c)*c*(1.0d0-2.0d0*e)*(-2.0d0*s)/4.0d0
        d2fdes(19) = +(1.0d0+c)*c*(1.0d0+2.0d0*e)*(-2.0d0*s)/4.0d0
        d2fdes(20) = -(1.0d0-c)*c*(1.0d0+2.0d0*e)*(-2.0d0*s)/4.0d0
        !Mid-face ξ = +-1
        ! d2fdes(22) = +(1.0d0+c)*c*(-2.0d0*e)*(-2.0d0*s)/2.0d0
        d2fdes(23) = +(1.0d0+c)*c*(-2.0d0*e)*(-2.0d0*s)/2.0d0  !Salome, ξ = +1
        ! d2fdes(24) = -(1.0d0-c)*c*(-2.0d0*e)*(-2.0d0*s)/2.0d0
        d2fdes(25) = -(1.0d0-c)*c*(-2.0d0*e)*(-2.0d0*s)/2.0d0  !Salome, ξ = -1
        !Mid-face η = +-1
        ! d2fdes(21) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-2.0d0*e)*(-2.0d0*s)/2.0d0
        d2fdes(22) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-2.0d0*e)*(-2.0d0*s)/2.0d0  !Salome, η = -1
        ! d2fdes(23) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+2.0d0*e)*(-2.0d0*s)/2.0d0
        d2fdes(24) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+2.0d0*e)*(-2.0d0*s)/2.0d0  !Salome, η = +1
        !Mid-face ζ = +-1
        ! d2fdes(25) = -(1.0d0+c)*(1.0d0-c)*(-2.0d0*e)*(1.0d0-2.0d0*s)/2.0d0
        d2fdes(21) = -(1.0d0+c)*(1.0d0-c)*(-2.0d0*e)*(1.0d0-2.0d0*s)/2.0d0  !Salome, ζ = -1
        d2fdes(26) = +(1.0d0+c)*(1.0d0-c)*(-2.0d0*e)*(1.0d0+2.0d0*s)/2.0d0
        !Central
        d2fdes(27) = +(1.0d0+c)*(1.0d0-c)*(-2.0d0*e)*(-2.0d0*s)
        
        !dφ/dζ
        !Corners
        dfds(1)  = -(1.0d0-c)*c*(1.0d0-e)*e*(1.0d0-2.0d0*s)/8.0d0
        dfds(2)  = +(1.0d0+c)*c*(1.0d0-e)*e*(1.0d0-2.0d0*s)/8.0d0
        dfds(3)  = -(1.0d0+c)*c*(1.0d0+e)*e*(1.0d0-2.0d0*s)/8.0d0
        dfds(4)  = +(1.0d0-c)*c*(1.0d0+e)*e*(1.0d0-2.0d0*s)/8.0d0
        dfds(5)  = +(1.0d0-c)*c*(1.0d0-e)*e*(1.0d0+2.0d0*s)/8.0d0
        dfds(6)  = -(1.0d0+c)*c*(1.0d0-e)*e*(1.0d0+2.0d0*s)/8.0d0
        dfds(7)  = +(1.0d0+c)*c*(1.0d0+e)*e*(1.0d0+2.0d0*s)/8.0d0
        dfds(8)  = -(1.0d0-c)*c*(1.0d0+e)*e*(1.0d0+2.0d0*s)/8.0d0
        !Mid-edge ξ = 0
        dfds(9)  = +(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e*(1.0d0-2.0d0*s)/4.0d0
        dfds(11) = -(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e*(1.0d0-2.0d0*s)/4.0d0
        dfds(13) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e*(1.0d0+2.0d0*s)/4.0d0
        dfds(15) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e*(1.0d0+2.0d0*s)/4.0d0
        !Mid-edge η = 0
        dfds(10) = -(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0-2.0d0*s)/4.0d0
        dfds(12) = +(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0-2.0d0*s)/4.0d0
        dfds(14) = +(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0+2.0d0*s)/4.0d0
        dfds(16) = -(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0+2.0d0*s)/4.0d0
        !Mid-edge ζ = 0
        dfds(17) = +(1.0d0-c)*c*(1.0d0-e)*e*(-2.0d0*s)/4.0d0
        dfds(18) = -(1.0d0+c)*c*(1.0d0-e)*e*(-2.0d0*s)/4.0d0
        dfds(19) = +(1.0d0+c)*c*(1.0d0+e)*e*(-2.0d0*s)/4.0d0
        dfds(20) = -(1.0d0-c)*c*(1.0d0+e)*e*(-2.0d0*s)/4.0d0
        !Mid-face ξ = +-1
        ! dfds(22) = +(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)*(-2.0d0*s)/2.0d0
        dfds(23) = +(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)*(-2.0d0*s)/2.0d0  !Salome, ξ = +1
        ! dfds(24) = -(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)*(-2.0d0*s)/2.0d0
        dfds(25) = -(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)*(-2.0d0*s)/2.0d0  !Salome, ξ = -1
        !Mid-face η = +-1
        ! dfds(21) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e*(-2.0d0*s)/2.0d0
        dfds(22) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e*(-2.0d0*s)/2.0d0  !Salome, η = -1
        ! dfds(23) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e*(-2.0d0*s)/2.0d0
        dfds(24) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e*(-2.0d0*s)/2.0d0  !Salome, η = +1
        !Mid-face ζ = +-1
        ! dfds(25) = -(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-2.0d0*s)/2.0d0
        dfds(21) = -(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-2.0d0*s)/2.0d0  !Salome, ζ = -1
        dfds(26) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+2.0d0*s)/2.0d0
        !Central
        dfds(27) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(-2.0d0*s)

        !d2φ/dζ2
        !Corners
        d2fds2(1)  = -(1.0d0-c)*c*(1.0d0-e)*e*(-2.0d0)/8.0d0
        d2fds2(2)  = +(1.0d0+c)*c*(1.0d0-e)*e*(-2.0d0)/8.0d0
        d2fds2(3)  = -(1.0d0+c)*c*(1.0d0+e)*e*(-2.0d0)/8.0d0
        d2fds2(4)  = +(1.0d0-c)*c*(1.0d0+e)*e*(-2.0d0)/8.0d0
        d2fds2(5)  = +(1.0d0-c)*c*(1.0d0-e)*e*(+2.0d0)/8.0d0
        d2fds2(6)  = -(1.0d0+c)*c*(1.0d0-e)*e*(+2.0d0)/8.0d0
        d2fds2(7)  = +(1.0d0+c)*c*(1.0d0+e)*e*(+2.0d0)/8.0d0
        d2fds2(8)  = -(1.0d0-c)*c*(1.0d0+e)*e*(+2.0d0)/8.0d0
        !Mid-edge ξ = 0
        d2fds2(9)  = +(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e*(-2.0d0)/4.0d0
        d2fds2(11) = -(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e*(-2.0d0)/4.0d0
        d2fds2(13) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e*(+2.0d0)/4.0d0
        d2fds2(15) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e*(+2.0d0)/4.0d0
        !Mid-edge η = 0
        d2fds2(10) = -(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)*(-2.0d0)/4.0d0
        d2fds2(12) = +(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)*(-2.0d0)/4.0d0
        d2fds2(14) = +(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)*(+2.0d0)/4.0d0
        d2fds2(16) = -(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)*(+2.0d0)/4.0d0
        !Mid-edge ζ = 0
        d2fds2(17) = +(1.0d0-c)*c*(1.0d0-e)*e*(-2.0d0)/4.0d0
        d2fds2(18) = -(1.0d0+c)*c*(1.0d0-e)*e*(-2.0d0)/4.0d0
        d2fds2(19) = +(1.0d0+c)*c*(1.0d0+e)*e*(-2.0d0)/4.0d0
        d2fds2(20) = -(1.0d0-c)*c*(1.0d0+e)*e*(-2.0d0)/4.0d0
        !Mid-face ξ = +-1
        ! d2fds2(22) = +(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)*(-2.0d0)/2.0d0
        d2fds2(23) = +(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)*(-2.0d0)/2.0d0  !Salome, ξ = +1
        ! d2fds2(24) = -(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)*(-2.0d0)/2.0d0
        d2fds2(25) = -(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)*(-2.0d0)/2.0d0  !Salome, ξ = -1
        !Mid-face η = +-1
        ! d2fds2(21) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e*(-2.0d0)/2.0d0
        d2fds2(22) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e*(-2.0d0)/2.0d0  !Salome, η = -1
        ! d2fds2(23) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e*(-2.0d0)/2.0d0
        d2fds2(24) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e*(-2.0d0)/2.0d0  !Salome, η = +1
        !Mid-face ζ = +-1
        ! d2fds2(25) = -(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(-2.0d0)/2.0d0
        d2fds2(21) = -(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(-2.0d0)/2.0d0  !Salome, ζ = -1
        d2fds2(26) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(+2.0d0)/2.0d0
        !Central
        d2fds2(27) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(-2.0d0)

    end subroutine setQuadraticHexahedron

    !----------------------------------------------------------------------- 

    subroutine setSerendipityHexahedron(bfn, dfdc, dfde, dfds, &
        d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes, c, e, s)

        implicit none

        PetscReal, intent(in) :: c, e, s
        PetscReal, dimension(:), intent(out) :: bfn, dfdc, dfde, dfds, &
            d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes

        !Basis function, φ
        !Corners
        bfn(1)  = -(1.0d0-c)*(1.0d0-e)*(1.0d0-s)*(2.0d0+c+e+s)/8.0d0
        bfn(2)  = -(1.0d0+c)*(1.0d0-e)*(1.0d0-s)*(2.0d0-c+e+s)/8.0d0
        bfn(3)  = -(1.0d0+c)*(1.0d0+e)*(1.0d0-s)*(2.0d0-c-e+s)/8.0d0
        bfn(4)  = -(1.0d0-c)*(1.0d0+e)*(1.0d0-s)*(2.0d0+c-e+s)/8.0d0
        bfn(5)  = -(1.0d0-c)*(1.0d0-e)*(1.0d0+s)*(2.0d0+c+e-s)/8.0d0
        bfn(6)  = -(1.0d0+c)*(1.0d0-e)*(1.0d0+s)*(2.0d0-c+e-s)/8.0d0
        bfn(7)  = -(1.0d0+c)*(1.0d0+e)*(1.0d0+s)*(2.0d0-c-e-s)/8.0d0
        bfn(8)  = -(1.0d0-c)*(1.0d0+e)*(1.0d0+s)*(2.0d0+c-e-s)/8.0d0
        !Mid-edge ξ = 0
        bfn(9)  = +(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*(1.0d0-s)/4.0d0
        bfn(11) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-s)/4.0d0
        bfn(13) = +(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*(1.0d0+s)/4.0d0
        bfn(15) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0+s)/4.0d0
        !Mid-edge η = 0
        bfn(10) = +(1.0d0+c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)/4.0d0
        bfn(12) = +(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)/4.0d0
        bfn(14) = +(1.0d0+c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)/4.0d0
        bfn(16) = +(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)/4.0d0
        !Mid-edge ζ = 0
        bfn(17) = +(1.0d0-c)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        bfn(18) = +(1.0d0+c)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        bfn(19) = +(1.0d0+c)*(1.0d0+e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        bfn(20) = +(1.0d0-c)*(1.0d0+e)*(1.0d0+s)*(1.0d0-s)/4.0d0

        !dφ/dξ
        !Corners
        dfdc(1)  = -((-1.0d0)*(2.0d0+c+e+s)+(1.0d0-c)*(+1.0d0))*(1.0d0-e)*(1.0d0-s)/8.0d0
        dfdc(2)  = -((+1.0d0)*(2.0d0-c+e+s)+(1.0d0+c)*(-1.0d0))*(1.0d0-e)*(1.0d0-s)/8.0d0
        dfdc(3)  = -((+1.0d0)*(2.0d0-c-e+s)+(1.0d0+c)*(-1.0d0))*(1.0d0+e)*(1.0d0-s)/8.0d0
        dfdc(4)  = -((-1.0d0)*(2.0d0+c-e+s)+(1.0d0-c)*(+1.0d0))*(1.0d0+e)*(1.0d0-s)/8.0d0
        dfdc(5)  = -((-1.0d0)*(2.0d0+c+e-s)+(1.0d0-c)*(+1.0d0))*(1.0d0-e)*(1.0d0+s)/8.0d0
        dfdc(6)  = -((+1.0d0)*(2.0d0-c+e-s)+(1.0d0+c)*(-1.0d0))*(1.0d0-e)*(1.0d0+s)/8.0d0
        dfdc(7)  = -((+1.0d0)*(2.0d0-c-e-s)+(1.0d0+c)*(-1.0d0))*(1.0d0+e)*(1.0d0+s)/8.0d0
        dfdc(8)  = -((-1.0d0)*(2.0d0+c-e-s)+(1.0d0-c)*(+1.0d0))*(1.0d0+e)*(1.0d0+s)/8.0d0
        !Mid-edge ξ = 0
        dfdc(9)  = +(-2.0d0*c)*(1.0d0-e)*(1.0d0-s)/4.0d0
        dfdc(11) = +(-2.0d0*c)*(1.0d0+e)*(1.0d0-s)/4.0d0
        dfdc(13) = +(-2.0d0*c)*(1.0d0-e)*(1.0d0+s)/4.0d0
        dfdc(15) = +(-2.0d0*c)*(1.0d0+e)*(1.0d0+s)/4.0d0
        !Mid-edge η = 0
        dfdc(10) = +(+1.0d0)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)/4.0d0
        dfdc(12) = +(-1.0d0)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)/4.0d0
        dfdc(14) = +(+1.0d0)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)/4.0d0
        dfdc(16) = +(-1.0d0)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)/4.0d0
        !Mid-edge ζ = 0
        dfdc(17) = +(-1.0d0)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        dfdc(18) = +(+1.0d0)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        dfdc(19) = +(+1.0d0)*(1.0d0+e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        dfdc(20) = +(-1.0d0)*(1.0d0+e)*(1.0d0+s)*(1.0d0-s)/4.0d0

        !d2φ/dξ2
        !Corners
        d2fdc2(1)  = -((-1.0d0)*(+1.0d0)+(-1.0d0)*(+1.0d0))*(1.0d0-e)*(1.0d0-s)/8.0d0
        d2fdc2(2)  = -((+1.0d0)*(-1.0d0)+(+1.0d0)*(-1.0d0))*(1.0d0-e)*(1.0d0-s)/8.0d0
        d2fdc2(3)  = -((+1.0d0)*(-1.0d0)+(+1.0d0)*(-1.0d0))*(1.0d0+e)*(1.0d0-s)/8.0d0
        d2fdc2(4)  = -((-1.0d0)*(+1.0d0)+(-1.0d0)*(+1.0d0))*(1.0d0+e)*(1.0d0-s)/8.0d0
        d2fdc2(5)  = -((-1.0d0)*(+1.0d0)+(-1.0d0)*(+1.0d0))*(1.0d0-e)*(1.0d0+s)/8.0d0
        d2fdc2(6)  = -((+1.0d0)*(-1.0d0)+(+1.0d0)*(-1.0d0))*(1.0d0-e)*(1.0d0+s)/8.0d0
        d2fdc2(7)  = -((+1.0d0)*(-1.0d0)+(+1.0d0)*(-1.0d0))*(1.0d0+e)*(1.0d0+s)/8.0d0
        d2fdc2(8)  = -((-1.0d0)*(+1.0d0)+(-1.0d0)*(+1.0d0))*(1.0d0+e)*(1.0d0+s)/8.0d0
        !Mid-edge ξ = 0
        d2fdc2(9)  = +(-2.0d0)*(1.0d0-e)*(1.0d0-s)/4.0d0
        d2fdc2(11) = +(-2.0d0)*(1.0d0+e)*(1.0d0-s)/4.0d0
        d2fdc2(13) = +(-2.0d0)*(1.0d0-e)*(1.0d0+s)/4.0d0
        d2fdc2(15) = +(-2.0d0)*(1.0d0+e)*(1.0d0+s)/4.0d0
        !Mid-edge η = 0
        d2fdc2(10) = +0.0d0
        d2fdc2(12) = +0.0d0
        d2fdc2(14) = +0.0d0
        d2fdc2(16) = +0.0d0
        !Mid-edge ζ = 0
        d2fdc2(17) = +0.0d0
        d2fdc2(18) = +0.0d0
        d2fdc2(19) = +0.0d0
        d2fdc2(20) = +0.0d0

        !d2φ/dξdη
        !Corners
        d2fdce(1)  = -( ((-1.0d0)*(+1.0d0)+0.0d0)*(1.0d0-e) + ((-1.0d0)*(2.0d0+c+e+s)+(1.0d0-c)*(+1.0d0))*(-1.0d0) )*(1.0d0-s)/8.0d0
        d2fdce(2)  = -( ((+1.0d0)*(+1.0d0)+0.0d0)*(1.0d0-e) + ((+1.0d0)*(2.0d0-c+e+s)+(1.0d0+c)*(-1.0d0))*(-1.0d0) )*(1.0d0-s)/8.0d0
        d2fdce(3)  = -( ((+1.0d0)*(-1.0d0)+0.0d0)*(1.0d0+e) + ((+1.0d0)*(2.0d0-c-e+s)+(1.0d0+c)*(-1.0d0))*(+1.0d0) )*(1.0d0-s)/8.0d0
        d2fdce(4)  = -( ((-1.0d0)*(-1.0d0)+0.0d0)*(1.0d0+e) + ((-1.0d0)*(2.0d0+c-e+s)+(1.0d0-c)*(+1.0d0))*(+1.0d0) )*(1.0d0-s)/8.0d0
        d2fdce(5)  = -( ((-1.0d0)*(+1.0d0)+0.0d0)*(1.0d0-e) + ((-1.0d0)*(2.0d0+c+e-s)+(1.0d0-c)*(+1.0d0))*(-1.0d0) )*(1.0d0+s)/8.0d0
        d2fdce(6)  = -( ((+1.0d0)*(+1.0d0)+0.0d0)*(1.0d0-e) + ((+1.0d0)*(2.0d0-c+e-s)+(1.0d0+c)*(-1.0d0))*(-1.0d0) )*(1.0d0+s)/8.0d0
        d2fdce(7)  = -( ((+1.0d0)*(-1.0d0)+0.0d0)*(1.0d0+e) + ((+1.0d0)*(2.0d0-c-e-s)+(1.0d0+c)*(-1.0d0))*(+1.0d0) )*(1.0d0+s)/8.0d0
        d2fdce(8)  = -( ((-1.0d0)*(-1.0d0)+0.0d0)*(1.0d0+e) + ((-1.0d0)*(2.0d0+c-e-s)+(1.0d0-c)*(+1.0d0))*(+1.0d0) )*(1.0d0+s)/8.0d0
        !Mid-edge ξ = 0
        d2fdce(9)  = +(-2.0d0*c)*(-1.0d0)*(1.0d0-s)/4.0d0
        d2fdce(11) = +(-2.0d0*c)*(+1.0d0)*(1.0d0-s)/4.0d0
        d2fdce(13) = +(-2.0d0*c)*(-1.0d0)*(1.0d0+s)/4.0d0
        d2fdce(15) = +(-2.0d0*c)*(+1.0d0)*(1.0d0+s)/4.0d0
        !Mid-edge η = 0
        d2fdce(10) = +(+1.0d0)*(-2.0d0*e)*(1.0d0-s)/4.0d0
        d2fdce(12) = +(-1.0d0)*(-2.0d0*e)*(1.0d0-s)/4.0d0
        d2fdce(14) = +(+1.0d0)*(-2.0d0*e)*(1.0d0+s)/4.0d0
        d2fdce(16) = +(-1.0d0)*(-2.0d0*e)*(1.0d0+s)/4.0d0
        !Mid-edge ζ = 0
        d2fdce(17) = +(-1.0d0)*(-1.0d0)*(1.0d0+s)*(1.0d0-s)/4.0d0
        d2fdce(18) = +(+1.0d0)*(-1.0d0)*(1.0d0+s)*(1.0d0-s)/4.0d0
        d2fdce(19) = +(+1.0d0)*(+1.0d0)*(1.0d0+s)*(1.0d0-s)/4.0d0
        d2fdce(20) = +(-1.0d0)*(+1.0d0)*(1.0d0+s)*(1.0d0-s)/4.0d0

        !d2φ/dξdζ
        !Corners
        d2fdcs(1)  = -( ((-1.0d0)*(+1.0d0)+0.0d0)*(1.0d0-s)+((-1.0d0)*(2.0d0+c+e+s)+(1.0d0-c)*(+1.0d0))*(-1.0d0) )*(1.0d0-e)/8.0d0
        d2fdcs(2)  = -( ((+1.0d0)*(+1.0d0)+0.0d0)*(1.0d0-s)+((+1.0d0)*(2.0d0-c+e+s)+(1.0d0+c)*(-1.0d0))*(-1.0d0) )*(1.0d0-e)/8.0d0
        d2fdcs(3)  = -( ((+1.0d0)*(+1.0d0)+0.0d0)*(1.0d0-s)+((+1.0d0)*(2.0d0-c-e+s)+(1.0d0+c)*(-1.0d0))*(-1.0d0) )*(1.0d0+e)/8.0d0
        d2fdcs(4)  = -( ((-1.0d0)*(+1.0d0)+0.0d0)*(1.0d0-s)+((-1.0d0)*(2.0d0+c-e+s)+(1.0d0-c)*(+1.0d0))*(-1.0d0) )*(1.0d0+e)/8.0d0
        d2fdcs(5)  = -( ((-1.0d0)*(-1.0d0)+0.0d0)*(1.0d0+s)+((-1.0d0)*(2.0d0+c+e-s)+(1.0d0-c)*(+1.0d0))*(+1.0d0) )*(1.0d0-e)/8.0d0
        d2fdcs(6)  = -( ((+1.0d0)*(-1.0d0)+0.0d0)*(1.0d0+s)+((+1.0d0)*(2.0d0-c+e-s)+(1.0d0+c)*(-1.0d0))*(+1.0d0) )*(1.0d0-e)/8.0d0
        d2fdcs(7)  = -( ((+1.0d0)*(-1.0d0)+0.0d0)*(1.0d0+s)+((+1.0d0)*(2.0d0-c-e-s)+(1.0d0+c)*(-1.0d0))*(+1.0d0) )*(1.0d0+e)/8.0d0
        d2fdcs(8)  = -( ((-1.0d0)*(-1.0d0)+0.0d0)*(1.0d0+s)+((-1.0d0)*(2.0d0+c-e-s)+(1.0d0-c)*(+1.0d0))*(+1.0d0) )*(1.0d0+e)/8.0d0
        !Mid-edge ξ = 0
        d2fdcs(9)  = +(-2.0d0*c)*(1.0d0-e)*(-1.0d0)/4.0d0
        d2fdcs(11) = +(-2.0d0*c)*(1.0d0+e)*(-1.0d0)/4.0d0
        d2fdcs(13) = +(-2.0d0*c)*(1.0d0-e)*(+1.0d0)/4.0d0
        d2fdcs(15) = +(-2.0d0*c)*(1.0d0+e)*(+1.0d0)/4.0d0
        !Mid-edge η = 0
        d2fdcs(10) = +(+1.0d0)*(1.0d0+e)*(1.0d0-e)*(-1.0d0)/4.0d0
        d2fdcs(12) = +(-1.0d0)*(1.0d0+e)*(1.0d0-e)*(-1.0d0)/4.0d0
        d2fdcs(14) = +(+1.0d0)*(1.0d0+e)*(1.0d0-e)*(+1.0d0)/4.0d0
        d2fdcs(16) = +(-1.0d0)*(1.0d0+e)*(1.0d0-e)*(+1.0d0)/4.0d0
        !Mid-edge ζ = 0
        d2fdcs(17) = +(-1.0d0)*(1.0d0-e)*(-2.0d0*s)/4.0d0
        d2fdcs(18) = +(+1.0d0)*(1.0d0-e)*(-2.0d0*s)/4.0d0
        d2fdcs(19) = +(+1.0d0)*(1.0d0+e)*(-2.0d0*s)/4.0d0
        d2fdcs(20) = +(-1.0d0)*(1.0d0+e)*(-2.0d0*s)/4.0d0

        !dφ/dη
        !Corners
        dfde(1)  = -(1.0d0-c)*( (-1.0d0)*(2.0d0+c+e+s)+(1.0d0-e)*(+1.0d0) )*(1.0d0-s)/8.0d0
        dfde(2)  = -(1.0d0+c)*( (-1.0d0)*(2.0d0-c+e+s)+(1.0d0-e)*(+1.0d0) )*(1.0d0-s)/8.0d0
        dfde(3)  = -(1.0d0+c)*( (+1.0d0)*(2.0d0-c-e+s)+(1.0d0+e)*(-1.0d0) )*(1.0d0-s)/8.0d0
        dfde(4)  = -(1.0d0-c)*( (+1.0d0)*(2.0d0+c-e+s)+(1.0d0+e)*(-1.0d0) )*(1.0d0-s)/8.0d0
        dfde(5)  = -(1.0d0-c)*( (-1.0d0)*(2.0d0+c+e-s)+(1.0d0-e)*(+1.0d0) )*(1.0d0+s)/8.0d0
        dfde(6)  = -(1.0d0+c)*( (-1.0d0)*(2.0d0-c+e-s)+(1.0d0-e)*(+1.0d0) )*(1.0d0+s)/8.0d0
        dfde(7)  = -(1.0d0+c)*( (+1.0d0)*(2.0d0-c-e-s)+(1.0d0+e)*(-1.0d0) )*(1.0d0+s)/8.0d0
        dfde(8)  = -(1.0d0-c)*( (+1.0d0)*(2.0d0+c-e-s)+(1.0d0+e)*(-1.0d0) )*(1.0d0+s)/8.0d0
        !Mid-edge ξ = 0
        dfde(9)  = +(1.0d0+c)*(1.0d0-c)*(-1.0d0)*(1.0d0-s)/4.0d0
        dfde(11) = +(1.0d0+c)*(1.0d0-c)*(+1.0d0)*(1.0d0-s)/4.0d0
        dfde(13) = +(1.0d0+c)*(1.0d0-c)*(-1.0d0)*(1.0d0+s)/4.0d0
        dfde(15) = +(1.0d0+c)*(1.0d0-c)*(+1.0d0)*(1.0d0+s)/4.0d0
        !Mid-edge η = 0
        dfde(10) = +(1.0d0+c)*(-2.0d0*e)*(1.0d0-s)/4.0d0
        dfde(12) = +(1.0d0-c)*(-2.0d0*e)*(1.0d0-s)/4.0d0
        dfde(14) = +(1.0d0+c)*(-2.0d0*e)*(1.0d0+s)/4.0d0
        dfde(16) = +(1.0d0-c)*(-2.0d0*e)*(1.0d0+s)/4.0d0
        !Mid-edge ζ = 0
        dfde(17) = +(1.0d0-c)*(-1.0d0)*(1.0d0+s)*(1.0d0-s)/4.0d0
        dfde(18) = +(1.0d0+c)*(-1.0d0)*(1.0d0+s)*(1.0d0-s)/4.0d0
        dfde(19) = +(1.0d0+c)*(+1.0d0)*(1.0d0+s)*(1.0d0-s)/4.0d0
        dfde(20) = +(1.0d0-c)*(+1.0d0)*(1.0d0+s)*(1.0d0-s)/4.0d0

        !d2φ/dη2
        !Corners
        d2fde2(1)  = -(1.0d0-c)*( (-1.0d0)*(+1.0d0)+(-1.0d0)*(+1.0d0) )*(1.0d0-s)/8.0d0
        d2fde2(2)  = -(1.0d0+c)*( (-1.0d0)*(+1.0d0)+(-1.0d0)*(+1.0d0) )*(1.0d0-s)/8.0d0
        d2fde2(3)  = -(1.0d0+c)*( (+1.0d0)*(-1.0d0)+(+1.0d0)*(-1.0d0) )*(1.0d0-s)/8.0d0
        d2fde2(4)  = -(1.0d0-c)*( (+1.0d0)*(-1.0d0)+(+1.0d0)*(-1.0d0) )*(1.0d0-s)/8.0d0
        d2fde2(5)  = -(1.0d0-c)*( (-1.0d0)*(+1.0d0)+(-1.0d0)*(+1.0d0) )*(1.0d0+s)/8.0d0
        d2fde2(6)  = -(1.0d0+c)*( (-1.0d0)*(+1.0d0)+(-1.0d0)*(+1.0d0) )*(1.0d0+s)/8.0d0
        d2fde2(7)  = -(1.0d0+c)*( (+1.0d0)*(-1.0d0)+(+1.0d0)*(-1.0d0) )*(1.0d0+s)/8.0d0
        d2fde2(8)  = -(1.0d0-c)*( (+1.0d0)*(-1.0d0)+(+1.0d0)*(-1.0d0) )*(1.0d0+s)/8.0d0
        !Mid-edge ξ = 0
        d2fde2(9)  = +0.0d0
        d2fde2(11) = +0.0d0
        d2fde2(13) = +0.0d0
        d2fde2(15) = +0.0d0
        !Mid-edge η = 0
        d2fde2(10) = +(1.0d0+c)*(-2.0d0)*(1.0d0-s)/4.0d0
        d2fde2(12) = +(1.0d0-c)*(-2.0d0)*(1.0d0-s)/4.0d0
        d2fde2(14) = +(1.0d0+c)*(-2.0d0)*(1.0d0+s)/4.0d0
        d2fde2(16) = +(1.0d0-c)*(-2.0d0)*(1.0d0+s)/4.0d0
        !Mid-edge ζ = 0
        d2fde2(17) = +0.0d0
        d2fde2(18) = +0.0d0
        d2fde2(19) = +0.0d0
        d2fde2(20) = +0.0d0

        !d2φ/dηζ
        !Corners
        d2fdes(1)  = -(1.0d0-c)*( ((-1.0d0)*(+1.0d0)+0.0d0)*(1.0d0-s) + ((-1.0d0)*(2.0d0+c+e+s)+(1.0d0-e)*(+1.0d0))*(-1.0d0) )/8.0d0
        d2fdes(2)  = -(1.0d0+c)*( ((-1.0d0)*(+1.0d0)+0.0d0)*(1.0d0-s) + ((-1.0d0)*(2.0d0-c+e+s)+(1.0d0-e)*(+1.0d0))*(-1.0d0) )/8.0d0
        d2fdes(3)  = -(1.0d0+c)*( ((+1.0d0)*(+1.0d0)+0.0d0)*(1.0d0-s) + ((+1.0d0)*(2.0d0-c-e+s)+(1.0d0+e)*(-1.0d0))*(-1.0d0) )/8.0d0
        d2fdes(4)  = -(1.0d0-c)*( ((+1.0d0)*(+1.0d0)+0.0d0)*(1.0d0-s) + ((+1.0d0)*(2.0d0+c-e+s)+(1.0d0+e)*(-1.0d0))*(-1.0d0) )/8.0d0
        d2fdes(5)  = -(1.0d0-c)*( ((-1.0d0)*(-1.0d0)+0.0d0)*(1.0d0+s) + ((-1.0d0)*(2.0d0+c+e-s)+(1.0d0-e)*(+1.0d0))*(+1.0d0) )/8.0d0
        d2fdes(6)  = -(1.0d0+c)*( ((-1.0d0)*(-1.0d0)+0.0d0)*(1.0d0+s) + ((-1.0d0)*(2.0d0-c+e-s)+(1.0d0-e)*(+1.0d0))*(+1.0d0) )/8.0d0
        d2fdes(7)  = -(1.0d0+c)*( ((+1.0d0)*(-1.0d0)+0.0d0)*(1.0d0+s) + ((+1.0d0)*(2.0d0-c-e-s)+(1.0d0+e)*(-1.0d0))*(+1.0d0) )/8.0d0
        d2fdes(8)  = -(1.0d0-c)*( ((+1.0d0)*(-1.0d0)+0.0d0)*(1.0d0+s) + ((+1.0d0)*(2.0d0+c-e-s)+(1.0d0+e)*(-1.0d0))*(+1.0d0) )/8.0d0

        !Mid-edge ξ = 0
        d2fdes(9)  = +(1.0d0+c)*(1.0d0-c)*(-1.0d0)*(-1.0d0)/4.0d0
        d2fdes(11) = +(1.0d0+c)*(1.0d0-c)*(+1.0d0)*(-1.0d0)/4.0d0
        d2fdes(13) = +(1.0d0+c)*(1.0d0-c)*(-1.0d0)*(+1.0d0)/4.0d0
        d2fdes(15) = +(1.0d0+c)*(1.0d0-c)*(+1.0d0)*(+1.0d0)/4.0d0
        !Mid-edge η = 0
        d2fdes(10) = +(1.0d0+c)*(-2.0d0*e)*(-1.0d0)/4.0d0
        d2fdes(12) = +(1.0d0-c)*(-2.0d0*e)*(-1.0d0)/4.0d0
        d2fdes(14) = +(1.0d0+c)*(-2.0d0*e)*(+1.0d0)/4.0d0
        d2fdes(16) = +(1.0d0-c)*(-2.0d0*e)*(+1.0d0)/4.0d0
        !Mid-edge ζ = 0
        d2fdes(17) = +(1.0d0-c)*(-1.0d0)*(-2.0d0*s)/4.0d0
        d2fdes(18) = +(1.0d0+c)*(-1.0d0)*(-2.0d0*s)/4.0d0
        d2fdes(19) = +(1.0d0+c)*(+1.0d0)*(-2.0d0*s)/4.0d0
        d2fdes(20) = +(1.0d0-c)*(+1.0d0)*(-2.0d0*s)/4.0d0
        
        !dφ/dζ
        !Corners
        dfds(1)  = -(1.0d0-c)*(1.0d0-e)*( (-1.0d0)*(2.0d0+c+e+s)+(1.0d0-s)*(+1.0d0) )/8.0d0
        dfds(2)  = -(1.0d0+c)*(1.0d0-e)*( (-1.0d0)*(2.0d0-c+e+s)+(1.0d0-s)*(+1.0d0) )/8.0d0
        dfds(3)  = -(1.0d0+c)*(1.0d0+e)*( (-1.0d0)*(2.0d0-c-e+s)+(1.0d0-s)*(+1.0d0) )/8.0d0
        dfds(4)  = -(1.0d0-c)*(1.0d0+e)*( (-1.0d0)*(2.0d0+c-e+s)+(1.0d0-s)*(+1.0d0) )/8.0d0
        dfds(5)  = -(1.0d0-c)*(1.0d0-e)*( (+1.0d0)*(2.0d0+c+e-s)+(1.0d0+s)*(-1.0d0) )/8.0d0
        dfds(6)  = -(1.0d0+c)*(1.0d0-e)*( (+1.0d0)*(2.0d0-c+e-s)+(1.0d0+s)*(-1.0d0) )/8.0d0
        dfds(7)  = -(1.0d0+c)*(1.0d0+e)*( (+1.0d0)*(2.0d0-c-e-s)+(1.0d0+s)*(-1.0d0) )/8.0d0
        dfds(8)  = -(1.0d0-c)*(1.0d0+e)*( (+1.0d0)*(2.0d0+c-e-s)+(1.0d0+s)*(-1.0d0) )/8.0d0
        !Mid-edge ξ = 0
        dfds(9)  = +(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*(-1.0d0)/4.0d0
        dfds(11) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(-1.0d0)/4.0d0
        dfds(13) = +(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*(+1.0d0)/4.0d0
        dfds(15) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(+1.0d0)/4.0d0
        !Mid-edge η = 0
        dfds(10) = +(1.0d0+c)*(1.0d0+e)*(1.0d0-e)*(-1.0d0)/4.0d0
        dfds(12) = +(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(-1.0d0)/4.0d0
        dfds(14) = +(1.0d0+c)*(1.0d0+e)*(1.0d0-e)*(+1.0d0)/4.0d0
        dfds(16) = +(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(+1.0d0)/4.0d0
        !Mid-edge ζ = 0
        dfds(17) = +(1.0d0-c)*(1.0d0-e)*(-2.0d0*s)/4.0d0
        dfds(18) = +(1.0d0+c)*(1.0d0-e)*(-2.0d0*s)/4.0d0
        dfds(19) = +(1.0d0+c)*(1.0d0+e)*(-2.0d0*s)/4.0d0
        dfds(20) = +(1.0d0-c)*(1.0d0+e)*(-2.0d0*s)/4.0d0

        !d2φ/dζ2
        !Corners
        d2fds2(1)  = -(1.0d0-c)*(1.0d0-e)*( (-1.0d0)*(+1.0d0)+(-1.0d0)*(+1.0d0) )/8.0d0
        d2fds2(2)  = -(1.0d0+c)*(1.0d0-e)*( (-1.0d0)*(+1.0d0)+(-1.0d0)*(+1.0d0) )/8.0d0
        d2fds2(3)  = -(1.0d0+c)*(1.0d0+e)*( (-1.0d0)*(+1.0d0)+(-1.0d0)*(+1.0d0) )/8.0d0
        d2fds2(4)  = -(1.0d0-c)*(1.0d0+e)*( (-1.0d0)*(+1.0d0)+(-1.0d0)*(+1.0d0) )/8.0d0
        d2fds2(5)  = -(1.0d0-c)*(1.0d0-e)*( (+1.0d0)*(-1.0d0)+(+1.0d0)*(-1.0d0) )/8.0d0
        d2fds2(6)  = -(1.0d0+c)*(1.0d0-e)*( (+1.0d0)*(-1.0d0)+(+1.0d0)*(-1.0d0) )/8.0d0
        d2fds2(7)  = -(1.0d0+c)*(1.0d0+e)*( (+1.0d0)*(-1.0d0)+(+1.0d0)*(-1.0d0) )/8.0d0
        d2fds2(8)  = -(1.0d0-c)*(1.0d0+e)*( (+1.0d0)*(-1.0d0)+(+1.0d0)*(-1.0d0) )/8.0d0
        !Mid-edge ξ = 0
        d2fds2(9)  = +0.0d0
        d2fds2(11) = +0.0d0
        d2fds2(13) = +0.0d0
        d2fds2(15) = +0.0d0
        !Mid-edge η = 0
        d2fds2(10) = +0.0d0
        d2fds2(12) = +0.0d0
        d2fds2(14) = +0.0d0
        d2fds2(16) = +0.0d0
        !Mid-edge ζ = 0
        d2fds2(17) = +(1.0d0-c)*(1.0d0-e)*(-2.0d0)/4.0d0
        d2fds2(18) = +(1.0d0+c)*(1.0d0-e)*(-2.0d0)/4.0d0
        d2fds2(19) = +(1.0d0+c)*(1.0d0+e)*(-2.0d0)/4.0d0
        d2fds2(20) = +(1.0d0-c)*(1.0d0+e)*(-2.0d0)/4.0d0

    end subroutine setSerendipityHexahedron

    !----------------------------------------------------------------------- 

    subroutine setLinearHexahedron_Bfn(xi, bfn)

        implicit none

        PetscReal, dimension(:), intent(in)  :: xi
        PetscReal, dimension(:), intent(out) :: bfn

        PetscReal :: c, e, s

        c = xi(1) ; e = xi(2) ; s = xi(3)

        !Basis function, φ
        bfn(1)  = (1.0d0-c)*(1.0d0-e)*(1.0d0-s)/8.0d0
        bfn(2)  = (1.0d0+c)*(1.0d0-e)*(1.0d0-s)/8.0d0
        bfn(3)  = (1.0d0+c)*(1.0d0+e)*(1.0d0-s)/8.0d0
        bfn(4)  = (1.0d0-c)*(1.0d0+e)*(1.0d0-s)/8.0d0
        bfn(5)  = (1.0d0-c)*(1.0d0-e)*(1.0d0+s)/8.0d0
        bfn(6)  = (1.0d0+c)*(1.0d0-e)*(1.0d0+s)/8.0d0
        bfn(7)  = (1.0d0+c)*(1.0d0+e)*(1.0d0+s)/8.0d0
        bfn(8)  = (1.0d0-c)*(1.0d0+e)*(1.0d0+s)/8.0d0

    end subroutine setLinearHexahedron_Bfn

    !----------------------------------------------------------------------- 

    subroutine setQuadraticHexahedron_Bfn(xi, bfn)

        implicit none

        PetscReal, dimension(:), intent(in)  :: xi
        PetscReal, dimension(:), intent(out) :: bfn

        PetscReal :: c, e, s

        c = xi(1) ; e = xi(2) ; s = xi(3)

        !Corners
        bfn(1)  = -(1.0d0-c)*c*(1.0d0-e)*e*(1.0d0-s)*s/8.0d0
        bfn(2)  = +(1.0d0+c)*c*(1.0d0-e)*e*(1.0d0-s)*s/8.0d0
        bfn(3)  = -(1.0d0+c)*c*(1.0d0+e)*e*(1.0d0-s)*s/8.0d0
        bfn(4)  = +(1.0d0-c)*c*(1.0d0+e)*e*(1.0d0-s)*s/8.0d0
        bfn(5)  = +(1.0d0-c)*c*(1.0d0-e)*e*(1.0d0+s)*s/8.0d0
        bfn(6)  = -(1.0d0+c)*c*(1.0d0-e)*e*(1.0d0+s)*s/8.0d0
        bfn(7)  = +(1.0d0+c)*c*(1.0d0+e)*e*(1.0d0+s)*s/8.0d0
        bfn(8)  = -(1.0d0-c)*c*(1.0d0+e)*e*(1.0d0+s)*s/8.0d0
        !Mid-edge ξ = 0
        bfn(9)  = +(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e*(1.0d0-s)*s/4.0d0
        bfn(11) = -(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e*(1.0d0-s)*s/4.0d0
        bfn(13) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e*(1.0d0+s)*s/4.0d0
        bfn(15) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e*(1.0d0+s)*s/4.0d0
        !Mid-edge η = 0
        bfn(10) = -(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)*s/4.0d0
        bfn(12) = +(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)*s/4.0d0
        bfn(14) = +(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*s/4.0d0
        bfn(16) = -(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*s/4.0d0
        !Mid-edge ζ = 0
        bfn(17) = +(1.0d0-c)*c*(1.0d0-e)*e*(1.0d0+s)*(1.0d0-s)/4.0d0
        bfn(18) = -(1.0d0+c)*c*(1.0d0-e)*e*(1.0d0+s)*(1.0d0-s)/4.0d0
        bfn(19) = +(1.0d0+c)*c*(1.0d0+e)*e*(1.0d0+s)*(1.0d0-s)/4.0d0
        bfn(20) = -(1.0d0-c)*c*(1.0d0+e)*e*(1.0d0+s)*(1.0d0-s)/4.0d0
        !Mid-face ξ = +-1
        bfn(22) = +(1.0d0+c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/2.0d0
        bfn(24) = -(1.0d0-c)*c*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/2.0d0
        !Mid-face η = +-1
        bfn(21) = -(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*e*(1.0d0+s)*(1.0d0-s)/2.0d0
        bfn(23) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*e*(1.0d0+s)*(1.0d0-s)/2.0d0
        !Mid-face ζ = +-1
        bfn(25) = -(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)*s/2.0d0
        bfn(26) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*s/2.0d0
        !Central
        bfn(27) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)

    end subroutine setQuadraticHexahedron_Bfn

    !----------------------------------------------------------------------- 

    subroutine setSerendipityHexahedron_Bfn(xi, bfn)

        implicit none

        PetscReal, dimension(:), intent(in)  :: xi
        PetscReal, dimension(:), intent(out) :: bfn

        PetscReal :: c, e, s

        c = xi(1) ; e = xi(2) ; s = xi(3)

        !Corners
        bfn(1)  = -(1.0d0-c)*(1.0d0-e)*(1.0d0-s)*(2.0d0+c+e+s)/8.0d0
        bfn(2)  = -(1.0d0+c)*(1.0d0-e)*(1.0d0-s)*(2.0d0-c+e+s)/8.0d0
        bfn(3)  = -(1.0d0+c)*(1.0d0+e)*(1.0d0-s)*(2.0d0-c-e+s)/8.0d0
        bfn(4)  = -(1.0d0-c)*(1.0d0+e)*(1.0d0-s)*(2.0d0+c-e+s)/8.0d0
        bfn(5)  = -(1.0d0-c)*(1.0d0-e)*(1.0d0+s)*(2.0d0+c+e-s)/8.0d0
        bfn(6)  = -(1.0d0+c)*(1.0d0-e)*(1.0d0+s)*(2.0d0-c+e-s)/8.0d0
        bfn(7)  = -(1.0d0+c)*(1.0d0+e)*(1.0d0+s)*(2.0d0-c-e-s)/8.0d0
        bfn(8)  = -(1.0d0-c)*(1.0d0+e)*(1.0d0+s)*(2.0d0+c-e-s)/8.0d0
        !Mid-edge ξ = 0
        bfn(9)  = +(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*(1.0d0-s)/4.0d0
        bfn(11) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0-s)/4.0d0
        bfn(13) = +(1.0d0+c)*(1.0d0-c)*(1.0d0-e)*(1.0d0+s)/4.0d0
        bfn(15) = +(1.0d0+c)*(1.0d0-c)*(1.0d0+e)*(1.0d0+s)/4.0d0
        !Mid-edge η = 0
        bfn(10) = +(1.0d0+c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)/4.0d0
        bfn(12) = +(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(1.0d0-s)/4.0d0
        bfn(14) = +(1.0d0+c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)/4.0d0
        bfn(16) = +(1.0d0-c)*(1.0d0+e)*(1.0d0-e)*(1.0d0+s)/4.0d0
        !Mid-edge ζ = 0
        bfn(17) = +(1.0d0-c)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        bfn(18) = +(1.0d0+c)*(1.0d0-e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        bfn(19) = +(1.0d0+c)*(1.0d0+e)*(1.0d0+s)*(1.0d0-s)/4.0d0
        bfn(20) = +(1.0d0-c)*(1.0d0+e)*(1.0d0+s)*(1.0d0-s)/4.0d0

    end subroutine setSerendipityHexahedron_Bfn

    !-----------------------------------------------------------------------

    subroutine setHexahedron_dBfn(f_sub, xi, dfdc, dfde, dfds)

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

    end subroutine setHexahedron_dBfn

    !-----------------------------------------------------------------------

    subroutine setHexahedron_d2Bfn(f_sub, xi, d2fdc2, d2fde2, d2fds2, d2fdce, d2fdcs, d2fdes)

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

    end subroutine setHexahedron_d2Bfn

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

end module GaussHexahedron_mod
