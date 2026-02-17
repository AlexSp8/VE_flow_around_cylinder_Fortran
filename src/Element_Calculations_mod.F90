
#include <petsc/finclude/petscksp.h>

module ElementCalculations_mod

    use petscksp
    use MPIParameters_mod, only: Rank
    use PhysicalParameters_mod, only: ProblemParameters
    use FEMParameters_mod, only: FEParameters
    use MeshParameters_mod, only: MeshParameters
    use GaussParameters_mod, only: GaussIntegration
    use SolutionVariables_mod, only: SolutionArraysType
    use ElementVariables_mod, only: GaussPointQuantities, NodeArrays

    implicit none

    interface
        subroutine basis(idir, Nbf, Elem, GsPt)
            use ElementVariables_mod, only: GaussPointQuantities, NodeArrays
            implicit none
            PetscInt, dimension(3), intent(in) :: idir
            PetscInt, intent(in) :: Nbf
            type(NodeArrays), intent(inout) :: Elem
            type(GaussPointQuantities), intent(inout) :: GsPt
        end subroutine basis
    end interface

    procedure(basis), pointer :: basisInflow_p, basisMain_p

    contains

    subroutine allocateLocalArrays(Problem, Problem_Stability, FE, Elem)

        use Tools_mod, only: allocateArray

        implicit none

        type(ProblemParameters), intent(in) :: Problem, Problem_Stability
        type(FEParameters), intent(in) :: FE
        type(NodeArrays), intent(inout) :: Elem

        PetscInt :: Nbf, Neq, Nex, Neq_proj, Neq_S, Nex_S

        Nbf = FE%Nbf
        Neq = Problem%Neq
        Nex = Problem%Nex
        Neq_proj = Problem%Neq_proj

        call allocateArray(Elem%TEMP_TL, Nbf, Neq)
        call allocateArray(Elem%TEMP_TLo, Nbf, Neq)
        call allocateArray(Elem%TEMP_TLb, Nbf, Neq)
        call allocateArray(Elem%TEMP_EU, Nex)
        ! call allocateArray(Elem%Jump, Neq)
        call allocateArray(Elem%dFdXi, 3, Nbf)
        call allocateArray(Elem%d2FdXi2, 3, Nbf)
        call allocateArray(Elem%Xi_loc, Nbf, 3)

        call allocateArray(Elem%TEMP_TL_proj, Nbf, Neq_proj)

        Neq_S = Problem_Stability%Neq
        Nex_S = Problem_Stability%Nex
        
        call allocateArray(Elem%TEMP_TL_b, Nbf, Neq_S)
        call allocateArray(Elem%TEMP_TL_d, Nbf, Neq_S)
        call allocateArray(Elem%TEMP_EU_b, Nex_S)
        call allocateArray(Elem%TEMP_EU_d, Nex_S)

    end subroutine allocateLocalArrays

    !---------------------------------------------------------------------

    subroutine setBasisPointers(Ndim_main, Ndim_inflow)

        implicit none

        PetscInt, intent(in) :: Ndim_main, Ndim_inflow

        select case (Ndim_main)
        case (1)
            basisMain_p => basis1D
        case (2)
            basisMain_p => basis2D
        case (3)
            basisMain_p => basis3D
        case default
            write(*,'(a)') 'Wrong Ndim_main in setBasisPointers!'
            stop
        end select

        select case (Ndim_inflow)
        case (1)
            basisInflow_p => basis1D
        case (2)
            basisInflow_p => basis2D
        case (3)
            basisInflow_p => basis3D
        case default
            write(*,'(a)') 'Wrong Ndim_inflow in setBasisPointers!'
            stop
        end select

    end subroutine setBasisPointers

    !---------------------------------------------------------------------

    subroutine copyToLocal(Nbf, Xi_Mesh_Rank, Sol, iel_rank, Elem)

        use Tools_mod, only: getRankNode

        implicit none

        PetscInt, intent(in) :: Nbf
        PetscReal, dimension(:,:), intent(in) :: Xi_Mesh_Rank
        type(SolutionArraysType), intent(in) :: Sol
        PetscInt, intent(in) :: iel_rank
        type(NodeArrays), intent(inout) :: Elem

        PetscInt :: rnod, inod

        do inod = 1, Nbf

            rnod = getRankNode(iel_rank,inod,Nbf)

            Elem%TEMP_TL(inod,:) = Sol%TL_Rank(rnod,:)
            Elem%TEMP_TLo(inod,:) = Sol%TLo_Rank(rnod,:)
            Elem%TEMP_TLb(inod,:) = Sol%TLb_Rank(rnod,:)

            Elem%Xi_loc(inod,:) = Xi_Mesh_Rank(rnod,:)

            Elem%TEMP_TL_proj(inod,:) = Sol%TL_proj_Rank(rnod,:)

            Elem%TEMP_TL_b(inod,:) = Sol%TL_b_Rank(rnod,:)
            Elem%TEMP_TL_d(inod,:) = Sol%TL_d_Rank(rnod,:)

        end do

        Elem%TEMP_EU(:) = Sol%EU(:)

        Elem%TEMP_EU_b(:) = Sol%EU_b(:)
        Elem%TEMP_EU_d(:) = Sol%EU_d(:)

    end subroutine copyToLocal

    !---------------------------------------------------------------------

    subroutine copyToLocalGlobal(Nbf, Mesh, Sol, iel, Elem)

        implicit none

        PetscInt, intent(in) :: Nbf
        type(MeshParameters), intent(in) :: Mesh
        type(SolutionArraysType), intent(in) :: Sol
        PetscInt, intent(in) :: iel
        type(NodeArrays), intent(inout) :: Elem

        PetscInt :: gnod, inod

        do inod = 1, Nbf
            gnod = Mesh%Connectivity(iel,inod)
            Elem%TEMP_TL(inod,:) = Sol%TL(gnod,:)
            Elem%TEMP_TLo(inod,:) = Sol%TLo(gnod,:)
            Elem%TEMP_TLb(inod,:) = Sol%TLb(gnod,:)
            Elem%Xi_loc(inod,:) = Mesh%Xi_Mesh(gnod,:)
        end do

        Elem%TEMP_EU(:) = Sol%EU(:)

    end subroutine copyToLocalGlobal

    !---------------------------------------------------------------------

    subroutine setElementBasisFunctionsEdge(ig, ied, GaussInt, Elem)

        implicit none

        integer, intent(in) :: ig, ied
        type(GaussIntegration), intent(in), target :: GaussInt
        type(NodeArrays), intent(inout) :: Elem

        Elem%bfn_p => GaussInt%Bfn_E (:,ig,ied)
        Elem%dfdc_p => GaussInt%dfdc_E(:,ig,ied)
        Elem%dfde_p => GaussInt%dfde_E(:,ig,ied)
        Elem%dfds_p => GaussInt%dfds_E(:,ig,ied)
        Elem%d2fdc2_p => GaussInt%d2fdc2_E(:,ig,ied)
        Elem%d2fde2_p => GaussInt%d2fde2_E(:,ig,ied)
        Elem%d2fds2_p => GaussInt%d2fds2_E(:,ig,ied)
        Elem%d2fdce_p => GaussInt%d2fdce_E(:,ig,ied)
        Elem%d2fdcs_p => GaussInt%d2fdcs_E(:,ig,ied)
        Elem%d2fdes_p => GaussInt%d2fdes_E(:,ig,ied)

        ! Elem%dfdci_p => GaussInt%dfdci_E(:,:,ig,ied)

    end subroutine setElementBasisFunctionsEdge

    !---------------------------------------------------------------------

    subroutine setElementBasisFunctionsFace(ig, ifc, GaussInt, Elem)

        implicit none

        integer, intent(in) :: ig, ifc
        type(GaussIntegration), intent(in), target :: GaussInt
        type(NodeArrays), intent(inout) :: Elem

        Elem%bfn_p => GaussInt%Bfn_F (:,ig,ifc)
        Elem%dfdc_p => GaussInt%dfdc_F(:,ig,ifc)
        Elem%dfde_p => GaussInt%dfde_F(:,ig,ifc)
        Elem%dfds_p => GaussInt%dfds_F(:,ig,ifc)
        Elem%d2fdc2_p => GaussInt%d2fdc2_F(:,ig,ifc)
        Elem%d2fde2_p => GaussInt%d2fde2_F(:,ig,ifc)
        Elem%d2fds2_p => GaussInt%d2fds2_F(:,ig,ifc)
        Elem%d2fdce_p => GaussInt%d2fdce_F(:,ig,ifc)
        Elem%d2fdcs_p => GaussInt%d2fdcs_F(:,ig,ifc)
        Elem%d2fdes_p => GaussInt%d2fdes_F(:,ig,ifc)

        ! Elem%dfdci_p => GaussInt%dfdci_F(:,:,ig,ifc)

    end subroutine setElementBasisFunctionsFace

    !---------------------------------------------------------------------

    subroutine setElementBasisFunctions(ig, GaussInt, Elem)

        implicit none

        integer, intent(in) :: ig
        type(GaussIntegration), intent(in), target :: GaussInt
        type(NodeArrays), intent(inout) :: Elem

        Elem%bfn_p => GaussInt%Bfn(:,ig)
        Elem%dfdc_p => GaussInt%dfdc(:,ig)
        Elem%dfde_p => GaussInt%dfde(:,ig)
        Elem%dfds_p => GaussInt%dfds(:,ig)
        Elem%d2fdc2_p => GaussInt%d2fdc2(:,ig)
        Elem%d2fde2_p => GaussInt%d2fde2(:,ig)
        Elem%d2fds2_p => GaussInt%d2fds2(:,ig)
        Elem%d2fdce_p => GaussInt%d2fdce(:,ig)
        Elem%d2fdcs_p => GaussInt%d2fdcs(:,ig)
        Elem%d2fdes_p => GaussInt%d2fdes(:,ig)

        ! Elem%dfdci_p => GaussInt%dfdci(:,:,ig)

    end subroutine setElementBasisFunctions

    !---------------------------------------------------------------------

    subroutine setFlowQuantities(Problem, Nbf, Elem, GsPt)

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nbf
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(inout) :: GsPt

        call setFlowUnknowns(Problem, Elem, GsPt)

        call strongFormEquations(Problem, Elem%TEMP_EU, GsPt)

        call setStabilizationTerms(Nbf, Elem%dFdXi, GsPt)

    end subroutine setFlowQuantities

    !---------------------------------------------------------------------

    subroutine basis1D(idir, Nbf, Elem, GsPt)

        implicit none

        PetscInt, dimension(1), intent(in) :: idir
        PetscInt, intent(in) :: Nbf
        type(NodeArrays), intent(inout) :: Elem
        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt :: inod, j
        PetscReal :: dcdx

        GsPt%Xi(:) = 0.0d0
        GsPt%dXidci(:,:) = 0.0d0
        do inod = 1, Nbf
            do j = 1, 3
                GsPt%Xi(j) = GsPt%Xi(j) + Elem%bfn_p(inod)*Elem%Xi_loc(inod,j)
                GsPt%dXidci(1,j) = GsPt%dXidci(1,j) + Elem%dfdc_p(inod)*Elem%Xi_loc(inod,j)
                ! GsPt%dXidci(:,j) = GsPt%dXidci(:,j) + Elem%dfdci_p(:,inod)*Elem%Xi_loc(inod,j)
            end do
        end do

        associate( dir1 => idir(1), dxdc => GsPt%dXidci(1,dir1), &
                dfdc => Elem%dfdc_p, d2fdc2 => Elem%d2fdc2_p )

            GsPt%detJ = dxdc

            dcdx = 1.0d0/dxdc

            Elem%dFdXi(:,:) = 0.0d0 ; Elem%d2FdXi2(:,:) = 0.0d0
            do inod = 1, Nbf
                Elem%dFdXi(dir1,inod) = dcdx*dfdc(inod)
                ! Elem%dFdXi(dir1,inod) = dcdx*Elem%dfdci_p(1,inod)
                Elem%d2FdXi2(dir1,inod) = (dcdx**2)*d2fdc2(inod)
            end do

        end associate

    end subroutine basis1D

    !---------------------------------------------------------------------

    subroutine basis2D(idir, Nbf, Elem, GsPt)

        implicit none

        PetscInt, dimension(2), intent(in) :: idir
        PetscInt, intent(in) :: Nbf
        type(NodeArrays), intent(inout) :: Elem
        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt :: inod, j
        PetscReal :: dcdx, dedx, dcdy, dedy

        GsPt%Xi(:) = 0.0d0
        GsPt%dXidci(:,:) = 0.0d0
        do inod = 1, Nbf
            do j = 1, 3
                GsPt%Xi(j) = GsPt%Xi(j) + Elem%bfn_p(inod)*Elem%Xi_loc(inod,j)
                GsPt%dXidci(1,j) = GsPt%dXidci(1,j) + Elem%dfdc_p(inod)*Elem%Xi_loc(inod,j)
                GsPt%dXidci(2,j) = GsPt%dXidci(2,j) + Elem%dfde_p(inod)*Elem%Xi_loc(inod,j)
                ! GsPt%dXidci(:,j) = GsPt%dXidci(:,j) + Elem%dfdci_p(:,inod)*Elem%Xi_loc(inod,j)
            end do
        end do

        associate( dir1 => idir(1), dir2 => idir(2), &
                dxdc => GsPt%dXidci(1,dir1), dydc => GsPt%dXidci(1,dir2), &
                dxde => GsPt%dXidci(2,dir1), dyde => GsPt%dXidci(2,dir2), &
                dfdc => Elem%dfdc_p, dfde => Elem%dfde_p, &
                d2fdc2 => Elem%d2fdc2_p, d2fdce => Elem%d2fdce_p, d2fde2 => Elem%d2fde2_p )

            GsPt%detJ = dxdc*dyde-dxde*dydc

            dcdx =  dyde/(GsPt%detJ)
            dedx = -dydc/(GsPt%detJ)
            dcdy = -dxde/(GsPt%detJ)
            dedy =  dxdc/(GsPt%detJ)

            Elem%dFdXi(:,:) = 0.0d0 ; Elem%d2FdXi2(:,:) = 0.0d0
            do inod = 1, Nbf

                Elem%dFdXi(dir1,inod) = (dcdx*dfdc(inod)+dedx*dfde(inod))
                Elem%dFdXi(dir2,inod) = (dcdy*dfdc(inod)+dedy*dfde(inod))
                ! Elem%dFdXi(dir1,inod) = (dcdx*Elem%dfdci_p(1,inod)+dedx*Elem%dfdci_p(2,inod))
                ! Elem%dFdXi(dir2,inod) = (dcdy*Elem%dfdci_p(1,inod)+dedy*Elem%dfdci_p(2,inod))

                Elem%d2FdXi2(dir1,inod) = (dcdx**2)*d2fdc2(inod) + (dedx**2)*d2fde2(inod) &
                    + 2.0d0*dcdx*dedx*d2fdce(inod)
                Elem%d2FdXi2(dir2,inod) = (dcdy**2)*d2fdc2(inod) + (dedy**2)*d2fde2(inod) &
                    + 2.0d0*dcdy*dedy*d2fdce(inod)

            end do

        end associate

    end subroutine basis2D

    !---------------------------------------------------------------------

    subroutine basis3D(idir, Nbf, Elem, GsPt)

        implicit none

        PetscInt, dimension(3), intent(in) :: idir
        PetscInt, intent(in) :: Nbf
        type(NodeArrays), intent(inout) :: Elem
        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt :: inod, j
        PetscReal :: dcdx, dcdy, dcdz, dedx, dedy, dedz, dsdx, dsdy, dsdz

        GsPt%Xi(:) = 0.0d0
        GsPt%dXidci(:,:) = 0.0d0
        do inod = 1, Nbf
            do j = 1, 3
                GsPt%Xi(j) = GsPt%Xi(j) + Elem%bfn_p(inod)*Elem%Xi_loc(inod,j)
                GsPt%dXidci(1,j) = GsPt%dXidci(1,j) + Elem%dfdc_p(inod)*Elem%Xi_loc(inod,j)
                GsPt%dXidci(2,j) = GsPt%dXidci(2,j) + Elem%dfde_p(inod)*Elem%Xi_loc(inod,j)
                GsPt%dXidci(3,j) = GsPt%dXidci(3,j) + Elem%dfds_p(inod)*Elem%Xi_loc(inod,j)
                ! GsPt%dXidci(:,j) = GsPt%dXidci(:,j) + Elem%dfdci_p(:,inod)*Elem%Xi_loc(inod,j)
            end do
        end do

        associate( dxdc => GsPt%dXidci(1,1), dydc => GsPt%dXidci(1,2), dzdc => GsPt%dXidci(1,3), &
                    dxde => GsPt%dXidci(2,1), dyde => GsPt%dXidci(2,2), dzde => GsPt%dXidci(2,3), &
                    dxds => GsPt%dXidci(3,1), dyds => GsPt%dXidci(3,2), dzds => GsPt%dXidci(3,3), &
                    dir1 => idir(1), dir2 => idir(2), dir3 => idir(3), &
                    dfdc => Elem%dfdc_p, dfde => Elem%dfde_p, dfds => Elem%dfds_p, &
                    d2fdc2 => Elem%d2fdc2_p, d2fde2 => Elem%d2fde2_p, d2fds2 => Elem%d2fds2_p, &
                    d2fdce => Elem%d2fdce_p, d2fdcs => Elem%d2fdcs_p, d2fdes => Elem%d2fdes_p )

            GsPt%detJ = dxdc*(dyde*dzds-dyds*dzde) - &
                       dxde*(dydc*dzds-dyds*dzdc) + &
                       dxds*(dydc*dzde-dyde*dzdc)

            dcdx =  (dyde*dzds-dyds*dzde)/(GsPt%detJ)
            dcdy = -(dxde*dzds-dzde*dxds)/(GsPt%detJ)
            dcdz =  (dxde*dyds-dyde*dxds)/(GsPt%detJ)

            dedx = -(dydc*dzds-dzdc*dyds)/(GsPt%detJ)
            dedy =  (dxdc*dzds-dzdc*dxds)/(GsPt%detJ)
            dedz = -(dxdc*dyds-dydc*dxds)/(GsPt%detJ)

            dsdx =  (dydc*dzde-dzdc*dyde)/(GsPt%detJ)
            dsdy = -(dxdc*dzde-dzdc*dxde)/(GsPt%detJ)
            dsdz =  (dxdc*dyde-dydc*dxde)/(GsPt%detJ)

            Elem%dFdXi(:,:) = 0.0d0 ; Elem%d2FdXi2(:,:) = 0.0d0
            do inod = 1, Nbf

                Elem%dFdXi(dir1,inod) = dfdc(inod)*dcdx + dfde(inod)*dedx + dfds(inod)*dsdx
                Elem%dFdXi(dir2,inod) = dfdc(inod)*dcdy + dfde(inod)*dedy + dfds(inod)*dsdy
                Elem%dFdXi(dir3,inod) = dfdc(inod)*dcdz + dfde(inod)*dedz + dfds(inod)*dsdz

                ! Elem%dFdXi(dir1,inod) = Elem%dfdci_p(1,inod)*dcdx &
                !     + Elem%dfdci_p(2,inod)*dedx + Elem%dfdci_p(3,inod)*dsdx

                ! Elem%dFdXi(dir2,inod) = Elem%dfdci_p(1,inod)*dcdy &
                !     + Elem%dfdci_p(2,inod)*dedy + Elem%dfdci_p(3,inod)*dsdy

                ! Elem%dFdXi(dir3,inod) = Elem%dfdci_p(1,inod)*dcdz &
                !     + Elem%dfdci_p(2,inod)*dedz + Elem%dfdci_p(3,inod)*dsdz

                Elem%d2FdXi2(dir1,inod) = (dcdx**2)*d2fdc2(inod) + (dedx**2)*d2fde2(inod) &
                    + (dsdx**2)*d2fds2(inod) + 2.0d0*dcdx*dedx*d2fdce(inod) &
                    + 2.0d0*dcdx*dsdx*d2fdcs(inod) + 2.0d0*dedx*dsdx*d2fdes(inod)

                Elem%d2FdXi2(dir2,inod) = (dcdy**2)*d2fdc2(inod) + (dedy**2)*d2fde2(inod) &
                    + (dsdy**2)*d2fds2(inod) + 2.0d0*dcdy*dedy*d2fdce(inod) &
                    + 2.0d0*dcdy*dsdy*d2fdcs(inod) + 2.0d0*dedy*dsdy*d2fdes(inod)
                    
                Elem%d2FdXi2(dir3,inod) = (dcdz**2)*d2fdc2(inod) + (dedz**2)*d2fde2(inod) &
                    + (dsdz**2)*d2fds2(inod) + 2.0d0*dcdz*dedz*d2fdce(inod) &
                    + 2.0d0*dcdz*dsdz*d2fdcs(inod) + 2.0d0*dedz*dsdz*d2fdes(inod)

            end do

        end associate

    end subroutine basis3D

    !---------------------------------------------------------------------

    subroutine setFlowUnknowns(Problem, Elem, GsPt)

        use Tools_mod, only: getStressComponent, getTensorComponent

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt :: i, j, k, idim, istart, Neq_f, Neq_s
        
        Neq_f = Problem%Neq_f

        GsPt%U_f(:) = 0.0d0
        GsPt%U_fo(:) = 0.0d0
        GsPt%U_fb(:) = 0.0d0
        GsPt%GU(:,:) = 0.0d0
        GsPt%G2U(:,:) = 0.0d0

        GsPt%P = 0.0d0
        GsPt%dPdXi(:) = 0.0d0
        
        do i = 1, Neq_f

            GsPt%U_f(i) = dot_product( Elem%TEMP_TL(:,i),Elem%bfn_p(:) )

            GsPt%U_fo(i) = dot_product( Elem%TEMP_TLo(:,i),Elem%bfn_p(:) )
            GsPt%U_fb(i) = dot_product( Elem%TEMP_TLb(:,i),Elem%bfn_p(:) )

            do j = 1, 3
                GsPt%GU(j,i)  = dot_product( Elem%TEMP_TL(:,i),Elem%dFdXi(j,:) )
                GsPt%G2U(j,i) = dot_product( Elem%TEMP_TL(:,i),Elem%d2FdXi2(j,:) )
            end do

            GsPt%dPdXi(i) = dot_product( Elem%TEMP_TL(:,Neq_f+1),Elem%dFdXi(i,:) )

        end do

        GsPt%P = dot_product( Elem%TEMP_TL(:,Neq_f+1),Elem%bfn_p(:) )

        GsPt%S(:,:) = 0.0d0
        GsPt%So(:,:) = 0.0d0
        GsPt%Sb(:,:) = 0.0d0
        GsPt%dSdXi(:,:,:) = 0.0d0

        Neq_s = Problem%Neq_s
        istart = size(Elem%TEMP_TL,2) - Neq_s + 1
        do i = istart, size(Elem%TEMP_TL,2)

            call getStressComponent(i,istart,j,k)

            GsPt%S(j,k) = dot_product( Elem%TEMP_TL(:,i),Elem%bfn_p(:) )
            GsPt%S(k,j) = GsPt%S(j,k)
            do idim = 1, 3
                GsPt%dSdXi(j,k,idim) = dot_product( Elem%TEMP_TL(:,i),Elem%dFdXi(idim,:) )
            end do
            GsPt%dSdXi(k,j,:) = GsPt%dSdXi(j,k,:)

            GsPt%So(j,k) = dot_product( Elem%TEMP_TLo(:,i),Elem%bfn_p(:) )
            GsPt%So(k,j) = GsPt%So(j,k)

            GsPt%Sb(j,k) = dot_product( Elem%TEMP_TLb(:,i),Elem%bfn_p(:) )
            GsPt%Sb(k,j) = GsPt%Sb(j,k)

        end do

        if (Problem%name == 'Inflow_1D' .or. Problem%name == 'Inflow_2D') then
            ! GsPt%U_f(2)  = 0.0d0  !Uy
            ! GsPt%U_f(3)  = 0.0d0  !Uz
            GsPt%GU(1,:) = 0.0d0  !d/dx
            GsPt%GU(:,2) = 0.0d0  !duy
            GsPt%GU(:,3) = 0.0d0  !duz
            GsPt%P = 0.0d0
            GsPt%dPdXi(:) = 0.0d0
            GsPt%dSdXi(:,:,:) = 0.0d0  !dSjk/dxi
        end if

        GsPt%GU_proj(:,:) = 0.0d0
        do i = 1, size(Elem%TEMP_TL_proj,2)
            call getTensorComponent(i,j,k)
            GsPt%GU_proj(j,k) = dot_product(Elem%TEMP_TL_proj(:,i),Elem%bfn_p(:))
        end do

    end subroutine setFlowUnknowns

    !---------------------------------------------------------------------

    subroutine strongFormEquations(Problem, TEMP_EU, GsPt)

        use ContinuationParameters_mod, only: Continuation_Method
        use PhysicalParameters_mod, only: Stress_Reform
        use ConstitutiveModels_mod, only: strongFormConEQ
        use ConstitutiveModelsSQRT_mod, only: strongFormConEQ_SQRT

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscReal, dimension(:), intent(in) :: TEMP_EU
        type(GaussPointQuantities), intent(inout) :: GsPt

        if (Continuation_Method == 'Arclength') then
            call setContinuationParameter(Problem, TEMP_EU)
        end if

        call setFlowKinematics(Problem%Ndim,GsPt)

        GsPt%dUfdt(:) = 0.0d0
        GsPt%dSdt(:,:) = 0.0d0
        if (Continuation_Method == 'Transient') then
            call setTimeDerivatives(GsPt)
        end if

        select case (Stress_Reform)
        case ('SQRT')
            call strongFormConEQ_SQRT(GsPt)
        case default
            call strongFormConEQ(GsPt)
        end select

        call strongFormMomentum(GsPt)

        call strongFormContinuity(GsPt)

    end subroutine strongFormEquations

    !---------------------------------------------------------------------

    subroutine setContinuationParameter(Problem, TEMP_EU)

        use ContinuationVariables_mod, only: Cvar1
        use SolutionVariables_mod, only: Sol_Main

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscReal, dimension(:), intent(in) :: TEMP_EU

        PetscReal :: ArN
        PetscInt :: Nex

        Nex = size(Sol_Main%EU)

        if (Problem%name == 'Inflow_1D' .or. Problem%name == 'Inflow_2D') then
            ArN = Sol_Main%EU(Nex)
        else
            ArN = TEMP_EU(Nex)
        end if

        ! WiN = WiN0*ArN
        Cvar1%p = (Cvar1%var0)*ArN

    end subroutine setContinuationParameter

    !---------------------------------------------------------------------

    subroutine setFlowKinematics(N,GsPt)

        use Tools_mod, only: traceTensor

        implicit none

        PetscInt, intent(in) :: N
        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscReal, dimension(3,3), parameter :: &
            I_tensor = reshape( [1,0,0,0,1,0,0,0,1], shape = shape(I_tensor) )

        GsPt%GUT(:,:) = transpose(GsPt%GU)

        GsPt%G_dot(:,:) = GsPt%GU(:,:) + GsPt%GUT(:,:)

        !Traceless
        GsPt%GU_tr(:,:) = GsPt%GU(:,:) !- traceTensor(GsPt%GU)*I_tensor(:,:)/3.0d0
        GsPt%GU_tr(1:N,1:N) = GsPt%GU(1:N,1:N) - traceTensor(GsPt%GU)*I_tensor(1:N,1:N)/N

        GsPt%GUT_tr(:,:) = transpose(GsPt%GU_tr)
        GsPt%G_dot_tr(:,:) = GsPt%GU_tr(:,:) + GsPt%GUT_tr(:,:)

        ! GsPt%Jump(:) = 0.0d0!Ux  - Uxo
        ! do inod = 1, Nbf
        !   rnod = getRankNode(iel_rank,inod,Nbf)
        !   GsPt%Jump(:) = GsPt%Jump(:) + (TL_Rank(rnod,:) - TLo_Rank(rnod,:)
        ! end do

    end subroutine setFlowKinematics

    !---------------------------------------------------------------------

    subroutine setTimeDerivatives(GsPt)

        use ContinuationVariables_mod, only: dt, Increment

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        if (Increment == 1) then
            GsPt%dUfdt(:) = (GsPt%U_f(:) - GsPt%U_fo(:))/dt
            GsPt%dSdt(:,:) = (GsPt%S(:,:) - GsPt%So(:,:))/dt
        else
            GsPt%dUfdt(:) = (3.0d0*GsPt%U_f(:) - 4.0d0*GsPt%U_fo(:) + GsPt%U_fb(:))/(2.0d0*dt)
            GsPt%dSdt(:,:) = (3.0d0*GsPt%S(:,:) - 4.0d0*GsPt%So(:,:) + GsPt%Sb(:,:))/(2.0d0*dt)
        end if

    end subroutine setTimeDerivatives

    !---------------------------------------------------------------------

    subroutine strongFormMomentum(GsPt)

        use PhysicalParameters_mod, only: ReN, BvN

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt :: i, j
        PetscReal, dimension(3) :: div_Tve, div_Ptot

        div_Tve(:) = 0.0d0
        do i = 1, 3
            do j = 1, 3
                div_Tve(i) = div_Tve(i) + GsPt%dTdXi(j,i,j)
            end do
        end do

        div_Ptot(:) = -GsPt%dPdXi(:) + div_Tve(:) &
                    + BvN*(GsPt%G2U(1,:)+GsPt%G2U(2,:)+GsPt%G2U(3,:))

        GsPt%Convection(:) = GsPt%dUfdt(:) + matmul(GsPt%U_f,GsPt%GU)

        GsPt%Momentum(:) = ReN*(GsPt%Convection(:)) - div_Ptot(:)

    end subroutine strongFormMomentum

    !---------------------------------------------------------------------

    subroutine strongFormContinuity(GsPt)

        use Tools_mod, only: traceTensor

        implicit none

        type(GaussPointQuantities), intent(inout) :: GsPt

        GsPt%Continuity = traceTensor(GsPt%GU)

    end subroutine strongFormContinuity

    !---------------------------------------------------------------------

    subroutine setStabilizationTerms(Nbf, dFdXi, GsPt)

        use Tools_mod, only: doubleDotProduct
        use PhysicalParameters_mod, only: BvN, ReN, WiN, Stress_Reform, DfN
        use ContinuationVariables_mod, only: dt

        implicit none

        PetscInt, intent(in) :: Nbf
        PetscReal, dimension(3,Nbf) :: dFdXi
        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt :: inod
        PetscReal :: Helem
        PetscScalar :: tsugn, Ha
        PetscScalar :: double_dot_T, double_dot_G, double_dot_S, double_dot_C

        Helem = 0.0d0
        tsugn = 0.0d0
        do inod = 1, Nbf
            ! Helem = Helem + sqrt(dFdXi(1,inod)**2 + dFdXi(2,inod)**2 + dFdXi(3,inod)**2) !2D Quad
            Helem = Helem + dFdXi(1,inod)**2 + dFdXi(2,inod)**2 + dFdXi(3,inod)**2 !3D
            tsugn = tsugn + abs(dot_product(GsPt%U_f,dFdXi(:,inod)))
        end do

        ! Helem = 1.0d0/Helem
        Helem = 1.0d0/sqrt(Helem)

        ! Helem = 0.0d0
        ! do ig = 1, GaussInt%NGauss_b
        !     call basis_p(Problem%idir, Nbf, Elem, GsPt)
        !     Helem = Helem + (GaussInt%Weights_b(ig))*abs(GsPt%detJ)
        ! end do
        ! Helem = sqrt(Helem)/4.0d0

        double_dot_T = doubleDotProduct(GsPt%T_ve,GsPt%T_ve)
        double_dot_G = doubleDotProduct(GsPt%G_dot_tr,GsPt%G_dot_tr)
        Ha = sqrt( ((1.0d0-BvN)/WiN)**2 + 0.5d0*double_dot_T )
        Ha = Ha/sqrt(1.0d0+0.5d0*double_dot_G)

        GsPt%tlsme = (ReN/dt)**2 + (ReN*tsugn)**2 + (BvN/(Helem**2))**2 + (Ha/(Helem**2))**2
        GsPt%tlsme = 1.0d0/sqrt(GsPt%tlsme)

        GsPt%tlsic = (Helem**2)/(GsPt%tlsme)

        double_dot_S = doubleDotProduct(GsPt%S,GsPt%S)

        GsPt%tdcs = 0.5d0*double_dot_S
        select case (Stress_Reform)
        case ('SQRT')
            double_dot_C = doubleDotProduct(GsPt%C_Tensor,GsPt%C_Tensor)
        case default
            GsPt%tdcs = (GsPt%tdcs) + ((1.0d0-BvN)/WiN)**2
            double_dot_C = double_dot_T
        end select

        GsPt%tdcs = (Helem**2)/sqrt(GsPt%tdcs)
        if (DfN > 1.0d-12) then
            GsPt%tdcs = 0.0d0
        end if

        GsPt%mag_con_EQ = doubleDotProduct(GsPt%Con_EQ,GsPt%Con_EQ)
        GsPt%mag_con_EQ = sqrt(0.5d0*GsPt%mag_con_EQ)
        
        GsPt%mag_con_EQ_time = doubleDotProduct(GsPt%Con_EQ_time,GsPt%Con_EQ_time)

        GsPt%tlsce = (WiN/dt)**2 + (WiN*tsugn)**2 +(WiN*sqrt(0.5d0*double_dot_G))**2 &
                    + (GsPt%f(2))**2 + (GsPt%f(3)*sqrt(0.5d0*double_dot_C))**2 &
                    + ((GsPt%tdcs)*(GsPt%mag_con_EQ)/(Helem**2))**2 &
                    + (DfN/(Helem**2))**2
        GsPt%tlsce = 1.0d0/sqrt(GsPt%tlsce)

    end subroutine setStabilizationTerms

    !---------------------------------------------------------------------

    subroutine setStabilizationTermsT_original(Problem, Nbf, Elem, GsPt)

        use Tools_mod, only: doubleDotProduct
        use PhysicalParameters_mod, only: BvN, ReN, WiN

        implicit none

        type(ProblemParameters), intent(in) :: Problem
        PetscInt, intent(in) :: Nbf
        type(NodeArrays), intent(in) :: Elem
        type(GaussPointQuantities), intent(inout) :: GsPt

        PetscInt  :: inod, j, Neq_f
        PetscReal :: tsugn, term, Helem, Uelem, Hugn, Ha
        PetscReal :: double_dot_T, double_dot_G, double_dot_S

        double_dot_T = doubleDotProduct(GsPt%T_ve,GsPt%T_ve)
        double_dot_G = doubleDotProduct(GsPt%G_dot,GsPt%G_dot)
        ! double_dot_G = doubleDotProduct(GsPt%G_dot_tr,GsPt%G_dot_tr)
        double_dot_S = doubleDotProduct(GsPt%S,GsPt%S)
        GsPt%mag_con_EQ = doubleDotProduct(GsPt%Con_EQ,GsPt%Con_EQ)
        GsPt%mag_con_EQ_time = doubleDotProduct(GsPt%Con_EQ_time,GsPt%Con_EQ_time)

        Neq_f = Problem%Neq_f

        Helem = 0.0d0
        Uelem = 0.0d0
        tsugn = 0.0d0
        do inod = 1, Nbf

            term = 0.0d0
            do j = 1, Neq_f
                term = term + Elem%TEMP_TL(inod,j)**2
            end do
            Uelem = Uelem + sqrt(term)

            Helem = Helem + sqrt(Elem%dFdXi(1,inod)**2 + Elem%dFdXi(2,inod)**2 + Elem%dFdXi(3,inod)**2) !2D
            ! Helem = Helem + (Elem%dFdXi(1,inod)**2 + Elem%dFdXi(2,inod)**2 + Elem%dFdXi(3,inod)**2) !3D

            tsugn = tsugn + abs( dot_product( GsPt%U_f,Elem%dFdXi(:,inod) ) )

        end do
        
        Uelem = Uelem/Nbf

        Helem = 1.0d0/Helem !2D
        ! Helem = 1.0d0/sqrt(Helem) !3D

        if (tsugn < 1.0d-8) tsugn = 1.0d-8
        Hugn = 0.0d0
        do j = 1, 3
            Hugn = Hugn + GsPt%U_f(j)**2
        end do
        Hugn = sqrt(Hugn)/tsugn
        if (Hugn < 1.0d-8) Hugn = 1.0d-8

        Ha = (((1.0d0-BvN)/WiN)**2)+0.5d0*double_dot_T
        Ha = sqrt(Ha/(1.0d0+0.5d0*double_dot_G))

        GsPt%tlsme = (ReN*Uelem/Hugn)**2 +(BvN/(Helem**2))**2 + (Ha/(Helem**2))**2
        GsPt%tlsme = 1.0d0/sqrt(GsPt%tlsme)

        GsPt%tlsic = (Helem**2)/GsPt%tlsme

        GsPt%tdcs = sqrt(((1.0d0-BvN)/WiN)**2+0.5d0*double_dot_S)
        GsPt%tdcs = (Helem**2)/GsPt%tdcs
        ! GsPt%tdcs = sqrt(0.5d0*GsPt%mag_con_EQ)*(Helem**2)/GsPt%tdcs

        GsPt%tlsce = (WiN*Uelem/Hugn)**2 +(WiN*sqrt(0.5d0*double_dot_G))**2 + GsPt%f_T**2 + 0.0d0
        GsPt%tlsce = 1.0d0/sqrt(GsPt%tlsce)

    end subroutine setStabilizationTermsT_original

    !---------------------------------------------------------------------

    subroutine setEdgeQuantities(ned, idir, FE_name, GsPt)

        implicit none

        PetscInt, intent(in) :: ned
        PetscInt, dimension(2), intent(in) :: idir
        character(*), intent(in) :: FE_name
        type(GaussPointQuantities), intent(inout) :: GsPt

        GsPt%normal(:) = 0.0d0
        GsPt%tangent(:) = 0.0d0
        GsPt%tangent2(:) = 0.0d0
        associate( dir1 => idir(1), dir2 => idir(2), &
                    dxdc => GsPt%dXidci(1,dir1), dydc => GsPt%dXidci(1,dir2), &
                    dxde => GsPt%dXidci(2,dir1), dyde => GsPt%dXidci(2,dir2) )

            select case (FE_name)

            case ('Triangle')

                select case(ned)
                case(1) !η=0
                    GsPt%normal(1) = dydc ; GsPt%normal(2) = -dxdc
                    ! GsPt%normal(:) = -GsPt%normal(:)
                case(2) !ξ=0
                    GsPt%normal(1) = dyde ; GsPt%normal(2) = -dxde
                    GsPt%normal(:) = -GsPt%normal(:)
                case(3) !1-ξ
                    GsPt%normal(1) = -dydc+dyde ; GsPt%normal(2) = dxdc-dxde
                    ! GsPt%normal(:) = -GsPt%normal(:)
                case default
                    write(*,'(a)') "Wrong ned choice in setEdgeQuantities!"
                    stop
                end select

            case ('Quadrangle')

                select case(ned)
                case(1) !η=-1
                    GsPt%normal(1) = dydc ; GsPt%normal(2) = -dxdc
                    ! GsPt%normal(:) = -GsPt%normal(:)
                case(2) !ξ=+1
                    GsPt%normal(1) = dyde ; GsPt%normal(2) = -dxde
                    ! GsPt%normal(:) = -GsPt%normal(:)
                case(3) !η=+1
                    GsPt%normal(1) = -dydc ; GsPt%normal(2) = dxdc
                    ! GsPt%normal(:) = -GsPt%normal(:)
                case(4) !ξ=-1
                    GsPt%normal(1) = -dyde ; GsPt%normal(2) = dxde
                    ! GsPt%normal(:) = -GsPt%normal(:)
                case default
                    write(*,'(a)') "Wrong ned choice in setEdgeQuantities!"
                    stop
                end select

            case default
                write(*,'(a)') "Wrong FE%name in setEdgeQuantities!"
                stop
            end select

        end associate

        GsPt%dL = sqrt(GsPt%normal(1)**2+GsPt%normal(2)**2)

        GsPt%normal(:) = GsPt%normal(:)/(GsPt%dL)
        ! GsPt%normal(:) = -GsPt%normal(:)
        ! print*, GsPt%normal(:)

        GsPt%tangent(1) = -GsPt%normal(2) ; GsPt%tangent(2) = +GsPt%normal(1)
        ! GsPt%tangent(:) = -GsPt%tangent(:)
        ! print*, GsPt%tangent(:)

    end subroutine setEdgeQuantities

    !---------------------------------------------------------------------

    subroutine setFaceQuantities(ifc, FE_name, GsPt)

        use Tools_mod, only: crossProduct3D

        implicit none

        PetscInt, intent(in) :: ifc
        character(*), intent(in) :: FE_name
        type(GaussPointQuantities), intent(inout) :: GsPt

        GsPt%normal(:) = 0.0d0
        GsPt%tangent(:) = 0.0d0
        GsPt%tangent2(:) = 0.0d0
        associate(dxdc => GsPt%dXidci(1,1), dydc => GsPt%dXidci(1,2), dzdc => GsPt%dXidci(1,3), &
                dxde => GsPt%dXidci(2,1), dyde => GsPt%dXidci(2,2), dzde => GsPt%dXidci(2,3), &
                dxds => GsPt%dXidci(3,1), dyds => GsPt%dXidci(3,2), dzds => GsPt%dXidci(3,3))

            select case (FE_name)

            case ('Tetrahedron')

                select case(ifc)
                case(1) !ξ=0
                    GsPt%tangent(1) = dxde
                    GsPt%tangent(2) = dyde
                    GsPt%tangent(3) = dzde

                    GsPt%tangent2(1) = dxds
                    GsPt%tangent2(2) = dyds
                    GsPt%tangent2(3) = dzds

                    GsPt%normal(1) = -(dyde*dzds-dzde*dyds)
                    GsPt%normal(2) = -(dzde*dxds-dxde*dzds)
                    GsPt%normal(3) = -(dxde*dyds-dyde*dxds)
                case(2) !η=0
                    GsPt%tangent(1) = dxds
                    GsPt%tangent(2) = dyds
                    GsPt%tangent(3) = dzds

                    GsPt%tangent2(1) = dxdc
                    GsPt%tangent2(2) = dydc
                    GsPt%tangent2(3) = dzdc

                    GsPt%normal(1) = -(dzdc*dyds-dydc*dzds)
                    GsPt%normal(2) = -(dxdc*dzds-dzdc*dxds)
                    GsPt%normal(3) = -(dydc*dxds-dxdc*dyds)
                case(3) !ζ=0
                    GsPt%tangent(1) = dxdc
                    GsPt%tangent(2) = dydc
                    GsPt%tangent(3) = dzdc

                    GsPt%tangent2(1) = dxde
                    GsPt%tangent2(2) = dyde
                    GsPt%tangent2(3) = dzde

                    GsPt%normal(1) = -(dydc*dzde-dzdc*dyde)
                    GsPt%normal(2) = -(dzdc*dxde-dxdc*dzde)
                    GsPt%normal(3) = -(dxdc*dyde-dydc*dxde)
                case(4) !1-ξ-η
                    GsPt%tangent(1) = dxdc-dxds
                    GsPt%tangent(2) = dydc-dyds
                    GsPt%tangent(3) = dzdc-dzds

                    GsPt%tangent2(1) = dxde-dxds
                    GsPt%tangent2(2) = dyde-dyds
                    GsPt%tangent2(3) = dzde-dzds

                    GsPt%normal(1) = (dydc-dyds)*(dzde-dzds)-(dzdc-dzds)*(dyde-dyds)
                    GsPt%normal(2) = (dzdc-dzds)*(dxde-dxds)-(dxdc-dxds)*(dzde-dzds)
                    GsPt%normal(3) = (dxdc-dxds)*(dyde-dyds)-(dydc-dyds)*(dxde-dxds)
                case default
                    write(*,'(a)') "Wrong ifc choice in Tetrahedron setFaceQuantities"
                    stop
                end select

            case ('Hexahedron')

                select case(ifc)
                case(1) !ζ=-1
                    GsPt%tangent(1) = dxdc
                    GsPt%tangent(2) = dydc
                    GsPt%tangent(3) = dzdc

                    GsPt%tangent2(1) = dxde
                    GsPt%tangent2(2) = dyde
                    GsPt%tangent2(3) = dzde

                    GsPt%normal(1) = -(dydc*dzde-dzdc*dyde)
                    GsPt%normal(2) = -(dzdc*dxde-dxdc*dzde)
                    GsPt%normal(3) = -(dxdc*dyde-dydc*dxde)
                case(2) !ζ=+1
                    GsPt%tangent(1) = dxdc
                    GsPt%tangent(2) = dydc
                    GsPt%tangent(3) = dzdc

                    GsPt%tangent2(1) = dxde
                    GsPt%tangent2(2) = dyde
                    GsPt%tangent2(3) = dzde

                    GsPt%normal(1) = dydc*dzde-dzdc*dyde
                    GsPt%normal(2) = dzdc*dxde-dxdc*dzde
                    GsPt%normal(3) = dxdc*dyde-dydc*dxde
                case(3) !η=-1
                    GsPt%tangent(1) = dxdc
                    GsPt%tangent(2) = dydc
                    GsPt%tangent(3) = dzdc

                    GsPt%tangent2(1) = dxds
                    GsPt%tangent2(2) = dyds
                    GsPt%tangent2(3) = dzds

                    GsPt%normal(1) = dydc*dzds-dzdc*dyds
                    GsPt%normal(2) = dzdc*dxds-dxdc*dzds
                    GsPt%normal(3) = dxdc*dyds-dydc*dxds
                case(4) !η=+1
                    GsPt%tangent(1) = dxdc
                    GsPt%tangent(2) = dydc
                    GsPt%tangent(3) = dzdc

                    GsPt%tangent2(1) = dxds
                    GsPt%tangent2(2) = dyds
                    GsPt%tangent2(3) = dzds

                    GsPt%normal(1) = -(dydc*dzds-dzdc*dyds)
                    GsPt%normal(2) = -(dzdc*dxds-dxdc*dzds)
                    GsPt%normal(3) = -(dxdc*dyds-dydc*dxds)
                case(5) !ξ=-1
                    GsPt%tangent(1) = dxde
                    GsPt%tangent(2) = dyde
                    GsPt%tangent(3) = dzde

                    GsPt%tangent2(1) = dxds
                    GsPt%tangent2(2) = dyds
                    GsPt%tangent2(3) = dzds

                    GsPt%normal(1) = -(dyde*dzds-dzde*dyds)
                    GsPt%normal(2) = -(dzde*dxds-dxde*dzds)
                    GsPt%normal(3) = -(dxde*dyds-dyde*dxds)
                case(6) !ξ=+1
                    GsPt%tangent(1) = dxde
                    GsPt%tangent(2) = dyde
                    GsPt%tangent(3) = dzde

                    GsPt%tangent2(1) = dxds
                    GsPt%tangent2(2) = dyds
                    GsPt%tangent2(3) = dzds

                    GsPt%normal(1) = dyde*dzds-dzde*dyds
                    GsPt%normal(2) = dzde*dxds-dxde*dzds
                    GsPt%normal(3) = dxde*dyds-dyde*dxds
                case default
                    write(*,'(a)') "Wrong ifc choice in Hexahedron setFaceQuantities!"
                    stop
                end select
            case default
                write(*,'(a)') "Wrong FE%name choice in setFaceQuantities!"
                stop
            end select

        end associate

        ! GsPt%normal(:) = crossProduct3D(GsPt%tangent1,GsPt%tangent2)

        GsPt%dL = sqrt(GsPt%normal(1)**2+GsPt%normal(2)**2+GsPt%normal(3)**2)

        GsPt%normal(:) = GsPt%normal(:)/(GsPt%dL)
        GsPt%normal(:) = -GsPt%normal(:) !for OBC
        ! write(*,'(i3,*(f12.4))') ifc, GsPt%normal
        ! pause

    end subroutine setFaceQuantities

end module ElementCalculations_mod
