
#include <petsc/finclude/petscksp.h>

module Gauss_mod

    contains

    subroutine setGaussParameters(FE, GaussInt)

        use FEMParameters_mod, only: FEParameters
        use GaussParameters_mod, only: GaussIntegration
        use GaussLine_mod, only: setGaussParametersLine
        use GaussTriangle_mod, only: setGaussParametersTriangle
        use GaussQuadrangle_mod, only: setGaussParametersQuadrangle
        use GaussTetrahedron_mod, only: setGaussParametersTetrahedron
        use GaussHexahedron_mod, only: setGaussParametersHexahedron

        implicit none

        type(FEParameters), intent(in) :: FE
        type(GaussIntegration), intent(out) :: GaussInt

        GaussInt%order = FE%order

        GaussInt%Nbf = FE%Nbf
        GaussInt%Ned = FE%Ned
        GaussInt%Nfc = FE%Nfc

        select case(FE%name)
        case ('Line')
            call setGaussParametersLine(GaussInt)
        case ('Triangle')
            call setGaussParametersTriangle(GaussInt)
        case ('Quadrangle')
            call setGaussParametersQuadrangle(GaussInt)
        case ('Tetrahedron')
            call setGaussParametersTetrahedron(GaussInt)
        case ('Hexahedron')
            call setGaussParametersHexahedron(GaussInt)
        case default
            return
        end select
        
    end subroutine setGaussParameters

end module Gauss_mod
