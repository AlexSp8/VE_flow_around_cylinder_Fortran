
#include <petsc/finclude/petscksp.h>

module Tools_mod

    implicit none

    interface traceTensor
        module procedure traceTensorReal
        module procedure traceTensorComplex
    end interface traceTensor

    interface doubleDotProduct
        module procedure doubleDotProductReal
        module procedure doubleDotProductRealComplex
        module procedure doubleDotProductComplexReal
        module procedure doubleDotProductComplex
    end interface doubleDotProduct

    interface gaussElimination
        module procedure gaussEliminationReal
        module procedure gaussEliminationRealComplex
    end interface gaussElimination

    interface perturbVariable
        module procedure perturbVariableReal
        module procedure perturbVariableComplex
    end interface perturbVariable

    interface allocateArray
        module procedure allocateIntArray1D
        module procedure allocateIntArray2D
        module procedure allocateIntArray3D
        module procedure allocateIntArray4D
        module procedure allocateRealArray1D
        module procedure allocateRealArray2D
        module procedure allocateRealArray3D
        module procedure allocateRealArray4D
        module procedure allocateScalarArray1D
        module procedure allocateScalarArray2D
        module procedure allocateScalarArray3D
        module procedure allocateScalarArray4D
        module procedure allocateCharacterArray1D
        module procedure allocateCharacterArray2D
        module procedure allocateCharacterArray3D
        module procedure allocateCharacterArray4D
    end interface allocateArray

    interface sortArray
        module procedure sortIntArray1D
        module procedure sortRealArray1D
        module procedure sortIntArray2D
        module procedure sortRealArray2D
    end interface sortArray

    interface sortFile
        module procedure sortFile1D
        module procedure sortFile2D
    end interface sortFile

    contains

    function crossProduct3D(a,b)

        implicit none

        PetscReal, dimension(3), intent(in) :: a, b
        PetscReal, dimension(3) :: crossProduct3D

        crossProduct3D(1) = a(2)*b(3) - a(3)*b(2)
        crossProduct3D(2) = a(3)*b(1) - a(1)*b(3)
        crossProduct3D(3) = a(1)*b(2) - a(2)*b(1)

    end function crossProduct3D

    !-----------------------------------------------------------------------  

    function determinant2D(A)

        implicit none

        PetscReal, dimension(2,2), intent(in) :: A
        PetscReal :: determinant2D

        determinant2D = A(1,1)*A(2,2) - A(1,2)*A(2,1)

    end function determinant2D

    !--------------------------------------------------------------------------

    function determinant3D(A)

        implicit none

        PetscReal, dimension(3,3), intent(in) :: A
        PetscReal :: determinant3D

        PetscReal, dimension(2,2) :: temp_A
        PetscReal :: term1, term2, term3

        temp_A(:,:) = A(2:3,2:3)
        term1 = A(1,1)*determinant2D(temp_A)

        temp_A(:,:) = A(2:3,1:3:2)
        term2 = A(1,2)*determinant2D(temp_A)

        temp_A(:,:) = A(2:3,1:2)
        term3 = A(1,3)*determinant2D(temp_A)

        determinant3D = term1-term2+term3

    end function determinant3D

    !-----------------------------------------------------------------------  

    subroutine inverseTensor1D(A, AI)

        implicit none

        PetscReal, dimension(1,1), intent(in) :: A
        PetscReal, dimension(1,1), intent(out) :: AI

        AI(1,1) =  1.0d0/A(1,1)
        
    end subroutine inverseTensor1D

    !-----------------------------------------------------------------------  

    subroutine inverseTensor2D(A, AI)

        implicit none

        PetscReal, dimension(:,:), intent(in) :: A
        PetscReal, dimension(:,:), intent(out) :: AI

        PetscReal :: detA

        detA = determinant2D(A)

        AI(:,:) = 0.0d0
        AI(1,1) =  A(2,2) ; AI(1,2) = -A(1,2)
        AI(2,1) = -A(2,1) ; AI(2,2) =  A(1,1)

        AI(:,:) = AI(:,:)/detA

    end subroutine inverseTensor2D

    !-----------------------------------------------------------------------  

    subroutine inverseTensor3D(A, AI)

        implicit none

        PetscReal, dimension(:,:), intent(in) :: A
        PetscReal, dimension(:,:), intent(out) :: AI

        ! PetscReal :: detA
        ! PetscReal, dimension(2,2) :: temp_A
        ! PetscInt :: i, j
        ! PetscInt, dimension(3,3), parameter :: &
        !     signs = reshape( [+1,-1,+1,-1,+1,-1,+1,-1,+1], shape(signs) )
        ! PetscInt, dimension(2,3), parameter :: &
        !     irows = reshape( [2,3, 1,3, 1,2], shape = shape(irows) )

        PetscReal :: detS, det11, det12, det13, det21, det22, det23, det31, det32, det33

        ! do i = 1, 3
        !     do j = 1, 3
        !         temp_A(:,:) = A(irows(:,i),irows(:,j))
        !         AI(i,j) = signs(i,j)*determinant2D(temp_A)
        !     end do
        ! end do

        ! AI(:,:) = transpose(AI)

        ! detA = determinant3D(A)

        ! AI(:,:) = AI(:,:)/detA


        det11 = A(2,2)*A(3,3) - A(2,3)*A(3,2)
        det12 = A(2,1)*A(3,3) - A(2,3)*A(3,1)
        det13 = A(2,1)*A(3,2) - A(2,2)*A(3,1)

        det21 = A(1,2)*A(3,3) - A(1,3)*A(3,2)
        det22 = A(1,1)*A(3,3) - A(1,3)*A(3,1)
        det23 = A(1,1)*A(3,2) - A(1,2)*A(3,1)

        det31 = A(1,2)*A(2,3) - A(1,3)*A(2,2)
        det32 = A(1,1)*A(2,3) - A(1,3)*A(2,1)
        det33 = A(1,1)*A(2,2) - A(1,2)*A(2,1)

        detS = A(1,1)*det11 - A(1,2)*det12 + A(1,3)*det13

        AI(:,:) = 0.0d0
        AI(1,1) =  det11 ; AI(1,2) = -det12 ; AI(1,3) =  det13
        AI(2,1) = -det21 ; AI(2,2) =  det22 ; AI(2,3) = -det23
        AI(3,1) =  det31 ; AI(3,2) = -det32 ; AI(3,3) =  det33

        AI = transpose(AI)/detS

    end subroutine inverseTensor3D

    !-----------------------------------------------------------------------  

    function traceTensorReal(A)

        implicit none

        PetscReal, dimension(:,:), intent(in) :: A
        PetscReal :: traceTensorReal

        PetscInt :: i

        traceTensorReal = 0.0d0
        do i = 1, size(A,1)
            traceTensorReal = traceTensorReal + A(i,i)
        end do

    end function traceTensorReal

    !-----------------------------------------------------------------------  

    function traceTensorComplex(A)

        implicit none

        PetscScalar, dimension(:,:), intent(in) :: A
        PetscScalar :: traceTensorComplex

        PetscInt :: i

        traceTensorComplex = 0.0d0
        do i = 1, size(A,1)
            traceTensorComplex = traceTensorComplex + A(i,i)
        end do

    end function traceTensorComplex

    !----------------------------------------------------------------------- 

    function vectorMagnitude(v)

        implicit none

        PetscReal, dimension(:), intent(in) :: v
        PetscReal :: vectorMagnitude

        PetscInt :: i

        vectorMagnitude = 0.0d0
        do i = 1, size(v)
            vectorMagnitude = vectorMagnitude + v(i)**2
        end do

        vectorMagnitude = sqrt(vectorMagnitude)

    end function vectorMagnitude

    !-----------------------------------------------------------------------  

    function doubleDotProductReal(A, B)

        implicit none

        PetscReal, dimension(:,:), intent(in)  :: A, B
        PetscReal :: doubleDotProductReal

        PetscInt :: i, j

        doubleDotProductReal = 0.0d0
        do i = 1, size(A,1)
            do j = 1, size(A,2)
                doubleDotProductReal = doubleDotProductReal + A(i,j)*B(j,i)
            end do
        end do

    end function doubleDotProductReal

    !----------------------------------------------------------------------- 

    function doubleDotProductRealComplex(A, B)

        implicit none

        PetscReal, dimension(:,:), intent(in)  :: A
        PetscScalar, dimension(:,:), intent(in)  :: B
        PetscScalar :: doubleDotProductRealComplex

        PetscInt :: i, j

        doubleDotProductRealComplex = 0.0d0
        do i = 1, size(A,1)
            do j = 1, size(A,2)
                doubleDotProductRealComplex = doubleDotProductRealComplex + A(i,j)*B(j,i)
            end do
        end do

    end function doubleDotProductRealComplex

    !----------------------------------------------------------------------- 

    function doubleDotProductComplexReal(A, B)

        implicit none

        PetscScalar, dimension(:,:), intent(in)  :: A
        PetscReal, dimension(:,:), intent(in)  :: B
        PetscScalar :: doubleDotProductComplexReal

        PetscInt :: i, j

        doubleDotProductComplexReal = 0.0d0
        do i = 1, size(A,1)
            do j = 1, size(A,2)
                doubleDotProductComplexReal = doubleDotProductComplexReal + A(i,j)*B(j,i)
            end do
        end do

    end function doubleDotProductComplexReal

    !----------------------------------------------------------------------- 

    function doubleDotProductComplex(A, B)

        implicit none

        PetscScalar, dimension(:,:), intent(in)  :: A, B
        PetscScalar :: doubleDotProductComplex

        PetscInt :: i, j

        doubleDotProductComplex = 0.0d0
        do i = 1, size(A,1)
            do j = 1, size(A,2)
                doubleDotProductComplex = doubleDotProductComplex + A(i,j)*B(j,i)
            end do
        end do

    end function doubleDotProductComplex

    !----------------------------------------------------------------------- 

    function tensorMagnitude(A) result(mag)

        implicit none

        PetscReal, dimension(:,:), intent(in) :: A
        PetscReal :: mag

        PetscInt :: i, j

        mag = 0.0d0
        do i = 1, size(A,1)
            do j = 1, size(A,2)
                mag = mag + A(i,j)**2
            end do
        end do

        mag = sqrt(0.5d0*mag)

    end function tensorMagnitude

    !----------------------------------------------------------------------- 

    function normFrobenius(A)

        implicit none

        PetscReal, dimension(:,:), intent(in) :: A
        PetscReal :: normFrobenius

        PetscReal :: double_dot_A

        double_dot_A = doubleDotProduct(A,A)

        normFrobenius = sqrt(double_dot_A)

    end function normFrobenius

    !----------------------------------------------------------------------- 

    function dyadicProduct(a, b) result(T)

        implicit none

        PetscReal, dimension(:), intent(in) :: a, b
        PetscReal, dimension(size(a), size(b)) :: T

        PetscInt :: i, j

        do i = 1, size(a)
            do j = 1, size(b)
                T(i,j) = a(i)*b(j)
            end do
        end do

    end function dyadicProduct

    !--------------------------------------------------------------------------

    subroutine extrapolationLagrange(xp, xb, xo, x, Lb, Lo, L)

        implicit none

        PetscReal, intent(in)  :: xp, xb, xo, x
        PetscReal, intent(out) :: Lb, Lo, L

        Lb = (xp-xo)*(xp-x) /((xb-xo)*(xb-x))
        Lo = (xp-xb)*(xp-x) /((xo-xb)*(xo-x))
        L  = (xp-xb)*(xp-xo)/((x -xb)*(x -xo))

    end subroutine extrapolationLagrange

    !--------------------------------------------------------------------------

    function interpolationLinear(x1, x2, y1, y2, y0) result(x_value)

        implicit none

        PetscReal, intent(in)  :: x1, x2, y1, y2, y0
        PetscReal :: x_value

        PetscReal :: slope

        slope = (y2-y1)/(x2-x1)

        x_value = x1 + (y0-y1)/slope

    end function interpolationLinear

    !---------------------------------------------------------------------

    function gaussEliminationReal(A0, b0) result(x)

        implicit none

        PetscReal, dimension(:,:), intent(in) :: A0
        PetscReal, dimension(:), intent(in) :: b0
        PetscReal, dimension(size(b0)) :: x

        PetscReal, dimension(size(x)) :: s, b
        PetscReal, dimension(size(x),size(x)) :: A
        PetscInt :: i, j, k, l, n
        PetscReal :: factor, sum, pivot, store, tol

        b(:) = b0(:)
        A(:,:) = A0(:,:)

        n = size(x)
        tol = 1.d-8

        k_loop:do k = 1, n-1

            !scaling: s(i) stores the largest element of row i
            do i = k, n             !rows loop
                s(i) = 0.d0
                do j = k, n           !columns of ith row
                    s(i) = max(s(i), abs(A(i,j)))
                end do
            end do

            !pivoting: change row with the one that has the largest element below
            pivot = abs(A(k,k)/s(k))
            l = k

            do j = k+1, n                             !for every element below kth column
                if(abs(A(j,k)/s(j)) > pivot) then      !if an element below is larger
                    pivot = abs(A(j,k)/s(j))             !change pivot value
                    l = j                                 !change position index of largest value of column
                end if
            end do

            if(pivot .lt. tol)then
                print*, "the matrix is singular!"
                return
            end if

            if (l /= k) then                          !if the pivot has changed
                do j = k, n                             !interchange the rows
                    store  = A(k,j)
                    A(k,j) = A(l,j)
                    A(l,j) = store
                end do
                store = b(k)                          !interchange the rhs
                b(k)  = b(l)
                b(l)  = store
            end if

            !forward elimination
            do i = k+1, n
                factor = A(i,k)/A(k,k)
                do j = k+1, n
                    A(i,j) = A(i,j)-factor*A(k,j)
                end do
                b(i) = b(i)-factor*b(k)
            end do

        end do k_loop

        !back substitution
        x(n) = b(n)/A(n,n)
        do i = n-1, 1,-1
            sum = b(i)
            do j = i+1, n
                sum = sum-A(i,j)*x(j)
            end do
            x(i) = sum/A(i,i)
        end do

    end function gaussEliminationReal

    !---------------------------------------------------------------------

    function gaussEliminationRealComplex(A0, b0) result(x)

        implicit none

        PetscReal, dimension(:,:), intent(in) :: A0
        PetscScalar, dimension(:), intent(in) :: b0
        PetscScalar, dimension(size(b0)) :: x

        PetscReal, dimension(size(x)) :: s
        PetscScalar, dimension(size(x)) :: b
        PetscReal, dimension(size(x),size(x)) :: A
        PetscInt :: i, j, k, l, n
        PetscReal :: factor, sum, pivot, store, tol

        b(:) = b0(:)
        A(:,:) = A0(:,:)

        n = size(x)
        tol = 1.d-8

        k_loop:do k = 1, n-1

            !scaling: s(i) stores the largest element of row i
            do i = k, n             !rows loop
                s(i) = 0.d0
                do j = k, n           !columns of ith row
                    s(i) = max(s(i), abs(A(i,j)))
                end do
            end do

            !pivoting: change row with the one that has the largest element below
            pivot = abs(A(k,k)/s(k))
            l = k

            do j = k+1, n                             !for every element below kth column
                if(abs(A(j,k)/s(j)) > pivot) then      !if an element below is larger
                    pivot = abs(A(j,k)/s(j))             !change pivot value
                    l = j                                 !change position index of largest value of column
                end if
            end do

            if(pivot .lt. tol)then
                print*, "the matrix is singular!"
                return
            end if

            if (l /= k) then                          !if the pivot has changed
                do j = k, n                             !interchange the rows
                    store  = A(k,j)
                    A(k,j) = A(l,j)
                    A(l,j) = store
                end do
                store = b(k)                          !interchange the rhs
                b(k)  = b(l)
                b(l)  = store
            end if

            !forward elimination
            do i = k+1, n
                factor = A(i,k)/A(k,k)
                do j = k+1, n
                    A(i,j) = A(i,j)-factor*A(k,j)
                end do
                b(i) = b(i)-factor*b(k)
            end do

        end do k_loop

        !back substitution
        x(n) = b(n)/A(n,n)
        do i = n-1, 1,-1
            sum = b(i)
            do j = i+1, n
                sum = sum-A(i,j)*x(j)
            end do
            x(i) = sum/A(i,i)
        end do

    end function gaussEliminationRealComplex

    !---------------------------------------------------------------------

    function dim(var)
        
        implicit none

        character(len=*), intent(in) :: var

        PetscInt :: dim

        select case(var)

        case('X')
            dim = 1
        case('Y')
            dim = 2
        case('Z')
            dim = 3
        case default
            write(*,'(A)') 'Wrong dimension in dim'
            stop
        end select

    end function dim

    !---------------------------------------------------------------------

    function var2D(var)
        
        implicit none

        character(len=*), intent(in) :: var

        PetscInt :: var2D

        select case(var)

        case('Ux')
            var2D = 1
        case('Uy')
            var2D = 2
        case('P')
            var2D = 3
        case('Txx')
            var2D = 4
        case('Txy')
            var2D = 5
        case('Tyy')
            var2D = 6
        case default
            write(*,'(A)') 'Wrong variable in var2D'
            stop
        end select

    end function var2D

    !---------------------------------------------------------------------

    function var3D(var)
        
        implicit none

        character(len=*), intent(in) :: var

        PetscInt :: var3D

        select case(var)

        case('Ux')
            var3D = 1
        case('Uy')
            var3D = 2
        case('Uz')
            var3D = 3
        case('P')
            var3D = 4
        case('Txx')
            var3D = 5
        case('Txy')
            var3D = 6
        case('Tyy')
            var3D = 7
        case('Txz')
            var3D = 8
        case('Tyz')
            var3D = 9
        case('Tzz')
            var3D = 10
        case default
            write(*,'(A)') 'Wrong variable in var3D'
            stop
        end select

    end function var3D

    !---------------------------------------------------------------------

    function extraEQ(var)
        
        implicit none

        character(len=*), intent(in) :: var

        PetscInt :: extraEQ

        select case(var)
        case('dPL')
            extraEQ = 1
        case('ArN')
            extraEQ = 2
        case default
            write(*,'(A)') 'Wrong EQ in extraEQ'
            stop
        end select

    end function extraEQ

    !---------------------------------------------------------------------

    function extraEQInflow(var)
        
        implicit none

        character(len=*), intent(in) :: var

        PetscInt :: extraEQInflow

        select case(var)

        case('dPL')
            extraEQInflow = 1
        case default
            write(*,'(A)') 'Wrong EQ in extraEQInflow'
            stop
        end select

    end function extraEQInflow

    !-----------------------------------------------------------------------  

    function getRows(nodes,eqs,Neq)

        implicit none

        PetscInt, intent(in) :: Neq
        PetscInt, dimension(:), intent(in) :: nodes, eqs
        PetscInt, dimension(size(nodes)*size(eqs)) :: getRows
        
        PetscInt :: i, lb, ub

        getRows = 0
        do i = 1, size(nodes)
            lb = size(eqs)*(i-1) + 1
            ub = size(eqs)*i
            getRows(lb:ub) = (nodes(i)-1)*Neq + eqs(:)
        end do

    end function getRows

    !---------------------------------------------------------------------

    function getRankNode(iel_rank, inod, Nbf) result(rnod)
        
        implicit none

        PetscInt, intent(in) :: iel_rank, inod, Nbf
        PetscInt :: rnod

        rnod = Nbf*(iel_rank-1)+inod

    end function getRankNode

    !---------------------------------------------------------------------

    subroutine getStressComponent(i, istart, j, k)
        
        implicit none

        PetscInt, intent(in) :: i, istart
        PetscInt, intent(out) :: j, k

        if (i == istart) then
            j = 1 ; k = 1
        else if (i == istart+1) then
            j = 1 ; k = 2
        else if (i == istart+2) then
            j = 2 ; k = 2
        else if (i == istart+3) then
            j = 1 ; k = 3
        else if (i == istart+4) then
            j = 2 ; k = 3
        else if (i == istart+5) then
            j = 3 ; k = 3
        else
            write(*,'(A)') 'Out of range getStressComponent!'
            stop
        end if

    end subroutine getStressComponent

    !---------------------------------------------------------------------

    subroutine getTensorComponent(i, j, k)
        
        implicit none

        PetscInt, intent(in) :: i
        PetscInt, intent(out) :: j, k

        if (i == 1) then
            j = 1 ; k = 1
        else if (i == 2) then
            j = 1 ; k = 2
        else if (i == 3) then
            j = 2 ; k = 1
        else if (i == 4) then
            j = 2 ; k = 2
        else if (i == 5) then
            j = 1 ; k = 3
        else if (i == 6) then
            j = 3 ; k = 1
        else if (i == 7) then
            j = 2 ; k = 3
        else if (i == 8) then
            j = 3 ; k = 2
        else if (i == 9) then
            j = 3 ; k = 3
        else
            write(*,'(A)') 'Out of range getStressComponent!'
            stop
        end if

    end subroutine getTensorComponent
    
    !-----------------------------------------------------------------------

    function perturbVariableReal(x) result(eps)

        implicit none

        PetscReal, intent(in) :: x
        PetscReal :: eps

        eps = 1.0d-8!*DMAX1(1.0d0,DABS(x))*SIGN(1.0d0,x)

    end function perturbVariableReal

    !--------------------------------------------------------------------------

    function perturbVariableComplex(x) result(eps)

        implicit none

        PetscScalar, intent(in) :: x
        PetscScalar :: eps

        complex(8) :: i = (0.0, 1.0)

        eps = 1.0d-8 + 1.0d-8*i

    end function perturbVariableComplex

    !---------------------------------------------------------------------

    function indexInArray(array,value) result(index)

        implicit none

        PetscInt, dimension(:), intent(in) :: array
        PetscInt, intent(in) :: value
        PetscInt :: index

        PetscInt :: i

        index = 0
        do i = 1, size(array)
            if (array(i) == value) then
                index = i
                exit
            end if
        end do

    end function indexInArray

    !---------------------------------------------------------------------

    function arrayClosestValue(array, value) result(x)

        implicit none

        PetscReal, dimension(:), intent(in) :: array
        PetscReal, intent(in) :: value
        PetscReal :: x

        PetscInt :: i
        PetscReal :: min_dif, dif

        x = array(1)
        min_dif = abs(array(1) - value)
        do i = 2, size(array)
            dif = abs(array(i) - value)
            if (dif < min_dif) then
                min_dif = dif
                x = array(i)
            end if
        end do

    end function arrayClosestValue

    !---------------------------------------------------------------------

    subroutine sortIntArray1D(x,y)

        implicit none

        PetscInt, dimension(:), intent(inout) :: x, y

        PetscInt :: imin, i, j
        PetscInt :: xmin, temp_x, temp_y

        do j = 1, size(x)-1
            xmin = x(j)
            imin = j
            do i = j+1, size(x)
                if (x(i) < xmin) then
                    xmin = x(i)
                    imin = i
                end if
            end do

            temp_x = x(j)
            x(j) = x(imin)
            x(imin) = temp_x

            temp_y = y(j)
            y(j) = y(imin)
            y(imin) = temp_y

        end do

    end subroutine sortIntArray1D

    !---------------------------------------------------------------------

    subroutine sortRealArray1D(x,y)

        implicit none

        PetscReal, dimension(:), intent(inout) :: x, y

        PetscInt :: imin, i, j
        PetscReal :: xmin, temp_x, temp_y

        do j = 1, size(x)-1
            xmin = x(j)
            imin = j
            do i = j+1, size(x)
                if (x(i) < xmin) then
                    xmin = x(i)
                    imin = i
                end if
            end do

            temp_x = x(j)
            x(j) = x(imin)
            x(imin) = temp_x

            temp_y = y(j)
            y(j) = y(imin)
            y(imin) = temp_y

        end do

    end subroutine sortRealArray1D

    !---------------------------------------------------------------------

    subroutine sortIntArray2D(xi)

        implicit none

        PetscInt, dimension(:,:), intent(inout) :: xi

        PetscInt :: imin, i, j, k
        PetscInt :: xmin, temp_x, temp_y

        do j = 1, size(xi,1)-1
            xmin = xi(j,1)
            imin = j
            do i = j+1, size(xi,1)
                if (xi(i,1) < xmin) then
                    xmin = xi(i,1)
                    imin = i
                end if
            end do

            temp_x = xi(j,1)
            xi(j,1) = xi(imin,1)
            xi(imin,1) = temp_x

            do k = 2, size(xi,2)
                temp_y = xi(j,k)
                xi(j,k) = xi(imin,k)
                xi(imin,k) = temp_y
            end do

        end do

    end subroutine sortIntArray2D

    !---------------------------------------------------------------------

    subroutine sortRealArray2D(xi)

        implicit none

        PetscReal, dimension(:,:), intent(inout) :: xi

        PetscInt :: imin, i, j, k
        PetscReal :: xmin, temp_x, temp_y

        do j = 1, size(xi,1)-1
            xmin = xi(j,1)
            imin = j
            do i = j+1, size(xi,1)
                if (xi(i,1) < xmin) then
                    xmin = xi(i,1)
                    imin = i
                end if
            end do

            temp_x = xi(j,1)
            xi(j,1) = xi(imin,1)
            xi(imin,1) = temp_x

            do k = 2, size(xi,2)
                temp_y = xi(j,k)
                xi(j,k) = xi(imin,k)
                xi(imin,k) = temp_y
            end do

        end do

    end subroutine sortRealArray2D

    !---------------------------------------------------------------------

    subroutine sortFile1D(fn)

        implicit none

        character(*), intent(in) :: fn

        PetscInt :: ioerr, n, i
        PetscReal :: dum1, dum2
        PetscReal, allocatable, dimension(:) :: x, y

        open(8,file=fn,action='read')

        n = 0
        do
            read(8,*,iostat=ioerr) dum1, dum2
            if (ioerr /= 0) exit
            n = n+1
        end do

        rewind(8)

        allocate(x(n)) ; allocate(y(n))

        do i = 1, size(x)
            read(8,*,iostat=ioerr) x(i), y(i)
        end do

        close(8)

        call sortArray(x,y)
        
        open(8,file=fn,action='write')

        do i = 1, size(x)
            write(8,'(*(es16.8))') x(i), y(i)
        end do

        close(8)

        deallocate(x)
        deallocate(y)

    end subroutine sortFile1D

    !---------------------------------------------------------------------

    subroutine sortFile2D(fn, ndim)

        implicit none

        character(*), intent(in) :: fn
        PetscInt, intent(in) :: ndim

        PetscInt :: ioerr, n, i
        PetscReal :: dum1, dum2
        PetscReal, allocatable, dimension(:,:) :: xi

        open(8,file=fn,action='read')

        n = 0
        do
            read(8,*,iostat=ioerr) dum1, dum2
            if (ioerr /= 0) exit
            n = n+1
        end do

        rewind(8)

        allocate(xi(n,ndim))
        n = 1
        do
            read(8,*,iostat=ioerr) ( xi(n,i), i = 1, ndim)
            if (ioerr /= 0) exit
            n = n+1
        end do

        close(8)

        call sortArray(xi)
        
        open(8,file=fn,action='write')

        do n = 1, size(xi,1)
            write(8,'(*(es16.8))') ( xi(n,i), i = 1, ndim)
        end do

        close(8)

        deallocate(xi)

    end subroutine sortFile2D

    !---------------------------------------------------------------------

    subroutine allocateRealArray1D(array,n)

        implicit none

        PetscReal, allocatable, intent(out) :: array(:)
        PetscInt, intent(in) :: n

        allocate(array(n))
        array = 0.0d0

    end subroutine allocateRealArray1D

    !---------------------------------------------------------------------

    subroutine allocateIntArray1D(array,n)

        implicit none

        PetscInt, allocatable, intent(out) :: array(:)
        PetscInt, intent(in) :: n

        allocate(array(n))
        array = 0

    end subroutine allocateIntArray1D

    !---------------------------------------------------------------------

    subroutine allocateCharacterArray1D(array,n)

        implicit none

        character(*), allocatable, intent(out) :: array(:)
        PetscInt, intent(in) :: n

        allocate(array(n))
        array = ''

    end subroutine allocateCharacterArray1D

    !---------------------------------------------------------------------

    subroutine allocateRealArray2D(array,n,m)

        implicit none

        PetscReal, allocatable, intent(out) :: array(:,:)
        PetscInt, intent(in) :: n,m

        allocate(array(n,m))
        array = 0.0d0

    end subroutine allocateRealArray2D

    !---------------------------------------------------------------------

    subroutine allocateIntArray2D(array,n,m)

        implicit none

        PetscInt, allocatable, intent(out) :: array(:,:)
        PetscInt, intent(in) :: n, m

        allocate(array(n,m))
        array = 0

    end subroutine allocateIntArray2D

    !---------------------------------------------------------------------

    subroutine allocateCharacterArray2D(array,n,m)

        implicit none

        character(*), allocatable, intent(out) :: array(:,:)
        PetscInt, intent(in) :: n,m

        allocate(array(n,m))
        array = ''

    end subroutine allocateCharacterArray2D

    !---------------------------------------------------------------------

    subroutine allocateRealArray3D(array,n,m,k)

        implicit none

        PetscReal, allocatable, intent(out) :: array(:,:,:)
        PetscInt, intent(in) :: n,m,k

        allocate(array(n,m,k))
        array = 0.0d0

    end subroutine allocateRealArray3D

    !---------------------------------------------------------------------

    subroutine allocateIntArray3D(array,n,m,k)

        implicit none

        PetscInt, allocatable, intent(out) :: array(:,:,:)
        PetscInt, intent(in) :: n,m,k

        allocate(array(n,m,k))
        array = 0

    end subroutine allocateIntArray3D

    !---------------------------------------------------------------------

    subroutine allocateCharacterArray3D(array,n,m,k)

        implicit none

        character(*), allocatable, intent(out) :: array(:,:,:)
        PetscInt, intent(in) :: n,m,k

        allocate(array(n,m,k))
        array = ''

    end subroutine allocateCharacterArray3D

    !---------------------------------------------------------------------

    subroutine allocateRealArray4D(array,n,m,k,l)

        implicit none

        PetscReal, allocatable, intent(out) :: array(:,:,:,:)
        PetscInt, intent(in) :: n,m,k,l

        allocate(array(n,m,k,l))
        array = 0.0d0

    end subroutine allocateRealArray4D

    !---------------------------------------------------------------------

    subroutine allocateIntArray4D(array,n,m,k,l)

        implicit none

        PetscInt, allocatable, intent(out) :: array(:,:,:,:)
        PetscInt, intent(in) :: n,m,k,l

        allocate(array(n,m,k,l))
        array = 0

    end subroutine allocateIntArray4D

    !---------------------------------------------------------------------

    subroutine allocateCharacterArray4D(array,n,m,k,l)

        implicit none

        character(*), allocatable, intent(out) :: array(:,:,:,:)
        PetscInt, intent(in) :: n,m,k,l

        allocate(array(n,m,k,l))
        array = ''

    end subroutine allocateCharacterArray4D

    !---------------------------------------------------------------------

    subroutine allocateScalarArray1D(array,n)

        implicit none

        PetscScalar, allocatable, intent(out) :: array(:)
        PetscInt, intent(in) :: n

        PetscScalar, parameter :: i = (0.0d0,1.0d0)
        
        allocate(array(n))
        array = 0.0d0 + 0.0d0*i

    end subroutine allocateScalarArray1D

    !---------------------------------------------------------------------

    subroutine allocateScalarArray2D(array,n,m)

        implicit none

        PetscScalar, allocatable, intent(out) :: array(:,:)
        PetscInt, intent(in) :: n,m

        PetscScalar, parameter :: i = (0.0d0,1.0d0)
        
        allocate(array(n,m))
        array = 0.0d0 + 0.0d0*i

    end subroutine allocateScalarArray2D

    !---------------------------------------------------------------------

    subroutine allocateScalarArray3D(array,n,m,k)

        implicit none

        PetscScalar, allocatable, intent(out) :: array(:,:,:)
        PetscInt, intent(in) :: n,m,k

        PetscScalar, parameter :: i = (0.0d0,1.0d0)
        
        allocate(array(n,m,k))
        array = 0.0d0 + 0.0d0*i

    end subroutine allocateScalarArray3D

    !---------------------------------------------------------------------

    subroutine allocateScalarArray4D(array,n,m,k,l)

        implicit none

        PetscScalar, allocatable, intent(out) :: array(:,:,:,:)
        PetscInt, intent(in) :: n,m,k,l

        PetscScalar, parameter :: i = (0.0d0,1.0d0)
        
        allocate(array(n,m,k,l))
        array = 0.0d0 + 0.0d0*i

    end subroutine allocateScalarArray4D

    !----------------------------------------------------------------------- 

    subroutine makeDirectories(continuation_method)

        implicit none

        character(len=*), intent(in) :: continuation_method

        call execute_command_line('mkdir -p Results/')
        call execute_command_line('cp -r delete_non_sol.sh Results/')
        call execute_command_line('mkdir -p Results/Base/')
        call execute_command_line('mkdir -p Results/Base/Info')
        if (continuation_method == 'Stability') then
            call execute_command_line('mkdir -p Results/Stability/')
        end if

    end subroutine makeDirectories

    !---------------------------------------------------------------------

    subroutine writeElapsedTime(start_time,end_time,str)

        implicit none

        PetscReal,        intent(in) :: start_time, end_time
        character(len=*), intent(in) :: str

        PetscReal :: total_time, seconds
        PetscInt  :: days, hours, minutes

        total_time = end_time-start_time
        days    = int(total_time/8.64d4)
        hours   = int((total_time-days*8.64d4)/3.6d3)
        minutes = int((total_time-days*8.64d4-hours*3.6d3)/60.0d0)
        seconds = total_time-days*8.64d4-hours*3.6d3-minutes*60.0d0

        write(*, 60) str//': ', days, hours, minutes, seconds
        write(*, '(a)') '--------------------------------------------------------------'

        60 format(A, T60, I2, " d, ", I2, " hr, ", I2, " min, ", F6.3, " s")

    end subroutine writeElapsedTime

    !----------------------------------------------------------------------- 

    subroutine memoryUsage()

        implicit none

        character(len=100) :: command, string
        PetscInt :: status

        command = "free -h >> memory_usage.DAT"
        call execute_command_line(command, cmdstat=status)
        if (status /= 0) then
            write(*,'(a)') "Error executing the command:", trim(command)
            stop
        endif

        command = 'echo " " >> memory_usage.DAT'
        call execute_command_line(command, cmdstat=status)

    end subroutine memoryUsage

    !----------------------------------------------------------------------- 

    function updateEMA(x_new, a, EMAo) result(EMA)

        implicit none

        PetscReal, intent(in)  :: x_new     !new value
        PetscReal, intent(in)  :: a         !smoothing factor (0 < a < 1)
        PetscReal, intent(in)  :: EMAo      !previous EMA value
        PetscReal :: EMA                    !updated EMA value

        EMA = a*max(0.0d0,x_new) + (1.0d0-a)*EMAo

    end function updateEMA

end module Tools_mod
