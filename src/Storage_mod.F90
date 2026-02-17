
#include <petsc/finclude/petscksp.h>

module Storage_mod

    use petscksp

    contains

    subroutine storeResidual(Res, grow1, petsc_vec)

        implicit none

        PetscScalar, dimension(:,:), intent(in) :: Res
        PetscInt, dimension(:), intent(in) :: grow1
        Vec, intent(inout) :: petsc_vec

        PetscInt :: inod, ieq, irow, Nbf, Neq
        PetscInt, dimension(size(Res,1)*size(Res,2)) :: grows
        PetscScalar, dimension(size(Res,1)*size(Res,2)) :: values
        PetscErrorCode :: ierr

        Nbf = size(Res,1)
        Neq = size(Res,2)

        irow = 0
        do inod = 1, Nbf
            do ieq = 1, Neq
                irow = irow+1
                grows(irow) = grow1(inod) + ieq - 1 - 1
                values(irow) = Res(inod,ieq)
            end do
        end do

        call VecSetValues(petsc_vec, size(grows), grows, values, ADD_VALUES, ierr)

    end subroutine storeResidual

    !----------------------------------------------------------------------

    subroutine storeJacobian(Jac, grow1, petsc_mat)

        implicit none

        PetscScalar, dimension(:,:,:,:), intent(in) :: Jac
        PetscInt, dimension(:), intent(in) :: grow1
        Mat, intent(inout) :: petsc_mat

        PetscInt :: inod, jnod, ieq, jeq, idim, jdim, kdim, ldim
        PetscInt :: irow, jcol, Jtot
        PetscInt, dimension(size(Jac,3)*size(Jac,4)) :: grows
        PetscInt, dimension(size(Jac,1)*size(Jac,2)) :: gcols
        PetscScalar, dimension(size(Jac,1)*size(Jac,2)*size(Jac,3)*size(Jac,4)) :: values
        PetscErrorCode :: ierr

        idim = size(Jac,1) ; jdim = size(Jac,2) ; kdim = size(Jac,3) ; ldim = size(Jac,4)

        irow = 0 ; Jtot = 0
        do jnod = 1, size(Jac,3)
            do jeq = 1, size(Jac,4)
                irow = irow+1 
                grows(irow) = grow1(jnod) + jeq - 1 - 1
                jcol = 0
                do inod = 1, size(Jac,1)
                    do ieq = 1, size(Jac,2)
                        jcol = jcol+1
                        gcols(jcol) = grow1(inod) + ieq - 1 - 1

                        Jtot = Jtot+1
                        values(Jtot) = Jac(inod,ieq,jnod,jeq)
                    end do
                end do
            end do
        end do
        
        call MatSetValues(petsc_mat, size(grows), grows, &
                          size(gcols), gcols, values, ADD_VALUES, ierr)

    end subroutine storeJacobian

    !----------------------------------------------------------------------

    subroutine storeJacobian_Ac(Ac, grow1, petsc_mat, jcol)

        implicit none

        PetscInt, dimension(:), intent(in) :: grow1
        PetscScalar, dimension(:,:), intent(in) :: Ac
        Mat, intent(inout) :: petsc_mat
        PetscInt, intent(in) :: jcol

        PetscInt :: inod, ieq, irow, Jtot, Nbf, Neq
        PetscInt, dimension(size(Ac,1)*size(Ac,2)) :: grows
        PetscScalar, dimension(size(Ac,1)*size(Ac,2)) :: values
        PetscErrorCode :: ierr

        Nbf = size(Ac,1)
        Neq = size(Ac,2)

        irow = 0
        do inod = 1, Nbf
            do ieq = 1, Neq

                irow = irow+1

                grows(irow) = grow1(inod) + ieq - 1 - 1

                values(irow) = Ac(inod,ieq)

            end do
        end do

        call MatSetValues(petsc_mat, size(grows), grows, 1, jcol, values, ADD_VALUES, ierr)

    end subroutine storeJacobian_Ac

    !----------------------------------------------------------------------

    subroutine storeResidualGlobal(Res_extra, petsc_vec, irow)

        implicit none

        PetscScalar, intent(in) :: Res_extra
        PetscInt, intent(in) :: irow
        Vec, intent(inout) :: petsc_vec

        PetscScalar :: value
        PetscErrorCode :: ierr

        value = Res_extra
        call VecSetValues(petsc_vec, 1, irow, value, ADD_VALUES, ierr)

    end subroutine storeResidualGlobal

    !----------------------------------------------------------------------

    subroutine storeJacobian_Ar(grow1, Ar, petsc_mat, irow)

        implicit none

        PetscInt, intent(in) :: irow
        PetscInt, dimension(:), intent(in) :: grow1
        PetscScalar, dimension(:,:), intent(in) :: Ar
        Mat, intent(inout) :: petsc_mat

        PetscInt :: ieq, inod, jcol, Jtot, Nbf, Neq
        PetscInt, dimension(size(Ar,1)*size(Ar,2)) :: gcols
        PetscScalar, dimension(size(Ar,1)*size(Ar,2)) :: values
        PetscErrorCode :: ierr

        Nbf = size(Ar,1)
        Neq = size(Ar,2)

        jcol = 0
        do inod = 1, Nbf
            do ieq = 1, Neq

                jcol = jcol+1

                gcols(jcol) = grow1(inod) + ieq - 1 - 1

                values(jcol) = Ar(inod,ieq)

            end do
        end do

        call MatSetValues(petsc_mat, 1, irow, size(gcols), gcols, values, ADD_VALUES, ierr)

    end subroutine storeJacobian_Ar

    !----------------------------------------------------------------------

    subroutine storeJacobian_Ah(Ah, petsc_mat, inex, irow)

        implicit none

        PetscScalar, dimension(:), intent(in) :: Ah
        Mat, intent(inout) :: petsc_mat
        PetscInt, intent(in) :: inex, irow

        PetscInt :: ieq, jcol, nex
        PetscInt, dimension(size(Ah)) :: global_cols
        PetscScalar, dimension(size(Ah)) :: values
        PetscErrorCode :: ierr

        Nex = size(Ah)
        jcol = 0
        do ieq = 1, Nex

            jcol = jcol+1

            global_cols(jcol) = irow - inex + ieq

            values(jcol) = Ah(ieq)

        end do

        call MatSetValues(petsc_mat, 1, irow, size(global_cols), global_cols, values, ADD_VALUES, ierr)

    end subroutine storeJacobian_Ah

    !----------------------------------------------------------------------

    subroutine storeJacobianPeriodic(Jac, grow1, grow1_periodic, petsc_mat)

        implicit none

        PetscScalar, dimension(:,:,:,:), intent(in) :: Jac
        PetscInt, dimension(:), intent(in) :: grow1
        PetscInt, dimension(:), intent(in) :: grow1_periodic
        Mat, intent(inout) :: petsc_mat

        PetscInt :: inod, jnod, ieq, jeq, idim, jdim, kdim, ldim
        PetscInt :: irow, jcol, Jtot
        PetscInt, dimension(size(Jac,3)*size(Jac,4)) :: grows
        PetscInt, dimension(size(Jac,1)*size(Jac,2)) :: gcols
        PetscScalar, dimension(size(Jac,1)*size(Jac,2)*size(Jac,3)*size(Jac,4)) :: values
        PetscErrorCode :: ierr

        idim = size(Jac,1) ; jdim = size(Jac,2) ; kdim = size(Jac,3) ; ldim = size(Jac,4)

        irow = 0 ; Jtot = 0

        do jnod = 1, size(Jac,3)
            do jeq = 1, size(Jac,4)
                irow = irow+1 
                grows(irow) = grow1_periodic(jnod) + jeq - 1 - 1
                jcol = 0
                do inod = 1, size(Jac,1)
                    do ieq = 1, size(Jac,2)
                        jcol = jcol+1
                        gcols(jcol) = grow1(inod) + ieq - 1 - 1

                        Jtot = Jtot+1
                        values(Jtot) = Jac(inod,ieq,jnod,jeq)
                    end do
                end do
            end do
        end do
        
        call MatSetValues(petsc_mat, size(grows), grows, &
                          size(gcols), gcols, values, ADD_VALUES, ierr)

    end subroutine storeJacobianPeriodic

end module Storage_mod
