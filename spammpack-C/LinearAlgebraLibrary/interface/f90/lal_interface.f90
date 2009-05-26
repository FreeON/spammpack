module lal

  use lal_types

  interface

    function f90_lal_get (i, j, A)
      use lal_types
      integer, intent(in)           :: i, j
      type(lal_matrix), intent(in)  :: A
      real(double)                  :: f90_lal_get
    end function f90_lal_get

    subroutine f90_lal_set (i, j, Aij, A)
      use lal_types
      integer, intent(in)           :: i, j
      real(double), intent(in)      :: Aij
      type(lal_matrix), intent(in)  :: A
    end subroutine f90_lal_set

    function f90_lal_allocate (M, N, A)
      use lal_types
      integer, intent(in)             :: M, N
      type(lal_matrix), intent(inout) :: A
      integer                         :: f90_lal_allocate
    end function f90_lal_allocate

    subroutine f90_lal_free (A)
      use lal_types
      type(lal_matrix) :: A
    end subroutine f90_lal_free

    subroutine f90_lal_zero (A)
      use lal_types
      type(lal_matrix), intent(inout) :: A
    end subroutine f90_lal_zero

    subroutine f90_lal_rand (A)
      use lal_types
      type(lal_matrix), intent(inout) :: A
    end subroutine f90_lal_rand

    function f90_lal_equals (A, B, tolerance)
      use lal_types
      type(lal_matrix), intent(in) :: A, B
      real(double)                 :: tolerance
      integer                      :: f90_lal_equals
    end function f90_lal_equals

    subroutine f90_lal_dgemm (transA, transB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)
      use lal_types
      character        :: transA, transB
      integer          :: M, N, K, lda, ldb, ldc
      real(double)     :: alpha, beta
      type(lal_matrix) :: A, B, C
    end subroutine f90_lal_dgemm

    subroutine f90_lal_print (A)
      use lal_types
      type(lal_matrix) :: A
    end subroutine f90_lal_print

  end interface

end module lal
