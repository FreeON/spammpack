!> @file

#include "config_fortran.h"

MODULE SpAMMPACK

  use, intrinsic :: iso_c_binding
  use spammpack_types
  USE SpAMMPACK_ALGEBRA

  implicit none

  interface spamm_new

    !> Interface to spamm_new().
    subroutine spamm_new (ndim, N, chunk_tier, use_linear_tree, A)
      use, intrinsic :: iso_c_binding
      integer(c_int) :: ndim
      integer(c_int), dimension(:) :: N
      integer(c_int) :: chunk_tier
      logical :: use_linear_tree
      type(c_ptr), intent(out) :: A
    end subroutine spamm_new

  end interface spamm_new

  interface

    !> Interface for spamm_convert_dense_to_spamm().
    subroutine spamm_convert_dense_to_spamm (ndim, N, chunk_tier, use_linear_tree, A_dense, A) &
        bind(C, name = "spamm_convert_dense_to_spamm_interface")
      use, intrinsic :: iso_c_binding
      integer(c_int), intent(in) :: ndim
      integer(c_int), dimension(ndim), intent(in) :: N
      integer(c_int), intent(in) :: chunk_tier
      integer(c_int), intent(in) :: use_linear_tree
      real(kind = c_float), dimension(*) :: A_dense
      type(c_ptr), intent(out) :: A
    end subroutine spamm_convert_dense_to_spamm

  end interface

contains

END MODULE SPAMMPACK
