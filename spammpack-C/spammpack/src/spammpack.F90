!> @file

#include "config_fortran.h"

MODULE SpAMMPACK

  USE, INTRINSIC :: iso_c_binding
  USE SpAMMPACK_ALGEBRA
  USE SpAMMPACK_CHUNK
  USE SpAMMPACK_MANAGEMENT
  USE SpAMMPACK_TYPES

  IMPLICIT NONE

  INTERFACE

    !> Interface for spamm_convert_dense_to_spamm().
    SUBROUTINE spamm_convert_dense_to_spamm (ndim, N, chunk_tier, use_linear_tree, A_dense, A) &
        BIND(C, NAME = "spamm_convert_dense_to_spamm_interface")
      USE, INTRINSIC :: iso_c_binding
      INTEGER(c_int), INTENT(IN) :: ndim
      INTEGER(c_int), DIMENSION(ndim), INTENT(IN) :: N
      INTEGER(c_int), INTENT(IN) :: chunk_tier
      INTEGER(c_int), INTENT(IN) :: use_linear_tree
      REAL(KIND = c_float), DIMENSION(*) :: A_dense
      TYPE(c_ptr), INTENT(OUT) :: A
    END subroutine spamm_convert_dense_to_spamm

  END INTERFACE

CONTAINS

END MODULE SPAMMPACK
