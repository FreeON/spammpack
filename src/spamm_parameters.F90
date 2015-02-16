module spamm_parameters

  use spamm_realkind

  implicit none

  INTEGER,          PARAMETER :: SpAMM_BLOCK_SIZE = 2

  real(SPAMM_KIND), parameter :: SpAMM_Zero  = 0D0
  real(SPAMM_KIND), parameter :: SpAMM_Half  = 5D-1
  real(SPAMM_KIND), parameter :: SpAMM_One   = 1D0
  real(SPAMM_KIND), parameter :: SpAMM_Two   = 2D0
  real(SPAMM_KIND), parameter :: SpAMM_Three = 3D0
  real(SPAMM_KIND), parameter :: SpAMM_Four  = 4D0
  real(SPAMM_KIND), parameter :: SpAMM_Eight = 8D0

  real(SPAMM_KIND), parameter :: SpAMM_BIG_DBL = huge(SpAMM_One)
  integer,          parameter :: SpAMM_BIG_INT = 2**28

contains

end module spamm_parameters
