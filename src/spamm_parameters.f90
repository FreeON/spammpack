module spamm_numbers
  use spamm_real_precision
  implicit none
  !> Define the number zero.
  !real(SPAMM_KIND), parameter :: SpAMM_Zero = 0D0
  !> Define the number 1/2.
  !real(SPAMM_KIND), parameter :: SpAMM_Half = 5D-1
  !> Define the number 1.
  real(SPAMM_KIND), parameter :: SpAMM_One = 1D0
  !> Define the number 2.
  real(SPAMM_KIND), parameter :: SpAMM_Two = 2D0
  !> Define the number 4.
  real(SPAMM_KIND), parameter :: SpAMM_Four = 4D0
  !> Define the number 8.
  real(SPAMM_KIND), parameter :: SpAMM_Eight = 8D0
  !> Bigest machine double for ONX_KIND
  !real(SPAMM_KIND), parameter :: SpAMM_BIG_DBL = huge(SpAMM_One)
  !> Bigest machine int for int*4
  integer, parameter :: SpAMM_BIG_INT = 2**28
contains
end module spamm_numbers
