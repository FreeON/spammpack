module spamm_parameters

  implicit none

  integer, parameter :: SPAMM_KIND = SPAMM_KIND_EXPRESSION

  integer, parameter :: SBS = SPAMM_BLOCK_SIZE
  integer, parameter :: SBS2 = SPAMM_BLOCK_SIZE**2
  integer, parameter :: SBS3 = SPAMM_BLOCK_SIZE**3

  real(SPAMM_KIND), parameter :: SpAMM_normclean  = 1e-12_SPAMM_KIND
  real(SPAMM_KIND), parameter :: SpAMM_init  = 123456789e10_SPAMM_KIND

  !real(SPAMM_KIND), parameter :: SpAMM_Zero = 0e0_SPAMM_KIND
  real(SPAMM_KIND), parameter :: SpAMM_Half = 5e-1_SPAMM_KIND
  real(SPAMM_KIND), parameter :: SpAMM_One = 1e0_SPAMM_KIND
  real(SPAMM_KIND), parameter :: SpAMM_Two = 2e0_SPAMM_KIND
  real(SPAMM_KIND), parameter :: SpAMM_Three = 3e0_SPAMM_KIND
  real(SPAMM_KIND), parameter :: SpAMM_Four = 4e0_SPAMM_KIND
  real(SPAMM_KIND), parameter :: SpAMM_Eight = 8e0_SPAMM_KIND

  real(SPAMM_KIND), parameter :: SpAMM_BIG_DBL = huge(SpAMM_One)
  integer, parameter :: SpAMM_BIG_INT = 2**28

end module spamm_parameters
