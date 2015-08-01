module spamm_parameters

  use spamm_realkind

  implicit none

  ! h2o_350_6-311G**, blk = 8, tau= 0.05 
  ! h2o_350_6-311G**, blk = 64, tau= 0.1 

!  water, tau=1d-2/1d-3 and tau_y=1d-4/1d-5
!  INTEGER,          PARAMETER :: SpAMM_BLOCK_SIZE=2

!  tubes, tau=1d-2/1d-3 and tau_y=1d-4/1d-5
  INTEGER,          PARAMETER :: SpAMM_BLOCK_SIZE=8

!  INTEGER,          PARAMETER :: SpAMM_BLOCK_SIZE = 512

  INTEGER,          PARAMETER :: SBS              = SpAMM_BLOCK_SIZE
  INTEGER,          PARAMETER :: SBS2             = SpAMM_BLOCK_SIZE**2
  INTEGER,          PARAMETER :: SBS3             = SpAMM_BLOCK_SIZE**3
  !
  REAL(SpAMM_KIND), parameter :: SpAMM_normclean  = 1d-12
  REAL(SpAMM_KIND), parameter :: SpAMM_init  = 123456789d10

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

!to 0.13522%,    1.83456%,    0.35217%
! z 0.17314%,    2.80286%,    0.47922%
! h 0.16605%,    2.59818%,    0.47009%
