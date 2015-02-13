module spamm_decoration
  use spamm_real_precision
  implicit none

  type :: decoration_1d
     !> The integer length the vector 
     integer                           :: Dimen
     !> The axis-aligned bounding box.
     integer,  pointer, dimension(0:1) :: BndBx
     !> The square of the Frobenius norm.
     real(SPAMM_KIND)                  :: Norm2 = -1
     !> The number of non-zero elements.
     real(kind(0d0))                   :: FlOps = -1
     !> The number of non-zero elements.
     real(kind(0d0))                   :: Non0s = -1
  end type decoration_1d

  type :: decoration_2d
     !> The integer dimension of the native matrix
     integer,           dimension(1:2)     :: Dimem
     !> The axis-aligned bounding box.
     integer,  pointer, dimension(0:1,1:2) :: BndBx
     !> The square of the Frobenius norm.
     real(SPAMM_KIND)                      :: Norm2 = -1
     !> The number of non-zero elements.
     real(kind(0d0))                       :: FlOps = -1
     !> The number of non-zero elements.
     real(kind(0d0))                       :: Non0s = -1
  end type decoration_2d

contains
end module spamm_decoration
