!----------------------------------------------------------------------------------
! Structures for nested SpAMM solvers
!
MODULE SpAMMsand_structures

  USE  spammpack

  implicit none

  ! SpAMM sandwich 4 matrix functions: f(|a>) = |x_0>.|x_1> ... |x_m>
  type :: SpAMMsand_tree_2d_symm

     ! number of matrices in the nest ...
     integer                                             :: mats
     ! the SpAMM squared threshold for each algebra
     integer,                  dimension(:), allocatable :: tau2
     ! the preconditioner and residuals ...
     type(SpAMM_tree_2d_symm), dimension(:), pointer     :: trix

  end type SpAMMsand_tree_2d_symm

  ! --
contains

  ! --
end module spammsand_structures ! ... and we're out ...
