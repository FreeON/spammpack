! cmake -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS="-O0 -g -fbounds-check -Wall -fbacktrace -finit-real=nan"
module spammpack

  use spamm_parameters
  use spamm_structures
  use spamm_decoration
  use spamm_xstructors
  use spamm_conversion
  use spamm_elementals
  use spamm_nbdyalgbra

contains

end module spammpack
