! cmake -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_Fortran_FLAGS="-O0 -g -fbounds-check -Wall -fbacktrace -finit-real=nan"

! cmake -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_Fortran_FLAGS="-DLAPACK_FOUND -stand f08 -O0 -g -extend-source -debug all -check all -warn unused -traceback"

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
