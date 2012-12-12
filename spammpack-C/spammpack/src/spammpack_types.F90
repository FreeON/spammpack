MODULE SpAMMPACK_TYPES

  USE iso_c_binding

  TYPE SpAMM
    type(c_ptr) :: matrix
  END TYPE SpAMM

END MODULE SpAMMPACK_TYPES
