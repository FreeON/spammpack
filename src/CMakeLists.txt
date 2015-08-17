add_library(spamm
  spamm.F90
  spamm_types.F90)

install(TARGETS spamm
  INCLUDES DESTINATION include
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)