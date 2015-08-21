add_executable(spammsand_invsqrt
  spammsand_inverse_squareroot.F90
  spammsand_rqi_extremals.F90
  spammsand_structures.F90
  spammsand.F90
  test_utilities.F90
  mmio.f)
set_target_properties(spammsand_invsqrt
  PROPERTIES
  COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS}")

include_directories(${CMAKE_BINARY_DIR}/src)
target_link_libraries(spammsand_invsqrt spamm)

install(TARGETS spammsand_invsqrt
  RUNTIME DESTINATION bin)