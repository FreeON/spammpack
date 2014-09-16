add_library( utilities
  test_utilities.F90
  mmio.f
  )

add_executable( spamm_multiply spamm_multiply.F90)
target_link_libraries( spamm_multiply utilities spammpack_serial_shared )
target_include_directories( spamm_multiply PRIVATE ${CMAKE_BINARY_DIR}/src )

if( OPENMP_FOUND )
  add_executable( spamm_multiply-threaded spamm_multiply.F90)
  set_target_properties( spamm_multiply-threaded
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS}
    )
  target_link_libraries( spamm_multiply-threaded utilities spammpack_threaded_shared )
  target_include_directories( spamm_multiply-threaded PRIVATE ${CMAKE_BINARY_DIR}/src )
endif()

# More unit-like tests.
set( UNIT_TESTS
  add_identify_to_2nd_order.F90
  allocate_2nd_order.F90
  copy_2nd_order.F90
  delete_2nd_order.F90
  lowdin.F90
  multiply_1st_order_x_scalar.F90
  multiply_1st_order_x_scalar_2.F90
  multiply_2nd_order_x_1st_order.F90
  multiply_2nd_order_x_2nd_order.F90
  multiply_2nd_order_x_2nd_order_2.F90
  multiply_2nd_order_x_2nd_order_3.F90
  multiply_2nd_order_x_scalar.F90
  norm_2nd_order.F90
  trace_2nd_order.F90
  trace_2nd_order_2.F90
  )

if( LAPACK_FOUND )
  set( UNIT_LAPACK "--lapack-linker-flags ${LAPACK_LINKER_FLAGS}" )
else()
  set( UNIT_LAPACK "" )
endif()

execute_process( COMMAND
  ${CMAKE_CURRENT_SOURCE_DIR}/generate_unit_tests.py
  --spammpack-lib spammpack_serial_shared
  --output ${CMAKE_CURRENT_BINARY_DIR}/unit-tests.cmake
  ${UNIT_LAPACK}
  ${UNIT_TESTS}
  )

include( ${CMAKE_CURRENT_BINARY_DIR}/unit-tests.cmake )
