add_library( utilities
  test_utilities.F90
  mmio.f
  )

add_executable( spamm-multiply spamm-multiply.F90 )
add_dependencies( spamm-multiply spammpack-serial-shared )
target_link_libraries( spamm-multiply utilities spammpack-serial-shared )

add_executable( spamm-lowdin spamm-lowdin.F90 )
add_dependencies( spamm-lowdin spammpack-serial-shared )
target_link_libraries( spamm-lowdin utilities spammpack-serial-shared )

include_directories( ${CMAKE_CURRENT_BINARY_DIR}/../src/include-serial )

# More unit-like tests.
set( UNIT_TESTS
  add_2nd_order_to_2nd_order.F90
  add_identity_to_2nd_order.F90
  allocate_2nd_order.F90
  copy_2nd_order.F90
  copy_matrix_column_to_vector.F90
  copy_order_1_to_order_1.F90
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
  set_2nd_order.F90
  trace_2nd_order.F90
  trace_2nd_order_2.F90
  )

execute_process( COMMAND ${PYTHON_EXECUTABLE}
  ${CMAKE_CURRENT_SOURCE_DIR}/generate_unit_tests.py
  --spammpack-lib spammpack-serial-shared
  --output ${CMAKE_CURRENT_BINARY_DIR}/unit-tests.cmake
  --lapack-libraries ${LAPACK_LIBRARIES}
  --
  ${UNIT_TESTS}
  )

include( ${CMAKE_CURRENT_BINARY_DIR}/unit-tests.cmake )
