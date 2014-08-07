add_library( spammtests
  test_utilities.F90
  spammtests.F90
  )
target_include_directories( spammtests PRIVATE ${CMAKE_BINARY_DIR}/src )

add_executable( spamm_multiply spamm_multiply.F90)
target_link_libraries( spamm_multiply spammpack_serial_shared spammtests )
target_include_directories( spamm_multiply PRIVATE ${CMAKE_BINARY_DIR}/src )

if( OPENMP_FOUND )
  add_executable( spamm_multiply-threaded spamm_multiply.F90)
  set_target_properties( spamm_multiply-threaded
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS}
    )
  target_link_libraries( spamm_multiply-threaded spammpack_threaded_shared spammtests )
  target_include_directories( spamm_multiply-threaded PRIVATE ${CMAKE_BINARY_DIR}/src )
endif()

add_executable( multiply_2nd_order_x_2nd_order multiply_2nd_order_x_2nd_order.F90)
target_link_libraries( multiply_2nd_order_x_2nd_order spammpack_serial_shared spammtests )
target_include_directories( multiply_2nd_order_x_2nd_order PRIVATE ${CMAKE_BINARY_DIR}/src )
add_test( multiply_2nd_order_x_2nd_order multiply_2nd_order_x_2nd_order )

add_executable( allocate_2nd_order allocate_2nd_order.F90)
target_link_libraries( allocate_2nd_order spammpack_serial_shared spammtests )
target_include_directories( allocate_2nd_order PRIVATE ${CMAKE_BINARY_DIR}/src )
add_test( allocate_2nd_order allocate_2nd_order )
