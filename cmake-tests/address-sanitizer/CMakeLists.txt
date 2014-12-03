project( address-sanitizer Fortran )
cmake_minimum_required( VERSION 2.8 )
add_executable( dummy dummy.F90 )
set_target_properties( dummy
  PROPERTIES
  COMPILE_FLAGS "-fsanitize=address"
  LINK_FLAGS "-fsanitize=address"
  )
