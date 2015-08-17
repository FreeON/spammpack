set(tests
  01-create)

include_directories(${CMAKE_BINARY_DIR}/src)

foreach(test ${tests})
  message(STATUS "Creating test ${test}")
  add_executable(${test} ${test}.F90)
  add_test(NAME ${test} COMMAND ${test})
endforeach()
