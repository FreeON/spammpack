
add_library( utilities
  test_utilities.F90
  mmio.f
  )

set_directory_properties(
  PROPERTIES
  COMPILE_DEFINITIONS "${SPAMM_EXTRA_DEFINES}"
  )


#add_executable( spamm-lowdin spamm-lowdin.F90 )
#add_dependencies( spamm-lowdin spammpack-serial-shared )
#target_link_libraries( spamm-lowdin
#  utilities
#  spammpack-serial-shared
#  ${LAPACK_LIBRARIES}
#  )

add_executable( spammsand_invsqrt 
                spammsand_inverse_squareroot.F90 
		spammsand_rqi_extremals.F90
                spammsand_structures.F90 
                spammsand.F90 
                test_utilities.F90  )

add_dependencies( spammsand_invsqrt spammpack-serial-shared )

target_link_libraries( spammsand_invsqrt
  spammpack-serial-shared
  ${LAPACK_LIBRARIES}
  )


#execute_process( COMMAND ${PYTHON_EXECUTABLE}
#  ${CMAKE_CURRENT_SOURCE_DIR}/generate_unit_tests.py
#  --spammpack-lib spammpack-serial-shared
#  --output ${CMAKE_CURRENT_BINARY_DIR}/unit-tests.cmake
#  --lapack-libraries ${LAPACK_LIBRARIES}
#  --
#  ${UNIT_TESTS}
#  )

#include( ${CMAKE_CURRENT_BINARY_DIR}/unit-tests.cmake )

include_directories( ${CMAKE_CURRENT_BINARY_DIR}/../src/include-serial )










#add_library( utilities
#  test_utilities.F90
#  mmio.f
#  )

#set_directory_properties(
#  PROPERTIES
#  COMPILE_DEFINITIONS "${SPAMM_EXTRA_DEFINES}"
#  )

#add_executable( spammsand_invsqrt 
#                spammsand_inverse_squareroot.F90 
#                spammsand_structures.F90 
#                spammsand.F90 
#                test_utilities.F90  )

#add_dependencies( spammsand_invsqrt spammpack-serial-shared )

#target_link_libraries( spammsand_invsqrt
#  utilities
#  spammpack-serial-shared
#  ${LAPACK_LIBRARIES}
#  )

#include_directories( ../src/include-serial ../src/ ../utilities )

