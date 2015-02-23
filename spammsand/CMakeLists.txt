# Copyright (c) 2014, Los Alamos National Laboratory
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
                spammsand_structures.F90 
                spammsand.F90 
                test_utilities.F90  )

add_dependencies( spammsand_invsqrt spammpack-serial-shared )

target_link_libraries( spammsand_invsqrt
  utilities
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

