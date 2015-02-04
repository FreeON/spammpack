# Copyright (c) 2015, Los Alamos National Laboratory
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

set_directory_properties(
  PROPERTIES
  COMPILE_DEFINITIONS "${SPAMM_EXTRA_DEFINES}"
  )

# More unit-like tests.
set( UNIT_TESTS
  add_2nd_order_to_2nd_order.F90
  add_identity_to_2nd_order.F90
  allocate_2nd_order.F90
  allocate_tree_1d.F90
  allocate_tree_2d.F90
  copy_2nd_order.F90
  copy_matrix_column_to_vector.F90
  copy_order_1_to_order_1.F90
  copy_tree_1d.F90
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

include_directories(
  ${CMAKE_BINARY_DIR}/src/include-serial
  ${CMAKE_BINARY_DIR}/utilities
  )
