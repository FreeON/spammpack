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

set(spammpack-sources
  spamm_realkind.F90
  spamm_parameters.F90
  spamm_structures.F90
  spamm_conversion.F90
  spamm_xstructors.F90
  spamm_decoration.F90
  spamm_elementals.F90
  spamm_nbdyalgbra.F90
  spamm_nbdyalgbra_times.F90
  spamm_nbdyalgbra_plus.F90
  spamm_nbdyalgbra_trace.F90
  spammpack.F90)

set(MODULE_FILES
  spamm_realkind.mod
  spamm_parameters.mod
  spamm_structures.mod
  spamm_conversion.mod
  spamm_xstructors.mod
  spamm_decoration.mod
  spamm_elementals.mod
  spamm_nbdyalgbra.mod
  spamm_nbdyalgbra_times.mod
  spamm_nbdyalgbra_plus.mod
  spamm_nbdyalgbra_trace.mod
  spammpack.mod)

set(LIBRARY_BASENAME "spammpack")

include(${CMAKE_SOURCE_DIR}/cmake-scripts/FindFortranKind.cmake)

set(KIND_EXPRESSION "kind(0d0)")
find_fortran_kind(${KIND_EXPRESSION})

#configure_file(spamm_real_precision.F90.in spamm_real_precision.F90)

set_directory_properties(
  PROPERTIES
  COMPILE_DEFINITIONS "${SPAMM_EXTRA_DEFINES}")

include(${CMAKE_SOURCE_DIR}/cmake-scripts/BuildLibrary.cmake)

build_library(FALSE "")
