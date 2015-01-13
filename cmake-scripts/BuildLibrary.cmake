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

# Defines the function build_library.
#
# The variable MODULE_FILES has to be set externally before calling this
# function. It should contain the names of the generated module files. The
# variable HEADER_FILES should contain the header files from the source tree
# to install.
function( build_library use_OpenMP COMPILE_FLAGS )

  if( use_OpenMP )
    set( extension "threaded")
  else()
    set( extension "serial" )
  endif()

  add_library( spammpack-${extension}-shared SHARED ${spammpack-sources} )
  set_target_properties( spammpack-${extension}-shared
    PROPERTIES
    POSITION_INDEPENDENT_CODE TRUE
    COMPILE_FLAGS "${COMPILE_FLAGS}"
    SOVERSION 0
    VERSION 0.0.0
    OUTPUT_NAME ${LIBRARY_BASENAME}_${extension}
    )

  add_library( spammpack-${extension}-static STATIC ${spammpack-sources} )
  # Without these additional (chaining) dependencies, a parallel build
  # gets all confused because of the way CMake builds Fortran sources:
  # The compilation is done in the src dir, but uses -o to write the
  # .o file into a subdirectory named after the target. CMake does not
  # use -J or something similar to place the generated .mod file in
  # that directory as well. Instead it copies the .mod file from the
  # src directory into the target directory.  Parallel builds are
  # still possible within a library, but not across unfortunately.
  add_dependencies( spammpack-${extension}-static spammpack-${extension}-shared )
  set_target_properties( spammpack-${extension}-static
    PROPERTIES
    COMPILE_FLAGS "${COMPILE_FLAGS}"
    OUTPUT_NAME ${LIBRARY_BASENAME}_${extension}
    )

  install(
    TARGETS spammpack-${extension}-static spammpack-${extension}-shared
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

  file( MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include-${extension}" )
  add_custom_command( TARGET spammpack-${extension}-static
    POST_BUILD
    COMMAND cp ${MODULE_FILES} "include-${extension}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
  add_custom_command( TARGET spammpack-${extension}-static
    POST_BUILD
    COMMAND cp ${HEADER_FILES} "${CMAKE_CURRENT_BINARY_DIR}/include-${extension}"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

  unset( ${extension}-modfiles )
  foreach( mod ${MODULE_FILES} )
    list( APPEND ${extension}-modfiles
      ${CMAKE_CURRENT_BINARY_DIR}/include-${extension}/${mod}
      )
  endforeach()

  install(
    FILES ${${extension}-modfiles}
    DESTINATION include-${extension}
    )

endfunction( build_library )
