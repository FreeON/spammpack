add_executable( chunk-multiply
  chunk_multiply.c
  )

set_target_properties( chunk-multiply
  PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  )