add_library( libchunk
  chunk.c
  chunk_block.c
  chunk_tiled.c
  chunk_tree.c
  chunk_mic.c
  chunk_mic.h
  lapack_interface.h
  matrix_4x4.S
  )
include_directories( ${CMAKE_BINARY_DIR} )
set_target_properties( libchunk
  PROPERTIES
  COMPILE_FLAGS "${STDC} ${OpenMP_C_FLAGS} ${MIC_CFLAGS}"
  LINK_FLAGS "${MIC_LDFLAGS}" )

set( headers
	chunk.h
	chunk_block.h
	chunk_tiled.h
	chunk_tree.h
	)

install( TARGETS libchunk
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib )
install( FILES ${headers}
  DESTINATION include )
