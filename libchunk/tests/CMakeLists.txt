add_executable( chunk_add chunk_add.c )
add_executable( chunk_add_identity chunk_add_identity.c )
add_executable( chunk_scale chunk_scale.c )
add_executable( chunk_trace chunk_trace.c )
add_executable( chunk_tree_alloc chunk_tree_alloc.c )
add_executable( chunk_multiply chunk_multiply.c bcsr.h bcsr.c )
add_executable( chunk_SP2 chunk_SP2.c bcsr.h bcsr.c )

include_directories( ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/src )
set_target_properties( chunk_add
  chunk_add_identity
  chunk_scale
  chunk_trace
  chunk_tree_alloc
  chunk_multiply
  chunk_SP2
  PROPERTIES
  COMPILE_FLAGS "${STDC} ${OpenMP_C_FLAGS} ${MIC_CFLAGS}"
  LINK_FLAGS "${OpenMP_C_FLAGS} ${MIC_LDFLAGS}" )
target_link_libraries( chunk_add libchunk ${MATHLIB} )
target_link_libraries( chunk_add_identity libchunk ${MATHLIB} )
target_link_libraries( chunk_scale libchunk ${MATHLIB} )
target_link_libraries( chunk_trace libchunk ${MATHLIB} )
target_link_libraries( chunk_tree_alloc libchunk ${MATHLIB} )
target_link_libraries( chunk_multiply libchunk ${MATHLIB} )
target_link_libraries( chunk_SP2 libchunk ${MATHLIB} )

add_test( chunk_add
  ${CMAKE_CURRENT_BINARY_DIR}/chunk_add )
add_test( chunk_add_identity
  ${CMAKE_CURRENT_BINARY_DIR}/chunk_add_identity )
add_test( chunk_scale
  ${CMAKE_CURRENT_BINARY_DIR}/chunk_scale )
add_test( chunk_trace
  ${CMAKE_CURRENT_BINARY_DIR}/chunk_trace )
add_test( chunk_tree_alloc
  ${CMAKE_CURRENT_BINARY_DIR}/chunk_tree_alloc )
