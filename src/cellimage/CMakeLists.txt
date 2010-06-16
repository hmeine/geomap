INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR})

VIGRA_ADD_NUMPY_MODULE(cellimage
  SOURCES
	cellimage_module.cxx
	cellimage_nodes.cxx
	cellimage_edges.cxx
	cellimage_faces.cxx
	cellimage_pyramid.cxx
	cellimage_dart+info.cxx
	foureightsegmentation.cxx
	cellstatistics.cxx
  LIBRARIES
    vigranumpy_core
)

ADD_DEPENDENCIES(vigranumpy cellimage)
