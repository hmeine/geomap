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

IF(PYTHON_PLATFORM MATCHES "^windows$")
    INSTALL(TARGETS cellimage RUNTIME DESTINATION ${VIGRANUMPY_INSTALL_DIR})
ELSE()
    INSTALL(TARGETS cellimage LIBRARY DESTINATION ${VIGRANUMPY_INSTALL_DIR})
ENDIF()

ADD_DEPENDENCIES(vigranumpy cellimage)
