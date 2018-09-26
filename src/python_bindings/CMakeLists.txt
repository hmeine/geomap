VIGRA_FIND_PACKAGE( Boost 1.40.0 COMPONENTS signals )

ADD_DEFINITIONS(-DNDEBUG)
ADD_DEFINITIONS(-DAVOID_DUPLICATE_DEFINITION_FROM_VIGRANUMPYCORE)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ../geomap)

VIGRA_ADD_NUMPY_MODULE(geomap
  SOURCES
    vectorconv.cxx
    pythondiff2d.cxx
    polygon.cxx
    statistics.cxx
    cppmapmodule_utils.cxx
    cppmapmodule_stats.cxx
    cppmapmodule.cxx
    geomapmodule.cxx
    subpixel_watersheds.cxx
    dsl.cxx
  LIBRARIES
    libgeomap
    vigranumpy_core
    ${Boost_SIGNALS_LIBRARIES}
)

IF(PYTHON_PLATFORM MATCHES "^windows$")
    INSTALL(TARGETS geomap RUNTIME DESTINATION ${VIGRANUMPY_INSTALL_DIR})
ELSE()
    INSTALL(TARGETS geomap LIBRARY DESTINATION ${VIGRANUMPY_INSTALL_DIR})
ENDIF()

ADD_DEPENDENCIES(vigranumpy geomap)
