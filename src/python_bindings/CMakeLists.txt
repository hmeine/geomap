include(FindPkgConfig)
PKG_CHECK_MODULES(sigcpp REQUIRED sigc++-2.0)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${sigcpp_INCLUDE_DIRS}
  ../cellimage)

VIGRA_ADD_NUMPY_MODULE(geomap
  SOURCES
    vectorconv.cxx
    pythondiff2d.cxx
    shapecontext.cxx
    geomapmodule.cxx
    polygon.cxx
    statistics.cxx
    features.cxx
    ellipse.cxx
    cppmapmodule_utils.cxx
    cppmapmodule_stats.cxx
    cppmapmodule.cxx
    cppmap.cxx
    cppmap_utils.cxx
    crackedgemap.cxx
    dsl.cxx
  LIBRARIES
    vigranumpy_core
    ${sigcpp_LIBRARIES}
)

ADD_DEPENDENCIES(vigranumpy geomap)
