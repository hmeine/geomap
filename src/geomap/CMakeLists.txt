VIGRA_FIND_PACKAGE( Boost 1.40.0 COMPONENTS signals2 )

ADD_DEFINITIONS(-DNDEBUG)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${Boost_INCLUDE_DIRS})

ADD_LIBRARY(libgeomap ${LIBTYPE}
    cppmap.cxx
    cppmap_utils.cxx
    crackedgemap.cxx
)

SET(libgeomap_OUTPUT_NAME "libgeomap" CACHE STRING "geomap core library output name")

SET_TARGET_PROPERTIES(libgeomap PROPERTIES
    OUTPUT_NAME ${libgeomap_OUTPUT_NAME})

INSTALL(TARGETS libgeomap
#        EXPORT vigra-targets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib${LIB_SUFFIX}
        ARCHIVE DESTINATION lib${LIB_SUFFIX})
