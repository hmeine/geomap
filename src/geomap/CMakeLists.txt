include(FindPkgConfig)
PKG_CHECK_MODULES(sigcpp REQUIRED sigc++-2.0)

ADD_DEFINITIONS(-DNDEBUG)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${sigcpp_INCLUDE_DIRS})

ADD_LIBRARY(libgeomap ${LIBTYPE}
    cppmap.cxx
    cppmap_utils.cxx
    crackedgemap.cxx
)

TARGET_LINK_LIBRARIES(libgeomap ${sigcpp_LIBRARIES})

INSTALL(TARGETS libgeomap
#        EXPORT vigra-targets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib${LIB_SUFFIX}
        ARCHIVE DESTINATION lib${LIB_SUFFIX})
