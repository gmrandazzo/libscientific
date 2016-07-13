# Install script for directory: /Users/marco/projects/libscientific-0.7.4/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/lib/libscientific.dylib")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/lib" TYPE SHARED_LIBRARY FILES "/Users/marco/projects/libscientific-0.7.4/build/src/libscientific.dylib")
  if(EXISTS "$ENV{DESTDIR}/usr/local/lib/libscientific.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/usr/local/lib/libscientific.dylib")
    execute_process(COMMAND "/usr/bin/install_name_tool"
      -id "libscientific.dylib"
      "$ENV{DESTDIR}/usr/local/lib/libscientific.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "$ENV{DESTDIR}/usr/local/lib/libscientific.dylib")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/include/scientific/scientificinfo.h;/usr/local/include/scientific/optimization.h;/usr/local/include/scientific/vector.h;/usr/local/include/scientific/matrix.h;/usr/local/include/scientific/array.h;/usr/local/include/scientific/numeric.h;/usr/local/include/scientific/graphs.h;/usr/local/include/scientific/algebra.h;/usr/local/include/scientific/statistic.h;/usr/local/include/scientific/metricspace.h;/usr/local/include/scientific/pca.h;/usr/local/include/scientific/pls.h;/usr/local/include/scientific/upca.h;/usr/local/include/scientific/upls.h;/usr/local/include/scientific/lda.h;/usr/local/include/scientific/clustering.h;/usr/local/include/scientific/mlr.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/include/scientific" TYPE FILE FILES
    "/Users/marco/projects/libscientific-0.7.4/src/scientificinfo.h"
    "/Users/marco/projects/libscientific-0.7.4/src/optimization.h"
    "/Users/marco/projects/libscientific-0.7.4/src/vector.h"
    "/Users/marco/projects/libscientific-0.7.4/src/matrix.h"
    "/Users/marco/projects/libscientific-0.7.4/src/array.h"
    "/Users/marco/projects/libscientific-0.7.4/src/numeric.h"
    "/Users/marco/projects/libscientific-0.7.4/src/graphs.h"
    "/Users/marco/projects/libscientific-0.7.4/src/algebra.h"
    "/Users/marco/projects/libscientific-0.7.4/src/statistic.h"
    "/Users/marco/projects/libscientific-0.7.4/src/metricspace.h"
    "/Users/marco/projects/libscientific-0.7.4/src/pca.h"
    "/Users/marco/projects/libscientific-0.7.4/src/pls.h"
    "/Users/marco/projects/libscientific-0.7.4/src/upca.h"
    "/Users/marco/projects/libscientific-0.7.4/src/upls.h"
    "/Users/marco/projects/libscientific-0.7.4/src/lda.h"
    "/Users/marco/projects/libscientific-0.7.4/src/clustering.h"
    "/Users/marco/projects/libscientific-0.7.4/src/mlr.h"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/include/scientific.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/usr/local/include" TYPE FILE FILES "/Users/marco/projects/libscientific-0.7.4/src/scientific.h")
endif()

