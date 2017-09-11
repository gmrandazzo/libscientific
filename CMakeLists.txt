project(libscientific C Fortran)
cmake_minimum_required(VERSION 2.6)
set(AUTHOR "Giuseppe Marco Randazzo")
set(VERSION "0.6.0")


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(NOT CMAKE_INSTALL_PREFIX)
  message(SEND_ERROR "Please specify the installation prefix -DCMAKE_INSTALL_PREFIX=...")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
set(CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
set(CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")

if(WIN32)
  message(STATUS "Build Release: WIN32")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O2 -pipe -mwindows -pedantic -Wall -static -DWIN32")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libgcc -static-libgfortran -s -Wl,-O1")
elseif(APPLE)
  if(cmake_build_type_tolower STREQUAL "debug")
    message(STATUS "Build Debug: APPLE")
    set(CMAKE_C_FLAGS "-std=c99 -O1 -pedantic -Wall -g")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static -static-libgcc -pthread")
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O1 -g -fsanitize=address -fno-omit-frame-pointer -pipe -Wall")
#     set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -g -fsanitize=address")
  else()
    message(STATUS "Build Release: APPLE")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O2 -pipe -Wall")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}  -static -static-libgcc -pthread")
  endif()
else(NOT WIN32)
  if(cmake_build_type_tolower STREQUAL "debug")
    set(CMAKE_BUILD_TYPE "Debug")
    message(STATUS "Build Debug: Linux")
    #set(CMAKE_C_FLAGS "-O0 -pedantic -Werror -g -DDEBUG")
    set(CMAKE_C_FLAGS "-std=c99 -O1 -pedantic -Wall -g -pthread")
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
  else()
    message(STATUS "Build Release: Linux")
    set(CMAKE_BUILD_TYPE "Release")
      if(NOT CMAKE_C_FLAGS)
        set(CMAKE_C_FLAGS "-std=c99 -O3 -Wall -g -pthread")
      endif()
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-O1")
  endif()
endif(WIN32)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")

message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

add_subdirectory(src)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake" IMMEDIATE @ONLY)
add_custom_target(uninstall "${CMAKE_COMMAND}" -P
        "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")
