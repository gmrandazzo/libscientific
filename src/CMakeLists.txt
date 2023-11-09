set(Scientific_C_SRCS memwrapper.c
                   scientificinfo.c
                   optimization.c
                   vector.c
                   list.c
                   matrix.c
                   tensor.c
                   numeric.c
                   graphs.c
                   algebra.c
                   statistic.c
                   metricspace.c
                   interpolate.c
                   preprocessing.c
                   ica.c
                   pca.c
                   cpca.c
                   pls.c
                   epls.c
                   upca.c
                   upls.c
                   lda.c
                   clustering.c
                   mlr.c
                   datasets.c
                   modelvalidation.c
                   io.c)

set(Scientific_C_H
                scientificinfo.h
                optimization.h
                vector.h
                list.h
                matrix.h
                tensor.h
                numeric.h
                graphs.h
                algebra.h
                statistic.h
                metricspace.h
                interpolate.h
                preprocessing.h
                ica.h
                pca.h
                pls.h
                epls.h
                cpca.h
                upca.h
                upls.h
                lda.h
                clustering.h
                mlr.h
                datasets.h
                modelvalidation.h
                io.h)

if(WIN32)
  message(STATUS "Set mingw options")
  set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS_RELEASE)
  set(CMAKE_EXECUTABLE_SUFFIX ".exe")
  set(CMAKE_STATIC_LIBRARY_SUFFIX ".dll")
endif()

add_library(libscientific SHARED ${Scientific_C_SRCS} ${Scientific_C_H})
set_target_properties(libscientific PROPERTIES VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH} SOVERSION ${VERSION_MAJOR} )
set_target_properties(libscientific PROPERTIES OUTPUT_NAME "scientific")
target_link_libraries(libscientific PUBLIC m pthread ${LAPACK_LIBRARIES} ${SQLite3_LIBRARIES})

if(WIN32)
  install(TARGETS libscientific ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib COMPONENT library)
  install(FILES ${Scientific_C_H} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/scientific)
  install(FILES scientific.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
else()
  install(TARGETS libscientific LIBRARY DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR})
  install(FILES ${Scientific_C_H} DESTINATION include/scientific)
  install(FILES scientific.h DESTINATION include)
endif()

if(PORTABLE_PYTHON_PACKAGE)
  # Copy the compiled library with the dependencies 
  # into the python_binding library
  function(copy_library TARGET_NAME DESTINATION_DIR LIBRARY_SRC)
        get_filename_component(libname ${LIBRARY_SRC} NAME)
        add_custom_command(
            TARGET ${TARGET_NAME}
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
                "${LIBRARY_SRC}"
                "${DESTINATION_DIR}/${libname}"
            COMMENT "Copying ${libname} to ${DESTINATION_DIR}"
        )
  endfunction()
  set(dst_dir "${CMAKE_SOURCE_DIR}/src/python_bindings/libscientific/")
  copy_library(libscientific ${dst_dir} "${CMAKE_BINARY_DIR}/src/libscientific${CMAKE_SHARED_LIBRARY_SUFFIX}")
  if(WIN32)
    get_filename_component(COMPILER_DIRECTORY ${CMAKE_C_COMPILER} DIRECTORY)
    message("PORTABLE_PYTHON_PACKAGE: Supported only with MinWG64")
    find_file(LIBGCCSSEHDLL NAMES "libgcc_s_seh-1.dll" PATHS ${COMPILER_DIRECTORY} NO_DEFAULT_PATH REQUIRED)
    copy_library(libscientific ${dst_dir} ${LIBGCCSSEHDLL})
    find_file(QUADMATHDLL NAMES "libquadmath-0.dll" PATHS ${COMPILER_DIRECTORY} NO_DEFAULT_PATH REQUIRED)
    copy_library(libscientific ${dst_dir} ${QUADMATHDLL})
    find_file(GFORTRANDLL NAMES "libgfortran-5.dll" PATHS ${COMPILER_DIRECTORY} NO_DEFAULT_PATH REQUIRED)
    copy_library(libscientific ${dst_dir} ${GFORTRANDLL})
    find_file(BLASDLL NAMES "libblas.dll" PATHS ${COMPILER_DIRECTORY} NO_DEFAULT_PATH REQUIRED)
    copy_library(libscientific ${dst_dir} ${BLASDLL})
    find_file(LAPACKDLL NAMES "liblapack.dll" PATHS ${COMPILER_DIRECTORY} NO_DEFAULT_PATH REQUIRED)
    copy_library(libscientific ${dst_dir} ${LAPACKDLL})
    find_file(WINPTHREADDLL NAMES "libwinpthread-1.dll" PATHS ${COMPILER_DIRECTORY} NO_DEFAULT_PATH REQUIRED)
    copy_library(libscientific ${dst_dir} ${WINPTHREADDLL})
    find_file(SQLITE3DLL NAMES "libsqlite3-0.dll" PATHS ${COMPILER_DIRECTORY} NO_DEFAULT_PATH REQUIRED)
    copy_library(libscientific ${dst_dir} ${SQLITE3DLL})
  elseif(APPLE)
    message("PORTABLE_PYTHON_PACKAGE: Supported")
  else()
    message("PORTABLE_PYTHON_PACKAGE: Supported")
    set(lib_paths /lib
      /lib64
      /usr/lib
      /usr/local/lib
      /usr/lib64
      /usr/local/lib64 
      /usr/lib/x86_64-linux-gnu
      /usr/lib/gcc/x86_64-linux-gnu
      /usr/lib/gcc/x86_64-redhat-linux
    )
    copy_library(libscientific ${dst_dir} ${LAPACK_LIBRARIES})
    copy_library(libscientific ${dst_dir} ${BLAS_LIBRARIES})
    find_file(GFORTRAN_SO NAMES libgfortran.so.3 libgfortran.so.4 libgfortran.so.5 PATHS ${lib_paths} NO_DEFAULT_PATH REQUIRED)    
    copy_library(libscientific ${dst_dir} ${GFORTRAN_SO})
    find_file(QUADMATH_SO NAMES libquadmath.so.0.0.0 PATHS ${lib_paths} NO_DEFAULT_PATH REQUIRED)
    copy_library(libscientific ${dst_dir} ${QUADMATH_SO})
  endif()
endif()
