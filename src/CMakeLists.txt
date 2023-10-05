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
target_link_libraries(libscientific pthread ${LAPACK_LIBRARIES} ${SQLite3_LIBRARIES})

if(WIN32)
  install(TARGETS libscientific ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib COMPONENT library)
  install(FILES ${Scientific_C_H} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/scientific)
  install(FILES scientific.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
else()
  install(TARGETS libscientific LIBRARY DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR})
  install(FILES ${Scientific_C_H} DESTINATION include/scientific)
  install(FILES scientific.h DESTINATION include)
endif()

# Copy the compiled library into the python_binding library
# in order to create a python library with a static library inside.
set(src_dir "${CMAKE_BINARY_DIR}/src")
set(dst_dir "${CMAKE_SOURCE_DIR}/src/python_bindings/libscientific")

set(libname "libscientific")
# Define the custom command to copy the library
add_custom_command(
    TARGET libscientific
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        "${src_dir}/${libname}${CMAKE_SHARED_LIBRARY_SUFFIX}"
        "${dst_dir}/${libname}${CMAKE_SHARED_LIBRARY_SUFFIX}"
    COMMENT "Copying ${libname} to ${dst_dir}"
)

# Optionally, create a custom target to trigger the copy as part of the build
add_custom_target(copy_${libname}
    DEPENDS ${libname}
)
