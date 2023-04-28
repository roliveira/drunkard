cmake_minimum_required(VERSION 3.5)
project(dlib_builder CXX)
include(ExternalProject)
find_package(Git REQUIRED)

ExternalProject_Add(
    dlib
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/dlib
    GIT_REPOSITORY "https://github.com/davisking/dlib"
    # SOURCE_DIR ${CMAKE_SOURCE_DIR}/externals/dlib
    GIT_TAG master
    TIMEOUT 10
    # DOWNLOAD_COMMAND ""
    UPDATE_COMMAND ${GIT_EXECUTABLE} pull
    CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        # -DLIB_INSTALL_DIR=${CMAKE_INSTALL_LIBDIR}
    LOG_DOWNLOAD ON
)

set(DLIB_MODULE_PATH 
    "${CMAKE_CURRENT_BINARY_DIR}"         
    CACHE INTERNAL 
    "${PROJECT_NAME}: Include Directories" 
    FORCE
)

set(DLIB_INCLUDES    
    "${CMAKE_CURRENT_BINARY_DIR}/include" 
    CACHE INTERNAL 
    "${PROJECT_NAME}: Include Directories" 
    FORCE
)

set(DLIB_LIB         
    "${CMAKE_CURRENT_BINARY_DIR}/lib"     
    CACHE INTERNAL 
    "${PROJECT_NAME}: Include Directories" 
    FORCE
)

set(dlib_LIBRARIES
    "${DLIB_LIB}/libdlib.a"
    CACHE INTERNAL
    "${PROJECT_NAME}: Library file" 
    FORCE
)