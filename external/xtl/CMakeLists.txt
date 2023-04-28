cmake_minimum_required(VERSION 3.5)
project(xtl_builder CXX)
include(ExternalProject)
find_package(Git REQUIRED)

ExternalProject_Add(
    xtl
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/xtl
    GIT_REPOSITORY https://github.com/QuantStack/xtl.git
    TIMEOUT 10
    UPDATE_COMMAND ${GIT_EXECUTABLE} pull
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}
    LOG_DOWNLOAD ON
)

set(XTL_MODULE_PATH 
    "${CMAKE_CURRENT_BINARY_DIR}"         
    CACHE INTERNAL 
    "${PROJECT_NAME}: Include Directories" 
    FORCE
)

set(XTL_INCLUDES    
    "${CMAKE_CURRENT_BINARY_DIR}/include" 
    CACHE INTERNAL 
    "${PROJECT_NAME}: Include Directories" 
    FORCE
)

set(XTL_LIB         
    "${CMAKE_CURRENT_BINARY_DIR}/lib"     
    CACHE INTERNAL 
    "${PROJECT_NAME}: Include Directories" 
    FORCE
)
