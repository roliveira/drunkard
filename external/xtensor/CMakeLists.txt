cmake_minimum_required(VERSION 3.5)
project(xtensor_builder CXX)
include(ExternalProject)
find_package(Git REQUIRED)

ExternalProject_Add(
    xtensor
    DEPENDS xtl
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/xtensor
    GIT_REPOSITORY https://github.com/QuantStack/xtensor.git
    TIMEOUT 10
    UPDATE_COMMAND ${GIT_EXECUTABLE} pull
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR} -Dxtl_DIR=${XTL_LIB}/cmake/xtl
    LOG_DOWNLOAD ON
)

set(XTENSOR_MODULE_PATH 
    "${CMAKE_CURRENT_BINARY_DIR}"         
    CACHE INTERNAL 
    "${PROJECT_NAME}: Include Directories" 
    FORCE
)

set(XTENSOR_INCLUDES    
    "${CMAKE_CURRENT_BINARY_DIR}/include" 
    CACHE INTERNAL 
    "${PROJECT_NAME}: Include Directories" 
    FORCE
)

set(XTENSOR_LIB         
    "${CMAKE_CURRENT_BINARY_DIR}/lib"     
    CACHE INTERNAL 
    "${PROJECT_NAME}: Include Directories" 
    FORCE
)
