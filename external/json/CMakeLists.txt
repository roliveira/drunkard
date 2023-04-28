cmake_minimum_required(VERSION 3.5)
project(json_builder CXX)
include(ExternalProject)
find_package(Git REQUIRED)

ExternalProject_Add(
    json
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    TIMEOUT 10
    UPDATE_COMMAND ${GIT_EXECUTABLE} pull
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR} -DJSON_BuildTests=OFF
    LOG_DOWNLOAD ON
)

set(JSON_MODULE_PATH 
    "${CMAKE_CURRENT_BINARY_DIR}"         
    CACHE INTERNAL 
    "${PROJECT_NAME}: Include Directories" 
    FORCE
)

set(JSON_INCLUDES    
    "${CMAKE_CURRENT_BINARY_DIR}/include" 
    CACHE INTERNAL 
    "${PROJECT_NAME}: Include Directories" 
    FORCE
)

set(JSON_LIB         
    "${CMAKE_CURRENT_BINARY_DIR}/lib"     
    CACHE INTERNAL 
    "${PROJECT_NAME}: Include Directories" 
    FORCE
)

