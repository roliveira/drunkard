
cmake_minimum_required(VERSION 3.5)
project(drunkard)

# Set c++ version
set(CMAKE_CXX_STANDARD 14)

# Set the cmake module path of the project
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Set the project directory names
set(DRUNKARD_DOCS_DIR     ${CMAKE_SOURCE_DIR}/docs    )
set(DRUNKARD_EXTERNAL_DIR ${CMAKE_SOURCE_DIR}/external)
set(DRUNKARD_INCLUDE_DIR  ${CMAKE_SOURCE_DIR}/include )
set(DRUNKARD_MODULES_DIR  ${CMAKE_SOURCE_DIR}/modules )
set(DRUNKARD_TESTS_DIR    ${CMAKE_SOURCE_DIR}/tests   )

# Version parsing
file(READ "${DRUNKARD_INCLUDE_DIR}/config.hpp" ver)

string(REGEX MATCH "DRUNKARD_VERSION_MAJOR ([0-9]*)" _ ${ver})
set(DRUNKARD_VERSION_MAJOR ${CMAKE_MATCH_1})
string(REGEX MATCH "DRUNKARD_VERSION_MINOR ([0-9]*)" _ ${ver})
set(DRUNKARD_VERSION_MINOR ${CMAKE_MATCH_1})
string(REGEX MATCH "DRUNKARD_VERSION_PATCH ([0-9]*)" _ ${ver})
set(DRUNKARD_VERSION_PATCH ${CMAKE_MATCH_1})

set(${PROJECT_NAME}_VERSION
    ${DRUNKARD_VERSION_MAJOR}.${DRUNKARD_VERSION_MINOR}.${DRUNKARD_VERSION_PATCH}
)

message(STATUS "Building ${PROJECT_NAME} v${${PROJECT_NAME}_VERSION}")

# Define build components for drunkard
option(BUILD_DOCS     "Build documentation."        OFF)
option(BUILD_EXTERNAL "Build external dependecies." ON )
option(BUILD_MODULES  "Build modules."              ON )
option(BUILD_TESTS    "Build tests."                OFF)

find_package(Boost     REQUIRED)
find_package(VTK       REQUIRED)
find_package(PETSc     REQUIRED)
find_package(MPI       REQUIRED)
find_package(Reaktoro  REQUIRED)
find_package(ThermoFun REQUIRED)
find_package(xtl       REQUIRED)
find_package(xtensor   REQUIRED)
find_package(nlohmann_json REQUIRED)

# Set the include directories
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${BOOST_INCLUDE_DIR})
include_directories(${PETSC_INCLUDES})
include_directories(${DRUNKARD_INCLUDE_DIR})
include_directories(${MPI_CXX_INCLUDE_PATH})

add_definitions(${PETSC_DEFINITIONS})


# Set include files
include(${VTK_USE_FILE})

# Add library 
add_subdirectory(include)

# Build the external requirements
add_subdirectory(external)
include_directories(
    ${JSON_INCLUDES}
    ${xtl_INCLUDE_DIRS}
    ${xtensor_INCLUDE_DIRS}
    ${DLIB_INCLUDES}
)


# Build the applications
if(BUILD_MODULES)
    add_subdirectory(modules)
else()
    add_subdirectory(modules EXCLUDE_FROM_ALL)
endif()


# Add target "modules" for manual building of modules, as `make modules`, if BUILD_MODULES is OFF
add_custom_target(modules
    COMMAND ${CMAKE_MAKE_PROGRAM}
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
)

