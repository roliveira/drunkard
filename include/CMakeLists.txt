cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(
	methods
	utils
	io
	geometry
	particle
)

set(CONFIG
	config.hpp
)

set(IO
    io/echo.hpp
	io/logger.hpp
	io/reader.hpp
    io/writer.hpp
    # io/parser.hpp
)

set(GEOMETRY
    geometry/geometry.hpp
)

set(PARTICLE
	particle/particle.hpp
)

set(METHODS
	methods/transport.hpp
	methods/continuous.hpp
	methods/flow.hpp
	methods/darcy.hpp
)

set(UTILS
	utils/index.hpp
	utils/interpolation.hpp
	utils/operations.hpp
	utils/quadrature.hpp
	utils/root.hpp
)

enable_testing()

# add_library(
# 	${PROJECT_NAME} 
# 	${CONFIG}
# 	${IO} 
# 	${GEOMETRY}
# 	${PARTICLE}
# 	${METHODS} 
# 	${UTILS}
# )

add_library(${PROJECT_NAME} INTERFACE)

# Build the external requirements
# add_dependencies(
# 	${PROJECT_NAME} 
# 	json
# 	xtl
# 	xtensor
# )

target_link_libraries(
	${PROJECT_NAME} 
	INTERFACE 
	# ${Boost_LIBRARIES} 
	${VTK_LIBRARIES}
)

target_include_directories(
	${PROJECT_NAME} 
	INTERFACE 
	methods
	utils
	io
	geometry
	particle
)

# install(TARGETS ${PROJECT_NAME} DESTINATION lib COMPONENT libraries)
