
add_executable(ctrw ctrw.cpp)
add_dependencies(ctrw ${PROJECT_NAME})
target_link_libraries(
    ctrw 
    ${Boost_LIBRARIES} 
    ${VTK_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${PETSC_LIBRARIES}
    ${REAKTORO_LIBRARIES}
    Reaktoro::Reaktoro
)

add_executable(rctrw rctrw.cpp)
add_dependencies(rctrw ${PROJECT_NAME})
target_link_libraries(
    rctrw
    ${Boost_LIBRARIES} 
    ${VTK_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${PETSC_LIBRARIES}
    ${REAKTORO_LIBRARIES}
    Reaktoro::Reaktoro
)
