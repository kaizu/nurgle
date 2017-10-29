cmake_minimum_required(VERSION 2.6)
project(nurgle)

set(OUTPUTS
    ${CMAKE_CURRENT_BINARY_DIR}/outputs/metabolism.csv
    ${CMAKE_CURRENT_BINARY_DIR}/outputs/compounds.csv
    )

find_package(PythonInterp ${PY_VERSION})

add_custom_command(
    OUTPUT ${OUTPUTS}
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_SOURCE_DIR} ${PYTHON_EXECUTABLE} -m gen -o ${CMAKE_CURRENT_BINARY_DIR}/outputs
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gen/__main__.py
    )

add_custom_target(
    outputs ALL DEPENDS ${OUTPUTS})

install(FILES ${OUTPUTS} DESTINATION share/nurgle)

find_package(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -stdlib=libc++ -O3")
include_directories(.)

add_subdirectory(nurgle)
add_subdirectory(samples)
