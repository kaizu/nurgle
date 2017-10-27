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

add_subdirectory(src)
