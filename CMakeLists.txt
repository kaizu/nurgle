cmake_minimum_required(VERSION 2.6)
project(nurgle)

# set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)

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

if(DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -stdlib=libc++ -O2 -g -Wall")
else(DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -stdlib=libc++ -O3")
endif(DEBUG)

include_directories(.)

add_subdirectory(nurgle)
add_subdirectory(samples)
add_subdirectory(notebooks)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/outputs/timecourse.csv
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_BINARY_DIR}/outputs ${CMAKE_CURRENT_BINARY_DIR}/samples/main .
    DEPENDS main
    DEPENDS outputs
    )
add_custom_target(
    timecourse DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/outputs/timecourse.csv)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/outputs/timecourse.csv DESTINATION share/nurgle OPTIONAL)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/outputs/index.html
    COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_CURRENT_BINARY_DIR}/outputs jupyter nbconvert --config ${CMAKE_CURRENT_SOURCE_DIR}/notebooks/ipython_nbconvert_config.py --execute ${CMAKE_CURRENT_SOURCE_DIR}/notebooks/index.ipynb --output ${CMAKE_CURRENT_BINARY_DIR}/outputs/index.html
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/outputs/timecourse.csv
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/notebooks/index.ipynb
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/notebooks/ipython_nbconvert_config.py
    )
add_custom_target(
    html DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/outputs/index.html)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/outputs/index.html DESTINATION share/nurgle OPTIONAL)
