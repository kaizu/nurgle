file(GLOB HPP_FILES *.hpp)
file(GLOB CPP_FILES *.cpp)

add_library(nurgle STATIC ${CPP_FILES} ${HPP_FILES})
# add_library(nurgle SHARED ${CPP_FILES} ${HPP_FILES})
# install(TARGETS nurgle DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION "include/nurgle")
