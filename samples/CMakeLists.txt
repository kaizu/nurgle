add_executable(main main.cpp)
target_link_libraries(main nurgle)
install(TARGETS main DESTINATION bin)
