add_subdirectory(lib)
add_subdirectory(short_rate)
add_subdirectory(LIBOR)

target_include_directories(LibAnalytical PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)