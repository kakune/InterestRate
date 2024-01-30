add_library(LibProcess random.cpp short_rate.cpp)

target_link_libraries(LibProcess LibMath)
target_include_directories(LibProcess PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)