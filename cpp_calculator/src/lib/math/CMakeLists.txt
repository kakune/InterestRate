add_library(LibMath special_functions.cpp interpolate_1d.cpp integral_1d.cpp)

target_include_directories(LibMath PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)