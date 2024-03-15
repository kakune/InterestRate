add_library(LibMath special_functions.cpp interpolate_1d.cpp interpolate_multi.cpp integral_1d.cpp matrix.cpp)

target_include_directories(LibMath PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)