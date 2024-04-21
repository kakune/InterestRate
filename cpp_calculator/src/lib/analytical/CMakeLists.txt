add_library(LibAnalytical Black_Scholes.cpp Black76.cpp)

target_link_libraries(LibAnalytical LibMath)
target_include_directories(LibAnalytical PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)