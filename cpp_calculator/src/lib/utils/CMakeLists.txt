add_library(LibUtils parameters.cpp csv.cpp)

target_include_directories(LibUtils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)