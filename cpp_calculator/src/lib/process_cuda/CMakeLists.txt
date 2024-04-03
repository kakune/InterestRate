add_library(LibProcessCUDA asset.cu asset.cpp)

target_include_directories(LibProcessCUDA PUBLIC /usr/local/cuda/include)

target_link_libraries(LibProcessCUDA PRIVATE LibMath LibAnalytical)
target_include_directories(LibProcessCUDA PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)