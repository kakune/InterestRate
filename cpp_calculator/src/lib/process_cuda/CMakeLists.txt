# ライブラリの追加
add_library(LibProcessCUDA asset.cu asset.cpp)


enable_language(CUDA)
target_include_directories(LibProcessCUDA PUBLIC /usr/local/cuda/include)

set_target_properties(LibProcessCUDA PROPERTIES CUDA_ARCHITECTURES "75")
target_compile_features(LibProcessCUDA PUBLIC cuda_std_17)
target_link_libraries(LibProcessCUDA PRIVATE LibMath LibAnalytical)
target_include_directories(LibProcessCUDA PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)