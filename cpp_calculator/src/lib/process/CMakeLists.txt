add_library(LibProcess random.cpp market_data.cpp)

if(CUDA_ENABLED)
    enable_language(CUDA)
    target_include_directories(LibProcess PUBLIC ${CUDAToolkit_INCLUDE_DIRS})
    target_link_libraries(LibProcess LibMath LibAnalytical CUDA::cudart)
    set_target_properties(LibProcess PROPERTIES CUDA_ARCHITECTURES "75")
    target_compile_features(LibProcess PUBLIC cuda_std_17)
else()
    target_link_libraries(LibProcess LibMath LibAnalytical)
endif()

target_include_directories(LibProcess PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
