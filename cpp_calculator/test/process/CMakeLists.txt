# add_executable(test_random test_random.cpp)
# target_link_libraries(test_random LibProcess)

# add_executable(test_market test_market.cpp)
# target_link_libraries(test_market LibProcess gtest_main)


if(CUDA_ENABLED)
    add_executable(test_random_cuda test_random_cuda.cu)
    target_link_libraries(test_random_cuda LibProcess)
endif()
