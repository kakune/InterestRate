find_package(GTest REQUIRED)
include(GoogleTest)

add_executable(test_parameters test_parameters.cpp)
target_link_libraries(test_parameters LibUtils)

add_executable(test_findroot_1d test_findroot_1d.cpp)
target_link_libraries(test_findroot_1d LibUtils)

add_executable(test_Black_Scholes test_Black_Scholes.cpp)
target_link_libraries(test_Black_Scholes LibAnalytical)

add_executable(test_random test_random.cpp)
target_link_libraries(test_random LibProcess)

add_executable(test_interpolate_1d test_interpolate_1d.cpp)
target_link_libraries(test_interpolate_1d LibMath GTest::GTest GTest::Main)
gtest_add_tests(TARGET test_interpolate_1d)

add_executable(test_interpolate_multi test_interpolate_multi.cpp)
target_link_libraries(test_interpolate_multi LibMath GTest::GTest GTest::Main)
gtest_add_tests(TARGET test_interpolate_multi)

add_executable(test_matrix test_matrix.cpp)
target_link_libraries(test_matrix LibMath GTest::GTest GTest::Main)
gtest_add_tests(TARGET test_matrix)

add_executable(test_special_functions test_special_functions.cpp)
target_link_libraries(test_special_functions LibMath)

add_executable(test_asset test_asset.cpp)
target_link_libraries(test_asset LibProcess)

add_executable(test_market test_market.cpp)
target_link_libraries(test_market LibProcess GTest::GTest GTest::Main)
gtest_add_tests(TARGET test_market)

add_executable(test_ODE test_ODE.cpp)
target_link_libraries(test_ODE LibMath GTest::GTest GTest::Main)
gtest_add_tests(TARGET test_ODE)

add_subdirectory(short_rate)


if(CUDA_ENABLED)
    # add_executable(test_random_cuda test_random_cuda.cpp)
    # target_link_libraries(test_random_cuda LibProcess LibProcessCUDA)

    add_executable(test_asset_cuda test_asset_cuda.cpp)
    target_link_libraries(test_asset_cuda LibProcess LibProcessCUDA)
endif()

# target_link_libraries(test_executable main_executable)

# add_test(NAME Test COMMAND test_parameters)