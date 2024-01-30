add_executable(test_parameters test_parameters.cpp)
target_link_libraries(test_parameters LibUtils)

add_executable(test_findroot_1D test_findroot_1D.cpp)
target_link_libraries(test_findroot_1D LibUtils)

add_executable(test_Black_Scholes test_Black_Scholes.cpp)
target_link_libraries(test_Black_Scholes LibAnalytical)

add_executable(test_random test_random.cpp)
target_link_libraries(test_random LibProcess)

add_executable(test_interpolate_1D test_interpolate_1D.cpp)
target_link_libraries(test_interpolate_1D LibMath)

add_executable(test_special_functions test_special_functions.cpp)
target_link_libraries(test_special_functions LibMath)

add_executable(test_short_rate test_short_rate.cpp)
target_link_libraries(test_short_rate LibProcess)

add_executable(test_asset test_asset.cpp)
target_link_libraries(test_asset LibProcess)

# target_link_libraries(test_executable main_executable)

# add_test(NAME Test COMMAND test_parameters)