add_executable(test_parameters test_parameters.cpp)
target_link_libraries(test_parameters LibUtils)

add_executable(test_findroot_1D test_findroot_1D.cpp)
target_link_libraries(test_findroot_1D LibUtils)

add_executable(test_Black_Scholes test_Black_Scholes.cpp)
target_link_libraries(test_Black_Scholes LibUtils)

# target_link_libraries(test_executable main_executable)

# add_test(NAME Test COMMAND test_parameters)