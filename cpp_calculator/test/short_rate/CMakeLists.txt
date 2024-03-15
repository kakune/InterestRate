add_executable(test_constant test_constant.cpp)
target_link_libraries(test_constant LibProcess gtest_main)

add_executable(test_HoLee test_HoLee.cpp)
target_link_libraries(test_HoLee LibProcess gtest_main)

add_executable(test_Vasicek test_Vasicek.cpp)
target_link_libraries(test_Vasicek LibProcess gtest_main)

add_executable(test_GSR test_GSR.cpp)
target_link_libraries(test_GSR LibProcess gtest_main)

add_executable(test_constant_affine test_constant_affine.cpp)
target_link_libraries(test_constant_affine LibProcess gtest_main)
