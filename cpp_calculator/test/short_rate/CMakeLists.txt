find_package(GTest REQUIRED)
include(GoogleTest)

add_executable(test_constant test_constant.cpp)
target_link_libraries(test_constant LibProcess GTest::GTest GTest::Main)
gtest_add_tests(TARGET test_constant)

add_executable(test_HoLee test_HoLee.cpp)
target_link_libraries(test_HoLee LibProcess GTest::GTest GTest::Main)
gtest_add_tests(TARGET test_HoLee)

add_executable(test_Vasicek test_Vasicek.cpp)
target_link_libraries(test_Vasicek LibProcess GTest::GTest GTest::Main)
gtest_add_tests(TARGET test_Vasicek)

add_executable(test_GSR test_GSR.cpp)
target_link_libraries(test_GSR LibProcess GTest::GTest GTest::Main)
gtest_add_tests(TARGET test_GSR)

add_executable(test_constant_affine test_constant_affine.cpp)
target_link_libraries(test_constant_affine LibProcess GTest::GTest GTest::Main)
gtest_add_tests(TARGET test_constant_affine)