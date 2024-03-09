add_executable(test_constant test_constant.cpp)
target_link_libraries(test_constant LibProcess)

add_executable(test_HoLee test_HoLee.cpp)
target_link_libraries(test_HoLee LibProcess)

add_executable(test_Vasicek test_Vasicek.cpp)
target_link_libraries(test_Vasicek LibProcess)