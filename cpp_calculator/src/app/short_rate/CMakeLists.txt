add_executable(calc_with_param calc_with_param.cpp)
target_link_libraries(calc_with_param LibUtils LibProcess)

add_executable(calc_with_market calc_with_market.cpp)
target_link_libraries(calc_with_market LibUtils LibProcess)