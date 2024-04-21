add_executable(test_LFM test_LFM.cpp)
target_link_libraries(test_LFM LibLIBOR gtest_main)

add_executable(test_SV test_SV.cpp)
target_link_libraries(test_SV LibLIBOR gtest_main)