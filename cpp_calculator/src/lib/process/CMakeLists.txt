add_library(LibProcess random.cpp asset.cpp market_data.cpp model_data.cpp)
target_link_libraries(LibProcess LibMath LibAnalytical)
target_include_directories(LibProcess PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
