add_subdirectory(short_rate_MC_one)
add_subdirectory(short_rate_PDE)

add_library(LibProcess random.cpp asset.cpp market_data.cpp ${SHORT_RATE_MC_SOURCES} ${SHORT_RATE_PDE_SOURCES})
target_link_libraries(LibProcess LibMath LibAnalytical)
target_include_directories(LibProcess PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
