add_subdirectory(short_rate)

add_library(LibProcess random.cpp asset.cpp market.cpp ${SHORT_RATE_SOURCES})

target_link_libraries(LibProcess LibMath LibAnalytical)
target_include_directories(LibProcess PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
