add_subdirectory(one-factor)
add_subdirectory(multi-factor)
add_subdirectory(PDE)

add_library(LibShortRate ${ONE_FACTOR_SOURCES} ${MULTI_FACTOR_SOURCES} ${PDE_SOURCES})
target_link_libraries(LibShortRate LibProcess LibMath LibAnalytical)
target_include_directories(LibShortRate PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
