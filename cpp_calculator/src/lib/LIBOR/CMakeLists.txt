add_subdirectory(forward)

add_library(LibLIBOR ${FORWARD_SOURCES})
target_link_libraries(LibLIBOR LibProcess LibMath LibAnalytical)
target_include_directories(LibLIBOR PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
