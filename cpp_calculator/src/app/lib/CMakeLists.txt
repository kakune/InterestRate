add_library(LibApp terms_creator.cpp short_rate_creator.cpp LIBOR_forward_creator.cpp)
target_link_libraries(LibApp LibUtils LibMath LibProcess LibShortRate LibLIBOR)
target_include_directories(LibApp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})