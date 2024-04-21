add_executable(SR_param SR_param.cpp)
target_link_libraries(SR_param LibUtils LibShortRate LibApp)

add_executable(SR_market SR_market.cpp)
target_link_libraries(SR_market LibUtils LibShortRate LibApp)