add_library(LibMath special_functions.cpp interpolate_1d.cpp interpolate_multi.cpp integral_1d.cpp matrix.cpp)

if(MKL_FOUND)
    set(MKL_LIBRARIES -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl)
    target_include_directories(LibMath PRIVATE ${MKL_INCLUDE_DIRS})
    target_link_libraries(LibMath PRIVATE ${MKL_LIBRARIES})
else()
    set(LAPACK_LIBRARIES -lblas -llapack -llapacke)
    target_include_directories(LibMath PRIVATE ${LAPACK_INCLUDE_DIRS})
    target_link_libraries(LibMath PRIVATE ${LAPACK_LIBRARIES})
endif()

target_include_directories(LibMath PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)