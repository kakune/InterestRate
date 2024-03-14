find_package(MKL)

if(MKL_FOUND)
  add_definitions(-DUSE_MKL)
  include_directories(${MKL_INCLUDE_DIRS})
  link_libraries(${MKL_LIBRARIES})
else()
  find_package(LAPACK REQUIRED)
  if(LAPACK_FOUND)
    include_directories(${LAPACK_INCLUDE_DIRS})
    link_libraries(${LAPACK_LIBRARIES})
  endif()
endif()

add_library(LibMath special_functions.cpp interpolate_1d.cpp interpolate_multi.cpp integral_1d.cpp matrix.cpp)

target_include_directories(LibMath PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)