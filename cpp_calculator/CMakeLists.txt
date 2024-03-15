cmake_minimum_required(VERSION 3.18)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)



if(CUDA_ENABLED)
    project(InterestRate LANGUAGES CXX CUDA)
    set(CMAKE_CUDA_FLAGS_DEBUG "-G")
    set(CMAKE_CUDA_FLAGS_RELEASE "-O3")
else()
    project(InterestRate LANGUAGES CXX)
endif()

find_package(MKL QUIET)
if(MKL_FOUND)
  find_package(OpenMP REQUIRED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
  add_definitions(-DNUSE_MKL)
  find_package(LAPACK REQUIRED)
endif()

add_subdirectory(src)

enable_testing()
add_subdirectory(test)


