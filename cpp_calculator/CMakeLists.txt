cmake_minimum_required(VERSION 3.18)

if(CUDA_ENABLED)
    project(InterestRate LANGUAGES CXX CUDA)
    set(CMAKE_CUDA_FLAGS_DEBUG "-G")
    set(CMAKE_CUDA_FLAGS_RELEASE "-O3")

else()
    project(InterestRate LANGUAGES CXX)

endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(ICPX_USED)
    set(CMAKE_CXX_COMPILER "icpx")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG -fp-model=fast -fiopenmp")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -traceback")
endif()


add_subdirectory(src)

enable_testing()
add_subdirectory(test)
