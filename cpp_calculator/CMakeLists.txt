cmake_minimum_required(VERSION 3.10)

project(InterestRate)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# find_package(GSL REQUIRED)
set(CMAKE_CXX_COMPILER "icpx")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG -fp-model=fast -fiopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -traceback")

add_subdirectory(src)

enable_testing()
add_subdirectory(test)

