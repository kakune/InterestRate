include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)

add_subdirectory(analytical)
add_subdirectory(math)
# add_subdirectory(process)
add_subdirectory(short_rate)
add_subdirectory(LIBOR)

add_executable(test_asset test_asset.cpp)
target_link_libraries(test_asset LibProcess)

if(CUDA_ENABLED)
    add_executable(test_asset_cuda test_asset_cuda.cpp)
    target_link_libraries(test_asset_cuda LibProcess LibProcessCUDA)
endif()