add_subdirectory(utils)
add_subdirectory(process)
add_subdirectory(analytical)
add_subdirectory(math)

if(CUDA_ENABLED)
    add_subdirectory(process_cuda)
endif()
