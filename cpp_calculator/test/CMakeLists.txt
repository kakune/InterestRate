include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)

add_subdirectory(analytical)
add_subdirectory(math)
add_subdirectory(process)
add_subdirectory(short_rate)
add_subdirectory(LIBOR)