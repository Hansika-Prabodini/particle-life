cmake_minimum_required(VERSION 3.16)
project(particle_life LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_library(particle_sim src/simulation.cpp)
target_include_directories(particle_sim PUBLIC src)
target_compile_options(particle_sim PRIVATE -Wall -Wextra -Wpedantic -O3 -ffp-contract=off)

add_executable(particle_life src/main.cpp)
target_link_libraries(particle_life PRIVATE particle_sim)
target_compile_options(particle_life PRIVATE -Wall -Wextra -Wpedantic -O3 -ffp-contract=off)

enable_testing()
add_executable(particle_life_tests tests/simulation_test.cpp)
target_link_libraries(particle_life_tests PRIVATE particle_sim)
target_compile_options(particle_life_tests PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME simulation COMMAND particle_life_tests)
