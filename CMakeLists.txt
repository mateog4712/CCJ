cmake_minimum_required(VERSION 3.10)
project(CCJ)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Directory for parameters
set(PARAMS_DIR "${CMAKE_SOURCE_DIR}/params")  # or CMAKE_INSTALL_PREFIX based path

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native") # conside -NDEBUG if I want to disable asserts for performance gain

SET(SOURCE
  src/W_final.cc 
  src/pseudo_loop.cc 
  src/pseudo_loop_gen.cc
  src/CCJ.cc
  src/cmdline.cc
  src/s_energy_matrix.cc
  src/part_func.cc
  src/dot_plot.cc
  src/stoch_backtrack.cc
)

set(vienna_SOURCE
    src/ViennaRNA/params/default.cc
    src/ViennaRNA/params/params.cc
    src/ViennaRNA/params/io.cc
    src/ViennaRNA/model.cc 
    src/ViennaRNA/utils.cc
)

add_library(RNA STATIC ${vienna_SOURCE})
# configure_file(RNA.pc.in RNA.pc @ONLY)

target_include_directories(RNA PRIVATE .)

add_executable(CCJ ${SOURCE})

target_link_libraries(CCJ PRIVATE RNA)

# Pass it as a compile definition
target_compile_definitions(CCJ PRIVATE 
    PARAMS_DIR="${PARAMS_DIR}"
)

# Install HFold into the bin directory
install(TARGETS CCJ RUNTIME DESTINATION bin)

include_directories(src)
