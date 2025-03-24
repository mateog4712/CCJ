# Ian Wark July 6 2017

SET(SOURCE
  W_final.cpp 
  pseudo_loop.cpp 
  V_final.cpp 
  VM_final.cpp
  h_common.cpp 
  trace_arrow.cpp
  candidate_list.cpp
  CCJ.cpp 
)

# pass SIMFOLD_HOME to binary
add_definitions("-DSIMFOLD_HOME=\"${SIMFOLD_HOME}\"")

# set path for executable to be created
# set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# create executable
add_executable(Sparse_CCJ ${SOURCE})
set_target_properties(Sparse_CCJ PROPERTIES OUTPUT_NAME CCJ)

# Install Sparse_CCJ into the bin directory (Needed for conda)
install(TARGETS Sparse_CCJ RUNTIME DESTINATION bin)

# enable required C++11 features
target_compile_features(Sparse_CCJ PRIVATE cxx_nullptr)

# link to simfold
target_link_libraries(Sparse_CCJ LINK_PUBLIC simfold)
