# Ian Wark July 6 2017

SET(SOURCE
  W_final.cpp 
  pseudo_loop.cpp 
  V_final.cpp 
  VM_final.cpp
  h_common.cpp 
  CCJ.cpp
  param_path.cc
)

# pass SIMFOLD_HOME to binary
add_definitions("-DSIMFOLD_HOME=\"${SIMFOLD_HOME}\"")

# set path for executable to be created
# set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# create executable
add_executable(Original_CCJ ${SOURCE})
set_target_properties(Original_CCJ PROPERTIES OUTPUT_NAME CCJ)

# Install Original_CCJ into the bin directory (Needed for conda)
install(TARGETS Original_CCJ RUNTIME DESTINATION bin)

# enable required C++11 features
target_compile_features(Original_CCJ PRIVATE cxx_nullptr)

# link to simfold
target_link_libraries(Original_CCJ LINK_PUBLIC simfold)
