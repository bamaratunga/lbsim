# Definition of the minimum required cmake Version
cmake_minimum_required(VERSION 3.0)
# Definition of the Project
# Later you can access the project variable like ${CFDLAB_SOURCE_DIR}
project(LBSIM VERSION 1.0 LANGUAGES CXX CUDA)

# Definition of the C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)

set(CMAKE_CUDA_ARCHITECTURES off)

# Find packages
# find_package(MPI)
# Find a package with different components e.g. BOOST
# find_package(Boost COMPONENTS filesystem REQUIRED)
# Require a package
# VTK Library
find_package(VTK REQUIRED)
if (NOT VTK_FOUND)
  message("VTK not found")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
include(${VTK_USE_FILE})

# Filesystem

set(gpp9 True)

if(gpp9)
    message("g++ Version is 9 or higher")
    add_compile_definitions(gpp9)
else()
    message("g++ Version is lower than Version 9")
    if (NOT APPLE)
      # Link Filesystem Libraries
      link_libraries(stdc++fs)
    endif()
endif()

# Add include directory
include_directories(include)

if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX /usr/local )
endif()

# Creating the executable of our project and the required dependencies
# the executable is called lbsim
add_executable(lbsim_cpp src/Case.cpp src/Boundary.cpp src/Initialize.cpp src/Processes.cpp)
add_executable(lbsim_cuda src/Lbmcuda.cu)

## Compiler options
# Turn on warnings
target_compile_options(lbsim_cpp PRIVATE -Wall -Wextra)
target_compile_options(lbsim_cuda PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-ftz=false -prec-div=true -prec-sqrt=true>)
target_compile_features(lbsim_cuda PUBLIC cxx_std_11)

# if you use external libraries you have to link them like
# target_link_libraries(lbsim_cpp PRIVATE MPI::MPI_CXX)
target_link_libraries(lbsim_cpp PRIVATE ${VTK_LIBRARIES})
target_link_libraries(lbsim_cuda PRIVATE ${VTK_LIBRARIES})

set_target_properties(lbsim_cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

install(TARGETS lbsim_cpp DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
install(TARGETS lbsim_cuda DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
# If you write tests, you can include your subdirectory (in this case tests) as done here
# Testing
