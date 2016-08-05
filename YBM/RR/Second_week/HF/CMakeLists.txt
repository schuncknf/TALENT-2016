# CMakeLists.txt, JGC, 03/08/2016

# If you install cmake this code will create a makefile etc for your code.
# The best way to use it is to make a "build" directory, change to this
# directory then type "cmake .." this command will call the CMakeLists.txt
# file in the parent directory, this will create a Makefile to compile your
# code as well as copy all necessary files to run the code in your build 
# directory.

# Using this we don't need to recompile each time someone pushes and destroys
# the executable on your machine - we can easily run the code in our own 
# build directories

#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required (VERSION 2.6)
project (Esercizio)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
file(GLOB sources ${PROJECT_SOURCE_DIR}../*.C)
file(GLOB headers ${PROJECT_SOURCE_DIR}../*.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the main, source and header files.
#
add_executable(Esercizio Esercizio.C ${sources} ${headers})

#----------------------------------------------------------------------------
# Copy all files necessary to run the code to your build directory.
#
set(FILES
  Oxy_density_proton.txt
  Oxy_density_neutron.txt
  )

foreach(_script ${FILES})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
#----------------------------------------------------------------------------
# End of CMakeLists.txt file
#
