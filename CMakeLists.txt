cmake_minimum_required(VERSION 2.8 FATAL_ERROR)



#=================== Set the example path ===========================
# Check the example directory to find more examples. By changing the 
# address below to the desired example, you set the compiler for it.
 
set(example_path ../examples/2d/von_karman_vortex/main.f90)
#====================================================================

#=================== Set project name ===============================
# The project executable file is made by this name.

set (project_name vortex)
#====================================================================


# The project name is applied.
project(${project_name})

# Fortran language is activated.
enable_language(Fortran)

# You can set the build type in terminal as: 
#$ cmake -DCMAKE_BUILD_TYPE=Debug ..
# Otherwise "Release" is the default build type.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Gnu compiler flags
# -w                             switchs off all the warnings.
# -wall                          switchs on  all the warnings.
# -march=native                  optimization for the local machine
# -O0, -O1, -O2, -O3, -Ofast     overall optimization levels
# -g0,g1,g2,g3					 levels of debug information

# This flags are applied to all bluid types
set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -w")
# Debug flags
set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
# Release flags
set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native")

# Fortran module files are placed here.
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# Copying lib_vtk module to the modules of the bindary directory.
# Make sure lib_vtk_io.mod file is placed in 
# "${CMAKE_SOURCE_DIR}/extlib" directory (see the guide file). 
configure_file(extlib/lib_vtk_io.mod modules/lib_vtk_io.mod COPYONLY)

# Create vtk directory needed to collect VTK files.
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/vtk)

# All the subdirectories having CMakeLists.txt are listed.
add_subdirectory(tools)
add_subdirectory(collision/pseudopotential)
add_subdirectory(collision/boundaryconditions)
add_subdirectory(collision/base)
add_subdirectory(imgproc)
add_subdirectory(initialization)
add_subdirectory(lattice)
add_subdirectory(lists/base)
add_subdirectory(lists/kinds)
add_subdirectory(lists/parts)
add_subdirectory(stream)
add_subdirectory(writeout)


# Adding library dependencies
# Because "lattice" is dependent to "tools" at the first line,
# lists_base is dependent to tools as well in the second line.
add_dependencies(lattice tools)
add_dependencies(lists_base lattice)
add_dependencies(stream lists_base)
add_dependencies(initialization lists_base)
add_dependencies(lists_kinds stream collision_base)
add_dependencies(lists_parts lists_kinds)
add_dependencies(writeout lists_parts)
add_dependencies(imgproc lists_parts)
add_dependencies(collision_base lists_base)
add_dependencies(collision_BC collision_base)
add_dependencies(collision_PP collision_BC lists_parts)


# The executable is built with the name of the project from the
# example path provided above.
add_executable(
    ${project_name}
    ${example_path}
    )

# Libraries needed by the executable are listed. Make sure lib_vtk
# library files are placed in "${CMAKE_SOURCE_DIR}/extlib", see the
# guide file. 
target_link_libraries(
	${project_name} 
	collision_PP
	collision_BC
	collision_base
	imgproc
	initialization
	lattice
	lists_base
	lists_parts
	lists_kinds
	stream
	writeout
	tools
	${CMAKE_SOURCE_DIR}/extlib/lib_vtk_io.o
	${CMAKE_SOURCE_DIR}/extlib/ir_precision.o
	${CMAKE_SOURCE_DIR}/extlib/lib_base64.o
)

# Display general status
message("---------------Status--------------")
message("Project name = ${project_name}")
message("Example path = ${example_path}")
message("Compiler = ${CMAKE_Fortran_COMPILER_ID}")
message("Build type = ${CMAKE_BUILD_TYPE}")
message("-----------------------------------")


