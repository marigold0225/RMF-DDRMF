# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/marigold/Documents/workspace/RMF&DDRMF/RMF-Y"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/marigold/Documents/workspace/RMF&DDRMF/RMF-Y/build"

# Include any dependencies generated for this target.
include CMakeFiles/run.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/run.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/run.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/run.dir/flags.make

CMakeFiles/run.dir/src/main.f90.o: CMakeFiles/run.dir/flags.make
CMakeFiles/run.dir/src/main.f90.o: /home/marigold/Documents/workspace/RMF&DDRMF/RMF-Y/src/main.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/marigold/Documents/workspace/RMF&DDRMF/RMF-Y/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/run.dir/src/main.f90.o"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c "/home/marigold/Documents/workspace/RMF&DDRMF/RMF-Y/src/main.f90" -o CMakeFiles/run.dir/src/main.f90.o

CMakeFiles/run.dir/src/main.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/run.dir/src/main.f90.i"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E "/home/marigold/Documents/workspace/RMF&DDRMF/RMF-Y/src/main.f90" > CMakeFiles/run.dir/src/main.f90.i

CMakeFiles/run.dir/src/main.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/run.dir/src/main.f90.s"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S "/home/marigold/Documents/workspace/RMF&DDRMF/RMF-Y/src/main.f90" -o CMakeFiles/run.dir/src/main.f90.s

# Object files for target run
run_OBJECTS = \
"CMakeFiles/run.dir/src/main.f90.o"

# External object files for target run
run_EXTERNAL_OBJECTS =

bin/run: CMakeFiles/run.dir/src/main.f90.o
bin/run: CMakeFiles/run.dir/build.make
bin/run: lib/libmylib.so
bin/run: CMakeFiles/run.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/marigold/Documents/workspace/RMF&DDRMF/RMF-Y/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable bin/run"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/run.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/run.dir/build: bin/run
.PHONY : CMakeFiles/run.dir/build

CMakeFiles/run.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/run.dir/cmake_clean.cmake
.PHONY : CMakeFiles/run.dir/clean

CMakeFiles/run.dir/depend:
	cd "/home/marigold/Documents/workspace/RMF&DDRMF/RMF-Y/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/marigold/Documents/workspace/RMF&DDRMF/RMF-Y" "/home/marigold/Documents/workspace/RMF&DDRMF/RMF-Y" "/home/marigold/Documents/workspace/RMF&DDRMF/RMF-Y/build" "/home/marigold/Documents/workspace/RMF&DDRMF/RMF-Y/build" "/home/marigold/Documents/workspace/RMF&DDRMF/RMF-Y/build/CMakeFiles/run.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/run.dir/depend

