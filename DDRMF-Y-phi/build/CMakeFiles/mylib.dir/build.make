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
CMAKE_SOURCE_DIR = "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/build"

# Include any dependencies generated for this target.
include CMakeFiles/mylib.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/mylib.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/mylib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mylib.dir/flags.make

CMakeFiles/mylib.dir/lib/cmd-progress/progress.f90.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/lib/cmd-progress/progress.f90.o: /home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/cmd-progress/progress.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/mylib.dir/lib/cmd-progress/progress.f90.o"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/cmd-progress/progress.f90" -o CMakeFiles/mylib.dir/lib/cmd-progress/progress.f90.o

CMakeFiles/mylib.dir/lib/cmd-progress/progress.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/mylib.dir/lib/cmd-progress/progress.f90.i"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/cmd-progress/progress.f90" > CMakeFiles/mylib.dir/lib/cmd-progress/progress.f90.i

CMakeFiles/mylib.dir/lib/cmd-progress/progress.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/mylib.dir/lib/cmd-progress/progress.f90.s"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/cmd-progress/progress.f90" -o CMakeFiles/mylib.dir/lib/cmd-progress/progress.f90.s

CMakeFiles/mylib.dir/lib/gobal-constants/gobal_constants.f90.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/lib/gobal-constants/gobal_constants.f90.o: /home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/gobal-constants/gobal_constants.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object CMakeFiles/mylib.dir/lib/gobal-constants/gobal_constants.f90.o"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/gobal-constants/gobal_constants.f90" -o CMakeFiles/mylib.dir/lib/gobal-constants/gobal_constants.f90.o

CMakeFiles/mylib.dir/lib/gobal-constants/gobal_constants.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/mylib.dir/lib/gobal-constants/gobal_constants.f90.i"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/gobal-constants/gobal_constants.f90" > CMakeFiles/mylib.dir/lib/gobal-constants/gobal_constants.f90.i

CMakeFiles/mylib.dir/lib/gobal-constants/gobal_constants.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/mylib.dir/lib/gobal-constants/gobal_constants.f90.s"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/gobal-constants/gobal_constants.f90" -o CMakeFiles/mylib.dir/lib/gobal-constants/gobal_constants.f90.s

CMakeFiles/mylib.dir/lib/mean-field-theory/MFT.f90.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/lib/mean-field-theory/MFT.f90.o: /home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/mean-field-theory/MFT.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building Fortran object CMakeFiles/mylib.dir/lib/mean-field-theory/MFT.f90.o"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/mean-field-theory/MFT.f90" -o CMakeFiles/mylib.dir/lib/mean-field-theory/MFT.f90.o

CMakeFiles/mylib.dir/lib/mean-field-theory/MFT.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/mylib.dir/lib/mean-field-theory/MFT.f90.i"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/mean-field-theory/MFT.f90" > CMakeFiles/mylib.dir/lib/mean-field-theory/MFT.f90.i

CMakeFiles/mylib.dir/lib/mean-field-theory/MFT.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/mylib.dir/lib/mean-field-theory/MFT.f90.s"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/mean-field-theory/MFT.f90" -o CMakeFiles/mylib.dir/lib/mean-field-theory/MFT.f90.s

CMakeFiles/mylib.dir/lib/model/choose_module.f90.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/lib/model/choose_module.f90.o: /home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/model/choose_module.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Building Fortran object CMakeFiles/mylib.dir/lib/model/choose_module.f90.o"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/model/choose_module.f90" -o CMakeFiles/mylib.dir/lib/model/choose_module.f90.o

CMakeFiles/mylib.dir/lib/model/choose_module.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/mylib.dir/lib/model/choose_module.f90.i"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/model/choose_module.f90" > CMakeFiles/mylib.dir/lib/model/choose_module.f90.i

CMakeFiles/mylib.dir/lib/model/choose_module.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/mylib.dir/lib/model/choose_module.f90.s"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/model/choose_module.f90" -o CMakeFiles/mylib.dir/lib/model/choose_module.f90.s

CMakeFiles/mylib.dir/lib/tov/TOV.f90.o: CMakeFiles/mylib.dir/flags.make
CMakeFiles/mylib.dir/lib/tov/TOV.f90.o: /home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/tov/TOV.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_5) "Building Fortran object CMakeFiles/mylib.dir/lib/tov/TOV.f90.o"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/tov/TOV.f90" -o CMakeFiles/mylib.dir/lib/tov/TOV.f90.o

CMakeFiles/mylib.dir/lib/tov/TOV.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/mylib.dir/lib/tov/TOV.f90.i"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/tov/TOV.f90" > CMakeFiles/mylib.dir/lib/tov/TOV.f90.i

CMakeFiles/mylib.dir/lib/tov/TOV.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/mylib.dir/lib/tov/TOV.f90.s"
	/usr/bin/f95 $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/lib/tov/TOV.f90" -o CMakeFiles/mylib.dir/lib/tov/TOV.f90.s

# Object files for target mylib
mylib_OBJECTS = \
"CMakeFiles/mylib.dir/lib/cmd-progress/progress.f90.o" \
"CMakeFiles/mylib.dir/lib/gobal-constants/gobal_constants.f90.o" \
"CMakeFiles/mylib.dir/lib/mean-field-theory/MFT.f90.o" \
"CMakeFiles/mylib.dir/lib/model/choose_module.f90.o" \
"CMakeFiles/mylib.dir/lib/tov/TOV.f90.o"

# External object files for target mylib
mylib_EXTERNAL_OBJECTS =

lib/libmylib.so: CMakeFiles/mylib.dir/lib/cmd-progress/progress.f90.o
lib/libmylib.so: CMakeFiles/mylib.dir/lib/gobal-constants/gobal_constants.f90.o
lib/libmylib.so: CMakeFiles/mylib.dir/lib/mean-field-theory/MFT.f90.o
lib/libmylib.so: CMakeFiles/mylib.dir/lib/model/choose_module.f90.o
lib/libmylib.so: CMakeFiles/mylib.dir/lib/tov/TOV.f90.o
lib/libmylib.so: CMakeFiles/mylib.dir/build.make
lib/libmylib.so: CMakeFiles/mylib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/build/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_6) "Linking Fortran shared library lib/libmylib.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mylib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mylib.dir/build: lib/libmylib.so
.PHONY : CMakeFiles/mylib.dir/build

CMakeFiles/mylib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mylib.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mylib.dir/clean

CMakeFiles/mylib.dir/depend:
	cd "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/build" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi" "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi" "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/build" "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/build" "/home/marigold/Documents/workspace/RMF&DDRMF/DDRMF-Y-phi/build/CMakeFiles/mylib.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/mylib.dir/depend

