# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.5.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.5.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/marco/projects/libscientific-0.7.4

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/marco/projects/libscientific-0.7.4/build

# Include any dependencies generated for this target.
include src/CMakeFiles/testoptimization.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/testoptimization.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/testoptimization.dir/flags.make

src/CMakeFiles/testoptimization.dir/testoptimization.c.o: src/CMakeFiles/testoptimization.dir/flags.make
src/CMakeFiles/testoptimization.dir/testoptimization.c.o: ../src/testoptimization.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/marco/projects/libscientific-0.7.4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/CMakeFiles/testoptimization.dir/testoptimization.c.o"
	cd /Users/marco/projects/libscientific-0.7.4/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/testoptimization.dir/testoptimization.c.o   -c /Users/marco/projects/libscientific-0.7.4/src/testoptimization.c

src/CMakeFiles/testoptimization.dir/testoptimization.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/testoptimization.dir/testoptimization.c.i"
	cd /Users/marco/projects/libscientific-0.7.4/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/marco/projects/libscientific-0.7.4/src/testoptimization.c > CMakeFiles/testoptimization.dir/testoptimization.c.i

src/CMakeFiles/testoptimization.dir/testoptimization.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/testoptimization.dir/testoptimization.c.s"
	cd /Users/marco/projects/libscientific-0.7.4/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/marco/projects/libscientific-0.7.4/src/testoptimization.c -o CMakeFiles/testoptimization.dir/testoptimization.c.s

src/CMakeFiles/testoptimization.dir/testoptimization.c.o.requires:

.PHONY : src/CMakeFiles/testoptimization.dir/testoptimization.c.o.requires

src/CMakeFiles/testoptimization.dir/testoptimization.c.o.provides: src/CMakeFiles/testoptimization.dir/testoptimization.c.o.requires
	$(MAKE) -f src/CMakeFiles/testoptimization.dir/build.make src/CMakeFiles/testoptimization.dir/testoptimization.c.o.provides.build
.PHONY : src/CMakeFiles/testoptimization.dir/testoptimization.c.o.provides

src/CMakeFiles/testoptimization.dir/testoptimization.c.o.provides.build: src/CMakeFiles/testoptimization.dir/testoptimization.c.o


# Object files for target testoptimization
testoptimization_OBJECTS = \
"CMakeFiles/testoptimization.dir/testoptimization.c.o"

# External object files for target testoptimization
testoptimization_EXTERNAL_OBJECTS =

src/testoptimization: src/CMakeFiles/testoptimization.dir/testoptimization.c.o
src/testoptimization: src/CMakeFiles/testoptimization.dir/build.make
src/testoptimization: src/libscientific.dylib
src/testoptimization: src/CMakeFiles/testoptimization.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/marco/projects/libscientific-0.7.4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable testoptimization"
	cd /Users/marco/projects/libscientific-0.7.4/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testoptimization.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/testoptimization.dir/build: src/testoptimization

.PHONY : src/CMakeFiles/testoptimization.dir/build

src/CMakeFiles/testoptimization.dir/requires: src/CMakeFiles/testoptimization.dir/testoptimization.c.o.requires

.PHONY : src/CMakeFiles/testoptimization.dir/requires

src/CMakeFiles/testoptimization.dir/clean:
	cd /Users/marco/projects/libscientific-0.7.4/build/src && $(CMAKE_COMMAND) -P CMakeFiles/testoptimization.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/testoptimization.dir/clean

src/CMakeFiles/testoptimization.dir/depend:
	cd /Users/marco/projects/libscientific-0.7.4/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/marco/projects/libscientific-0.7.4 /Users/marco/projects/libscientific-0.7.4/src /Users/marco/projects/libscientific-0.7.4/build /Users/marco/projects/libscientific-0.7.4/build/src /Users/marco/projects/libscientific-0.7.4/build/src/CMakeFiles/testoptimization.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/testoptimization.dir/depend

