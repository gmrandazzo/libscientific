# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.4

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.4.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.4.1/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/marco/projects/libscientific-0.7.4

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/marco/projects/libscientific-0.7.4/build

# Include any dependencies generated for this target.
include src/CMakeFiles/testmetricspace.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/testmetricspace.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/testmetricspace.dir/flags.make

src/CMakeFiles/testmetricspace.dir/testmetrics.c.o: src/CMakeFiles/testmetricspace.dir/flags.make
src/CMakeFiles/testmetricspace.dir/testmetrics.c.o: ../src/testmetrics.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/marco/projects/libscientific-0.7.4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object src/CMakeFiles/testmetricspace.dir/testmetrics.c.o"
	cd /Users/marco/projects/libscientific-0.7.4/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/testmetricspace.dir/testmetrics.c.o   -c /Users/marco/projects/libscientific-0.7.4/src/testmetrics.c

src/CMakeFiles/testmetricspace.dir/testmetrics.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/testmetricspace.dir/testmetrics.c.i"
	cd /Users/marco/projects/libscientific-0.7.4/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/marco/projects/libscientific-0.7.4/src/testmetrics.c > CMakeFiles/testmetricspace.dir/testmetrics.c.i

src/CMakeFiles/testmetricspace.dir/testmetrics.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/testmetricspace.dir/testmetrics.c.s"
	cd /Users/marco/projects/libscientific-0.7.4/build/src && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/marco/projects/libscientific-0.7.4/src/testmetrics.c -o CMakeFiles/testmetricspace.dir/testmetrics.c.s

src/CMakeFiles/testmetricspace.dir/testmetrics.c.o.requires:

.PHONY : src/CMakeFiles/testmetricspace.dir/testmetrics.c.o.requires

src/CMakeFiles/testmetricspace.dir/testmetrics.c.o.provides: src/CMakeFiles/testmetricspace.dir/testmetrics.c.o.requires
	$(MAKE) -f src/CMakeFiles/testmetricspace.dir/build.make src/CMakeFiles/testmetricspace.dir/testmetrics.c.o.provides.build
.PHONY : src/CMakeFiles/testmetricspace.dir/testmetrics.c.o.provides

src/CMakeFiles/testmetricspace.dir/testmetrics.c.o.provides.build: src/CMakeFiles/testmetricspace.dir/testmetrics.c.o


# Object files for target testmetricspace
testmetricspace_OBJECTS = \
"CMakeFiles/testmetricspace.dir/testmetrics.c.o"

# External object files for target testmetricspace
testmetricspace_EXTERNAL_OBJECTS =

src/testmetricspace: src/CMakeFiles/testmetricspace.dir/testmetrics.c.o
src/testmetricspace: src/CMakeFiles/testmetricspace.dir/build.make
src/testmetricspace: src/libscientific.dylib
src/testmetricspace: src/CMakeFiles/testmetricspace.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/marco/projects/libscientific-0.7.4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable testmetricspace"
	cd /Users/marco/projects/libscientific-0.7.4/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testmetricspace.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/testmetricspace.dir/build: src/testmetricspace

.PHONY : src/CMakeFiles/testmetricspace.dir/build

src/CMakeFiles/testmetricspace.dir/requires: src/CMakeFiles/testmetricspace.dir/testmetrics.c.o.requires

.PHONY : src/CMakeFiles/testmetricspace.dir/requires

src/CMakeFiles/testmetricspace.dir/clean:
	cd /Users/marco/projects/libscientific-0.7.4/build/src && $(CMAKE_COMMAND) -P CMakeFiles/testmetricspace.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/testmetricspace.dir/clean

src/CMakeFiles/testmetricspace.dir/depend:
	cd /Users/marco/projects/libscientific-0.7.4/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/marco/projects/libscientific-0.7.4 /Users/marco/projects/libscientific-0.7.4/src /Users/marco/projects/libscientific-0.7.4/build /Users/marco/projects/libscientific-0.7.4/build/src /Users/marco/projects/libscientific-0.7.4/build/src/CMakeFiles/testmetricspace.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/testmetricspace.dir/depend

