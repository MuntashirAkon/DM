# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Volumes/Fallout/Projects/DM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Volumes/Fallout/Projects/DM

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /Volumes/Fallout/Projects/DM/CMakeFiles /Volumes/Fallout/Projects/DM/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /Volumes/Fallout/Projects/DM/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named dm

# Build rule for target.
dm: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 dm
.PHONY : dm

# fast build rule for target.
dm/fast:
	$(MAKE) -f CMakeFiles/dm.dir/build.make CMakeFiles/dm.dir/build
.PHONY : dm/fast

Process.o: Process.cpp.o

.PHONY : Process.o

# target to build an object file
Process.cpp.o:
	$(MAKE) -f CMakeFiles/dm.dir/build.make CMakeFiles/dm.dir/Process.cpp.o
.PHONY : Process.cpp.o

Process.i: Process.cpp.i

.PHONY : Process.i

# target to preprocess a source file
Process.cpp.i:
	$(MAKE) -f CMakeFiles/dm.dir/build.make CMakeFiles/dm.dir/Process.cpp.i
.PHONY : Process.cpp.i

Process.s: Process.cpp.s

.PHONY : Process.s

# target to generate assembly for a file
Process.cpp.s:
	$(MAKE) -f CMakeFiles/dm.dir/build.make CMakeFiles/dm.dir/Process.cpp.s
.PHONY : Process.cpp.s

Set.o: Set.cpp.o

.PHONY : Set.o

# target to build an object file
Set.cpp.o:
	$(MAKE) -f CMakeFiles/dm.dir/build.make CMakeFiles/dm.dir/Set.cpp.o
.PHONY : Set.cpp.o

Set.i: Set.cpp.i

.PHONY : Set.i

# target to preprocess a source file
Set.cpp.i:
	$(MAKE) -f CMakeFiles/dm.dir/build.make CMakeFiles/dm.dir/Set.cpp.i
.PHONY : Set.cpp.i

Set.s: Set.cpp.s

.PHONY : Set.s

# target to generate assembly for a file
Set.cpp.s:
	$(MAKE) -f CMakeFiles/dm.dir/build.make CMakeFiles/dm.dir/Set.cpp.s
.PHONY : Set.cpp.s

main.o: main.cpp.o

.PHONY : main.o

# target to build an object file
main.cpp.o:
	$(MAKE) -f CMakeFiles/dm.dir/build.make CMakeFiles/dm.dir/main.cpp.o
.PHONY : main.cpp.o

main.i: main.cpp.i

.PHONY : main.i

# target to preprocess a source file
main.cpp.i:
	$(MAKE) -f CMakeFiles/dm.dir/build.make CMakeFiles/dm.dir/main.cpp.i
.PHONY : main.cpp.i

main.s: main.cpp.s

.PHONY : main.s

# target to generate assembly for a file
main.cpp.s:
	$(MAKE) -f CMakeFiles/dm.dir/build.make CMakeFiles/dm.dir/main.cpp.s
.PHONY : main.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... dm"
	@echo "... Process.o"
	@echo "... Process.i"
	@echo "... Process.s"
	@echo "... Set.o"
	@echo "... Set.i"
	@echo "... Set.s"
	@echo "... main.o"
	@echo "... main.i"
	@echo "... main.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

