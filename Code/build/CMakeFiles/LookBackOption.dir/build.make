# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

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
CMAKE_SOURCE_DIR = /mnt/c/Users/mouad/Documents/EK/PAGES/Cpp/LookBackOption

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/mouad/Documents/EK/PAGES/Cpp/LookBackOption/build

# Include any dependencies generated for this target.
include CMakeFiles/LookBackOption.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/LookBackOption.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/LookBackOption.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/LookBackOption.dir/flags.make

CMakeFiles/LookBackOption.dir/main.cpp.o: CMakeFiles/LookBackOption.dir/flags.make
CMakeFiles/LookBackOption.dir/main.cpp.o: ../main.cpp
CMakeFiles/LookBackOption.dir/main.cpp.o: CMakeFiles/LookBackOption.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/mouad/Documents/EK/PAGES/Cpp/LookBackOption/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/LookBackOption.dir/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/LookBackOption.dir/main.cpp.o -MF CMakeFiles/LookBackOption.dir/main.cpp.o.d -o CMakeFiles/LookBackOption.dir/main.cpp.o -c /mnt/c/Users/mouad/Documents/EK/PAGES/Cpp/LookBackOption/main.cpp

CMakeFiles/LookBackOption.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LookBackOption.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/mouad/Documents/EK/PAGES/Cpp/LookBackOption/main.cpp > CMakeFiles/LookBackOption.dir/main.cpp.i

CMakeFiles/LookBackOption.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LookBackOption.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/mouad/Documents/EK/PAGES/Cpp/LookBackOption/main.cpp -o CMakeFiles/LookBackOption.dir/main.cpp.s

# Object files for target LookBackOption
LookBackOption_OBJECTS = \
"CMakeFiles/LookBackOption.dir/main.cpp.o"

# External object files for target LookBackOption
LookBackOption_EXTERNAL_OBJECTS =

LookBackOption: CMakeFiles/LookBackOption.dir/main.cpp.o
LookBackOption: CMakeFiles/LookBackOption.dir/build.make
LookBackOption: /usr/lib/x86_64-linux-gnu/libpython3.8.so
LookBackOption: CMakeFiles/LookBackOption.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/mouad/Documents/EK/PAGES/Cpp/LookBackOption/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable LookBackOption"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/LookBackOption.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/LookBackOption.dir/build: LookBackOption
.PHONY : CMakeFiles/LookBackOption.dir/build

CMakeFiles/LookBackOption.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/LookBackOption.dir/cmake_clean.cmake
.PHONY : CMakeFiles/LookBackOption.dir/clean

CMakeFiles/LookBackOption.dir/depend:
	cd /mnt/c/Users/mouad/Documents/EK/PAGES/Cpp/LookBackOption/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/mouad/Documents/EK/PAGES/Cpp/LookBackOption /mnt/c/Users/mouad/Documents/EK/PAGES/Cpp/LookBackOption /mnt/c/Users/mouad/Documents/EK/PAGES/Cpp/LookBackOption/build /mnt/c/Users/mouad/Documents/EK/PAGES/Cpp/LookBackOption/build /mnt/c/Users/mouad/Documents/EK/PAGES/Cpp/LookBackOption/build/CMakeFiles/LookBackOption.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/LookBackOption.dir/depend

