# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

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
CMAKE_SOURCE_DIR = /home/local1/my-B1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/local1/my-B1/build

# Include any dependencies generated for this target.
include CMakeFiles/exampleB1.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/exampleB1.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/exampleB1.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/exampleB1.dir/flags.make

CMakeFiles/exampleB1.dir/exampleB1.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/exampleB1.cc.o: /home/local1/my-B1/exampleB1.cc
CMakeFiles/exampleB1.dir/exampleB1.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/local1/my-B1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/exampleB1.dir/exampleB1.cc.o"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/exampleB1.cc.o -MF CMakeFiles/exampleB1.dir/exampleB1.cc.o.d -o CMakeFiles/exampleB1.dir/exampleB1.cc.o -c /home/local1/my-B1/exampleB1.cc

CMakeFiles/exampleB1.dir/exampleB1.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/exampleB1.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/local1/my-B1/exampleB1.cc > CMakeFiles/exampleB1.dir/exampleB1.cc.i

CMakeFiles/exampleB1.dir/exampleB1.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/exampleB1.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/local1/my-B1/exampleB1.cc -o CMakeFiles/exampleB1.dir/exampleB1.cc.s

CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o: /home/local1/my-B1/src/ActionInitialization.cc
CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/local1/my-B1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o -MF CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o.d -o CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o -c /home/local1/my-B1/src/ActionInitialization.cc

CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/local1/my-B1/src/ActionInitialization.cc > CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.i

CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/local1/my-B1/src/ActionInitialization.cc -o CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.s

CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o: /home/local1/my-B1/src/DetectorConstruction.cc
CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/local1/my-B1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o -MF CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o.d -o CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o -c /home/local1/my-B1/src/DetectorConstruction.cc

CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/local1/my-B1/src/DetectorConstruction.cc > CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.i

CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/local1/my-B1/src/DetectorConstruction.cc -o CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.s

CMakeFiles/exampleB1.dir/src/EventAction.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/src/EventAction.cc.o: /home/local1/my-B1/src/EventAction.cc
CMakeFiles/exampleB1.dir/src/EventAction.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/local1/my-B1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/exampleB1.dir/src/EventAction.cc.o"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/src/EventAction.cc.o -MF CMakeFiles/exampleB1.dir/src/EventAction.cc.o.d -o CMakeFiles/exampleB1.dir/src/EventAction.cc.o -c /home/local1/my-B1/src/EventAction.cc

CMakeFiles/exampleB1.dir/src/EventAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/src/EventAction.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/local1/my-B1/src/EventAction.cc > CMakeFiles/exampleB1.dir/src/EventAction.cc.i

CMakeFiles/exampleB1.dir/src/EventAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/src/EventAction.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/local1/my-B1/src/EventAction.cc -o CMakeFiles/exampleB1.dir/src/EventAction.cc.s

CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o: /home/local1/my-B1/src/PrimaryGeneratorAction.cc
CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/local1/my-B1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o -MF CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o.d -o CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o -c /home/local1/my-B1/src/PrimaryGeneratorAction.cc

CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/local1/my-B1/src/PrimaryGeneratorAction.cc > CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.i

CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/local1/my-B1/src/PrimaryGeneratorAction.cc -o CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.s

CMakeFiles/exampleB1.dir/src/RunAction.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/src/RunAction.cc.o: /home/local1/my-B1/src/RunAction.cc
CMakeFiles/exampleB1.dir/src/RunAction.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/local1/my-B1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/exampleB1.dir/src/RunAction.cc.o"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/src/RunAction.cc.o -MF CMakeFiles/exampleB1.dir/src/RunAction.cc.o.d -o CMakeFiles/exampleB1.dir/src/RunAction.cc.o -c /home/local1/my-B1/src/RunAction.cc

CMakeFiles/exampleB1.dir/src/RunAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/src/RunAction.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/local1/my-B1/src/RunAction.cc > CMakeFiles/exampleB1.dir/src/RunAction.cc.i

CMakeFiles/exampleB1.dir/src/RunAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/src/RunAction.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/local1/my-B1/src/RunAction.cc -o CMakeFiles/exampleB1.dir/src/RunAction.cc.s

CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o: CMakeFiles/exampleB1.dir/flags.make
CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o: /home/local1/my-B1/src/SteppingAction.cc
CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o: CMakeFiles/exampleB1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/local1/my-B1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o -MF CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o.d -o CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o -c /home/local1/my-B1/src/SteppingAction.cc

CMakeFiles/exampleB1.dir/src/SteppingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB1.dir/src/SteppingAction.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/local1/my-B1/src/SteppingAction.cc > CMakeFiles/exampleB1.dir/src/SteppingAction.cc.i

CMakeFiles/exampleB1.dir/src/SteppingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB1.dir/src/SteppingAction.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/local1/my-B1/src/SteppingAction.cc -o CMakeFiles/exampleB1.dir/src/SteppingAction.cc.s

# Object files for target exampleB1
exampleB1_OBJECTS = \
"CMakeFiles/exampleB1.dir/exampleB1.cc.o" \
"CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o" \
"CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o" \
"CMakeFiles/exampleB1.dir/src/EventAction.cc.o" \
"CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o" \
"CMakeFiles/exampleB1.dir/src/RunAction.cc.o" \
"CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o"

# External object files for target exampleB1
exampleB1_EXTERNAL_OBJECTS =

exampleB1: CMakeFiles/exampleB1.dir/exampleB1.cc.o
exampleB1: CMakeFiles/exampleB1.dir/src/ActionInitialization.cc.o
exampleB1: CMakeFiles/exampleB1.dir/src/DetectorConstruction.cc.o
exampleB1: CMakeFiles/exampleB1.dir/src/EventAction.cc.o
exampleB1: CMakeFiles/exampleB1.dir/src/PrimaryGeneratorAction.cc.o
exampleB1: CMakeFiles/exampleB1.dir/src/RunAction.cc.o
exampleB1: CMakeFiles/exampleB1.dir/src/SteppingAction.cc.o
exampleB1: CMakeFiles/exampleB1.dir/build.make
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4Tree.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4FR.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4GMocren.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4visHepRep.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4RayTracer.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4VRML.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4ToolsSG.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4OpenGL.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4visQt3D.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4vis_management.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4modeling.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4interfaces.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4mctruth.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4geomtext.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4gdml.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4error_propagation.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4readout.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4physicslists.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4run.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4event.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4tracking.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4parmodels.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4processes.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4digits_hits.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4track.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4particles.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4geometry.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4materials.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4graphics_reps.so
exampleB1: /usr/lib64/libXm.so
exampleB1: /usr/lib64/libXmu.so
exampleB1: /usr/lib64/libXext.so
exampleB1: /usr/lib64/libXt.so
exampleB1: /usr/lib64/libICE.so
exampleB1: /usr/lib64/libSM.so
exampleB1: /usr/lib64/libX11.so
exampleB1: /usr/lib64/libQt6Widgets.so.6.6.2
exampleB1: /usr/lib64/libQt6OpenGL.so.6.6.2
exampleB1: /usr/lib64/libQt6Gui.so.6.6.2
exampleB1: /usr/lib64/libGL.so
exampleB1: /usr/lib64/libQt6Core.so.6.6.2
exampleB1: /usr/local/XercesC/3.2.5/lib/libxerces-c.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4analysis.so
exampleB1: /usr/lib64/libexpat.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4zlib.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4intercoms.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4global.so
exampleB1: /usr/local/geant4.11.2.1/lib64/libG4ptl.so.2.3.3
exampleB1: /usr/local/clhep/2.4.7.1/lib/libCLHEP-2.4.7.1.so
exampleB1: CMakeFiles/exampleB1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/local1/my-B1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX executable exampleB1"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/exampleB1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/exampleB1.dir/build: exampleB1
.PHONY : CMakeFiles/exampleB1.dir/build

CMakeFiles/exampleB1.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/exampleB1.dir/cmake_clean.cmake
.PHONY : CMakeFiles/exampleB1.dir/clean

CMakeFiles/exampleB1.dir/depend:
	cd /home/local1/my-B1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/local1/my-B1 /home/local1/my-B1 /home/local1/my-B1/build /home/local1/my-B1/build /home/local1/my-B1/build/CMakeFiles/exampleB1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/exampleB1.dir/depend

