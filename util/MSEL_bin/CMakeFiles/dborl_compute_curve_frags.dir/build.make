# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/guoy/workspace/MSEL/MSEL_src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/guoy/workspace/MSEL/MSEL_bin

# Include any dependencies generated for this target.
include CMakeFiles/dborl_compute_curve_frags.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/dborl_compute_curve_frags.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/dborl_compute_curve_frags.dir/flags.make

CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.o: CMakeFiles/dborl_compute_curve_frags.dir/flags.make
CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.o: /home/guoy/workspace/MSEL/MSEL_src/dborl_edge_det_link_main.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /home/guoy/workspace/MSEL/MSEL_bin/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.o -c /home/guoy/workspace/MSEL/MSEL_src/dborl_edge_det_link_main.cxx

CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/guoy/workspace/MSEL/MSEL_src/dborl_edge_det_link_main.cxx > CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.i

CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/guoy/workspace/MSEL/MSEL_src/dborl_edge_det_link_main.cxx -o CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.s

CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.o.requires:
.PHONY : CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.o.requires

CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.o.provides: CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.o.requires
	$(MAKE) -f CMakeFiles/dborl_compute_curve_frags.dir/build.make CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.o.provides.build
.PHONY : CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.o.provides

CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.o.provides.build: CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.o

# Object files for target dborl_compute_curve_frags
dborl_compute_curve_frags_OBJECTS = \
"CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.o"

# External object files for target dborl_compute_curve_frags
dborl_compute_curve_frags_EXTERNAL_OBJECTS =

dborl_compute_curve_frags: CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.o
dborl_compute_curve_frags: CMakeFiles/dborl_compute_curve_frags.dir/build.make
dborl_compute_curve_frags: MSEL_core/libMSEL_core.a
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libpdf1d.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libmbl.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvil_io.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libbrip.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libbil_algo.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvil_algo.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libgevd.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvil.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libgeotiff.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libopenjpeg2.so.2.0.0
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvtol.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvdgl.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libbsta.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvnl_io.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libbsol.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvsol.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvbl_io.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvgl_io.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvil1.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libtiff.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libjpeg.so
dborl_compute_curve_frags: /usr/lib/x86_64-linux-gnu/libpng.so
dborl_compute_curve_frags: /usr/lib/x86_64-linux-gnu/libz.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvpgl.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvbl.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvsl.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvgl_algo.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvgl.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvnl_algo.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvnl.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libnetlib.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libv3p_netlib.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvul.so
dborl_compute_curve_frags: /home/guoy/lemsvpe/vxl-bin-rel/lib/libvcl.so
dborl_compute_curve_frags: CMakeFiles/dborl_compute_curve_frags.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable dborl_compute_curve_frags"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dborl_compute_curve_frags.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/dborl_compute_curve_frags.dir/build: dborl_compute_curve_frags
.PHONY : CMakeFiles/dborl_compute_curve_frags.dir/build

CMakeFiles/dborl_compute_curve_frags.dir/requires: CMakeFiles/dborl_compute_curve_frags.dir/dborl_edge_det_link_main.cxx.o.requires
.PHONY : CMakeFiles/dborl_compute_curve_frags.dir/requires

CMakeFiles/dborl_compute_curve_frags.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/dborl_compute_curve_frags.dir/cmake_clean.cmake
.PHONY : CMakeFiles/dborl_compute_curve_frags.dir/clean

CMakeFiles/dborl_compute_curve_frags.dir/depend:
	cd /home/guoy/workspace/MSEL/MSEL_bin && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/guoy/workspace/MSEL/MSEL_src /home/guoy/workspace/MSEL/MSEL_src /home/guoy/workspace/MSEL/MSEL_bin /home/guoy/workspace/MSEL/MSEL_bin /home/guoy/workspace/MSEL/MSEL_bin/CMakeFiles/dborl_compute_curve_frags.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/dborl_compute_curve_frags.dir/depend

