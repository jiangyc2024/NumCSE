# Numerical Methods for CSE

This repository will host all the codes used in the lecture notes, and assignments.

**Additional links**

- [Debugger tutorial](https://gitlab.math.ethz.ch/tille/debugging-cpp-code-with-lldb)
- [Developer notes](Documentation/developer_notes.md)
- [Project specific cmake documentation](Documentation/cmake.md)

## Project Structure

- `LectureCodes` - all codes used in the lecture notes / script sorted by their subject
- `Documentation` - enhanced documentation for developers
- `MathGL` - documentation and example codes for MathGL
- `CMake` - macros, modules used by CMake

## How to use

In the LectureCodes you can find a folder for each code from the script.

There you can find different versions of the same code:

- C++
- Matlab (most of the codes)
- Python (most of the codes)

### Source Download and Compilation

	$ git clone git@gitlab.math.ethz.ch:NumCSE/NumCSE.git
	$ mkdir build && cd build
	$ cmake ..
	$ make

All binaries can then be found in the `bin` folder.

### Building C++-codes & usage with Eigen and MathGL/Figure:

All C++ codes should use the `add_executable_numcse` macro 
that works similar to cmake built in `add_executable` command, but 
automatically links with MathGL and places the binary in the correct folder.
For further information see [Documentation](Documentation/cmake.md)

**Example CMakeLists.txt file**

```
add_executable_numcse(main main.cpp)
```

### Third party libraries

Dependencies / Requirements

- C++ compiler (C++11 support required)
- git
- cmake
- mathgl
- eigen
- [gmp]
- [mpfr]