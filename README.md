# Numerical Methods for CSE

This repository will host all the codes used in the lecture notes, and assignments.

**Additional links**

- [Course webpage](https://www.sam.math.ethz.ch/~grsam/HS16/NumCSE/)
- [Course VVZ](http://www.vvz.ethz.ch/Vorlesungsverzeichnis/lerneinheitPre.do?lerneinheitId=109126&semkez=2016W&lang=de)
- [Git tutorial](https://gitlab.math.ethz.ch/tille/gitlab-introduction/blob/master/git/README.md)
- [Debugger tutorial](https://gitlab.math.ethz.ch/tille/debugging-cpp-code-with-lldb)

## Project Structure

- `LectureCodes` - all codes used in the lecture notes / script sorted by their subject
- `MathGL` - documentation and example codes for MathGL
- `CMake` - macros, modules used by CMake
- `Assigmnents` - all codes used in the assignment. The path for each problem is:
    - for templates: `Assigmnents/Codes/<Chapter>/<ProblemName>/templates_nolabels`
    - for solutions: `Assigmnents/Codes/<Chapter>/<ProblemName>/solutions_nolabels`
    - each one of these folder has an independent `CMake` file. Either within the cloned repository
      or using the `Download zip` button, you shoud be able to compile and execute the 
      codes for the problem using:

```
$ cmake .
$ make
$ ./executable_name
```

## How to use

In the [LectureCodes](LectureCodes/) you can find a folder for each code from the lecture notes script.

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

TIP: Using:

    $ make -j<number_ob_processes>
    
may (or may not) speed up the compilation time.

__Alternative download__ [zip](https://gitlab.math.ethz.ch/NumCSE/NumCSE/repository/archive.zip?ref=master)

### Third party libraries

Dependencies / Requirements

- C++ compiler (C++11 support required)
- git
- cmake
- mathgl
- eigen
- [gmp]
- [mpfr]
