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

	$ git clone https://gitlab.math.ethz.ch/NumCSE/NumCSE.git
	$ cd NumCSE
	$ mkdir build && cd build
	$ cmake ..
	$ make

All binaries can then be found in the `bin` folder.

In order to compile an individual executable, type:

    $ make <code-type>_<chapter>_<problem-name>_<executable-name>
    
Where:
- `<code-type>` is one of:
  - `assignment_codes`: for codes used in assignments
  - `lecture-codes`: for codes used in the lecture notes
- `<chapter>`: is a short name for the chapter
- `<problem-name>`: is a short name for the problem
- `<executable-name>`: is a short name for the executable
    
The corresponding executable will be located in:
- For assignments:


        $ ./bin/Assignments/Codes/<chapter>/<problem-name>/

- For lecture codes (most of them):

        $ ./bin/LectureCodes/<chapter>/<problem-name>/
        
    
    or
    
        $ ./bin/LectureCodes/<chapter>/<sub-chapter>/<problem-name>/

TIP: Using:

    $ make -j<number_ob_processes>
    
may (or may not) speed up the compilation time.

__Alternative download__ [zip](https://gitlab.math.ethz.ch/NumCSE/NumCSE/repository/archive.zip?ref=master)

### Third party libraries

Dependencies / Requirements

Required:
- C++ compiler (C++11 support required), tested only with gcc and clang
- git
- cmake
- mathgl (requires libpng, opengl and zlib)
- eigen
- boost
Optional:
- gmp
- mpfr
- MKL

### F.A.Q.

#### `Unable to find the required Boost libraries`
#### `Could NOT find ZLIB` or `Could not find PNG`

Some package may be missing on your machine.

- On a fresh install of *Ubuntu*:

    sudo apt-get install git cmake libpng++-dev freeglut3-dev libboost-dev-all
    
- On *Mac OS X*:

If you are missing `CMake` or `boost` on Mac OS X, the easiest way to obtain those packages is via [Homebrew](http://brew.sh/).
After Homebrew has been installed, you can install CMake and boost:

    brew install boost
    brew install cmake
    
If you are missing `zlib` (not tested):

    xcode-select --install

of (if above doesn't work):

    brew install zlib-devel

#### `fatal error: mgl2/mgl.h: No such file or directory`

MathGL installation missing (possibly because of previous point) or corrupted:
- make sure all dependencies are satisfied (in particular libpng, zlib, opengl), see previous point
- remove the possibly broken MathGL installation:


        cd build
        rm -R mathgl_install

- clone a fresh copy of the repository
