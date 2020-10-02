# Numerical Methods for CSE

This repository hosts all the codes used in the lecture notes and assignments.

**Additional links**

- [Course VVZ](http://www.vvz.ethz.ch/Vorlesungsverzeichnis/lerneinheit.view?semkez=2020W&ansicht=KATALOGDATEN&lerneinheitId=140998&lang=de)
- Documentation:
	- [Cppreference](https://en.cppreference.com/w/)
	- [Eigen](http://eigen.tuxfamily.org/dox/)
- Tutorials:
	- [Git tutorial](https://gitlab.math.ethz.ch/tille/gitlab-introduction/blob/master/git/README.md)
	- [Debugger tutorial](https://gitlab.math.ethz.ch/tille/debugging-cpp-code-with-lldb)

## Project Structure

- `LectureCodes` - all codes used in the lecture notes / script sorted by their subject; there are different versions of the same code: C++ / Matlab (for most of the codes) / Python (for most of the codes)
- `Assignments` - all codes used in the current assignments. The path for each problem is – note the CodeExpert subfolders:
    - for templates: `Assignments/Codes/<Chapter>/<ProblemName>/CodeExpert/template`
    - for solutions: `Assignments/Codes/<Chapter>/<ProblemName>/CodeExpert/solution`
- `MathGL` - legacy plotting; documentation and example codes for MathGL
- `CMake` - macros, modules used by CMake
- `CppTutorial` - tutorial for C++
- `Docker` - provides a Dockerfile for setup of a Docker container for the repository
- `FFT` - demos for FFT
- `MatplotlibC++` - contains the header for MatplotlibC++
- `OldExam` - code of old exam problems
- `Testing` - contains necessary files for independent testing
- `Utils` - folder for utility codes

## How to use

### Source Download and Compilation
	
	$ git clone https://gitlab.math.ethz.ch/NumCSE/NumCSE.git
	$ cd NumCSE
	$ mkdir build && cd build
	$ cmake ..
	$ make

All binaries can then be found in the `bin` folder.

In order to compile an individual executable, you can recover the list of all possible targets of `make`.
This can be done by typing:

    $ cmake --build . --target help

Then you can choose a target and run:

    $ make <target>

Targets will be in the format `<code-type>_<chapter>_<problem-name>_<executable-name>`, where:
- `<code-type>` is one of:
  - `assignment_codes`: for codes used in assignments
  - `lecture-codes`: for codes used in the lecture notes
- `<chapter>`: is a short name for the chapter
- `<problem-name>`: is a short name for the problem
- `<executable-name>`: is a short name for the executable

*OR* navigate to the assignment or lecture code folder inside `build`, then run 
`make` to compile all sources that correspond to that specific folder. To execute the binaries go to the specific subfolder in `build/bin` and run `./executable_name`.

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

### Working on homework problems
- All current homework problems are located in `Assignments/Codes/[...]/CodeExpert` subfolders! Please work only on these codes as the other ones are considered to be outdated and are currently not maintained.
- Every homework problem comes with a set of tests that can be found in `tests.cpp`. These tests are independent which means that you can fail a subproblem - let's say a - (or do not work on it at all) and call the function from this subproblem in a following subproblem - let's say b; the tests of a correct implementation of b will not fail because your solution will be linked against the mastersolution.
- To compile normal executables as well as tests, just go to your `build` folder, navigate to the `CodeExpert` subfolder of the specific exercise and then run `make`. Once compiled, you can find the executables in the specific `CodeExpert` subfolder of `build/bin`. To run them, call `./template` for your solution, `./solution` for the mastersolution and `./tests` for the tests.
- _Add-On:_ If you want to test multiple exercise at once or you do not want exact information provided by doctest-test-executables, you can invoke CTest in a subfolder with a makefile generated by CMake by writing `make test` (`make all test` if executables have not been compiled before). This will run all tests in this folder.

### Dependencies / Requirements

Required:
- C++ compiler (C++17 support required), tested only with gcc and clang
- Git
- Cmake
- Eigen (best to install it yourself, but will be installed during cmake-process if not found)
- Python
- Numpy
- Matplotlib

Optional:
- boost
- gmp
- mpfr
- MKL
- FFTW

### Visualization using `matplotlibcpp.h`

MatplotlibC++ provides a C++ interface to Python's matplotlib
plotting library. It is header-only and simple to use.
A documentation is available [here](https://matplotlib-cpp.readthedocs.io/en/latest/)

### Testing using `doctest.h`

Doctest is a fast, header-only unit testing framework used to test the homework codes.

### Known issues

#### `Unable to find the required Boost libraries`
#### `Could NOT find ZLIB` or `Could not find PNG`

Some package may be missing on your machine.

- On a fresh install of *Ubuntu*:

        sudo apt-get install git cmake libpng++-dev freeglut3-dev libboost-all-dev

- On *Mac OS X*:

If you are missing `CMake` or `boost` on Mac OS X, the easiest way to obtain those packages is via [Homebrew](http://brew.sh/).
After Homebrew has been installed, you can install CMake and boost:

    brew install boost
    brew install cmake

If you are missing `zlib` (not tested):

    xcode-select --install

of (if above doesn't work):

    brew install zlib-devel

If `libpng`is missing:

    brew install libpng

#### `fatal error: mgl2/mgl.h: No such file or directory`

MathGL installation missing (possibly because of previous point) or corrupted:
- make sure all dependencies are satisfied (in particular libpng, zlib, opengl), see previous point
- remove the possibly broken MathGL installation:


        cd build
        rm -R mathgl_install

- clone a fresh copy of the repository

### `dyld: Library not loaded: @rpath/libmgl.7.4.0.dylib`

Known issue (Mac OS X):
- run `cmake .` and `make` twice
- if this doesn't fix problem: try with clean `clone`
