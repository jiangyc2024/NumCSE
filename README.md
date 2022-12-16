# Numerical Methods for CSE

This repository hosts all the codes used in the lecture notes and assignments.

**Additional links**

- [Course VVZ](https://www.vorlesungen.ethz.ch/Vorlesungsverzeichnis/sucheLehrangebot.view?lang=de&search=on&semkez=2022W&lerneinheitstitel=Numerical+Methods+for+C*&famname=Hiptmair)
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
- `Docker` - provides files for setup of a Docker container for the repository
- `FFT` - demos for FFT
- `MatplotlibC++` - contains the header for MatplotlibC++
- `OldExam` - code of old exam problems
- `Testing` - contains necessary files for independent testing
- `Utils` - folder for utility codes

## How to use

### Download and Compilation with Docker (recommended)

First, please install [Docker](https://www.docker.com/products/docker-desktop/), if you have not done so already. For a more detailed explanation, see [the docker readme](Docker/README.md). Then run

    $ git clone https://gitlab.math.ethz.ch/NumCSE/NumCSE.git
    $ cd NumCSE
    $ Docker/enter.sh
    root@*:/build# make

### Download and Compilation without Docker

**Caution**: the requirements of this code base are rather complex. Proceeding without Docker may result in a rather large, continued effort to fulfill these requirements on your local machine. Please also note that fixing any issues encountered may be your own responsibility, since others cannot always reproduce them. As a starting point, you can read the [Dockerfile](Docker/base/Dockerfile) and use it as documentation to set up your own machine. Also, see our [known issues and dependencies](Docker/no_docker_notes.md). You may not need all of them for your specific task. After setup, you should be able to run
    
    $ git clone https://gitlab.math.ethz.ch/NumCSE/NumCSE.git
    $ cd NumCSE
    $ mkdir build && cd build
    $ cmake ..
    $ make

### Build process

All binaries can be found in the `bin` folder.

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

__Alternative download__ [zip](https://gitlab.math.ethz.ch/NumCSE/NumCSE/repository/archive.zip?ref=master)

### Working on homework problems
- All current homework problems are located in `Assignments/Codes/[...]/CodeExpert` subfolders! Please work only on these codes as the other ones are considered to be outdated and are currently not maintained.
- Every homework problem comes with a set of tests that can be found in `tests.cpp`. These tests are independent which means that you can fail a subproblem - let's say a - (or do not work on it at all) and call the function from this subproblem in a following subproblem - let's say b; the tests of a correct implementation of b will not fail because your solution will be linked against the mastersolution.
- To compile normal executables as well as tests, just go to your `build` folder, navigate to the `CodeExpert` subfolder of the specific exercise and then run `make`. Once compiled, you can find the executables in the specific `CodeExpert` subfolder of `build/bin`. To run them, call `./template` for your solution, `./solution` for the mastersolution and `./tests` for the tests.
- _Add-On:_ If you want to test multiple exercise at once or you do not want exact information provided by doctest-test-executables, you can invoke CTest in a subfolder with a makefile generated by CMake by writing `make test` (`make all test` if executables have not been compiled before). This will run all tests in this folder.

### Visualization using `matplotlibcpp.h`

MatplotlibC++ provides a C++ interface to Python's matplotlib
plotting library. It is header-only and simple to use.
A documentation is available [here](https://matplotlib-cpp.readthedocs.io/en/latest/)

### Testing using `doctest.h`

Doctest is a fast, header-only unit testing framework used to test the homework codes.