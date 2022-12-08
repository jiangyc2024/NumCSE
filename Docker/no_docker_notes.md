Note: this document is only relevant if you do not use the docker setup.

# Dependencies / Requirements

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

# Known issues

## `Unable to find the required Boost libraries`
## `Could NOT find ZLIB` or `Could not find PNG`

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

## `fatal error: mgl2/mgl.h: No such file or directory`

MathGL installation missing (possibly because of previous point) or corrupted:
- make sure all dependencies are satisfied (in particular libpng, zlib, opengl), see previous point
- remove the possibly broken MathGL installation:


        cd build
        rm -R mathgl_install

- clone a fresh copy of the repository

## `dyld: Library not loaded: @rpath/libmgl.7.4.0.dylib`

Known issue (Mac OS X):
- run `cmake .` and `make` twice
- if this doesn't fix problem: try with clean `clone`

## `'Python.h' file not found` or `Python3 not found.  Please install it to use any function of matplotlibcpp.h.`

Cmake failed to identify a proper Python3 distribution with matplotlib and numpy installed.
The easiest way is to install the anaconda distribution, activate conda envrionment using `conda activate base` or
`conda activate self_defined_env`. And then rerun the cmake configuration `cmake ..` in build folder.

## TIP
Using:

    $ make -j<number_of_processes>

may (or may not) speed up the compilation time 