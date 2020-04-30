# Setup of CodeExpert-compatible homework projects

This folder contains additional files for setting up a CodeExpert homework project.


## Instructions for setup in online environment
- Start from an empty project.
- Copy `scripts` and `conf.yml` from this directory to the solution tab.
- Copy all files from the `NumCSE/Testing` directory to the solution tab.
- Copy all files from the `CodeExpert/solution` directory that corresponds to the exercise to the solution tab.
- Copy all files from `CodeExpert/solution` of the exercise to the solution tab.
- Go to the template tab, change `conf.yml` by commenting the line `test: echo Testing not supported in solution.` and uncommenting the line `# test: /scripts/test.sh`. Then change the files that are to be written by the student to the ones in `CodeExpert/template` by removing code.
- Add all header files that the students should write again in the template tab and rename them to `solution.hpp` (if there is more than one `solution2.hpp`, etc.). `solution.hpp` should include all the other `solutionx.hpp` files. Set these files to **hidden**! 
- Change the visibility settings of all files. Recommended: `doctest.h`, `copy_and_tweak.py`, `conf.yml` and maybe `tests.cpp` hidden.
