# Problems

## CodeQuiz (9-2)
- some code style changes, `const`
- made test more like a unit test
- removed comments in `myfunction` in the template

## NewtonArctan (9-3)
- minimal changes in the test (test only in [1,2])

## QuadraticConvergence (9-4)
- some small comment changes
- renewed plot
- **TODO** problem with testing (testing script replaces a line which in this case is not desirable), i.e. *TESTING IS TURNED OFF FOR THIS PROBLEM*

## OrdConvIter (9-5)
- error in Latex (9-5.c): intervals in linspace must be **OPEN**
- ported to CodeExpert and matplotlibcpp
- refactored the code a bit
- replace Latex graphic

## RecursionOrder (9-6)
- changed function signature -> better to have a default upper bound of iterations than to hard-code it -> **ADJUST LATEX, IF OK** otherwise change it back
- Latex should give a hint on initial guesses
- fixed error that will give segfault – 10 iterates give different order than in latex -> change (11 iterates gives the written order)
