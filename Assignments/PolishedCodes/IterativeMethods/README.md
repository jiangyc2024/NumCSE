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

## NonLinElectr (9-7)
- ported to matplotlibcpp and CE
- **TODO**: create Udiode.eps plot for Latex (my eps creation did not work) and put in project root folder
- added Latex markers for code

## JuliaSet (9-8)
- removed `using namespace Eigen;`

## ModifiedNewton (9-9)
- no changes

## QuasiLinear (9-10)
- **typo in LATEX**: 9.10.3 iterate x^{(k+11)} is used; should be x^{(k+1)}
- remove TODO in subtask f in latex (already implemented)
- small cosmetic changes

## Radioactive (9-11) see LeastSquares problem folder
- completely reworked the code
- still have to compute solution.eps for Latex
- Latex seems to be outdated

## CircleAppr (9-12) see LeastSquares problem folder
- wrote this task from ground up (nearly), changed some interfaces, but Latex needs a lot of work
- comparison plot is quite heavy – maybe zoom in to see difference between circles
- TODO: generate eps plots

## LevelSetPlot (9-13)
- **TODO**

## SymRankOneApprox (9-14)
- **TODO**

