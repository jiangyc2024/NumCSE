# Problems

## Fruits (3-1)

- ported to CodeExpert

## BlockLSEPiv (3-2)

- ported to CodeExpert
- wrapped the tabulation task (3-2.f) into a function and added plotting using matplotlibcpp
- **TODO:** Generate a new eps file and include CPU specifications in LATEX.

## BlockLU (3-3)
- ported to CodeExpert

## CholeskyQR (3-4)
- ported to CodeExpert

## CircuitImpedance (3-5)
- has already been checked thoroughly last semester, so no real changes

## EfficientBandMult (3-6)
- ported to CodeExpert
- changed the initialization for Vector c, d in solvelseA to be a loop because of Eigen assertions (type is templated...); this also conforms to the template requirements
- **TODO:** There are some typos in the Latex

## Householder (3-7)
- added division by squaredNorm() as requested by task description
- 

## Lyapunov (3-8)
- small cosmetic changes

## PartitionedMatrix (3-9)
- small cosmetic changes

## RankOneInvit (3-10)
- small cosmetic changes

## SolvePermb (3-11)
- ported to CodeExpert

## StructuredLSE (3-12)
- ported to CodeExpert and added matplotlibcpp code
- buildA: set matrix to zero because some entries are not touched
- **TODO**: generate eps files
- changed sizes to go from 16 to 512 instead of 16 to 4096 because of runtime issues on CodeExpert

## TripletToCRS (3-13)
- removed default constructor of Triplet
- 