# Problems

## Fruits (2-1)

- ported to CodeExpert

## BlockLSEPiv (2-2)

- ported to CodeExpert
- wrapped the tabulation task (2-2.f) into a function and added plotting using matplotlibcpp
- **TODO**: Generate a new eps file

## BlockLU (2-3)
- ported to CodeExpert

## CholeskyQR (2-4)
- ported to CodeExpert

## CircuitImpedance (2-5)
- has already been checked thoroughly last semester, so no real changes

## EfficientBandMult (2-6)
- ported to CodeExpert
- changed the initialization for Vector c, d in solvelseA to be a loop because of Eigen assertions (type is templated...); this also conforms to the template requirements

## Householder (2-7)
- added division by squaredNorm() as requested by task description
-

## Lyapunov (2-8)
- small cosmetic changes

## PartitionedMatrix (2-9)
- small cosmetic changes

## RankOneInvit (2-10)
- small cosmetic changes

## SolvePermb (2-11)
- ported to CodeExpert

## StructuredLSE (2-12)
- ported to CodeExpert and added matplotlibcpp code
- buildA: set matrix to zero because some entries are not touched
- changed sizes to go from 16 to 512 instead of 16 to 4096 because of runtime issues on CodeExpert

## TripletToCRS (2-13)
- removed default constructor of Triplet
-
