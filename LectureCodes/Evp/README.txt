Translate MATLAB lecture codes into C++. 

Collect all of them in this folder, then we will eventually move them somewhere that makes more sense.

Follow the general structure given by the previous contributor, with some changes:

$matlab_title$
 | Eigen
   | CMakeLists.txt
   | main.cpp
   | $matlab_title$.hpp
   | $matlab_title$.cpp
 | MATLAB
   | $matlab_file$ // just copied from the given files


Use += type operators where possible.


TODO: 
gradit (in progress, anwillia)
algcg (in progress, anwillia)
plottheorate
pcgbase
trusssim
trussimdriver
sspowitrp2d
sspowitrp


Note: 
All files should have the following header

///////////////////////////////////////////////////////////////////////////
/// Demonstration code for lecture "Numerical Methods for CSE" @ ETH Zurich
/// Porting <MATLAB title>.m to C++/eIGEN.
/// (C) 2020 SAM, D-MATH
/// Author(s): <name>
/// Repository: https://gitlab.math.ethz.ch/NumCSE/NumCSE/
/// Do not remove this header.
//////////////////////////////////////////////////////////////////////////

