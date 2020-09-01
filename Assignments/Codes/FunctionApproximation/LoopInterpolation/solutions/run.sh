#! /bin/bash

if g++ -std=c++11 -I/usr/include/eigen3 main.cpp 
	then
	echo Compilation successful! Now Running...
	./a.out 
else
	echo Compilation error!
fi 
