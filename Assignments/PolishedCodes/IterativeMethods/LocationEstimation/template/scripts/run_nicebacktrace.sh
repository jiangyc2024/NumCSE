#!/bin/bash

echo "Compiling with debugging options enabled -g..."

mkdir -p bin

# Returns true if the compilation passes (even with warnings).
if g++ main.cpp -fdiagnostics-color=always -std=c++17 -g -DNICEBACKTRACE\
  -I/usr/include/eigen3/ \
  -ldl \
  -o bin/a.out
then
  echo "Compilation successful. Running the program."
  ./bin/a.out
else
  echo "Compilation failed"
fi
