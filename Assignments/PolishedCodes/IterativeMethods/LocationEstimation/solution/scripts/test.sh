#!/bin/bash

echo "Compiling tests..."

mkdir -p bin

g++ tests.cpp -fdiagnostics-color=always -std=c++17 \
-I/usr/include/eigen3/ \
-o bin/tests.out 

compile_result=$? # get the return code of the g++ command
                  # 0 means good, 1 or 2 means failure.
if [ $compile_result -eq 0 ]
  then
    echo "Compilation successful"
    echo "Running tests ..."
    bin/tests.out -npf -s
  else
    echo "Compilation failed. Try to debug with Compile and Run before testing."
fi
