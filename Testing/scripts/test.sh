#!/bin/bash

export PYTHONIOENCODING="UTF-8"

python3 copy_and_tweak.py

echo "Compiling tests ..."

mkdir -p bin

g++ tests.cpp -fdiagnostics-color=always -std=c++11 \
-I/usr/local/include/python3.7m \
-I/usr/local/lib/python3.7/site-packages/numpy/core/include \
-I/usr/include/eigen3/ \
-lpython3.7m \
-lpthread -lutil -ldl \
-Xlinker -export-dynamic \
-o bin/tests.out \
> /dev/null 2>&1 # ignore any messages output by g++

compile_result=$? # get the return code of the g++ command
                  # 0 means good, 1 or 2 means failure.
if [ $compile_result -eq 0 ]
  then
    echo "Compilation successful"
    echo "Running tests ..."
    timeout 5 bin/tests.out -npf -s
    rm copy.hpp
  else
    echo "Compilation failed. Try to debug with Compile and Run before testing."
fi
