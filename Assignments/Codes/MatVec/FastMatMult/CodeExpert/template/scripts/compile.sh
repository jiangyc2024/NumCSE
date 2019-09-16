#!/bin/bash

echo "Compiling ..."
# compile c code
#mkdir -p bin
#g++ test_strassen.cpp -I/usr/include/eigen3/  -fdiagnostics-color=always --pedantic -Wextra -Wall -std=c++14 \
#    -o bin/test_strassen.out 

#g++ timing.cpp -I/usr/include/eigen3/  -fdiagnostics-color=always --pedantic -Wextra -Wall -std=c++14 \
#    -o bin/timing.out 

mkdir -p bin

find . -iname '*.cpp' | sort | xargs \
  g++ -fdiagnostics-color=always -std=c++11 \
  -I/usr/local/include/python3.7m \
  -I/usr/local/lib/python3.7/site-packages/numpy/core/include \
  -I/usr/include/eigen3/ \
  -lpython3.7m \
  -lpthread -lutil -ldl \
  -Xlinker -export-dynamic \
  -o bin/a.out
  

echo "Compilation successful"