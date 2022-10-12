#!/bin/bash

echo "Compiling ..."

mkdir -p bin

# Returns true if the compilation passes (even with warnings).
if g++ main.cpp -fdiagnostics-color=always -std=c++17 \
  -I/usr/local/include/python3.7m \
  -I/usr/local/lib/python3.7/site-packages/numpy/core/include \
  -I/usr/include/eigen3/ \
  -lpython3.7m \
  -lpthread -lutil -ldl \
  -Xlinker -export-dynamic \
  -o bin/a.out
then
  echo "Compilation successful"
else
  echo "Compilation failed"
fi
