#!/bin/bash

export PYTHONIOENCODING="UTF-8"

python3 copy_and_tweak.py

echo "Compiling ..."

mkdir -p bin

if g++ tests.cpp -fdiagnostics-color=always -std=c++11 \
  -I/usr/local/include/python3.7m \
  -I/usr/local/lib/python3.7/site-packages/numpy/core/include \
  -I/usr/include/eigen3/ \
  -lpython3.7m \
  -lpthread -lutil -ldl \
  -Xlinker -export-dynamic \
  -o bin/tests.out
then
	echo "Compilation successful"
	echo "Running tests ..."
	timeout 5 bin/tests.out -npf -s
	rm copy.hpp
fi
