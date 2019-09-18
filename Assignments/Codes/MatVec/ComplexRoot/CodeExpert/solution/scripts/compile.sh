#!/bin/bash

echo "Compiling ..."
# compile c code
mkdir -p bin
find . -iname '*.cpp' | sort | xargs g++ -fdiagnostics-color=always --pedantic -Wextra -Wall -std=c++14 -o bin/a.out || exit 1

# compile java code
# javac -encoding UTF8 -cp "${WORKDIR}" -d java-bin Main.java || exit 1
echo "Compilation successful"