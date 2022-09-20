#!/bin/bash

g++ main.cpp -g -ldl -I /usr/include/eigen3 -std=c++17 -DNICEBACKTRACE

./a.out