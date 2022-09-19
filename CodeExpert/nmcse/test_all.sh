#! /bin/bash

for assignment in $(cat assignments.txt); do echo $assignment; ./full_run.sh $assignment "test"; done