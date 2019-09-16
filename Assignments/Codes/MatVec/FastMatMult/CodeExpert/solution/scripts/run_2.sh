#!/bin/bash


# I would like to delete this file, but I cannot.
export PYTHONIOENCODING="UTF-8"

# compile (call compile script)
bash "${WORKDIR}/scripts/compile.sh"

echo "Running code ..."
echo "                     "
echo "                     "
# run
echo "Task 1-4.b) 'Testing the correctness of the Strassen Alorithm' "
echo "                       "
bin/test_strassen.out

echo "                       "
echo "Task 1-4.c) 'Comparing the timings of A*B and the Strassen Algorithm' "
echo "                       "
bin/timing.out