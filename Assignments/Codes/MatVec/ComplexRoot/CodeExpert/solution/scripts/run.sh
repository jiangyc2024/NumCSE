#!/bin/bash

export PYTHONIOENCODING="UTF-8"


# compile (call compile script)
bash "${WORKDIR}/scripts/compile.sh"

# run

echo "Running script ..."
bin/a.out #| python3 main.py
#java -ea -cp java-bin Main