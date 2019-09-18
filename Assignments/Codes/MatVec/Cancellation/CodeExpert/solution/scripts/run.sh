#!/bin/bash

export PYTHONIOENCODING="UTF-8"

# compile (call compile script)
bash "${WORKDIR}/scripts/compile.sh"

# run
bin/a.out 