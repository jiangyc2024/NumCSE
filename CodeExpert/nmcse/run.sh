#! /bin/bash

# find absolute path of the directory of this script
parent_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd -P)"

# if no action is specified, use the "run" action
cx_action=${1-"run"}

# make an output directory, for the output files (plots, etc.)
mkdir -p "$parent_dir/projectfiles/cx_out"

# uncomment the /bin/bash in the end for debugging the container (this will skip the ACTION)
# increasing MAXFILESIZE was necessary for the test action, see https://docs.expert.ethz.ch/lecturers/#environment-variables
docker run -e ACTION=$cx_action -e MAXFILESIZE=20000 -e TIMEOUT=1000 -e CPUTIME=1000 -v "$parent_dir/projectfiles:/var/lib/cxrun/projectfiles" cx-nmcse #/bin/bash
