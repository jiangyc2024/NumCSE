#! /bin/bash

parent_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd -P)"

# necessary for the output files
mkdir -p "$parent_dir/projectfiles/cx_out"

# uncomment the /bin/bash in the end for debugging the container (this will skip the ACTION)
# increasing MAXFILESIZE was necessary for the test action
docker run -it -e ACTION=run -e MAXFILESIZE=100000 -v "$parent_dir/projectfiles:/var/lib/cxrun/projectfiles" cx-nmcse #/bin/bash
