#! /bin/bash

cd "$(dirname "$0")"

# uncomment the /bin/bash in the end for debugging the container (this will skip the ACTION)
# increasing MAXFILESIZE was necessary for the test action
docker run -it -e ACTION=test -e MAXFILESIZE=100000 -v "$(dirname "$0")/projectfiles:/var/lib/cxrun/projectfiles" cx-nmcse #/bin/bash
