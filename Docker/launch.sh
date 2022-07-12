#! /bin/bash

# find absolute path of the directory of this script
parent_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd -P)"

#path to the NumCSE repository
repo_root="$parent_dir/../.."

docker build -t nmcse-dev .
docker run -it -v "$repo_root:/numcse" nmcse-dev /bin/bash