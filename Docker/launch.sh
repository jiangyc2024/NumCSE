#! /bin/bash

# find absolute path of the directory of this script
parent_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd -P)"

# path to the NumCSE repository
repo_root="$parent_dir/.."

# build base image
docker build -t nmcse-dev-base "$parent_dir/base"

# build 'last-mile' top image locally
cd "$repo_root" && cp "$parent_dir/top/Dockerfile" . && docker build -t nmcse-dev-top . && rm Dockerfile
docker run -it -v "$repo_root:/numcse" nmcse-dev-top