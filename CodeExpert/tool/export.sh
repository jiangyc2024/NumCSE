#! /bin/bash
# This is a docker wrapper around export.js, you only need to specify the assignment and the output will be next to this file.
# Use this if you are a docker environment

parent_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd -P)"
src="node:16.15.0-buster" 
args=${@:1}

docker pull $src && \
	docker image tag $src cx-tool && \
	docker run -it -w /numcse/CodeExpert/tool -v "$parent_dir/../..:/numcse" cx-tool node export.js /numcse $args