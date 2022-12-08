#! /bin/bash
# 
# Example Usage: ./export.sh PolynomialInterpolation/EvaluatingDerivatives [<options>]
# 
# This is a docker wrapper around export.js, you only need to specify the assignment and the output will be next to this file.

parent_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd -P)"
src="node:16.15.0-buster" 
assignment_path="Assignments/PolishedCodes/$1"
run_cmd="node export.js /numcse $assignment_path ${@:2}"
wd="/numcse/CodeExpert/tool"

if [ "$numcse_is_docker" = true ]; then

	cd $wd && $run_cmd

else

	docker pull $src && \
		docker image tag $src cx-tool && \
		docker run -it -w $wd -v "$parent_dir/../..:/numcse" cx-tool $run_cmd

fi