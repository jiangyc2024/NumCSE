#! /bin/bash

container_name="numcse_dev_container"

# check if the container is not already running
if [ ! "$( docker ps | grep $container_name )" ]; then

	# find absolute path of the directory of this script
	parent_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")"; pwd -P)"

	# path to the NumCSE repository
	repo_root="$parent_dir/.."

	# build base image
	docker build -t numcse_dev_base "$parent_dir/base"

	# build 'last-mile' top image locally with altered build context
	cd "$repo_root" && cp "$parent_dir/top/Dockerfile" . && docker build -t numcse_dev_top . && rm Dockerfile

	#launch container and bind mount the numcse repository for developing codes
	docker run --name $container_name -it -v "$repo_root:/numcse" numcse_dev_top
	
else

	#launch a second shell in the existing container
	docker exec -it $container_name /bin/bash
fi
