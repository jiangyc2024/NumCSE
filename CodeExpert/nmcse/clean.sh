#! /bin/bash

# remove all unused containers and images

docker container prune -f
docker image prune -f