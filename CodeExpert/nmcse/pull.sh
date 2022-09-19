#! /bin/bash

src="cxhub.ethz.ch/cx/cxenv/nmcse:latest"
docker pull $src
docker image tag $src cx-nmcse