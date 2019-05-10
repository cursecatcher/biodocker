#!/bin/bash

REGISTRY=`cat REGISTRY`
if [ "AA${REGISTRY}" != "AA" ]; then
	REGISTRY="${REGISTRY}/"
fi

IMAGENAME=`cat IMAGENAME`

# Dump actual commit
git log  -n 1 > BUILD_COMMIT

echo "Going to build docker image for ${REGISTRY}/${IMAGENAME}"

docker build -t ${REGISTRY}${IMAGENAME} .
