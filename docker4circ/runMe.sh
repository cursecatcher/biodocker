#!/bin/bash

REGISTRY=`cat REGISTRY`
if [ "AA${REGISTRY}" != "AA" ]; then
        REGISTRY="${REGISTRY}/"
fi

IMAGENAME=`cat IMAGENAME`

echo "Going to run docker image ${REGISTRY}${IMAGENAME}"
docker run -it ${REGISTRY}${IMAGENAME}
