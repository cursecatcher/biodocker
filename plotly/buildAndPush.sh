#!/bin/bash

docker build -t "heatmap" .
docker tag heatmap:latest cursecatcher/heatmap
docker push cursecatcher/heatmap
