#!/bin/bash

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <path to Data and Results folder>"
  exit 1
fi

docker run -ti --volume $1/Data:/app/Data --volume $1/Results:/app/Results mave-db-model:latest