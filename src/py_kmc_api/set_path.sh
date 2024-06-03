#!/bin/bash

#run this way:
# source set_path.sh
SCRIPTPATH="$( cd "$(dirname "$BASH_SOURCE")" ; pwd -P )"
export PYTHONPATH="${PYTHONPATH}:${SCRIPTPATH}/../bin/"
