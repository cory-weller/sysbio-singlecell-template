#!/usr/bin/env bash

if [ ! "$CONDA_SRC" = "" ]; then
    source ${CONDA_SRC}
    conda activate ${CONDA_ENV}
fi

# arguments = modules to load
for _module in $@; do
    module load $_module || { echo "ERROR: module $_module not found.";  exit 1; }
done


python3 -u ${SCRIPT_PATH}