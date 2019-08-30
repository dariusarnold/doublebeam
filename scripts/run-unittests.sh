#!/usr/bin/env bash

# Run all unittests in the tests folder of this repo.

# exit immediately if a command fails
trap 'exit' ERR
# fix for conda activate not working in bash script 
# see https://github.com/conda/conda/issues/7980
eval "$(conda shell.bash hook)"
conda activate py37
#echo `"${0%/*}/.."`
cd "${0%/*}/.."
python -m unittest discover -s tests
