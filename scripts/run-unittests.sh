#!/usr/bin/env bash

# exit immediately if a command fails
trap 'exit' ERR

source activate py36
#echo `"${0%/*}/.."`
cd "${0%/*}/.."
python -m unittest discover -s tests
