#!/bin/bash
eval "$(conda shell.bash hook)"

IPOD_BASE=$1
IDR_BASE=$2

cd ${IPOD_BASE}
conda env create -f conda_environment.yaml
conda activate ipod2

cd ${IDR_BASE}
python setup.py install

conda deactivate

