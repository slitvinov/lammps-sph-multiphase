#! /bin/bash

set -e
set -u

# to visualize with paraview run
cd data
PYTHONPATH=/scratch/work/Pizza.py/src/ python ../../scripts/dump2ensight.py dump*.dat