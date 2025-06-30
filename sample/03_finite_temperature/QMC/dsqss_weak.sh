#!/bin/sh

#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 32
#SBATCH -t 00:30:00

set -e

source /home/issp/materiapps/oneapi_compiler_classic-2023.0.0--openmpi-4.1.5/dsqss/dsqssvars.sh
module list

sh exec_weak.sh
