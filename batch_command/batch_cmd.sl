#!/bin/sh

#SBATCH -C haswell
#SBATCH -q genepool_shared
#SBATCH -N 1
#SBATCH -t 02:00:00

ARGS=$(sed -n ${SLURM_ARRAY_TASK_ID}p batch_data.txt)
./batch_wrapper.sh ${ARGS}
