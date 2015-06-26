#! /usr/bin/env bash
#
# How to setup the job, starting in the hfsolver root directory:
#
# HFSOLVER_DIR=`pwd`
# mkdir -p /scratch/certik/hfsolver/run1
# cd /scratch/certik/hfsolver/run1
# cp $HFSOLVER_DIR/src/tests/fem/run.sh .
# msub run.sh
#
# modify the run.sh script (nodes, ppn, walltime).
#
#MSUB -l nodes=4:ppn=12     # number of nodes and cores (ppn)
#MSUB -l walltime=0:59:00   # max duration of the job
#MSUB -m be                 # send mail twice: begin and end
#MSUB -j oe                 # combine stdout and stderr
#MSUB -V                    # forward current environment variables to the job
# Pathname for output:
#MSUB -o myjob.out

##### These are shell commands

# Exit on error
set -e
# Echo each command
set -x

date
pwd

set +x
module list
module purge
module load gcc/4.8.2 openmpi
module list
set -x

# Number of MPI processes to run.
NUMPROC=${SLURM_NPROCS}
echo "Running $NUMPROC processes."

HF_DIR=/users/certik/repos/hfsolver/src/tests/fem

mpiexec -n $NUMPROC $HF_DIR/test_free_energy_bddc 4 8 8 9 4 4 3

date
echo 'Done'
