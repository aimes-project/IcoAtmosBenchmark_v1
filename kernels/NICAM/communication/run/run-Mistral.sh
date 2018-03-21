#! /bin/bash -x
#
# for DKRZ Mistral
#
#SBATCH --partition=compute
#SBATCH --account=ku0598
#SBATCH --job-name=IABknl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=00:10:00
set -e

export OMP_NUM_THREADS=1

HMDIR=`pwd`/../../../..

ln -svf ${HMDIR}/bin/communication.exe .
ln -svf ${HMDIR}/kernels/NICAM/communication/data/communication.cnf .

mpirun -np 2 ./communication.exe
