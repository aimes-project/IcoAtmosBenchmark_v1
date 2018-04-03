#! /bin/bash
#
# for DKRZ Mistral
#
#SBATCH --partition=compute
#SBATCH --account=ku0598
#SBATCH --job-name=IABknl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
set -e

export OMP_NUM_THREADS=1

HMDIR=`pwd`/../../../..

ln -sf ${HMDIR}/bin/comp_caldyn_vert.exe .
ln -sf ${HMDIR}/kernels/DYNAMICO/comp_caldyn_vert/data/snapshot.comp_caldyn_vert.pe000000

./comp_caldyn_vert.exe
