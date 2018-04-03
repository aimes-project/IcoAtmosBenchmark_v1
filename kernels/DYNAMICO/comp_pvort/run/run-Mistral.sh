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

ln -sf ${HMDIR}/bin/comp_pvort.exe .
ln -sf ${HMDIR}/kernels/DYNAMICO/comp_pvort/data/snapshot.comp_pvort.pe000000

./comp_pvort.exe
