#! /bin/bash
#
# for DKRZ Mistral
#
#SBATCH --partition=compute
#SBATCH --job-name=IABknl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:10:00
set -e

export OMP_NUM_THREADS=1

HMDIR=`pwd`/../../../..

ln -sf ${HMDIR}/bin/comp_geopot.exe .
ln -sf ${HMDIR}/kernels/DYNAMICO/comp_geopot/data/snapshot.comp_geopot.pe000000

./comp_geopot.exe
