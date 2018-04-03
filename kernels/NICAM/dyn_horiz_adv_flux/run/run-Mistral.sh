#! /bin/bash -x
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

ln -svf ${HMDIR}/bin/dyn_horiz_adv_flux.exe .
ln -svf ${HMDIR}/kernels/NICAM/dyn_horiz_adv_flux/data/snapshot.dyn_horiz_adv_flux.pe000000 .

./dyn_horiz_adv_flux.exe
