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

export OMP_NUM_THREADS=4
#export SCOREP_METRIC_PAPI=PAPI_FP_OPS

HMDIR=`pwd`/../../../..

ln -sf ${HMDIR}/bin/dyn_vi_rhow_solver.exe .
ln -sf ${HMDIR}/kernels/NICAM/dyn_vi_rhow_solver/data/vgrid40_600m_24km.dat .
ln -sf ${HMDIR}/kernels/NICAM/dyn_vi_rhow_solver/data/snapshot.dyn_vi_rhow_solver.pe000000 .

rm -rf ./epik_trace

scan -t -e epik_trace ./dyn_vi_rhow_solver.exe
