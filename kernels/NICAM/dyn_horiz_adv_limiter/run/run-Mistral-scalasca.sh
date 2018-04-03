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

export OMP_NUM_THREADS=4
#export SCOREP_METRIC_PAPI=PAPI_FP_OPS

HMDIR=`pwd`/../../../..

ln -svf ${HMDIR}/bin/dyn_horiz_adv_limiter.exe .
ln -svf ${HMDIR}/kernels/NICAM/dyn_horiz_adv_limiter/data/snapshot.dyn_horiz_adv_limiter.pe000000 .

rm -rf ./epik_trace

scan -t -e epik_trace ./dyn_horiz_adv_limiter.exe
