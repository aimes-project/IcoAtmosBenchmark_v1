#! /bin/bash -x
#PBS -q s
#PBS -l nodes=1:ppn=1
#PBS -N dyn_horiz_adv_limiter
#PBS -l walltime=00:10:00
#PBS -o OUT.log
#PBS -e ERR.log

export OMP_NUM_THREADS=1

cd $PBS_O_WORKDIR

HMDIR=`pwd`/../../../..

ln -svf ${HMDIR}/bin/dyn_horiz_adv_limiter.exe .
ln -svf ${HMDIR}/kernels/NICAM/dyn_horiz_adv_limiter/data/snapshot.dyn_horiz_adv_limiter.pe000000 .

./dyn_horiz_adv_limiter.exe
