#! /bin/bash -x
#PBS -q s
#PBS -l nodes=1:ppn=1
#PBS -N dyn_vi_rhow_solver
#PBS -l walltime=00:10:00
#PBS -o OUT.log
#PBS -e ERR.log

export OMP_NUM_THREADS=1

cd $PBS_O_WORKDIR

HMDIR=`pwd`/../../../..

ln -svf ${HMDIR}/bin/dyn_vi_rhow_solver.exe .
ln -svf ${HMDIR}/kernels/NICAM/dyn_vi_rhow_solver/data/vgrid40_600m_24km.dat .
ln -svf ${HMDIR}/kernels/NICAM/dyn_vi_rhow_solver/data/snapshot.dyn_vi_rhow_solver.pe000000 .

./dyn_vi_rhow_solver.exe
