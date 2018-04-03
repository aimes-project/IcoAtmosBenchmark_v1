#! /bin/bash -x
#PBS -q s
#PBS -l nodes=1:ppn=2
#PBS -N communication
#PBS -l walltime=00:10:00
#PBS -o OUT.log
#PBS -e ERR.log

export OMP_NUM_THREADS=1

cd $PBS_O_WORKDIR

HMDIR=`pwd`/../../../..

ln -svf ${HMDIR}/bin/communication.exe .
ln -svf ${HMDIR}/kernels/NICAM/communication/data/communication.cnf .

mpirun -np 2 ./communication.exe
