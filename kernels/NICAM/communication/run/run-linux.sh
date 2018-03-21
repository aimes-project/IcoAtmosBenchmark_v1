#! /bin/bash

export OMP_NUM_THREADS=1

HMDIR=`pwd`/../../../..

ln -sf ${HMDIR}/bin/communication.exe .
ln -sf ${HMDIR}/kernels/NICAM/communication/data/communication.cnf .

mpirun -np 2 ./communication.exe
