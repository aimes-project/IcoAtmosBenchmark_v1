#! /bin/bash

export OMP_NUM_THREADS=1

HMDIR=`pwd`/../../../..

ln -sf ${HMDIR}/bin/comp_pvort.exe .
ln -sf ${HMDIR}/kernels/DYNAMICO/comp_pvort/data/snapshot.comp_pvort.pe000000

./comp_pvort.exe
