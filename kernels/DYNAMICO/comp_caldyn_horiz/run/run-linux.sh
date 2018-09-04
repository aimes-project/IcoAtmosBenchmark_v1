#! /bin/bash

export OMP_NUM_THREADS=1

HMDIR=`pwd`/../../../..

ln -sf ${HMDIR}/bin/comp_caldyn_horiz.exe .
ln -sf ${HMDIR}/kernels/DYNAMICO/comp_caldyn_horiz/data/snapshot.comp_caldyn_horiz.pe000000

./comp_caldyn_horiz.exe
