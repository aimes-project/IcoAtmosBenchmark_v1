#! /bin/bash

export OMP_NUM_THREADS=1

HMDIR=`pwd`/../../../..

ln -sf ${HMDIR}/bin/dyn_horiz_adv_limiter.exe .
ln -sf ${HMDIR}/kernels/NICAM/dyn_horiz_adv_limiter/data/snapshot.dyn_horiz_adv_limiter.pe000000 .

./dyn_horiz_adv_limiter.exe
