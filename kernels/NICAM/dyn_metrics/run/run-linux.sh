#! /bin/bash

export OMP_NUM_THREADS=1

HMDIR=`pwd`/../../../..

ln -sf ${HMDIR}/bin/dyn_metrics.exe .
ln -sf ${HMDIR}/kernels/NICAM/dyn_metrics/data/snapshot.dyn_metrics.pe000000 .

./dyn_metrics.exe
