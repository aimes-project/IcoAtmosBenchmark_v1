#! /bin/bash
#
# for K computer
#
#PJM --rsc-list "rscgrp=micro"
#PJM --rsc-list "node=2"
#PJM --rsc-list "elapse=00:10:00"
#PJM -j
#PJM -s
#
. /work/system/Env_base
#
export PARALLEL=8
export OMP_NUM_THREADS=8
export XOS_MMM_L_ARENA_FREE=2

HMDIR=`pwd`/../../../..

ln -sf ${HMDIR}/bin/communication.exe .
ln -sf ${HMDIR}/kernels/NICAM/communication/data/communication.cnf .

rm -rf ./prof*

fapp -C -Ihwm -Hevent=Statistics -d prof -L 10 mpiexec ./communication.exe
