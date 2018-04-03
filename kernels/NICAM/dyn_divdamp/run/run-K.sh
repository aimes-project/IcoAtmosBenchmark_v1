#! /bin/bash -x
#
# for K computer
#
#PJM --rsc-list "rscgrp=micro"
#PJM --rsc-list "node=1"
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

ln -svf ${HMDIR}/bin/dyn_divdamp.exe .
ln -svf ${HMDIR}/kernels/NICAM/dyn_divdamp/data/vgrid40_600m_24km.dat .
ln -svf ${HMDIR}/kernels/NICAM/dyn_divdamp/data/snapshot.dyn_divdamp.pe000000 .

rm -rf ./prof*

fapp -C -Ihwm -Hevent=Statistics -d prof -L 10 ./dyn_divdamp.exe
