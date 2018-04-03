#! /bin/bash

export OMP_NUM_THREADS=1

HMDIR=`pwd`/../../../..

ln -sf ${HMDIR}/bin/dyn_divdamp.exe .
ln -sf ${HMDIR}/kernels/NICAM/dyn_divdamp/data/vgrid40_600m_24km.dat .
ln -sf ${HMDIR}/kernels/NICAM/dyn_divdamp/data/snapshot.dyn_divdamp.pe000000 .

./dyn_divdamp.exe
