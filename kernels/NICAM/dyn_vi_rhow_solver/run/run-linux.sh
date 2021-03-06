#! /bin/bash

export OMP_NUM_THREADS=1

HMDIR=`pwd`/../../../..

ln -sf ${HMDIR}/bin/dyn_vi_rhow_solver.exe .
ln -sf ${HMDIR}/kernels/NICAM/dyn_vi_rhow_solver/data/vgrid40_600m_24km.dat .
ln -sf ${HMDIR}/kernels/NICAM/dyn_vi_rhow_solver/data/snapshot.dyn_vi_rhow_solver.pe000000 .

./dyn_vi_rhow_solver.exe
