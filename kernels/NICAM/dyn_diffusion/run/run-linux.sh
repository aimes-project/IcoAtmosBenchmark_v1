#! /bin/bash

export OMP_NUM_THREADS=1

HMDIR=`pwd`/../../../..

ln -sf ${HMDIR}/bin/dyn_diffusion.exe .
ln -sf ${HMDIR}/kernels/NICAM/dyn_diffusion/data/snapshot.dyn_diffusion.pe000000 .

./dyn_diffusion.exe
