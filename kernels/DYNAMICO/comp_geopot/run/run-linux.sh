#! /bin/bash

export OMP_NUM_THREADS=1

HMDIR=`pwd`/../../../..

ln -sf ${HMDIR}/bin/comp_geopot.exe .
ln -sf ${HMDIR}/kernels/DYNAMICO/comp_geopot/data/snapshot.comp_geopot.pe000000

./comp_geopot.exe

