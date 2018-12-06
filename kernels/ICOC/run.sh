#!/bin/bash
DTYPE=float
BLOCK=64
TILE=3
echo "using $TILE-side tiled grid (Data type: $DTYPE, Block size: $BLOCK)"
mpicc *.c -o a.out -DGVAL=$DTYPE -DNBRS=$TILE -DBLKSIZE=$BLOCK -O3 -lnetcdf
mpiexec -n 1  --map-by node --display-map  ./a.out -G=1024 -m=0 -M=1000 -T=10 -R=input.cdf -W=output.cdf

