#ifndef MEMORY_H

#define MEMORY_H
#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
typedef enum { GRID_POS_CELL, GRID_POS_EDGE, GRID_POS_VERTEX } GRID_VAR_POSITION;
typedef enum { GRID_DIM_2D, GRID_DIM_3D } GRID_VAR_DIMENSION;
void *allocate_variable(GRID *, GRID_VAR_POSITION, GRID_VAR_DIMENSION, size_t);
void *deallocate_variable(void *, GRID *, GRID_VAR_POSITION, GRID_VAR_DIMENSION, size_t);
#endif
