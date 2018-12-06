#ifndef COMPONENT_H

#define COMPONENT_H
#include <stdlib.h>
#include <stdio.h>
#include "grid.h"
typedef struct {
    int loaded;
    void (*init) (GRID *);
    void (*compute) (GRID *);
    void (*io) (GRID *);
    double (*flops) (GRID *);
    double (*memory) (GRID *);
     uint64_t(*checksum) (GRID *);
    void (*cleanup) (GRID *);
} MODEL_COMPONENT;
#endif
