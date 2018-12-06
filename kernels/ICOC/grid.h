#ifndef GRID_H

#define GRID_H
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#define GVAL float
#define NBRS 3
//#define BLKSIZE 4
typedef struct {
    int i;
    int j;
} map_t;
typedef struct {
    int cellCount;
    int height;
    int edgeCount;
    int blkSize;
    int cBlkCnt;
    int eBlkCnt;
    int *restrict neighbor[NBRS];
    struct {
        char *name;
        int loc;
        int dim;
        union {
            int *restrict * restrict p2;
            int *restrict * restrict * restrict p3;
        } data_pointer;
    } *cNeighborIdx[NBRS];
    struct {
        char *name;
        int loc;
        int dim;
        union {
            int *restrict * restrict p2;
            int *restrict * restrict * restrict p3;
        } data_pointer;
    } *cNeighborBlk[NBRS];
    int *restrict cedge[NBRS];
    struct {
        char *name;
        int loc;
        int dim;
        union {
            int *restrict * restrict p2;
            int *restrict * restrict * restrict p3;
        } data_pointer;
    } *cEdgeIdx[NBRS];
    struct {
        char *name;
        int loc;
        int dim;
        union {
            int *restrict * restrict p2;
            int *restrict * restrict * restrict p3;
        } data_pointer;
    } *cEdgeBlk[NBRS];
    int *restrict ecell[2];
    struct {
        char *name;
        int loc;
        int dim;
        union {
            int *restrict * restrict p2;
            int *restrict * restrict * restrict p3;
        } data_pointer;
    } *eCellIdx[2];
    struct {
        char *name;
        int loc;
        int dim;
        union {
            int *restrict * restrict p2;
            int *restrict * restrict * restrict p3;
        } data_pointer;
    } *eCellBlk[2];
    int *restrict * restrict map1;
//in fact no need to keep
    map_t *restrict map2;
    GVAL edge_weights[NBRS][BLKSIZE];
    int mpi_world_size;
    int mpi_rank;
} GRID;
void init_grid(GRID * g, int cellCount, int height);
void get_indices_c(GRID * g, int blk, int *cb, int *ce);
void get_indices_e(GRID * g, int blk, int *eb, int *ee);
#endif
