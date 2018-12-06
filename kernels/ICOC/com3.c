#include <mpi.h>

extern int local_cell_blocks;
extern int local_edge_blocks;
#include "grid.h"
#include "memory.h"
#include "component.h"
#include "io.h"
#include <stdint.h>
void com3_init(GRID * g);
void com3_compute(GRID * g);
void com3_io(GRID * g);
double com3_flops(GRID * g);
double com3_memory(GRID * g);
uint64_t com3_checksum(GRID *);
void com3_cleanup(GRID * g);
void O2VertIntegration(GRID * g);
void O3Indirect2D(GRID * g);
void O4Indirect3D(GRID * g);
void O5Indirect3D(GRID * g);
struct {
    char *name;
    int loc;
    int dim;
    union {
        GVAL *restrict * restrict p2;
        GVAL *restrict * restrict * restrict p3;
    } data_pointer;
} *gv_vi;
struct {
    char *name;
    int loc;
    int dim;
    union {
        GVAL *restrict * restrict p2;
        GVAL *restrict * restrict * restrict p3;
    } data_pointer;
} *gv_ind2Dparam;
struct {
    char *name;
    int loc;
    int dim;
    union {
        GVAL *restrict * restrict p2;
        GVAL *restrict * restrict * restrict p3;
    } data_pointer;
} *gv_ind2Dvar;
struct {
    char *name;
    int loc;
    int dim;
    union {
        GVAL *restrict * restrict p2;
        GVAL *restrict * restrict * restrict p3;
    } data_pointer;
} *gv_ind3Dvar;
struct {
    char *name;
    int loc;
    int dim;
    union {
        int *restrict * restrict p2;
        int *restrict * restrict * restrict p3;
    } data_pointer;
} *t3DBlk;
struct {
    char *name;
    int loc;
    int dim;
    union {
        int *restrict * restrict p2;
        int *restrict * restrict * restrict p3;
    } data_pointer;
} *t3DIdx;
struct {
    char *name;
    int loc;
    int dim;
    union {
        int *restrict * restrict p2;
        int *restrict * restrict * restrict p3;
    } data_pointer;
} *t3DVer;
io_var_t io_gv_vi;
io_var_t io_gv_ind2Dvar;
io_var_t io_gv_ind3Dvar;
MODEL_COMPONENT com3 = { 0, com3_init, com3_compute, com3_io, com3_flops, com3_memory, com3_checksum, com3_cleanup };

extern struct {
    char *name;
    int loc;
    int dim;
    union {
        GVAL *restrict * restrict p2;
        GVAL *restrict * restrict * restrict p3;
    } data_pointer;
} *gv_grad;
void init_opO4(GRID * g)
{
    {
        int num_blocks = local_edge_blocks ? local_edge_blocks : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        t3DBlk = malloc(24);
        t3DBlk->name = "t3DBlk";
        t3DBlk->loc = 1;
        t3DBlk->dim = 3;
        t3DBlk->data_pointer.p3 = malloc((num_blocks * g->height * g->blkSize) * sizeof(int) + (num_blocks * g->height) * sizeof(char *) + (num_blocks) * sizeof(char *));
        char *pos = (char *) t3DBlk->data_pointer.p3 + num_blocks * sizeof(char *);
        char *pos2 = (char *) t3DBlk->data_pointer.p3 + num_blocks * sizeof(char *) + num_blocks * g->height * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            t3DBlk->data_pointer.p3[b] = (int * *) pos;
            pos += g->height * sizeof(char *);
            for (int k = 0; k < g->height; k++) {
                t3DBlk->data_pointer.p3[b][k] = (int *) pos2;
                pos2 += g->blkSize * sizeof(int);
                for (int e = 0; e < g->blkSize; e++) {
                    t3DBlk->data_pointer.p3[b][k][e] = (int) 0;
                }
            }
        }
    }
    {
        int num_blocks = local_edge_blocks ? local_edge_blocks : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        t3DIdx = malloc(24);
        t3DIdx->name = "t3DIdx";
        t3DIdx->loc = 1;
        t3DIdx->dim = 3;
        t3DIdx->data_pointer.p3 = malloc((num_blocks * g->height * g->blkSize) * sizeof(int) + (num_blocks * g->height) * sizeof(char *) + (num_blocks) * sizeof(char *));
        char *pos = (char *) t3DIdx->data_pointer.p3 + num_blocks * sizeof(char *);
        char *pos2 = (char *) t3DIdx->data_pointer.p3 + num_blocks * sizeof(char *) + num_blocks * g->height * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            t3DIdx->data_pointer.p3[b] = (int * *) pos;
            pos += g->height * sizeof(char *);
            for (int k = 0; k < g->height; k++) {
                t3DIdx->data_pointer.p3[b][k] = (int *) pos2;
                pos2 += g->blkSize * sizeof(int);
                for (int e = 0; e < g->blkSize; e++) {
                    t3DIdx->data_pointer.p3[b][k][e] = (int) 0;
                }
            }
        }
    }
    {
        size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t height_index = (0); height_index < (g->height); height_index++) {
                for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                    if (gv_grad->data_pointer.p3[(block_index)][(height_index)][(edge_index)] > 0) {
                        t3DBlk->data_pointer.p3[(block_index)][(height_index)][(edge_index)] = g->eCellBlk[0]->data_pointer.p2[(block_index)][(edge_index)];
                        t3DIdx->data_pointer.p3[(block_index)][(height_index)][(edge_index)] = g->eCellIdx[0]->data_pointer.p2[(block_index)][(edge_index)];
                    } else {
                        t3DBlk->data_pointer.p3[(block_index)][(height_index)][(edge_index)] = g->eCellBlk[1]->data_pointer.p2[(block_index)][(edge_index)];
                        t3DIdx->data_pointer.p3[(block_index)][(height_index)][(edge_index)] = g->eCellIdx[1]->data_pointer.p2[(block_index)][(edge_index)];
                    }
                }
            }
        }
    }
}

void init_opO5(GRID * g)
{
    {
        int num_blocks = local_edge_blocks ? local_edge_blocks : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        t3DVer = malloc(24);
        t3DVer->name = "t3DVer";
        t3DVer->loc = 1;
        t3DVer->dim = 3;
        t3DVer->data_pointer.p3 = malloc((num_blocks * g->height * g->blkSize) * sizeof(int) + (num_blocks * g->height) * sizeof(char *) + (num_blocks) * sizeof(char *));
        char *pos = (char *) t3DVer->data_pointer.p3 + num_blocks * sizeof(char *);
        char *pos2 = (char *) t3DVer->data_pointer.p3 + num_blocks * sizeof(char *) + num_blocks * g->height * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            t3DVer->data_pointer.p3[b] = (int * *) pos;
            pos += g->height * sizeof(char *);
            for (int k = 0; k < g->height; k++) {
                t3DVer->data_pointer.p3[b][k] = (int *) pos2;
                pos2 += g->blkSize * sizeof(int);
                for (int e = 0; e < g->blkSize; e++) {
                    t3DVer->data_pointer.p3[b][k][e] = (int) 0;
                }
            }
        }
    }
    {
        size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t height_index = (0); height_index < ((g->height - 2)); height_index++) {
                for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                    t3DVer->data_pointer.p3[(block_index)][(height_index)][(edge_index)] = height_index + 1;
                }
            }
        }
    }
    {
        size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            const size_t height_index = ((g->height - 1));
            {
                for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                    t3DVer->data_pointer.p3[(block_index)][(height_index)][(edge_index)] = (g->height - 1);
                }
            }
        }
    }
}

extern MODEL_COMPONENT com1;
void com3_init(GRID * g)
{
    com3.loaded = 1;
    if (!com1.loaded)
        com1.init(g);
    {
        int num_blocks = local_cell_blocks ? local_cell_blocks : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        gv_vi = malloc(24);
        gv_vi->name = "gv_vi";
        gv_vi->loc = 0;
        gv_vi->dim = 2;
        gv_vi->data_pointer.p2 = malloc((num_blocks * g->blkSize) * sizeof(GVAL) + (num_blocks) * sizeof(char *));
        char *pos = (char *) gv_vi->data_pointer.p2 + num_blocks * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            gv_vi->data_pointer.p2[b] = (GVAL *) pos;
            pos += g->blkSize * sizeof(GVAL);
            for (int c = 0; c < g->blkSize; c++) {
                gv_vi->data_pointer.p2[b][c] = (GVAL) 0;
            }
        }
    }
    {
        int num_blocks = local_edge_blocks ? local_edge_blocks : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        gv_ind2Dparam = malloc(24);
        gv_ind2Dparam->name = "gv_ind2Dparam";
        gv_ind2Dparam->loc = 1;
        gv_ind2Dparam->dim = 2;
        gv_ind2Dparam->data_pointer.p2 = malloc((num_blocks * g->blkSize) * sizeof(GVAL) + (num_blocks) * sizeof(char *));
        char *pos = (char *) gv_ind2Dparam->data_pointer.p2 + num_blocks * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            gv_ind2Dparam->data_pointer.p2[b] = (GVAL *) pos;
            pos += g->blkSize * sizeof(GVAL);
            for (int e = 0; e < g->blkSize; e++) {
                gv_ind2Dparam->data_pointer.p2[b][e] = (GVAL) 0;
            }
        }
    }
    {
        int num_blocks = local_edge_blocks ? local_edge_blocks : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        gv_ind2Dvar = malloc(24);
        gv_ind2Dvar->name = "gv_ind2Dvar";
        gv_ind2Dvar->loc = 1;
        gv_ind2Dvar->dim = 3;
        gv_ind2Dvar->data_pointer.p3 = malloc((num_blocks * g->height * g->blkSize) * sizeof(GVAL) + (num_blocks * g->height) * sizeof(char *) + (num_blocks) * sizeof(char *));
        char *pos = (char *) gv_ind2Dvar->data_pointer.p3 + num_blocks * sizeof(char *);
        char *pos2 = (char *) gv_ind2Dvar->data_pointer.p3 + num_blocks * sizeof(char *) + num_blocks * g->height * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            gv_ind2Dvar->data_pointer.p3[b] = (GVAL * *)pos;
            pos += g->height * sizeof(char *);
            for (int k = 0; k < g->height; k++) {
                gv_ind2Dvar->data_pointer.p3[b][k] = (GVAL *) pos2;
                pos2 += g->blkSize * sizeof(GVAL);
                for (int e = 0; e < g->blkSize; e++) {
                    gv_ind2Dvar->data_pointer.p3[b][k][e] = (GVAL) 0;
                }
            }
        }
    }
    {
        int num_blocks = local_edge_blocks ? local_edge_blocks : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        gv_ind3Dvar = malloc(24);
        gv_ind3Dvar->name = "gv_ind3Dvar";
        gv_ind3Dvar->loc = 1;
        gv_ind3Dvar->dim = 3;
        gv_ind3Dvar->data_pointer.p3 = malloc((num_blocks * g->height * g->blkSize) * sizeof(GVAL) + (num_blocks * g->height) * sizeof(char *) + (num_blocks) * sizeof(char *));
        char *pos = (char *) gv_ind3Dvar->data_pointer.p3 + num_blocks * sizeof(char *);
        char *pos2 = (char *) gv_ind3Dvar->data_pointer.p3 + num_blocks * sizeof(char *) + num_blocks * g->height * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            gv_ind3Dvar->data_pointer.p3[b] = (GVAL * *)pos;
            pos += g->height * sizeof(char *);
            for (int k = 0; k < g->height; k++) {
                gv_ind3Dvar->data_pointer.p3[b][k] = (GVAL *) pos2;
                pos2 += g->blkSize * sizeof(GVAL);
                for (int e = 0; e < g->blkSize; e++) {
                    gv_ind3Dvar->data_pointer.p3[b][k][e] = (GVAL) 0;
                }
            }
        }
    }
    init_opO4(g);
    init_opO5(g);
    io_read_register(g, "gv_ind2Dparam", (GVAL *) gv_ind2Dparam, FLOAT32, FLOAT32, GRID_POS_EDGE, GRID_DIM_2D);
    io_write_define(g, "gv_vi", (GVAL *) gv_vi, FLOAT32, GRID_POS_CELL, GRID_DIM_2D, &io_gv_vi);
    io_write_define(g, "gv_ind2Dvar", (GVAL *) gv_ind2Dvar, FLOAT32, GRID_POS_EDGE, GRID_DIM_3D, &io_gv_ind2Dvar);
    io_write_define(g, "gv_ind3Dvar", (GVAL *) gv_ind3Dvar, FLOAT32, GRID_POS_EDGE, GRID_DIM_3D, &io_gv_ind3Dvar);
}

void com3_compute(GRID * g)
{
    O2VertIntegration(g);
    O3Indirect2D(g);
    O4Indirect3D(g);
    O5Indirect3D(g);
}

void com3_io(GRID * g)
{
    io_write_announce(g, &io_gv_vi);
    io_write_announce(g, &io_gv_ind2Dvar);
    io_write_announce(g, &io_gv_ind3Dvar);
}

double com3_flops(GRID * g)
{
    double flop = (double) g->cellCount * (double) g->height + 3.0 * (double) g->edgeCount * (double) g->height + (double) g->edgeCount * (double) g->height;
    return flop;
}

double com3_memory(GRID * g)
{
    double mem = ((double) g->cellCount + (double) g->edgeCount + (double) g->edgeCount * (double) g->height + (double) g->edgeCount * (double) g->height) * (double) sizeof(GVAL) + (3.0 * (double) g->edgeCount * (double) g->height) * (double) sizeof(int);
    return mem / (1024 * 1024);
}

uint64_t com3_checksum(GRID * g)
{
    uint64_t ret = 0;
    {
        size_t min_block = g->mpi_rank == (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t cell_index = (0); cell_index < (g->blkSize); cell_index++) {
                ret += (uint64_t) gv_vi->data_pointer.p2[(block_index)][(cell_index)];
            }
        }
    }
    {
        size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t height_index = (0); height_index < (g->height); height_index++) {
                for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                    ret += (uint64_t) gv_ind2Dvar->data_pointer.p3[(block_index)][(height_index)][(edge_index)];
                    ret += (uint64_t) gv_ind3Dvar->data_pointer.p3[(block_index)][(height_index)][(edge_index)];
                }
            }
        }
    }
    return ret;
}

void com3_cleanup(GRID * g)
{
    com3.loaded = 0;
    free((void *) gv_vi->data_pointer.p2);
    free((void *) gv_ind2Dparam->data_pointer.p2);
    free((void *) gv_ind2Dvar->data_pointer.p3);
    free((void *) gv_ind3Dvar->data_pointer.p3);
    free((void *) t3DBlk->data_pointer.p3);
    free((void *) t3DIdx->data_pointer.p3);
    free((void *) t3DVer->data_pointer.p3);
}
