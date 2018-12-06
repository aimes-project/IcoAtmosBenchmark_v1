#include <mpi.h>

extern int local_cell_blocks;
extern int local_edge_blocks;
#include "grid.h"
#include "memory.h"
#include "component.h"
#include "io.h"
#include <stdint.h>
void com4_init(GRID * g);
void com4_compute(GRID * g);
void com4_io(GRID * g);
double com4_flops(GRID * g);
double com4_memory(GRID * g);
uint64_t com4_checksum(GRID *);
void com4_cleanup(GRID * g);
void grad(GRID * g);
void O6precIndxDb(GRID * g);
void O7precIndxNb(GRID * g);
void O8ind(GRID * g);
GVAL *restrict * restrict * restrict gv_precInd;
int *restrict t6Blk;
int *restrict t6Ver;
int *restrict t6Ind;
int *restrict t7Blk;
int *restrict t7Ind;
struct {
    char *name;
    int loc;
    int dim;
    union {
        GVAL *restrict * restrict p2;
        GVAL *restrict * restrict * restrict p3;
    } data_pointer;
} *gv_o8param[3];
struct {
    char *name;
    int loc;
    int dim;
    union {
        GVAL *restrict * restrict p2;
        GVAL *restrict * restrict * restrict p3;
    } data_pointer;
} *gv_o8par2;
struct {
    char *name;
    int loc;
    int dim;
    union {
        GVAL *restrict * restrict p2;
        GVAL *restrict * restrict * restrict p3;
    } data_pointer;
} *gv_o8var;
io_var_t io_gv_precInd;
io_var_t io_gv_o8var;
MODEL_COMPONENT com4 = { 0, com4_init, com4_compute, com4_io, com4_flops, com4_memory, com4_checksum, com4_cleanup };

void init_opO6(GRID * g)
{
    t6Blk = malloc((g->eBlkCnt * g->height * g->blkSize) * sizeof(int));
    t6Ver = malloc((g->eBlkCnt * g->height * g->blkSize) * sizeof(int));
    t6Ind = malloc((g->eBlkCnt * g->height * g->blkSize) * sizeof(int));
    {
        size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t height_index = (0); height_index < (g->height); height_index++) {
                for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                    t6Blk[block_index * g->height * g->blkSize + height_index * g->blkSize + edge_index] = block_index;
                    t6Ver[block_index * g->height * g->blkSize + height_index * g->blkSize + edge_index] = height_index;
                    t6Ind[block_index * g->height * g->blkSize + height_index * g->blkSize + edge_index] = edge_index;
                }
            }
        }
    }
    int eb, ee;
#pragma omp parallel for
    for (int b = 0; b < g->eBlkCnt; b++) {
        get_indices_e(g, b, &eb, &ee);
        for (int k = 0; k < g->height; k++) {
            for (int e = eb; e < ee; e++) {
                t6Blk[b * g->height * g->blkSize + k * ee + e] = b;
                t6Ver[b * g->height * g->blkSize + k * ee + e] = k;
                t6Ind[b * g->height * g->blkSize + k * ee + e] = e;
            }
        }
    }
}

void init_opO7(GRID * g)
{
    t7Blk = malloc((g->edgeCount) * sizeof(int));
    t7Ind = malloc((g->edgeCount) * sizeof(int));
    {
        size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                t7Blk[block_index * g->blkSize + edge_index] = block_index;
                t7Ind[block_index * g->blkSize + edge_index] = edge_index;
            }
        }
    }
}

extern MODEL_COMPONENT com1;
void com4_init(GRID * g)
{
    com4.loaded = 1;
    if (!com1.loaded)
        com1.init(g);
    gv_precInd = (GVAL * restrict * restrict * restrict) allocate_variable(g, GRID_POS_EDGE, GRID_DIM_3D, sizeof(GVAL));
    init_opO6(g);
    init_opO7(g);
    for (int i = 0; i < 3; i++) {
        int num_blocks = local_cell_blocks ? local_cell_blocks : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        gv_o8param[i] = malloc(24);
        gv_o8param[i]->name = "gv_o8param[ i]";
        gv_o8param[i]->loc = 0;
        gv_o8param[i]->dim = 2;
        gv_o8param[i]->data_pointer.p2 = malloc((num_blocks * g->blkSize) * sizeof(GVAL) + (num_blocks) * sizeof(char *));
        char *pos = (char *) gv_o8param[i]->data_pointer.p2 + num_blocks * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            gv_o8param[i]->data_pointer.p2[b] = (GVAL *) pos;
            pos += g->blkSize * sizeof(GVAL);
            for (int c = 0; c < g->blkSize; c++) {
                gv_o8param[i]->data_pointer.p2[b][c] = (GVAL) 0;
            }
        }
    }
    {
        int num_blocks = local_cell_blocks ? local_cell_blocks : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        gv_o8par2 = malloc(24);
        gv_o8par2->name = "gv_o8par2";
        gv_o8par2->loc = 0;
        gv_o8par2->dim = 3;
        gv_o8par2->data_pointer.p3 = malloc((num_blocks * g->height * g->blkSize) * sizeof(GVAL) + (num_blocks * g->height) * sizeof(char *) + (num_blocks) * sizeof(char *));
        char *pos = (char *) gv_o8par2->data_pointer.p3 + num_blocks * sizeof(char *);
        char *pos2 = (char *) gv_o8par2->data_pointer.p3 + num_blocks * sizeof(char *) + num_blocks * g->height * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            gv_o8par2->data_pointer.p3[b] = (GVAL * *)pos;
            pos += g->height * sizeof(char *);
            for (int k = 0; k < g->height; k++) {
                gv_o8par2->data_pointer.p3[b][k] = (GVAL *) pos2;
                pos2 += g->blkSize * sizeof(GVAL);
                for (int c = 0; c < g->blkSize; c++) {
                    gv_o8par2->data_pointer.p3[b][k][c] = (GVAL) 0;
                }
            }
        }
    }
    {
        int num_blocks = local_cell_blocks ? local_cell_blocks : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        gv_o8var = malloc(24);
        gv_o8var->name = "gv_o8var";
        gv_o8var->loc = 0;
        gv_o8var->dim = 3;
        gv_o8var->data_pointer.p3 = malloc((num_blocks * g->height * g->blkSize) * sizeof(GVAL) + (num_blocks * g->height) * sizeof(char *) + (num_blocks) * sizeof(char *));
        char *pos = (char *) gv_o8var->data_pointer.p3 + num_blocks * sizeof(char *);
        char *pos2 = (char *) gv_o8var->data_pointer.p3 + num_blocks * sizeof(char *) + num_blocks * g->height * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            gv_o8var->data_pointer.p3[b] = (GVAL * *)pos;
            pos += g->height * sizeof(char *);
            for (int k = 0; k < g->height; k++) {
                gv_o8var->data_pointer.p3[b][k] = (GVAL *) pos2;
                pos2 += g->blkSize * sizeof(GVAL);
                for (int c = 0; c < g->blkSize; c++) {
                    gv_o8var->data_pointer.p3[b][k][c] = (GVAL) 0;
                }
            }
        }
    }
    io_read_register(g, "gv_o8param0", (GVAL *) gv_o8param[0], FLOAT32, FLOAT32, GRID_POS_CELL, GRID_DIM_2D);
    io_read_register(g, "gv_o8param1", (GVAL *) gv_o8param[1], FLOAT32, FLOAT32, GRID_POS_CELL, GRID_DIM_2D);
    io_read_register(g, "gv_o8param2", (GVAL *) gv_o8param[2], FLOAT32, FLOAT32, GRID_POS_CELL, GRID_DIM_2D);
    io_read_register(g, "gv_o8par2", (GVAL *) gv_o8par2, FLOAT32, FLOAT32, GRID_POS_CELL, GRID_DIM_3D);
    io_write_define(g, "gv_precInd", (GVAL *) gv_precInd, FLOAT32, GRID_POS_EDGE, GRID_DIM_3D, &io_gv_precInd);
    io_write_define(g, "gv_o8var", (GVAL *) gv_o8var, FLOAT32, GRID_POS_CELL, GRID_DIM_3D, &io_gv_o8var);
}

void com4_compute(GRID * g)
{
    grad(g);
    O6precIndxDb(g);
    O7precIndxNb(g);
    O8ind(g);
}

void com4_io(GRID * g)
{
    io_write_announce(g, &io_gv_precInd);
    io_write_announce(g, &io_gv_o8var);
}

double com4_flops(GRID * g)
{
    double flop = (double) g->edgeCount * (double) g->height + (double) g->edgeCount * (double) g->height + (double) g->edgeCount * (double) g->height + 14.0 * (double) g->cellCount * (double) (g->height - 1) + 20.0 * (double) g->cellCount;
    return flop;
}

double com4_memory(GRID * g)
{
    double mem = ((double) g->edgeCount * (double) g->height + 3.0 * (double) g->cellCount + (double) g->cellCount * (double) g->height + (double) g->cellCount * (double) g->height) * (double) sizeof(GVAL) + (3.0 * (double) g->eBlkCnt * (double) g->height * (double) g->blkSize + 2.0 * (double) g->edgeCount) * (double) sizeof(int);
    return mem / (1024 * 1024);
}

uint64_t com4_checksum(GRID * g)
{
    uint64_t ret = 0;
    {
        size_t min_block = g->mpi_rank == (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t height_index = (0); height_index < (g->height); height_index++) {
                for (size_t cell_index = (0); cell_index < (g->blkSize); cell_index++) {
                    ret += (uint64_t) gv_o8var->data_pointer.p3[(block_index)][(height_index)][(cell_index)];
                }
            }
        }
    }
    return ret;
}

void com4_cleanup(GRID * g)
{
    com4.loaded = 0;
    deallocate_variable((void *) gv_precInd, g, GRID_POS_EDGE, GRID_DIM_3D, sizeof(GVAL));
    free((void *) t6Blk);
    free((void *) t6Ver);
    free((void *) t6Ind);
    free((void *) t7Blk);
    free((void *) t7Ind);
    for (int i = 0; i < 3; i++)
        free((void *) gv_o8param[i]->data_pointer.p2);
    free((void *) gv_o8par2->data_pointer.p3);
    free((void *) gv_o8var->data_pointer.p3);
}
