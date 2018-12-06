#include <mpi.h>

extern int local_cell_blocks;
extern int local_edge_blocks;
#include "grid.h"
#include "memory.h"
#include "component.h"
#include "io.h"
#include <stdint.h>
void com1_init(GRID * g);
void com1_compute(GRID * g);
void com1_io(GRID * g);
double com1_flops(GRID * g);
double com1_memory(GRID * g);
uint64_t com1_checksum(GRID *);
void com1_cleanup(GRID * g);
void grad(GRID * g);
void dvg(GRID * g);
void step(GRID * g);
struct {
    char *name;
    int loc;
    int dim;
    union {
        GVAL *restrict * restrict p2;
        GVAL *restrict * restrict * restrict p3;
    } data_pointer;
} *gv_temp;
struct {
    char *name;
    int loc;
    int dim;
    union {
        GVAL *restrict * restrict p2;
        GVAL *restrict * restrict * restrict p3;
    } data_pointer;
} *gv_grad;
struct {
    char *name;
    int loc;
    int dim;
    union {
        GVAL *restrict * restrict p2;
        GVAL *restrict * restrict * restrict p3;
    } data_pointer;
} *gv_dvg;
static io_var_t io_gv_temp;
static io_var_t io_gv_grad;
static io_var_t io_gv_dvg;
MODEL_COMPONENT com1 = { 0, com1_init, com1_compute, com1_io, com1_flops, com1_memory, com1_checksum, com1_cleanup };

void com1_init(GRID * g)
{
    com1.loaded = 1;
    {
        int num_blocks = local_cell_blocks ? local_cell_blocks : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        gv_temp = malloc(24);
        gv_temp->name = "gv_temp";
        gv_temp->loc = 0;
        gv_temp->dim = 3;
        gv_temp->data_pointer.p3 = malloc((num_blocks * g->height * g->blkSize) * sizeof(GVAL) + (num_blocks * g->height) * sizeof(char *) + (num_blocks) * sizeof(char *));
        char *pos = (char *) gv_temp->data_pointer.p3 + num_blocks * sizeof(char *);
        char *pos2 = (char *) gv_temp->data_pointer.p3 + num_blocks * sizeof(char *) + num_blocks * g->height * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            gv_temp->data_pointer.p3[b] = (GVAL * *)pos;
            pos += g->height * sizeof(char *);
            for (int k = 0; k < g->height; k++) {
                gv_temp->data_pointer.p3[b][k] = (GVAL *) pos2;
                pos2 += g->blkSize * sizeof(GVAL);
                for (int c = 0; c < g->blkSize; c++) {
                    gv_temp->data_pointer.p3[b][k][c] = (GVAL) 0;
                }
            }
        }
    }
    {
        int num_blocks = local_edge_blocks ? local_edge_blocks : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        gv_grad = malloc(24);
        gv_grad->name = "gv_grad";
        gv_grad->loc = 1;
        gv_grad->dim = 3;
        gv_grad->data_pointer.p3 = malloc((num_blocks * g->height * g->blkSize) * sizeof(GVAL) + (num_blocks * g->height) * sizeof(char *) + (num_blocks) * sizeof(char *));
        char *pos = (char *) gv_grad->data_pointer.p3 + num_blocks * sizeof(char *);
        char *pos2 = (char *) gv_grad->data_pointer.p3 + num_blocks * sizeof(char *) + num_blocks * g->height * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            gv_grad->data_pointer.p3[b] = (GVAL * *)pos;
            pos += g->height * sizeof(char *);
            for (int k = 0; k < g->height; k++) {
                gv_grad->data_pointer.p3[b][k] = (GVAL *) pos2;
                pos2 += g->blkSize * sizeof(GVAL);
                for (int e = 0; e < g->blkSize; e++) {
                    gv_grad->data_pointer.p3[b][k][e] = (GVAL) 0;
                }
            }
        }
    }
    {
        int num_blocks = local_cell_blocks ? local_cell_blocks : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        gv_dvg = malloc(24);
        gv_dvg->name = "gv_dvg";
        gv_dvg->loc = 0;
        gv_dvg->dim = 3;
        gv_dvg->data_pointer.p3 = malloc((num_blocks * g->height * g->blkSize) * sizeof(GVAL) + (num_blocks * g->height) * sizeof(char *) + (num_blocks) * sizeof(char *));
        char *pos = (char *) gv_dvg->data_pointer.p3 + num_blocks * sizeof(char *);
        char *pos2 = (char *) gv_dvg->data_pointer.p3 + num_blocks * sizeof(char *) + num_blocks * g->height * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            gv_dvg->data_pointer.p3[b] = (GVAL * *)pos;
            pos += g->height * sizeof(char *);
            for (int k = 0; k < g->height; k++) {
                gv_dvg->data_pointer.p3[b][k] = (GVAL *) pos2;
                pos2 += g->blkSize * sizeof(GVAL);
                for (int c = 0; c < g->blkSize; c++) {
                    gv_dvg->data_pointer.p3[b][k][c] = (GVAL) 0;
                }
            }
        }
    }
    io_read_register(g, "gv_temp", (GVAL *) gv_temp, FLOAT32, FLOAT32, GRID_POS_CELL, GRID_DIM_3D);
    io_write_define(g, "gv_temp", (GVAL *) gv_temp, FLOAT32, GRID_POS_CELL, GRID_DIM_3D, &io_gv_temp);
    io_write_define(g, "gv_grad", (GVAL *) gv_grad, FLOAT32, GRID_POS_EDGE, GRID_DIM_3D, &io_gv_grad);
    io_write_define(g, "gv_dvg", (GVAL *) gv_dvg, FLOAT32, GRID_POS_CELL, GRID_DIM_3D, &io_gv_dvg);
}

void com1_compute(GRID * g)
{
    grad(g);
    dvg(g);
    step(g);
}

void com1_io(GRID * g)
{
    io_write_announce(g, &io_gv_grad);
    io_write_announce(g, &io_gv_dvg);
}

double com1_flops(GRID * g)
{
    double flop = (double) g->edgeCount * (double) g->height + 2.0 * (double) NBRS * (double) g->cellCount * (double) g->height + 2.0 * (double) g->cellCount * (double) g->height;
    return flop;
}

double com1_memory(GRID * g)
{
    double mem = ((double) g->edgeCount * (double) g->height + 2.0 * (double) g->cellCount * (double) g->height) * sizeof(GVAL);
    return mem / (1024 * 1024);
}

uint64_t com1_checksum(GRID * g)
{
    uint64_t ret = 0;
    {
        size_t min_block = g->mpi_rank == (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t height_index = (0); height_index < (g->height); height_index++) {
                for (size_t cell_index = (0); cell_index < (g->blkSize); cell_index++) {
                    ret += (uint64_t) gv_temp->data_pointer.p3[(block_index)][(height_index)][(cell_index)];
                    ret += (uint64_t) gv_dvg->data_pointer.p3[(block_index)][(height_index)][(cell_index)];
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
                    ret += (uint64_t) gv_grad->data_pointer.p3[(block_index)][(height_index)][(edge_index)];
                }
            }
        }
    }
    return ret;
}

void com1_cleanup(GRID * g)
{
    com1.loaded = 0;
    free((void *) gv_temp->data_pointer.p3);
    free((void *) gv_grad->data_pointer.p3);
    free((void *) gv_dvg->data_pointer.p3);
}
