#include <mpi.h>

#include "grid.h"
extern struct {
    char *name;
    int loc;
    int dim;
    union {
        GVAL *restrict * restrict p2;
        GVAL *restrict * restrict * restrict p3;
    } data_pointer;
} *gv_temp;
extern struct {
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
void O4Indirect3D(GRID * g)
{
    {
        size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t height_index = (0); height_index < (g->height); height_index++) {
                for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                    gv_ind3Dvar->data_pointer.p3[(block_index)][(height_index)][(edge_index)] = 0;
                }
            }
        }
    }
}
