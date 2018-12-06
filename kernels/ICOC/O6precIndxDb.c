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
} *gv_grad;
extern GVAL *restrict * restrict * restrict gv_precInd;
extern int *restrict t6Blk;
extern int *restrict t6Ver;
extern int *restrict t6Ind;
void O6precIndxDb(GRID * g)
{
    {
        size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t height_index = (0); height_index < (g->height); height_index++) {
                for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                    int i = block_index * g->height * g->blkSize + height_index * g->blkSize + edge_index;
                    gv_precInd[t6Blk[i]][t6Ver[i]][t6Ind[i]] = 0.5 * gv_grad->data_pointer.p3[(block_index)][(height_index)][(edge_index)];
                }
            }
        }
    }
}
