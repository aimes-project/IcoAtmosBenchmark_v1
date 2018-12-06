#include <mpi.h>

int *cn_c;
int *ce_c;
int *ec_c;
int *cn_crem;
int *ce_crem;
int *ec_crem;
int *neighbor_map;
int *cedge_map;
int *ecell_map;
int *neighbor_maprem;
int *cedge_maprem;
int *ecell_maprem;
GVAL **neighbor_2Dbuf;
GVAL **neighbor_3Dbuf;
GVAL **cedge_2Dbuf;
GVAL **cedge_3Dbuf;
GVAL **ecell_2Dbuf;
GVAL **ecell_3Dbuf;
GVAL **neighbor_2Dbufrem;
GVAL **neighbor_3Dbufrem;
GVAL **cedge_2Dbufrem;
GVAL **cedge_3Dbufrem;
GVAL **ecell_2Dbufrem;
GVAL **ecell_3Dbufrem;
MPI_Request *mpi_send_requests;
MPI_Request *mpi_recv_requests;
int comm_tag;
int local_cell_blocks;
int local_edge_blocks;
#include "grid.h"
#include "memory.h"
int transform(int n, int x, int y)
{
    int rx, ry, s, d = 0, t;
    for (s = n / 2; s > 0; s /= 2) {
        rx = (x & s) > 0;
        ry = (y & s) > 0;
        d += s * s * ((3 * rx) ^ ry);
        if (ry == 0) {
            if (rx) {
                x = n - 1 - x;
                y = n - 1 - y;
            }
            t = x;
            x = y;
            y = t;
        }
    }
    return d;
}

void create_maps(GRID * g)
{
    int x = g->height, y = g->height;
    g->map1 = malloc(x * sizeof(int *));
    g->map2 = malloc(g->cellCount * sizeof(map_t));
    for (int i = 0; i < x; i++) {
        g->map1[i] = malloc(y * sizeof(int));
        for (int j = 0; j < y; j++) {
            int t = transform(x, i, j);
            g->map1[i][j] = t;
            g->map2[t].i = i;
            g->map2[t].j = j;
        }
    }
}

#ifndef NBRS
#error "please define NBRS"
#endif
#if NBRS==3
int calc_edge_count(GRID * g)
{
    return (g->cellCount * 3) / 2;
}

void tessellation(GRID * g)
{
    for (int i = 0; i < NBRS; i++) {
        g->neighbor[i] = malloc((g->cellCount) * sizeof(int));
        g->cedge[i] = malloc((g->cellCount) * sizeof(int));
    }
    for (int i = 0; i < 2; i++) {
        g->ecell[i] = malloc((g->edgeCount) * sizeof(int));
    }
    int x = g->height, y = g->height;
    for (int i = 0; i < x - 1; i++)
        for (int j = 0; j < y; j++)
            g->neighbor[0][g->map1[i][j]] = g->map1[i + 1][j];
    for (int j = 0; j < y; j++)
        g->neighbor[0][g->map1[x - 1][j]] = g->map1[0][j];
    for (int i = 1; i < x; i++)
        for (int j = 0; j < y; j++)
            g->neighbor[1][g->map1[i][j]] = g->map1[i - 1][j];
    for (int j = 0; j < y; j++)
        g->neighbor[1][g->map1[0][j]] = g->map1[x - 1][j];
    for (int i = 0; i < x; i += 2)
        g->neighbor[2][g->map1[i][0]] = g->map1[i][y - 1];
    for (int i = 0; i < x; i += 2)
        for (int j = 2; j < y; j += 2)
            g->neighbor[2][g->map1[i][j]] = g->map1[i][j - 1];
    for (int i = 1; i < x; i += 2)
        for (int j = 1; j < y; j += 2)
            g->neighbor[2][g->map1[i][j]] = g->map1[i][j - 1];
    for (int i = 1; i < x; i += 2)
        for (int j = 0; j < y - 1; j += 2)
            g->neighbor[2][g->map1[i][j]] = g->map1[i][j + 1];
    for (int i = 0; i < x; i += 2)
        for (int j = 1; j < y - 1; j += 2)
            g->neighbor[2][g->map1[i][j]] = g->map1[i][j + 1];
    for (int i = y % 2; i < x; i += 2)
        g->neighbor[2][g->map1[i][y - 1]] = g->map1[i][0];
    for (int c = 0; c < g->cellCount; c++) {
        g->cedge[0][c] = (c * 3) / 2;
        g->cedge[1][g->neighbor[0][c]] = g->cedge[0][c];
        g->ecell[0][g->cedge[0][c]] = g->neighbor[0][c];
        g->ecell[1][g->cedge[0][c]] = c;
    }
    for (int c = 0; c < g->cellCount; c += 2) {
        g->cedge[2][c] = (c * 3) / 2 + 2;
        g->cedge[2][g->neighbor[2][c]] = g->cedge[2][c];
        g->ecell[0][g->cedge[2][c]] = c;
        g->ecell[1][g->cedge[2][c]] = g->neighbor[2][c];
    }
}

void init_edge_weights(GRID * g)
{
    GVAL j = -1.0;
    for (int i = 0; i < BLKSIZE; i++) {
        g->edge_weights[0][i] = 1.0;
        g->edge_weights[1][i] = -1.0;
        g->edge_weights[2][i] = j;
        j = j * -1;
    }
}
#elif NBRS==4
int calc_edge_count(GRID * g)
{
    return (int) (((float) g->cellCount * NBRS / 2.0) + (2.0 * g->height));
}

void tessellation(GRID * g)
{
    for (int i = 0; i < NBRS; i++) {
        g->neighbor[i] = malloc((g->cellCount) * sizeof(int));
        g->cedge[i] = malloc((g->cellCount) * sizeof(int));
    }
    for (int i = 0; i < 2; i++) {
        g->ecell[i] = malloc((g->edgeCount) * sizeof(int));
    }
    int x = g->height, y = g->height;
    for (int i = 0; i < x; i++) {
//
        for (int j = 0; j < y; j++) {
            if (i < x - 1)
                g->neighbor[0][g->map1[i][j]] = g->map1[i + 1][j];
            else
                g->neighbor[0][g->map1[i][j]] = g->cellCount;
            if (j < y - 1)
                g->neighbor[1][g->map1[i][j]] = g->map1[i][j + 1];
            else
                g->neighbor[1][g->map1[i][j]] = g->cellCount;
            if (i > 0)
                g->neighbor[2][g->map1[i][j]] = g->map1[i - 1][j];
            else
                g->neighbor[2][g->map1[i][j]] = g->cellCount;
            if (j > 0)
                g->neighbor[3][g->map1[i][j]] = g->map1[i][j - 1];
            else
                g->neighbor[3][g->map1[i][j]] = g->cellCount;
            g->neighbor[0][g->cellCount] = g->neighbor[1][g->cellCount] = g->neighbor[2][g->cellCount] = g->neighbor[3][g->cellCount] = g->cellCount;
//
            g->cedge[0][g->map1[i][j]] = g->map1[i][j];
            g->cedge[1][g->map1[i][j]] = g->cellCount + g->map1[i][j];
            if (i > 0)
                g->cedge[2][g->map1[i][j]] = g->map1[i - 1][j];
            else
                g->cedge[2][g->map1[i][j]] = g->cellCount * 2 + j;
            if (j > 0)
                g->cedge[3][g->map1[i][j]] = g->cellCount + g->map1[i][j - 1];
            else
                g->cedge[3][g->map1[i][j]] = g->cellCount * 2 + y + i;
            g->ecell[1][g->map1[i][j]] = g->map1[i][j];
            g->ecell[1][g->cellCount + g->map1[i][j]] = g->map1[i][j];
            if (i > 0)
                g->ecell[0][g->map1[i - 1][j]] = g->map1[i][j];
            else
                g->ecell[0][g->cellCount * 2 + j] = g->map1[i][j];
            if (j > 0)
                g->ecell[0][g->cellCount + g->map1[i][j - 1]] = g->map1[i][j];
            else
                g->ecell[0][g->cellCount * 2 + y + i] = g->map1[i][j];
            g->ecell[0][g->map1[x - 1][j]] = g->cellCount;
//TODO: out of loop
            g->ecell[0][g->cellCount + g->map1[i][y - 1]] = g->cellCount;
//
            g->ecell[1][g->cellCount * 2 + j] = g->cellCount;
//TODO: out of loop
            g->ecell[1][g->cellCount * 2 + y + i] = g->cellCount;
        }
    }
}

//
void init_edge_weights(GRID * g)
{
    for (int i = 0; i < BLKSIZE; i++) {
        g->edge_weights[0][i] = 1.0;
        g->edge_weights[1][i] = 1.0;
        g->edge_weights[2][i] = -1.0;
        g->edge_weights[3][i] = -1.0;
    }
}
#elif NBRS==6
int calc_edge_count(GRID * g)
{
    return (int) (((float) g->cellCount * NBRS / 2.0) + (4.0 * g->height) - 1);
}

void tessellation(GRID * g)
{
    for (int i = 0; i < NBRS; i++) {
        g->neighbor[i] = malloc((g->cellCount) * sizeof(int));
        g->cedge[i] = malloc((g->cellCount) * sizeof(int));
    }
    for (int i = 0; i < 2; i++) {
        g->ecell[i] = malloc((g->edgeCount) * sizeof(int));
    }
    int x = g->height, y = g->height;
    for (int i = 0; i < x; i++)
        for (int j = 0; j < y; j += 2)
            g->neighbor[0][g->map1[i][j]] = g->map1[i][j + 1];
    for (int i = 0; i < x - 1; i++)
        for (int j = 1; j < y - 1; j += 2)
            g->neighbor[0][g->map1[i][j]] = g->map1[i + 1][j + 1];
    for (int j = 1; j < y - 1; j += 2)
        g->neighbor[0][g->map1[x - 1][j]] = g->cellCount;
    for (int i = 0; i < x; i++)
        g->neighbor[0][g->map1[i][y - 1]] = g->cellCount;
    for (int i = 0; i < x - 1; i++)
        for (int j = 0; j < y; j++)
            g->neighbor[1][g->map1[i][j]] = g->map1[i + 1][j];
    for (int j = 0; j < y; j++)
        g->neighbor[1][g->map1[x - 1][j]] = g->cellCount;
    for (int i = 0; i < x; i++)
        g->neighbor[2][g->map1[i][0]] = g->cellCount;
    for (int i = 0; i < x; i++)
        for (int j = 2; j < y; j += 2)
            g->neighbor[2][g->map1[i][j]] = g->map1[i][j - 1];
    for (int i = 0; i < x - 1; i++)
        for (int j = 1; j < y; j += 2)
            g->neighbor[2][g->map1[i][j]] = g->map1[i + 1][j - 1];
    for (int j = 1; j < y; j += 2)
        g->neighbor[2][g->map1[x - 1][j]] = g->cellCount;
    for (int i = 0; i < x; i++)
        g->neighbor[3][g->map1[i][0]] = g->cellCount;
    for (int i = 1; i < x; i++)
        for (int j = 2; j < y; j += 2)
            g->neighbor[3][g->map1[i][j]] = g->map1[i - 1][j - 1];
    for (int j = 2; j < y; j += 2)
        g->neighbor[3][g->map1[0][j]] = g->cellCount;
    for (int i = 0; i < x; i++)
        for (int j = 1; j < y; j += 2)
            g->neighbor[3][g->map1[i][j]] = g->map1[i][j - 1];
    for (int i = 1; i < x; i++)
        for (int j = 0; j < y; j++)
            g->neighbor[4][g->map1[i][j]] = g->map1[i - 1][j];
    for (int j = 0; j < y; j++)
        g->neighbor[4][g->map1[0][j]] = g->cellCount;
    for (int i = 1; i < x; i++)
        for (int j = 0; j < y; j += 2)
            g->neighbor[5][g->map1[i][j]] = g->map1[i - 1][j + 1];
    for (int j = 0; j < y; j += 2)
        g->neighbor[5][g->map1[0][j]] = g->cellCount;
    for (int i = 0; i < x; i++)
        for (int j = 1; j < y - 1; j += 2)
            g->neighbor[5][g->map1[i][j]] = g->map1[i][j + 1];
    for (int i = 0; i < x; i++)
        g->neighbor[5][g->map1[i][y - 1]] = g->cellCount;
    for (int i = 0; i < x; i++)
        for (int j = 0; j < y; j++) {
            g->cedge[0][g->map1[i][j]] = g->map1[i][j];
            g->cedge[1][g->map1[i][j]] = g->cellCount + g->map1[i][j];
            g->cedge[2][g->map1[i][j]] = g->cellCount * 2 + g->map1[i][j];
        }
    for (int i = 0; i < x; i++)
        for (int j = 0; j < y; j++) {
            if (j == 0 || ((i == 0) && (j % 2 == 0)))
                g->cedge[3][g->map1[i][j]] = g->cellCount * 3 + (j == 0 ? i : x + j / 2 - 1);
            else
                g->cedge[3][g->map1[i][j]] = g->cedge[0][g->neighbor[3][g->map1[i][j]]];
            if (i == 0)
                g->cedge[4][g->map1[i][j]] = g->cellCount * 3 + x + x / 2 + j - 1;
            else
                g->cedge[4][g->map1[i][j]] = g->cedge[1][g->neighbor[4][g->map1[i][j]]];
            if ((j == y - 1) || ((i == 0) && (j % 2 == 0)))
                g->cedge[5][g->map1[i][j]] = g->cellCount * 3 + x + x / 2 + y - 1 + (j == y - 1 ? y / 2 + i : j / 2);
            else
                g->cedge[5][g->map1[i][j]] = g->cedge[2][g->neighbor[5][g->map1[i][j]]];
        }
    for (int c = 0; c < g->cellCount; c++) {
        g->ecell[0][g->cedge[0][c]] = g->neighbor[0][c];
        g->ecell[1][g->cedge[0][c]] = c;
        g->ecell[0][g->cedge[1][c]] = g->neighbor[1][c];
        g->ecell[1][g->cedge[1][c]] = c;
        g->ecell[0][g->cedge[2][c]] = g->neighbor[2][c];
        g->ecell[1][g->cedge[2][c]] = c;
        g->ecell[0][g->cedge[3][c]] = c;
        g->ecell[1][g->cedge[3][c]] = g->neighbor[3][c];
        g->ecell[0][g->cedge[4][c]] = c;
        g->ecell[1][g->cedge[4][c]] = g->neighbor[4][c];
        g->ecell[0][g->cedge[5][c]] = c;
        g->ecell[1][g->cedge[5][c]] = g->neighbor[5][c];
    }
}

void init_edge_weights(GRID * g)
{
    for (int i = 0; i < BLKSIZE; i++) {
        g->edge_weights[0][i] = g->edge_weights[1][i] = g->edge_weights[2][i] = 1.0;
        g->edge_weights[3][i] = g->edge_weights[4][i] = g->edge_weights[5][i] = -1.0;
    }
}
#else
#error "supported shapes are traiangles,rectangles and hexagons"
#endif
void init_blocking(GRID * g)
{
    for (int i = 0; i < NBRS; i++) {
        {
            int num_blocks = local_cell_blocks ? local_cell_blocks : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
            g->cNeighborIdx[i] = malloc(24);
            g->cNeighborIdx[i]->name = "g->cNeighborIdx[ i]";
            g->cNeighborIdx[i]->loc = 0;
            g->cNeighborIdx[i]->dim = 2;
            g->cNeighborIdx[i]->data_pointer.p2 = malloc((num_blocks * g->blkSize) * sizeof(int) + (num_blocks) * sizeof(char *));
            char *pos = (char *) g->cNeighborIdx[i]->data_pointer.p2 + num_blocks * sizeof(char *);
            for (int b = 0; b < num_blocks; b++) {
                g->cNeighborIdx[i]->data_pointer.p2[b] = (int *) pos;
                pos += g->blkSize * sizeof(int);
                for (int c = 0; c < g->blkSize; c++) {
                    g->cNeighborIdx[i]->data_pointer.p2[b][c] = (int) 0;
                }
            }
        }
        {
            int num_blocks = local_cell_blocks ? local_cell_blocks : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
            g->cNeighborBlk[i] = malloc(24);
            g->cNeighborBlk[i]->name = "g->cNeighborBlk[ i]";
            g->cNeighborBlk[i]->loc = 0;
            g->cNeighborBlk[i]->dim = 2;
            g->cNeighborBlk[i]->data_pointer.p2 = malloc((num_blocks * g->blkSize) * sizeof(int) + (num_blocks) * sizeof(char *));
            char *pos = (char *) g->cNeighborBlk[i]->data_pointer.p2 + num_blocks * sizeof(char *);
            for (int b = 0; b < num_blocks; b++) {
                g->cNeighborBlk[i]->data_pointer.p2[b] = (int *) pos;
                pos += g->blkSize * sizeof(int);
                for (int c = 0; c < g->blkSize; c++) {
                    g->cNeighborBlk[i]->data_pointer.p2[b][c] = (int) 0;
                }
            }
        }
    }
    int first_cBlock = g->mpi_rank * ((g->cBlkCnt + g->mpi_world_size - 1) / g->mpi_world_size);
    int first_eBlock = g->mpi_rank * ((g->eBlkCnt + g->mpi_world_size - 1) / g->mpi_world_size);
    for (int i = 0; i < NBRS; i++) {
        {
            size_t min_block = g->mpi_rank == (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
            size_t max_block = g->mpi_rank < (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
            for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
                for (size_t cell_index = (0); cell_index < (g->blkSize); cell_index++) {
                    if ((g->mpi_rank == g->mpi_world_size - 1) && (block_index == (g->cBlkCnt - 1) % ((g->cBlkCnt + g->mpi_world_size - 1) / g->mpi_world_size)) && (cell_index > (g->cellCount - 1) % g->blkSize)) {
                        g->cNeighborIdx[i]->data_pointer.p2[(block_index)][(cell_index)] = cell_index;
                        g->cNeighborBlk[i]->data_pointer.p2[(block_index)][(cell_index)] = block_index;
                    } else {
                        g->cNeighborIdx[i]->data_pointer.p2[(block_index)][(cell_index)] = g->neighbor[i][(first_cBlock + block_index) * g->blkSize + cell_index] % g->blkSize;
                        g->cNeighborBlk[i]->data_pointer.p2[(block_index)][(cell_index)] = g->neighbor[i][(first_cBlock + block_index) * g->blkSize + cell_index] / g->blkSize;
                    }
                }
            }
        }
    }
    for (int i = 0; i < NBRS; i++) {
        {
            int num_blocks = local_cell_blocks ? local_cell_blocks : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
            g->cEdgeIdx[i] = malloc(24);
            g->cEdgeIdx[i]->name = "g->cEdgeIdx[ i]";
            g->cEdgeIdx[i]->loc = 0;
            g->cEdgeIdx[i]->dim = 2;
            g->cEdgeIdx[i]->data_pointer.p2 = malloc((num_blocks * g->blkSize) * sizeof(int) + (num_blocks) * sizeof(char *));
            char *pos = (char *) g->cEdgeIdx[i]->data_pointer.p2 + num_blocks * sizeof(char *);
            for (int b = 0; b < num_blocks; b++) {
                g->cEdgeIdx[i]->data_pointer.p2[b] = (int *) pos;
                pos += g->blkSize * sizeof(int);
                for (int c = 0; c < g->blkSize; c++) {
                    g->cEdgeIdx[i]->data_pointer.p2[b][c] = (int) 0;
                }
            }
        }
        {
            int num_blocks = local_cell_blocks ? local_cell_blocks : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
            g->cEdgeBlk[i] = malloc(24);
            g->cEdgeBlk[i]->name = "g->cEdgeBlk[ i]";
            g->cEdgeBlk[i]->loc = 0;
            g->cEdgeBlk[i]->dim = 2;
            g->cEdgeBlk[i]->data_pointer.p2 = malloc((num_blocks * g->blkSize) * sizeof(int) + (num_blocks) * sizeof(char *));
            char *pos = (char *) g->cEdgeBlk[i]->data_pointer.p2 + num_blocks * sizeof(char *);
            for (int b = 0; b < num_blocks; b++) {
                g->cEdgeBlk[i]->data_pointer.p2[b] = (int *) pos;
                pos += g->blkSize * sizeof(int);
                for (int c = 0; c < g->blkSize; c++) {
                    g->cEdgeBlk[i]->data_pointer.p2[b][c] = (int) 0;
                }
            }
        }
    }
    for (int i = 0; i < NBRS; i++) {
        {
            size_t min_block = g->mpi_rank == (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
            size_t max_block = g->mpi_rank < (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
            for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
                for (size_t cell_index = (0); cell_index < (g->blkSize); cell_index++) {
                    if ((g->mpi_rank == g->mpi_world_size - 1) && (block_index == (g->cBlkCnt - 1) % ((g->cBlkCnt + g->mpi_world_size - 1) / g->mpi_world_size)) && (cell_index > (g->cellCount - 1) % g->blkSize)) {
                        g->cEdgeIdx[i]->data_pointer.p2[(block_index)][(cell_index)] = 0;
                        g->cEdgeBlk[i]->data_pointer.p2[(block_index)][(cell_index)] = first_eBlock;
                    } else {
                        g->cEdgeIdx[i]->data_pointer.p2[(block_index)][(cell_index)] = g->cedge[i][(first_cBlock + block_index) * g->blkSize + cell_index] % g->blkSize;
                        g->cEdgeBlk[i]->data_pointer.p2[(block_index)][(cell_index)] = g->cedge[i][(first_cBlock + block_index) * g->blkSize + cell_index] / g->blkSize;
                    }
                }
            }
        }
    }
    for (int i = 0; i < 2; i++) {
        {
            int num_blocks = local_edge_blocks ? local_edge_blocks : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
            g->eCellIdx[i] = malloc(24);
            g->eCellIdx[i]->name = "g->eCellIdx[ i]";
            g->eCellIdx[i]->loc = 1;
            g->eCellIdx[i]->dim = 2;
            g->eCellIdx[i]->data_pointer.p2 = malloc((num_blocks * g->blkSize) * sizeof(int) + (num_blocks) * sizeof(char *));
            char *pos = (char *) g->eCellIdx[i]->data_pointer.p2 + num_blocks * sizeof(char *);
            for (int b = 0; b < num_blocks; b++) {
                g->eCellIdx[i]->data_pointer.p2[b] = (int *) pos;
                pos += g->blkSize * sizeof(int);
                for (int e = 0; e < g->blkSize; e++) {
                    g->eCellIdx[i]->data_pointer.p2[b][e] = (int) 0;
                }
            }
        }
        {
            int num_blocks = local_edge_blocks ? local_edge_blocks : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
            g->eCellBlk[i] = malloc(24);
            g->eCellBlk[i]->name = "g->eCellBlk[ i]";
            g->eCellBlk[i]->loc = 1;
            g->eCellBlk[i]->dim = 2;
            g->eCellBlk[i]->data_pointer.p2 = malloc((num_blocks * g->blkSize) * sizeof(int) + (num_blocks) * sizeof(char *));
            char *pos = (char *) g->eCellBlk[i]->data_pointer.p2 + num_blocks * sizeof(char *);
            for (int b = 0; b < num_blocks; b++) {
                g->eCellBlk[i]->data_pointer.p2[b] = (int *) pos;
                pos += g->blkSize * sizeof(int);
                for (int e = 0; e < g->blkSize; e++) {
                    g->eCellBlk[i]->data_pointer.p2[b][e] = (int) 0;
                }
            }
        }
    }
    for (int i = 0; i < 2; i++) {
        {
            size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
            size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
            for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
                for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                    if ((g->mpi_rank == g->mpi_world_size - 1) && (block_index == (g->eBlkCnt - 1) % ((g->eBlkCnt + g->mpi_world_size - 1) / g->mpi_world_size)) && (edge_index > (g->edgeCount - 1) % g->blkSize)) {
                        g->eCellIdx[i]->data_pointer.p2[(block_index)][(edge_index)] = 0;
                        g->eCellBlk[i]->data_pointer.p2[(block_index)][(edge_index)] = first_cBlock;
                    } else {
                        g->eCellIdx[i]->data_pointer.p2[(block_index)][(edge_index)] = g->ecell[i][(first_eBlock + block_index) * g->blkSize + edge_index] % g->blkSize;
                        g->eCellBlk[i]->data_pointer.p2[(block_index)][(edge_index)] = g->ecell[i][(first_eBlock + block_index) * g->blkSize + edge_index] / g->blkSize;
                    }
                }
            }
        }
    }
}

void init_grid(GRID * g, int cellCount, int height)
{
    g->cellCount = cellCount * cellCount;
//cellCount;
    g->height = cellCount;
    g->edgeCount = calc_edge_count(g);
    g->blkSize = BLKSIZE;
    g->cBlkCnt = (g->cellCount + g->blkSize - 1) / g->blkSize;
    g->eBlkCnt = (g->edgeCount + g->blkSize - 1) / g->blkSize;
    create_maps(g);
    tessellation(g);
    g->height = height;
    init_edge_weights(g);
    init_blocking(g);
    {
        int cell_min = 0;
        int cell_max = g->mpi_rank > (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        int edge_min = 0;
        int edge_max = g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        int *cn_H = malloc(g->mpi_world_size * sizeof(int) * 2);
        cn_c = malloc(g->mpi_world_size * sizeof(int));
        for (int i = 0; i < g->mpi_world_size; i++)
            cn_c[i] = 0;
        for (int b = cell_min; b < cell_max; b++) {
            for (int c = (0); c < (g->blkSize); c++) {
                for (int n = (0); n < 3; n++) {
                    if (g->cNeighborBlk[n]->data_pointer.p2[b][c] >= g->cBlkCnt)
                        continue;
                    if (g->cNeighborBlk[n]->data_pointer.p2[b][c] < g->mpi_rank * ((((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size)) + cell_min || g->cNeighborBlk[n]->data_pointer.p2[b][c] >= g->mpi_rank * ((((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size)) + cell_max) {
                        cn_c[g->cNeighborBlk[n]->data_pointer.p2[b][c] / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size)]++;
                    }
                }
            }
        }
        cn_H[0] = cn_c[0] + cell_max * g->blkSize;
        for (int i = 1; i < g->mpi_world_size; i++) {
            cn_H[2 * i] = cn_c[i] + cn_H[2 * i - 2];
        }
        int ml = 0;
        for (int i = 0; i < g->mpi_world_size; i++) {
            ml += cn_c[i];
        }
        neighbor_map = malloc(ml * sizeof(int) * 5);
        for (int i = 0; i < g->mpi_world_size; i++) {
            cn_H[2 * i + 1] = cn_H[2 * i] % g->blkSize;
            cn_H[2 * i] = cn_H[2 * i] / g->blkSize;
        }
        int *tp = malloc(g->mpi_world_size * sizeof(int) * 2);
        tp[0] = cell_max;
        tp[1] = 0;
        for (int i = 1; i < g->mpi_world_size; i++) {
            tp[i * 2] = cn_H[i * 2 - 2];
            tp[i * 2 + 1] = cn_H[i * 2 - 1];
        }
        int *mi = malloc(g->mpi_world_size * sizeof(int));
        mi[0] = 0;
        for (int i = 1; i < g->mpi_world_size; i++)
            mi[i] = 5 * cn_c[i - 1] + mi[i - 1];
        for (int b = cell_min; b < cell_max; b++) {
            for (int c = (0); c < (g->blkSize); c++) {
                for (int n = (0); n < 3; n++) {
                    if (g->cNeighborBlk[n]->data_pointer.p2[b][c] >= g->cBlkCnt || g->cNeighborBlk[n]->data_pointer.p2[b][c] < 0) {
                        g->cNeighborBlk[n]->data_pointer.p2[b][c] = -1;
                    } else if (g->cNeighborBlk[n]->data_pointer.p2[b][c] < g->mpi_rank * ((((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size)) + cell_min || g->cNeighborBlk[n]->data_pointer.p2[b][c] >= g->mpi_rank * ((((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size)) + cell_max) {
                        int pn = g->cNeighborBlk[n]->data_pointer.p2[b][c] / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
                        neighbor_map[mi[pn]++] = pn;
                        neighbor_map[mi[pn]++] = g->cNeighborBlk[n]->data_pointer.p2[b][c] % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
                        neighbor_map[mi[pn]++] = g->cNeighborIdx[n]->data_pointer.p2[b][c];
                        neighbor_map[mi[pn]++] = tp[pn * 2];
                        neighbor_map[mi[pn]++] = tp[pn * 2 + 1];
                        g->cNeighborBlk[n]->data_pointer.p2[b][c] = tp[pn * 2];
                        g->cNeighborIdx[n]->data_pointer.p2[b][c] = tp[pn * 2 + 1];
                        if (++tp[pn * 2 + 1] == g->blkSize) {
                            tp[pn * 2]++;
                            tp[pn * 2 + 1] = 0;
                        }
                    } else {
                        g->cNeighborBlk[n]->data_pointer.p2[b][c] = g->cNeighborBlk[n]->data_pointer.p2[b][c] % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
                    }
                }
            }
        }
        int *ce_H = malloc(g->mpi_world_size * sizeof(int) * 2);
        ce_c = malloc(g->mpi_world_size * sizeof(int));
        for (int i = 0; i < g->mpi_world_size; i++)
            ce_c[i] = 0;
        for (int b = cell_min; b < cell_max; b++) {
            for (int c = (0); c < (g->blkSize); c++) {
                for (int n = (0); n < 3; n++) {
                    if (g->cEdgeBlk[n]->data_pointer.p2[b][c] < g->mpi_rank * ((((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size)) + edge_min || g->cEdgeBlk[n]->data_pointer.p2[b][c] >= g->mpi_rank * ((((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size)) + edge_max) {
                        ce_c[g->cEdgeBlk[n]->data_pointer.p2[b][c] / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size)]++;
                    }
                }
            }
        }
        ce_H[0] = ce_c[0] + edge_max * g->blkSize;
        for (int i = 1; i < g->mpi_world_size; i++) {
            ce_H[2 * i] = ce_c[i] + ce_H[2 * i - 2];
        }
        ml = 0;
        for (int i = 0; i < g->mpi_world_size; i++) {
            ml += ce_c[i];
        }
        cedge_map = malloc(ml * sizeof(int) * 5);
        for (int i = 0; i < g->mpi_world_size; i++) {
            ce_H[2 * i + 1] = ce_H[2 * i] % g->blkSize;
            ce_H[2 * i] = ce_H[2 * i] / g->blkSize;
        }
        local_edge_blocks = ce_H[g->mpi_world_size * 2 - 2] + 1;
        tp[0] = edge_max;
        tp[1] = 0;
        for (int i = 1; i < g->mpi_world_size; i++) {
            tp[i * 2] = ce_H[i * 2 - 2];
            tp[i * 2 + 1] = ce_H[i * 2 - 1];
        }
        mi[0] = 0;
        for (int i = 1; i < g->mpi_world_size; i++)
            mi[i] = 5 * ce_c[i - 1] + mi[i - 1];
        for (int b = cell_min; b < cell_max; b++) {
            for (int c = (0); c < (g->blkSize); c++) {
                for (int n = (0); n < 3; n++) {
                    if (g->cEdgeBlk[n]->data_pointer.p2[b][c] < g->mpi_rank * ((((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size)) + edge_min || g->cEdgeBlk[n]->data_pointer.p2[b][c] >= g->mpi_rank * ((((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size)) + edge_max) {
                        int pn = g->cEdgeBlk[n]->data_pointer.p2[b][c] / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
                        cedge_map[mi[pn]++] = pn;
                        cedge_map[mi[pn]++] = g->cEdgeBlk[n]->data_pointer.p2[b][c] % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
                        cedge_map[mi[pn]++] = g->cEdgeIdx[n]->data_pointer.p2[b][c];
                        cedge_map[mi[pn]++] = tp[pn * 2];
                        cedge_map[mi[pn]++] = tp[pn * 2 + 1];
                        g->cEdgeBlk[n]->data_pointer.p2[b][c] = tp[pn * 2];
                        g->cEdgeIdx[n]->data_pointer.p2[b][c] = tp[pn * 2 + 1];
                        if (++tp[pn * 2 + 1] == g->blkSize) {
                            tp[pn * 2]++;
                            tp[pn * 2 + 1] = 0;
                        }
                    } else {
                        g->cEdgeBlk[n]->data_pointer.p2[b][c] = g->cEdgeBlk[n]->data_pointer.p2[b][c] % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
                    }
                }
            }
        }
        int *ec_H = malloc(g->mpi_world_size * sizeof(int) * 2);
        ec_c = malloc(g->mpi_world_size * sizeof(int));
        for (int i = 0; i < g->mpi_world_size; i++)
            ec_c[i] = 0;
        for (int b = edge_min; b < edge_max; b++) {
            for (int e = (0); e < (g->blkSize); e++) {
                for (int n = (0); n < 2; n++) {
                    if (g->eCellBlk[n]->data_pointer.p2[b][e] >= g->cBlkCnt)
                        continue;
                    if (g->eCellBlk[n]->data_pointer.p2[b][e] < g->mpi_rank * ((((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size)) + cell_min || g->eCellBlk[n]->data_pointer.p2[b][e] >= g->mpi_rank * ((((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size)) + cell_max) {
                        ec_c[g->eCellBlk[n]->data_pointer.p2[b][e] / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size)]++;
                    }
                }
            }
        }
        ec_H[0] = ec_c[0] + cn_H[g->mpi_world_size * 2 - 2] * g->blkSize + cn_H[g->mpi_world_size * 2 - 1];
        for (int i = 1; i < g->mpi_world_size; i++) {
            ec_H[2 * i] = ec_c[i] + ec_H[2 * i - 2];
        }
        ml = 0;
        for (int i = 0; i < g->mpi_world_size; i++) {
            ml += ec_c[i];
        }
        ecell_map = malloc(ml * sizeof(int) * 5);
        for (int i = 0; i < g->mpi_world_size; i++) {
            ec_H[2 * i + 1] = ec_H[2 * i] % g->blkSize;
            ec_H[2 * i] = ec_H[2 * i] / g->blkSize;
        }
        local_cell_blocks = ec_H[g->mpi_world_size * 2 - 2] + 1;
        tp[0] = cn_H[g->mpi_world_size * 2 - 2];
        tp[1] = cn_H[g->mpi_world_size * 2 - 1];
        for (int i = 1; i < g->mpi_world_size; i++) {
            tp[i * 2] = ec_H[i * 2 - 2];
            tp[i * 2 + 1] = ec_H[i * 2 - 1];
        }
        mi[0] = 0;
        for (int i = 1; i < g->mpi_world_size; i++)
            mi[i] = 5 * ec_c[i - 1] + mi[i - 1];
        for (int b = edge_min; b < edge_max; b++) {
            for (int e = (0); e < (g->blkSize); e++) {
                for (int n = (0); n < 2; n++) {
                    if (g->eCellBlk[n]->data_pointer.p2[b][e] >= g->cBlkCnt || g->eCellBlk[n]->data_pointer.p2[b][e] < 0) {
                        g->eCellBlk[n]->data_pointer.p2[b][e] = -1;
                    } else if (g->eCellBlk[n]->data_pointer.p2[b][e] < g->mpi_rank * ((((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size)) + cell_min || g->eCellBlk[n]->data_pointer.p2[b][e] >= g->mpi_rank * ((((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size)) + cell_max) {
                        int pn = g->eCellBlk[n]->data_pointer.p2[b][e] / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
                        ecell_map[mi[pn]++] = pn;
                        ecell_map[mi[pn]++] = g->eCellBlk[n]->data_pointer.p2[b][e] % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
                        ecell_map[mi[pn]++] = g->eCellIdx[n]->data_pointer.p2[b][e];
                        ecell_map[mi[pn]++] = tp[pn * 2];
                        ecell_map[mi[pn]++] = tp[pn * 2 + 1];
                        g->eCellBlk[n]->data_pointer.p2[b][e] = tp[pn * 2];
                        g->eCellIdx[n]->data_pointer.p2[b][e] = tp[pn * 2 + 1];
                        if (++tp[pn * 2 + 1] == g->blkSize) {
                            tp[pn * 2]++;
                            tp[pn * 2 + 1] = 0;
                        }
                    } else {
                        g->eCellBlk[n]->data_pointer.p2[b][e] = g->eCellBlk[n]->data_pointer.p2[b][e] % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
                    }
                }
            }
        }
        free(tp);
        mpi_send_requests = malloc(g->mpi_world_size * 2 * sizeof(MPI_Request));
        mpi_recv_requests = &mpi_send_requests[g->mpi_world_size];
        mpi_send_requests[g->mpi_rank] = MPI_REQUEST_NULL;
        mpi_recv_requests[g->mpi_rank] = MPI_REQUEST_NULL;
        comm_tag = 10;
        cn_crem = malloc(g->mpi_world_size * sizeof(int));
        ce_crem = malloc(g->mpi_world_size * sizeof(int));
        ec_crem = malloc(g->mpi_world_size * sizeof(int));
        for (int pn = 0; pn < g->mpi_world_size; pn++) {
            if (pn != g->mpi_rank) {
                MPI_Isend(&cn_c[pn], 1, MPI_INT, pn, 0, MPI_COMM_WORLD, &mpi_send_requests[pn]);
                MPI_Irecv(&cn_crem[pn], 1, MPI_INT, pn, 0, MPI_COMM_WORLD, &mpi_recv_requests[pn]);
                MPI_Isend(&ce_c[pn], 1, MPI_INT, pn, 1, MPI_COMM_WORLD, &mpi_send_requests[pn]);
                MPI_Irecv(&ce_crem[pn], 1, MPI_INT, pn, 1, MPI_COMM_WORLD, &mpi_recv_requests[pn]);
                MPI_Isend(&ec_c[pn], 1, MPI_INT, pn, 2, MPI_COMM_WORLD, &mpi_send_requests[pn]);
                MPI_Irecv(&ec_crem[pn], 1, MPI_INT, pn, 2, MPI_COMM_WORLD, &mpi_recv_requests[pn]);
            } else
                cn_c[g->mpi_rank] = ce_c[g->mpi_rank] = ec_c[g->mpi_rank] = cn_crem[g->mpi_rank] = ce_crem[g->mpi_rank] = ec_crem[g->mpi_rank] = 0;
        }
        MPI_Waitall(g->mpi_world_size * 2, mpi_send_requests, MPI_STATUSES_IGNORE);
        for (int i = 1; i < g->mpi_world_size; i++) {
            cn_c[i] += cn_c[i - 1];
            cn_crem[i] += cn_crem[i - 1];
        }
        neighbor_maprem = malloc(cn_crem[g->mpi_world_size - 1] * sizeof(int) * 2);
        for (int i = 1; i < g->mpi_world_size; i++) {
            ce_c[i] += ce_c[i - 1];
            ce_crem[i] += ce_crem[i - 1];
        }
        cedge_maprem = malloc(ce_crem[g->mpi_world_size - 1] * sizeof(int) * 2);
        for (int i = 1; i < g->mpi_world_size; i++) {
            ec_c[i] += ec_c[i - 1];
            ec_crem[i] += ec_crem[i - 1];
        }
        ecell_maprem = malloc(ec_crem[g->mpi_world_size - 1] * sizeof(int) * 2);
        for (int pn = 0; pn < g->mpi_world_size; pn++) {
            if (pn != g->mpi_rank) {
                int *buf = malloc((cn_c[pn] - (pn ? cn_c[pn - 1] : 0)) * sizeof(int) * 2);
                for (int i = 0; i < (cn_c[pn] - (pn ? cn_c[pn - 1] : 0)); i++) {
                    buf[2 * i] = neighbor_map[(pn ? cn_c[pn - 1] * 5 : 0) + 5 * i + 1];
                    buf[2 * i + 1] = neighbor_map[(pn ? cn_c[pn - 1] * 5 : 0) + 5 * i + 2];
                }
                MPI_Isend(buf, (cn_c[pn] - (pn ? cn_c[pn - 1] : 0)) * 2, MPI_INT, pn, 3, MPI_COMM_WORLD, &mpi_send_requests[pn]);
                MPI_Irecv(&neighbor_maprem[(pn ? cn_crem[pn - 1] * 2 : 0)], (cn_crem[pn] - (pn ? cn_crem[pn - 1] : 0)) * 2, MPI_INT, pn, 3, MPI_COMM_WORLD, &mpi_recv_requests[pn]);
                MPI_Wait(&mpi_recv_requests[pn], MPI_STATUS_IGNORE);
                MPI_Wait(&mpi_send_requests[pn], MPI_STATUS_IGNORE);
                free(buf);
                buf = malloc((ce_c[pn] - (pn ? ce_c[pn - 1] : 0)) * sizeof(int) * 2);
                for (int i = 0; i < (ce_c[pn] - (pn ? ce_c[pn - 1] : 0)); i++) {
                    buf[2 * i] = cedge_map[(pn ? ce_c[pn - 1] * 5 : 0) + 5 * i + 1];
                    buf[2 * i + 1] = cedge_map[(pn ? ce_c[pn - 1] * 5 : 0) + 5 * i + 2];
                }
                MPI_Isend(buf, (ce_c[pn] - (pn ? ce_c[pn - 1] : 0)) * 2, MPI_INT, pn, 4, MPI_COMM_WORLD, &mpi_send_requests[pn]);
                MPI_Irecv(&cedge_maprem[(pn ? ce_crem[pn - 1] * 2 : 0)], (ce_crem[pn] - (pn ? ce_crem[pn - 1] : 0)) * 2, MPI_INT, pn, 4, MPI_COMM_WORLD, &mpi_recv_requests[pn]);
                MPI_Wait(&mpi_recv_requests[pn], MPI_STATUS_IGNORE);
                MPI_Wait(&mpi_send_requests[pn], MPI_STATUS_IGNORE);
                free(buf);
                buf = malloc((ec_c[pn] - (pn ? ec_c[pn - 1] : 0)) * sizeof(int) * 2);
                for (int i = 0; i < (ec_c[pn] - (pn ? ec_c[pn - 1] : 0)); i++) {
                    buf[2 * i] = ecell_map[(pn ? ec_c[pn - 1] * 5 : 0) + 5 * i + 1];
                    buf[2 * i + 1] = ecell_map[(pn ? ec_c[pn - 1] * 5 : 0) + 5 * i + 2];
                }
                MPI_Isend(buf, (ec_c[pn] - (pn ? ec_c[pn - 1] : 0)) * 2, MPI_INT, pn, 5, MPI_COMM_WORLD, &mpi_send_requests[pn]);
                MPI_Irecv(&ecell_maprem[(pn ? ec_crem[pn - 1] * 2 : 0)], (ec_crem[pn] - (pn ? ec_crem[pn - 1] : 0)) * 2, MPI_INT, pn, 5, MPI_COMM_WORLD, &mpi_recv_requests[pn]);
                MPI_Wait(&mpi_recv_requests[pn], MPI_STATUS_IGNORE);
                MPI_Wait(&mpi_send_requests[pn], MPI_STATUS_IGNORE);
                free(buf);
            }
        }
        neighbor_2Dbufrem = malloc(g->mpi_world_size * sizeof(GVAL *));
        neighbor_3Dbufrem = malloc(g->mpi_world_size * sizeof(GVAL *));
        cedge_2Dbufrem = malloc(g->mpi_world_size * sizeof(GVAL *));
        cedge_3Dbufrem = malloc(g->mpi_world_size * sizeof(GVAL *));
        ecell_2Dbufrem = malloc(g->mpi_world_size * sizeof(GVAL *));
        ecell_3Dbufrem = malloc(g->mpi_world_size * sizeof(GVAL *));
        for (int pn = 0; pn < g->mpi_world_size; pn++) {
            if (pn != g->mpi_rank) {
                neighbor_2Dbufrem[pn] = malloc((cn_crem[pn] - (pn ? cn_crem[pn - 1] : 0)) * sizeof(GVAL));
                neighbor_3Dbufrem[pn] = malloc((cn_crem[pn] - (pn ? cn_crem[pn - 1] : 0)) * g->height * sizeof(GVAL));
                cedge_2Dbufrem[pn] = malloc((ce_crem[pn] - (pn ? ce_crem[pn - 1] : 0)) * sizeof(GVAL));
                cedge_3Dbufrem[pn] = malloc((ce_crem[pn] - (pn ? ce_crem[pn - 1] : 0)) * g->height * sizeof(GVAL));
                ecell_2Dbufrem[pn] = malloc((ec_crem[pn] - (pn ? ec_crem[pn - 1] : 0)) * sizeof(GVAL));
                ecell_3Dbufrem[pn] = malloc((ec_crem[pn] - (pn ? ec_crem[pn - 1] : 0)) * g->height * sizeof(GVAL));
            }
        }
        neighbor_2Dbuf = malloc(g->mpi_world_size * sizeof(GVAL *));
        neighbor_3Dbuf = malloc(g->mpi_world_size * sizeof(GVAL *));
        cedge_2Dbuf = malloc(g->mpi_world_size * sizeof(GVAL *));
        cedge_3Dbuf = malloc(g->mpi_world_size * sizeof(GVAL *));
        ecell_2Dbuf = malloc(g->mpi_world_size * sizeof(GVAL *));
        ecell_3Dbuf = malloc(g->mpi_world_size * sizeof(GVAL *));
        for (int pn = 0; pn < g->mpi_world_size; pn++) {
            if (pn != g->mpi_rank) {
                neighbor_2Dbuf[pn] = malloc((cn_c[pn] - (pn ? cn_c[pn - 1] : 0)) * sizeof(GVAL));
                neighbor_3Dbuf[pn] = malloc((cn_c[pn] - (pn ? cn_c[pn - 1] : 0)) * g->height * sizeof(GVAL));
                cedge_2Dbuf[pn] = malloc((ce_c[pn] - (pn ? ce_c[pn - 1] : 0)) * sizeof(GVAL));
                cedge_3Dbuf[pn] = malloc((ce_c[pn] - (pn ? ce_c[pn - 1] : 0)) * g->height * sizeof(GVAL));
                ecell_2Dbuf[pn] = malloc((ec_c[pn] - (pn ? ec_c[pn - 1] : 0)) * sizeof(GVAL));
                ecell_3Dbuf[pn] = malloc((ec_c[pn] - (pn ? ec_c[pn - 1] : 0)) * g->height * sizeof(GVAL));
            }
        }
    }
}

void get_indices_c(GRID * g, int blk, int *cb, int *ce)
{
    *cb = 0;
    *ce = blk == g->cBlkCnt - 1 ? g->cellCount % g->blkSize == 0 ? g->blkSize : g->cellCount % g->blkSize : g->blkSize;
}

void get_indices_e(GRID * g, int blk, int *eb, int *ee)
{
    *eb = 0;
    *ee = blk == g->eBlkCnt - 1 ? g->edgeCount % g->blkSize == 0 ? g->blkSize : g->edgeCount % g->blkSize : g->blkSize;
}
