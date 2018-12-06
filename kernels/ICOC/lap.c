#include <mpi.h>

extern int *cn_c;
extern int *ce_c;
extern int *ec_c;
extern int *cn_crem;
extern int *ce_crem;
extern int *ec_crem;
extern int *neighbor_map;
extern int *cedge_map;
extern int *ecell_map;
extern int *neighbor_maprem;
extern int *cedge_maprem;
extern int *ecell_maprem;
extern GVAL **neighbor_2Dbuf;
extern GVAL **neighbor_3Dbuf;
extern GVAL **cedge_2Dbuf;
extern GVAL **cedge_3Dbuf;
extern GVAL **ecell_2Dbuf;
extern GVAL **ecell_3Dbuf;
extern GVAL **neighbor_2Dbufrem;
extern GVAL **neighbor_3Dbufrem;
extern GVAL **cedge_2Dbufrem;
extern GVAL **cedge_3Dbufrem;
extern GVAL **ecell_2Dbufrem;
extern GVAL **ecell_3Dbufrem;
extern MPI_Request *mpi_send_requests;
extern MPI_Request *mpi_recv_requests;
extern int comm_tag;
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
} *gv_temp_alt;
void lap(GRID * g)
{
    {
        {
            comm_tag++;
            for (int pn = 0; pn < g->mpi_world_size; pn++) {
                if (pn != g->mpi_rank) {
                    for (int i = 0; i < (cn_crem[pn] - (pn ? cn_crem[pn - 1] : 0)); i++) {
                        for (int k = 0; k < g->height; k++)
                            neighbor_3Dbufrem[pn][g->height * i + k] = gv_temp->data_pointer.p3[neighbor_maprem[(pn ? cn_crem[pn - 1] * 2 : 0) + 2 * i]][k][neighbor_maprem[(pn ? cn_crem[pn - 1] * 2 : 0) + 2 * i + 1]];
                    }
                    MPI_Isend(neighbor_3Dbufrem[pn], (cn_crem[pn] - (pn ? cn_crem[pn - 1] : 0)) * g->height, MPI_FLOAT, pn, comm_tag, MPI_COMM_WORLD, &mpi_send_requests[pn]);
                    MPI_Irecv(neighbor_3Dbuf[pn], (cn_c[pn] - (pn ? cn_c[pn - 1] : 0)) * g->height, MPI_FLOAT, pn, comm_tag, MPI_COMM_WORLD, &mpi_recv_requests[pn]);
                }
            }
            MPI_Waitall(g->mpi_world_size * 2, mpi_send_requests, MPI_STATUSES_IGNORE);
            for (int pn = 0; pn < g->mpi_world_size; pn++) {
                if (pn != g->mpi_rank) {
                    for (int i = 0; i < (cn_c[pn] - (pn ? cn_c[pn - 1] : 0)); i++) {
                        for (int k = 0; k < g->height; k++)
                            gv_temp->data_pointer.p3[neighbor_map[(pn ? cn_c[pn - 1] * 5 : 0) + 5 * i + 3]][k][neighbor_map[(pn ? cn_c[pn - 1] * 5 : 0) + 5 * i + 4]] = neighbor_3Dbuf[pn][g->height * i + k];
                    }
                }
            }
        }
        size_t min_block = g->mpi_rank == (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t height_index = (0); height_index < (g->height); height_index++) {
                for (size_t cell_index = (0); cell_index < (g->blkSize); cell_index++) {
                    gv_temp_alt->data_pointer.p3[(block_index)][(height_index)][(cell_index)] = (gv_temp->data_pointer.p3[(g->cNeighborBlk[(0)]->data_pointer.p2[(block_index)][(cell_index)])][(height_index)][(g->cNeighborIdx[(0)]->data_pointer.p2[(block_index)][(cell_index)])] + gv_temp->data_pointer.p3[(g->cNeighborBlk[(1)]->data_pointer.p2[(block_index)][(cell_index)])][(height_index)][(g->cNeighborIdx[(1)]->data_pointer.p2[(block_index)][(cell_index)])] + gv_temp->data_pointer.p3[(g->cNeighborBlk[(2)]->data_pointer.p2[(block_index)][(cell_index)])][(height_index)][(g->cNeighborIdx[(2)]->data_pointer.p2[(block_index)][(cell_index)])]) / 3.0;
                }
            }
        }
    }
    struct {
        char *name;
        int loc;
        int dim;
        union {
            GVAL *restrict * restrict p2;
            GVAL *restrict * restrict * restrict p3;
        } data_pointer;
    } *x = gv_temp_alt;
    gv_temp_alt = gv_temp;
    gv_temp = x;
}
