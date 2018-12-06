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
} *gv_grad;
extern struct {
    char *name;
    int loc;
    int dim;
    union {
        GVAL *restrict * restrict p2;
        GVAL *restrict * restrict * restrict p3;
    } data_pointer;
} *gv_dvg;
void dvg(GRID * g)
{
    {
        {
            comm_tag++;
            for (int pn = 0; pn < g->mpi_world_size; pn++) {
                if (pn != g->mpi_rank) {
                    for (int i = 0; i < (ce_crem[pn] - (pn ? ce_crem[pn - 1] : 0)); i++) {
                        for (int k = 0; k < g->height; k++)
                            cedge_3Dbufrem[pn][g->height * i + k] = gv_grad->data_pointer.p3[cedge_maprem[(pn ? ce_crem[pn - 1] * 2 : 0) + 2 * i]][k][cedge_maprem[(pn ? ce_crem[pn - 1] * 2 : 0) + 2 * i + 1]];
                    }
                    MPI_Isend(cedge_3Dbufrem[pn], (ce_crem[pn] - (pn ? ce_crem[pn - 1] : 0)) * g->height, MPI_FLOAT, pn, comm_tag, MPI_COMM_WORLD, &mpi_send_requests[pn]);
                    MPI_Irecv(cedge_3Dbuf[pn], (ce_c[pn] - (pn ? ce_c[pn - 1] : 0)) * g->height, MPI_FLOAT, pn, comm_tag, MPI_COMM_WORLD, &mpi_recv_requests[pn]);
                }
            }
            MPI_Waitall(g->mpi_world_size * 2, mpi_send_requests, MPI_STATUSES_IGNORE);
            for (int pn = 0; pn < g->mpi_world_size; pn++) {
                if (pn != g->mpi_rank) {
                    for (int i = 0; i < (ce_c[pn] - (pn ? ce_c[pn - 1] : 0)); i++) {
                        for (int k = 0; k < g->height; k++)
                            gv_grad->data_pointer.p3[cedge_map[(pn ? ce_c[pn - 1] * 5 : 0) + 5 * i + 3]][k][cedge_map[(pn ? ce_c[pn - 1] * 5 : 0) + 5 * i + 4]] = cedge_3Dbuf[pn][g->height * i + k];
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
                    gv_dvg->data_pointer.p3[(block_index)][(height_index)][(cell_index)] = 0;
                    for (int n = 0; n < NBRS; n++) {
                        gv_dvg->data_pointer.p3[(block_index)][(height_index)][(cell_index)] += g->edge_weights[n][cell_index % BLKSIZE] * gv_grad->data_pointer.p3[(g->cEdgeBlk[(n)]->data_pointer.p2[(block_index)][(cell_index)])][(height_index)][(g->cEdgeIdx[(n)]->data_pointer.p2[(block_index)][(cell_index)])];
                    }
                }
            }
        }
    }
}
