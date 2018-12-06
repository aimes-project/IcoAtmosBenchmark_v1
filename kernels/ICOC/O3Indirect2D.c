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
} *gv_ind2Dparam;
extern struct {
    char *name;
    int loc;
    int dim;
    union {
        GVAL *restrict * restrict p2;
        GVAL *restrict * restrict * restrict p3;
    } data_pointer;
} *gv_ind2Dvar;
void O3Indirect2D(GRID * g)
{
    {
        {
            comm_tag++;
            for (int pn = 0; pn < g->mpi_world_size; pn++) {
                if (pn != g->mpi_rank) {
                    for (int i = 0; i < (ec_crem[pn] - (pn ? ec_crem[pn - 1] : 0)); i++) {
                        ecell_2Dbufrem[pn][i] = gv_ind2Dparam->data_pointer.p2[ecell_maprem[(pn ? ec_crem[pn - 1] * 2 : 0) + 2 * i]][ecell_maprem[(pn ? ec_crem[pn - 1] * 2 : 0) + 2 * i + 1]];
                    }
                    MPI_Isend(ecell_2Dbufrem[pn], (ec_crem[pn] - (pn ? ec_crem[pn - 1] : 0)), MPI_FLOAT, pn, comm_tag, MPI_COMM_WORLD, &mpi_send_requests[pn]);
                    MPI_Irecv(ecell_2Dbuf[pn], (ec_c[pn] - (pn ? ec_c[pn - 1] : 0)), MPI_FLOAT, pn, comm_tag, MPI_COMM_WORLD, &mpi_recv_requests[pn]);
                }
            }
            MPI_Waitall(g->mpi_world_size * 2, mpi_send_requests, MPI_STATUSES_IGNORE);
            for (int pn = 0; pn < g->mpi_world_size; pn++) {
                if (pn != g->mpi_rank) {
                    for (int i = 0; i < (ec_c[pn] - (pn ? ec_c[pn - 1] : 0)); i++) {
                        gv_ind2Dparam->data_pointer.p2[ecell_map[(pn ? ec_c[pn - 1] * 5 : 0) + 5 * i + 3]][ecell_map[(pn ? ec_c[pn - 1] * 5 : 0) + 5 * i + 4]] = ecell_2Dbuf[pn][i];
                    }
                }
            }
        }
        {
            comm_tag++;
            for (int pn = 0; pn < g->mpi_world_size; pn++) {
                if (pn != g->mpi_rank) {
                    for (int i = 0; i < (ec_crem[pn] - (pn ? ec_crem[pn - 1] : 0)); i++) {
                        for (int k = 0; k < g->height; k++)
                            ecell_3Dbufrem[pn][g->height * i + k] = gv_temp->data_pointer.p3[ecell_maprem[(pn ? ec_crem[pn - 1] * 2 : 0) + 2 * i]][k][ecell_maprem[(pn ? ec_crem[pn - 1] * 2 : 0) + 2 * i + 1]];
                    }
                    MPI_Isend(ecell_3Dbufrem[pn], (ec_crem[pn] - (pn ? ec_crem[pn - 1] : 0)) * g->height, MPI_FLOAT, pn, comm_tag, MPI_COMM_WORLD, &mpi_send_requests[pn]);
                    MPI_Irecv(ecell_3Dbuf[pn], (ec_c[pn] - (pn ? ec_c[pn - 1] : 0)) * g->height, MPI_FLOAT, pn, comm_tag, MPI_COMM_WORLD, &mpi_recv_requests[pn]);
                }
            }
            MPI_Waitall(g->mpi_world_size * 2, mpi_send_requests, MPI_STATUSES_IGNORE);
            for (int pn = 0; pn < g->mpi_world_size; pn++) {
                if (pn != g->mpi_rank) {
                    for (int i = 0; i < (ec_c[pn] - (pn ? ec_c[pn - 1] : 0)); i++) {
                        for (int k = 0; k < g->height; k++)
                            gv_temp->data_pointer.p3[ecell_map[(pn ? ec_c[pn - 1] * 5 : 0) + 5 * i + 3]][k][ecell_map[(pn ? ec_c[pn - 1] * 5 : 0) + 5 * i + 4]] = ecell_3Dbuf[pn][g->height * i + k];
                    }
                }
            }
        }
        size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t height_index = (0); height_index < (g->height); height_index++) {
                for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                    gv_ind2Dvar->data_pointer.p3[(block_index)][(height_index)][(edge_index)] = gv_ind2Dparam->data_pointer.p2[(g->eCellBlk[(0)]->data_pointer.p2[(block_index)][(edge_index)])][(g->eCellIdx[(0)]->data_pointer.p2[(block_index)][(edge_index)])] * gv_temp->data_pointer.p3[(g->eCellBlk[(0)]->data_pointer.p2[(block_index)][(edge_index)])][(height_index)][(g->eCellIdx[(0)]->data_pointer.p2[(block_index)][(edge_index)])] - gv_ind2Dparam->data_pointer.p2[(g->eCellBlk[(1)]->data_pointer.p2[(block_index)][(edge_index)])][(g->eCellIdx[(1)]->data_pointer.p2[(block_index)][(edge_index)])] * gv_temp->data_pointer.p3[(g->eCellBlk[(1)]->data_pointer.p2[(block_index)][(edge_index)])][(height_index)][(g->eCellIdx[(1)]->data_pointer.p2[(block_index)][(edge_index)])];
                }
            }
        }
    }
}
