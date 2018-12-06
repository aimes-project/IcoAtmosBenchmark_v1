#include "grid.h"

extern int local_cell_blocks;
extern int local_edge_blocks;
#include "memory.h"
#include "component.h"
#include "io.h"
#include <stdint.h>
void com5_init(GRID * g);
void com5_compute(GRID * g);
void com5_io(GRID * g);
double com5_flops(GRID * g);
double com5_memory(GRID * g);
uint64_t com5_checksum(GRID *);
void com5_cleanup(GRID * g);
void lap(GRID * g);
struct {
    char *name;
    int loc;
    int dim;
    union {
        GVAL *restrict * restrict p2;
        GVAL *restrict * restrict * restrict p3;
    } data_pointer;
} *gv_temp_alt;
MODEL_COMPONENT com5 = { 0, com5_init, com5_compute, com5_io, com5_flops, com5_memory, com5_checksum, com5_cleanup };

extern MODEL_COMPONENT com1;
void com5_init(GRID * g)
{
    com5.loaded = 1;
    if (!com1.loaded)
        com1.init(g);
    {
        int num_blocks = local_cell_blocks ? local_cell_blocks : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        gv_temp_alt = malloc(24);
        gv_temp_alt->name = "gv_temp_alt";
        gv_temp_alt->loc = 0;
        gv_temp_alt->dim = 3;
        gv_temp_alt->data_pointer.p3 = malloc((num_blocks * g->height * g->blkSize) * sizeof(GVAL) + (num_blocks * g->height) * sizeof(char *) + (num_blocks) * sizeof(char *));
        char *pos = (char *) gv_temp_alt->data_pointer.p3 + num_blocks * sizeof(char *);
        char *pos2 = (char *) gv_temp_alt->data_pointer.p3 + num_blocks * sizeof(char *) + num_blocks * g->height * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            gv_temp_alt->data_pointer.p3[b] = (GVAL * *)pos;
            pos += g->height * sizeof(char *);
            for (int k = 0; k < g->height; k++) {
                gv_temp_alt->data_pointer.p3[b][k] = (GVAL *) pos2;
                pos2 += g->blkSize * sizeof(GVAL);
                for (int c = 0; c < g->blkSize; c++) {
                    gv_temp_alt->data_pointer.p3[b][k][c] = (GVAL) 0;
                }
            }
        }
    }
}

void com5_io(GRID * g)
{
}

void com5_compute(GRID * g)
{
    lap(g);
}

double com5_flops(GRID * g)
{
    double flop = 59.0 * (double) g->cellCount * (double) (g->height - 2);
    return flop;
}

double com5_memory(GRID * g)
{
    return 0.0;
}

uint64_t com5_checksum(GRID * g)
{
    return 0;
}

void com5_cleanup(GRID * g)
{
    free((void *) gv_temp_alt->data_pointer.p3);
}
