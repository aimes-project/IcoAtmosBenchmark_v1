#include "memory.h"

void *allocate_variable(GRID * g, GRID_VAR_POSITION gp, GRID_VAR_DIMENSION gd, size_t sz)
{
    void *ret = 0;
    int cb, ce, eb, ee;
    if (gp == GRID_POS_CELL) {
        if (gd == GRID_DIM_2D) {
            void **ptr = malloc((g->cBlkCnt) * sizeof(void *));
            for (int b = 0; b < g->cBlkCnt; b++) {
                get_indices_c(g, b, &cb, &ce);
                ptr[b] = malloc(ce * sz);
            }
            ret = (void *) ptr;
        } else if (gd == GRID_DIM_3D) {
            void ***ptr = malloc((g->cBlkCnt) * sizeof(void * *));
            for (int b = 0; b < g->cBlkCnt; b++) {
                get_indices_c(g, b, &cb, &ce);
                ptr[b] = malloc(g->height * sizeof(void *));
                for (int k = 0; k < g->height; k++) {
                    ptr[b][k] = malloc(ce * sz);
                }
            }
            ret = (void *) ptr;
        }
    } else if (gp == GRID_POS_EDGE) {
        if (gd == GRID_DIM_2D) {
            void **ptr = malloc((g->eBlkCnt) * sizeof(void *));
            for (int b = 0; b < g->eBlkCnt; b++) {
                get_indices_e(g, b, &eb, &ee);
                ptr[b] = malloc(ee * sz);
            }
            ret = (void *) ptr;
        } else if (gd == GRID_DIM_3D) {
            void ***ptr = malloc((g->eBlkCnt) * sizeof(void * *));
            for (int b = 0; b < g->eBlkCnt; b++) {
                get_indices_e(g, b, &eb, &ee);
                ptr[b] = malloc(g->height * sizeof(void *));
                for (int k = 0; k < g->height; k++) {
                    ptr[b][k] = malloc(ee * sz);
                }
            }
            ret = (void *) ptr;
        }
    }
    return ret;
}

void *deallocate_variable(void *var, GRID * g, GRID_VAR_POSITION gp, GRID_VAR_DIMENSION gd, size_t sz)
{
    int cb, ce, eb, ee;
    if (gp == GRID_POS_CELL) {
        if (gd == GRID_DIM_2D) {
            void **ptr = (void * *) var;
            for (int b = 0; b < g->cBlkCnt; b++) {
                free(ptr[b]);
            }
            free(ptr);
        } else if (gd == GRID_DIM_3D) {
            void ***ptr = (void * * *) var;
            for (int b = 0; b < g->cBlkCnt; b++) {
                for (int k = 0; k < g->height; k++) {
                    free(ptr[b][k]);
                }
                free(ptr[b]);
            }
            free(ptr);
        }
    } else if (gp == GRID_POS_EDGE) {
        if (gd == GRID_DIM_2D) {
            void **ptr = (void * *) var;
            for (int b = 0; b < g->eBlkCnt; b++) {
                free(ptr[b]);
            }
            free(ptr);
        } else if (gd == GRID_DIM_3D) {
            void ***ptr = (void * * *) var;
            for (int b = 0; b < g->eBlkCnt; b++) {
                for (int k = 0; k < g->height; k++) {
                    free(ptr[b][k]);
                }
                free(ptr[b]);
            }
            free(ptr);
        }
    }
    return (void *) 0;
}
