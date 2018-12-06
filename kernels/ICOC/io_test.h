#include <mpi.h>

extern int local_cell_blocks;
extern int local_edge_blocks;
#include "io.h"
void TestIO(GRID * g)
{
// TODO: seed the rngfalse
    printf("Init write...\n");
    io_write_init(g, "temp_netcdf_test_output.cdf");
    struct {
        char *name;
        int loc;
        int dim;
        union {
            GVAL *restrict * restrict p2;
            GVAL *restrict * restrict * restrict p3;
        } data_pointer;
    } *test_output1;
    io_var_t io_test_output1;
    struct {
        char *name;
        int loc;
        int dim;
        union {
            GVAL *restrict * restrict p2;
            GVAL *restrict * restrict * restrict p3;
        } data_pointer;
    } *test_input1;
    io_var_t io_test_input1;
    struct {
        char *name;
        int loc;
        int dim;
        union {
            GVAL *restrict * restrict p2;
            GVAL *restrict * restrict * restrict p3;
        } data_pointer;
    } *test_output2;
    io_var_t io_test_output2;
    struct {
        char *name;
        int loc;
        int dim;
        union {
            GVAL *restrict * restrict p2;
            GVAL *restrict * restrict * restrict p3;
        } data_pointer;
    } *test_input2;
    io_var_t io_test_input2;
    struct {
        char *name;
        int loc;
        int dim;
        union {
            GVAL *restrict * restrict p2;
            GVAL *restrict * restrict * restrict p3;
        } data_pointer;
    } *test_output3;
    io_var_t io_test_output3;
    struct {
        char *name;
        int loc;
        int dim;
        union {
            GVAL *restrict * restrict p2;
            GVAL *restrict * restrict * restrict p3;
        } data_pointer;
    } *test_input3;
    io_var_t io_test_input3;
    struct {
        char *name;
        int loc;
        int dim;
        union {
            GVAL *restrict * restrict p2;
            GVAL *restrict * restrict * restrict p3;
        } data_pointer;
    } *test_output4;
    io_var_t io_test_output4;
    struct {
        char *name;
        int loc;
        int dim;
        union {
            GVAL *restrict * restrict p2;
            GVAL *restrict * restrict * restrict p3;
        } data_pointer;
    } *test_input4;
    io_var_t io_test_input4;
    {
        int num_blocks = local_cell_blocks ? local_cell_blocks : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        test_output1 = malloc(24);
        test_output1->name = "test_output1";
        test_output1->loc = 0;
        test_output1->dim = 3;
        test_output1->data_pointer.p3 = malloc((num_blocks * g->height * g->blkSize) * sizeof(GVAL) + (num_blocks * g->height) * sizeof(char *) + (num_blocks) * sizeof(char *));
        char *pos = (char *) test_output1->data_pointer.p3 + num_blocks * sizeof(char *);
        char *pos2 = (char *) test_output1->data_pointer.p3 + num_blocks * sizeof(char *) + num_blocks * g->height * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            test_output1->data_pointer.p3[b] = (GVAL * *)pos;
            pos += g->height * sizeof(char *);
            for (int k = 0; k < g->height; k++) {
                test_output1->data_pointer.p3[b][k] = (GVAL *) pos2;
                pos2 += g->blkSize * sizeof(GVAL);
                for (int c = 0; c < g->blkSize; c++) {
                    test_output1->data_pointer.p3[b][k][c] = (GVAL) 0;
                }
            }
        }
    }
    {
        int num_blocks = local_cell_blocks ? local_cell_blocks : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        test_input1 = malloc(24);
        test_input1->name = "test_input1";
        test_input1->loc = 0;
        test_input1->dim = 3;
        test_input1->data_pointer.p3 = malloc((num_blocks * g->height * g->blkSize) * sizeof(GVAL) + (num_blocks * g->height) * sizeof(char *) + (num_blocks) * sizeof(char *));
        char *pos = (char *) test_input1->data_pointer.p3 + num_blocks * sizeof(char *);
        char *pos2 = (char *) test_input1->data_pointer.p3 + num_blocks * sizeof(char *) + num_blocks * g->height * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            test_input1->data_pointer.p3[b] = (GVAL * *)pos;
            pos += g->height * sizeof(char *);
            for (int k = 0; k < g->height; k++) {
                test_input1->data_pointer.p3[b][k] = (GVAL *) pos2;
                pos2 += g->blkSize * sizeof(GVAL);
                for (int c = 0; c < g->blkSize; c++) {
                    test_input1->data_pointer.p3[b][k][c] = (GVAL) 0;
                }
            }
        }
    }
    {
        int num_blocks = local_edge_blocks ? local_edge_blocks : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        test_output2 = malloc(24);
        test_output2->name = "test_output2";
        test_output2->loc = 1;
        test_output2->dim = 3;
        test_output2->data_pointer.p3 = malloc((num_blocks * g->height * g->blkSize) * sizeof(GVAL) + (num_blocks * g->height) * sizeof(char *) + (num_blocks) * sizeof(char *));
        char *pos = (char *) test_output2->data_pointer.p3 + num_blocks * sizeof(char *);
        char *pos2 = (char *) test_output2->data_pointer.p3 + num_blocks * sizeof(char *) + num_blocks * g->height * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            test_output2->data_pointer.p3[b] = (GVAL * *)pos;
            pos += g->height * sizeof(char *);
            for (int k = 0; k < g->height; k++) {
                test_output2->data_pointer.p3[b][k] = (GVAL *) pos2;
                pos2 += g->blkSize * sizeof(GVAL);
                for (int e = 0; e < g->blkSize; e++) {
                    test_output2->data_pointer.p3[b][k][e] = (GVAL) 0;
                }
            }
        }
    }
    {
        int num_blocks = local_edge_blocks ? local_edge_blocks : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        test_input2 = malloc(24);
        test_input2->name = "test_input2";
        test_input2->loc = 1;
        test_input2->dim = 3;
        test_input2->data_pointer.p3 = malloc((num_blocks * g->height * g->blkSize) * sizeof(GVAL) + (num_blocks * g->height) * sizeof(char *) + (num_blocks) * sizeof(char *));
        char *pos = (char *) test_input2->data_pointer.p3 + num_blocks * sizeof(char *);
        char *pos2 = (char *) test_input2->data_pointer.p3 + num_blocks * sizeof(char *) + num_blocks * g->height * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            test_input2->data_pointer.p3[b] = (GVAL * *)pos;
            pos += g->height * sizeof(char *);
            for (int k = 0; k < g->height; k++) {
                test_input2->data_pointer.p3[b][k] = (GVAL *) pos2;
                pos2 += g->blkSize * sizeof(GVAL);
                for (int e = 0; e < g->blkSize; e++) {
                    test_input2->data_pointer.p3[b][k][e] = (GVAL) 0;
                }
            }
        }
    }
    {
        int num_blocks = local_cell_blocks ? local_cell_blocks : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        test_output3 = malloc(24);
        test_output3->name = "test_output3";
        test_output3->loc = 0;
        test_output3->dim = 2;
        test_output3->data_pointer.p2 = malloc((num_blocks * g->blkSize) * sizeof(GVAL) + (num_blocks) * sizeof(char *));
        char *pos = (char *) test_output3->data_pointer.p2 + num_blocks * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            test_output3->data_pointer.p2[b] = (GVAL *) pos;
            pos += g->blkSize * sizeof(GVAL);
            for (int c = 0; c < g->blkSize; c++) {
                test_output3->data_pointer.p2[b][c] = (GVAL) 0;
            }
        }
    }
    {
        int num_blocks = local_cell_blocks ? local_cell_blocks : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        test_input3 = malloc(24);
        test_input3->name = "test_input3";
        test_input3->loc = 0;
        test_input3->dim = 2;
        test_input3->data_pointer.p2 = malloc((num_blocks * g->blkSize) * sizeof(GVAL) + (num_blocks) * sizeof(char *));
        char *pos = (char *) test_input3->data_pointer.p2 + num_blocks * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            test_input3->data_pointer.p2[b] = (GVAL *) pos;
            pos += g->blkSize * sizeof(GVAL);
            for (int c = 0; c < g->blkSize; c++) {
                test_input3->data_pointer.p2[b][c] = (GVAL) 0;
            }
        }
    }
    {
        int num_blocks = local_edge_blocks ? local_edge_blocks : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        test_output4 = malloc(24);
        test_output4->name = "test_output4";
        test_output4->loc = 1;
        test_output4->dim = 2;
        test_output4->data_pointer.p2 = malloc((num_blocks * g->blkSize) * sizeof(GVAL) + (num_blocks) * sizeof(char *));
        char *pos = (char *) test_output4->data_pointer.p2 + num_blocks * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            test_output4->data_pointer.p2[b] = (GVAL *) pos;
            pos += g->blkSize * sizeof(GVAL);
            for (int e = 0; e < g->blkSize; e++) {
                test_output4->data_pointer.p2[b][e] = (GVAL) 0;
            }
        }
    }
    {
        int num_blocks = local_edge_blocks ? local_edge_blocks : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
        test_input4 = malloc(24);
        test_input4->name = "test_input4";
        test_input4->loc = 1;
        test_input4->dim = 2;
        test_input4->data_pointer.p2 = malloc((num_blocks * g->blkSize) * sizeof(GVAL) + (num_blocks) * sizeof(char *));
        char *pos = (char *) test_input4->data_pointer.p2 + num_blocks * sizeof(char *);
        for (int b = 0; b < num_blocks; b++) {
            test_input4->data_pointer.p2[b] = (GVAL *) pos;
            pos += g->blkSize * sizeof(GVAL);
            for (int e = 0; e < g->blkSize; e++) {
                test_input4->data_pointer.p2[b][e] = (GVAL) 0;
            }
        }
    }
    printf("Generate test values...\n");
    {
        size_t min_block = g->mpi_rank == (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t height_index = (0); height_index < (g->height); height_index++) {
                for (size_t cell_index = (0); cell_index < (g->blkSize); cell_index++) {
                    test_output1->data_pointer.p3[(block_index)][(height_index)][(cell_index)] = (float) rand();
                }
            }
        }
    }
    {
        size_t min_block = g->mpi_rank == (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t cell_index = (0); cell_index < (g->blkSize); cell_index++) {
                test_output3->data_pointer.p2[(block_index)][(cell_index)] = (float) rand();
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
                    test_output2->data_pointer.p3[(block_index)][(height_index)][(edge_index)] = (float) rand();
                }
            }
        }
    }
    {
        size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                test_output4->data_pointer.p2[(block_index)][(edge_index)] = (float) rand();
            }
        }
    }
    printf("Setup output variable...\n");
    io_write_define(g, "test_output1", (GVAL *) test_output1, FLOAT32, GRID_POS_CELL, GRID_DIM_3D, &io_test_output1);
    io_write_define(g, "test_output2", (GVAL *) test_output2, FLOAT32, GRID_POS_EDGE, GRID_DIM_3D, &io_test_output2);
    io_write_define(g, "test_output3", (GVAL *) test_output3, FLOAT32, GRID_POS_CELL, GRID_DIM_2D, &io_test_output3);
    io_write_define(g, "test_output4", (GVAL *) test_output4, FLOAT32, GRID_POS_EDGE, GRID_DIM_2D, &io_test_output4);
    io_write_registration_complete(g);
    io_write_announce(g, &io_test_output1);
    io_write_announce(g, &io_test_output2);
    io_write_announce(g, &io_test_output3);
    io_write_announce(g, &io_test_output4);
    printf("Writing to disk...\n");
    io_write_start(g);
    io_write_finalize(g);
    printf("Init read...\n");
    io_read_init(g, "temp_netcdf_test_output.cdf");
    io_read_register(g, "test_output1", (GVAL *) test_input1, FLOAT32, FLOAT32, GRID_POS_CELL, GRID_DIM_3D);
    io_read_register(g, "test_output2", (GVAL *) test_input2, FLOAT32, FLOAT32, GRID_POS_EDGE, GRID_DIM_3D);
    io_read_register(g, "test_output3", (GVAL *) test_input3, FLOAT32, FLOAT32, GRID_POS_CELL, GRID_DIM_2D);
    io_read_register(g, "test_output4", (GVAL *) test_input4, FLOAT32, FLOAT32, GRID_POS_EDGE, GRID_DIM_2D);
    printf("Reading from disk...\n");
    io_read_start();
    printf("Comparing initial values with disk values...\n");
    int success = 1;
    size_t total = 0;
    size_t count = 0;
    {
        size_t min_block = g->mpi_rank == (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t height_index = (0); height_index < (g->height); height_index++) {
                for (size_t cell_index = (0); cell_index < (g->blkSize); cell_index++) {
                    if (test_output1->data_pointer.p3[(block_index)][(height_index)][(cell_index)] != test_input1->data_pointer.p3[(block_index)][(height_index)][(cell_index)]) {
                        count++;
                        success = 0;
                    }
                }
            }
        }
    }
    printf("%d/%d errors for CELL 3D\n", count, g->height * g->cellCount);
    total += count;
    count = 0;
    {
        size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t height_index = (0); height_index < (g->height); height_index++) {
                for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                    if (test_output2->data_pointer.p3[(block_index)][(height_index)][(edge_index)] != test_input2->data_pointer.p3[(block_index)][(height_index)][(edge_index)]) {
                        count++;
                        success = 0;
                    }
                }
            }
        }
    }
    printf("%d/%d errors for EDGE 3D\n", count, g->height * g->edgeCount);
    total += count;
    count = 0;
    {
        size_t min_block = g->mpi_rank == (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t cell_index = (0); cell_index < (g->blkSize); cell_index++) {
                if (test_output3->data_pointer.p2[(block_index)][(cell_index)] != test_input3->data_pointer.p2[(block_index)][(cell_index)]) {
                    count++;
                    success = 0;
                }
            }
        }
    }
    printf("%d/%d errors for CELL 2D\n", count, g->cellCount);
    total += count;
    count = 0;
    {
        size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
        size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
        for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
            for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                if (test_output4->data_pointer.p2[(block_index)][(edge_index)] != test_input4->data_pointer.p2[(block_index)][(edge_index)]) {
                    count++;
                    success = 0;
                }
            }
        }
    }
    printf("%d/%d errors for EDGE 2D\n", count, g->edgeCount);
    total += count;
    count = 0;
    printf("%d/%d errors total\n", count, g->height * g->cellCount + g->height * g->edgeCount + g->cellCount + g->cellCount);
    if (success) {
        printf("\nnetcdf io test\x1B[32m succeded\x1B[0m\n");
    } else {
        printf("\nnetcdf io test\x1B[31m failed\x1B[0m\n");
    }
}
