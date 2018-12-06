#include <mpi.h>

#include <netcdf.h>
#include "grid.h"
#include "io.h"
#include <assert.h>
static int ncid_w = 0;
static int ncid_r = 0;
// NOTE: do we need extra dimensions for edge variables?
static int Dims_2D_Cell[2];
static int Dims_1D_Cell[1];
static int Dims_2D_Edge[2];
static int Dims_1D_Edge[1];
write_queue WriteQueue;
// NOTE: this is for buffering the IO
int IOBuffersInitialized = 0;
float **CellBuffer2D;
float **EdgeBuffer2D;
float *CellBuffer1D;
float *EdgeBuffer1D;
void InitIOBuffers(GRID * g)
{
    if (!IOBuffersInitialized) {
        CellBuffer2D = malloc((sizeof(float *) * g->height));
        for (int i = 0; i < g->height; i++)
            CellBuffer2D[i] = malloc((sizeof(float) * g->cellCount));
        EdgeBuffer2D = malloc((sizeof(float *) * g->height));
        for (int i = 0; i < g->height; i++)
            EdgeBuffer2D[i] = malloc((sizeof(float) * g->edgeCount));
        CellBuffer1D = malloc((sizeof(float) * g->cellCount));
        EdgeBuffer1D = malloc((sizeof(float) * g->edgeCount));
        IOBuffersInitialized = 1;
    }
}

void FreeIOBuffers(GRID * g)
{
    if (IOBuffersInitialized) {
        for (int i = 0; i < g->height; i++)
            free(CellBuffer2D[i]);
        for (int i = 0; i < g->height; i++)
            free(EdgeBuffer2D[i]);
        free(CellBuffer1D);
        free(EdgeBuffer1D);
        free(CellBuffer2D);
        free(EdgeBuffer2D);
        IOBuffersInitialized = 0;
    }
}

///////////////////////////////////////////
// THIS IS BROKEN FOR MULTIPLE PROCESSES //
///////////////////////////////////////////
/*
  we either need to implement pnetcdf or some kind of mpi sync
*/
// NOTE: just for testing
#include <sys/time.h>
double stime()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double r = 1000000.0 * tv.tv_sec + tv.tv_usec;
    return r / 1000000.0;
}

#define PrintErrorAndExit(error) if(Error){printf("NetCDF ERROR: Received errorcode %d in function %s (Line: %d).\n  - Error message: %s\n", Error, __func__ ,__LINE__, nc_strerror(Error));exit(1);}
void io_write_init(GRID * g, const char *filename)
{
    int ret, Error;
    Error = nc_create(filename, NC_CLOBBER, &ncid_w);
//| NC_NETCDF4
    PrintErrorAndExit(Error);
    int DimHeight, DimCell, DimEdge;
    Error = nc_def_dim(ncid_w, "height", g->height, &DimHeight);
    PrintErrorAndExit(Error);
    Error = nc_def_dim(ncid_w, "cell", g->cellCount, &DimCell);
    PrintErrorAndExit(Error);
    Error = nc_def_dim(ncid_w, "edge", g->edgeCount, &DimEdge);
    PrintErrorAndExit(Error);
// Create 2D dimensions
    Dims_2D_Cell[0] = DimHeight;
    Dims_2D_Cell[1] = DimCell;
    Dims_2D_Edge[0] = DimHeight;
    Dims_2D_Edge[1] = DimEdge;
// Create 1D dimensions
    Dims_1D_Cell[0] = DimCell;
    Dims_1D_Edge[0] = DimEdge;
    WriteQueue.Head = WriteQueue.Tail = WriteQueue.isLooped = 0;
    InitIOBuffers(g);
}

void io_read_init(GRID * g, const char *filename)
{
    int Error = nc_open(filename, NC_CLOBBER, &ncid_r);
// | NC_NETCDF4
    PrintErrorAndExit(Error);
    InitIOBuffers(g);
}

// TODO: time this function
void io_read_register(GRID * g, const char *variable_name, GVAL * out_buff, datatype_t file_type, datatype_t memory_type, GRID_VAR_POSITION var_pos, GRID_VAR_DIMENSION var_dim)
{
    int ret;
    int varid;
// Get the variable ID based on its name
    int Error = nc_inq_varid(ncid_r, variable_name, &varid);
    PrintErrorAndExit(Error);
// Check the position and dimesion
    int Pos, Dim;
    Error = nc_get_att_int(ncid_r, varid, "position", &Pos);
    PrintErrorAndExit(Error);
    Error = nc_get_att_int(ncid_r, varid, "dimension", &Dim);
    PrintErrorAndExit(Error);
    assert(Pos == var_pos);
    assert(Dim == var_dim);
    if (var_dim == GRID_DIM_3D) {
        if (var_pos == GRID_POS_CELL) {
// NOTE: read data... still need to find out why get_var does not work
            for (int h = 0; h < g->height; h++) {
                size_t Start[] = { h, 0 };
                size_t Count[] = { 1, g->cellCount };
                Error = nc_get_vara_float(ncid_r, varid, Start, Count, &CellBuffer2D[h][0]);
            }
            struct {
                char *name;
                int loc;
                int dim;
                union {
                    GVAL *restrict * restrict p2;
                    GVAL *restrict * restrict * restrict p3;
                } data_pointer;
            } *var = (struct {
                      char *name; int loc; int dim; union {
                      GVAL * restrict * restrict p2; GVAL * restrict * restrict * restrict p3;} data_pointer;} * *restrict * restrict * restrict) out_buff;
/* TODO: find out why memcpy does not work
	    FOREACH cell IN grid | cell{0..0}
	    {
		int cb,ce;
		get_indices_c(g,cell.block,&cb,&ce);
		memcpy(&var[cell],&CellBuffer2D[cell.height][cell.block*g->blkSize], ce);
	    }
	    */
            {
                size_t min_block = g->mpi_rank == (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
                size_t max_block = g->mpi_rank < (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
                for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
                    for (size_t height_index = (0); height_index < (g->height); height_index++) {
                        for (size_t cell_index = (0); cell_index < (g->blkSize); cell_index++) {
                            var->data_pointer.p3[(block_index)][(height_index)][(cell_index)] = CellBuffer2D[height_index][(block_index * g->blkSize) + cell_index];
                        }
                    }
                }
            }
        } else {
            for (int h = 0; h < g->height; h++) {
                size_t Start[] = { h, 0 };
                size_t Count[] = { 1, g->edgeCount };
                Error = nc_get_vara_float(ncid_r, varid, Start, Count, &EdgeBuffer2D[h][0]);
            }
            struct {
                char *name;
                int loc;
                int dim;
                union {
                    GVAL *restrict * restrict p2;
                    GVAL *restrict * restrict * restrict p3;
                } data_pointer;
            } *var = (struct {
                      char *name; int loc; int dim; union {
                      GVAL * restrict * restrict p2; GVAL * restrict * restrict * restrict p3;} data_pointer;} * *restrict * restrict * restrict) out_buff;
//Error = nc_get_var_float(ncid_r, varid, &EdgeBuffer2D[0][0]);
/*
	    FOREACH edge IN grid | edge{0..0}
	    {
		int eb,ee;
		get_indices_e(g,edge.block,&eb,&ee);
		memcpy(&var[edge],&EdgeBuffer2D[edge.height][edge.block*g->blkSize], ee);
	    }
	    */
            {
                size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
                size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
                for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
                    for (size_t height_index = (0); height_index < (g->height); height_index++) {
                        for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                            var->data_pointer.p3[(block_index)][(height_index)][(edge_index)] = EdgeBuffer2D[height_index][(block_index * g->blkSize) + edge_index];
                        }
                    }
                }
            }
        }
    } else {
        if (var_pos == GRID_POS_CELL) {
            Error = nc_get_var_float(ncid_r, varid, &CellBuffer1D[0]);
            struct {
                char *name;
                int loc;
                int dim;
                union {
                    GVAL *restrict * restrict p2;
                    GVAL *restrict * restrict * restrict p3;
                } data_pointer;
            } *var = (struct {
                      char *name; int loc; int dim; union {
                      GVAL * restrict * restrict p2; GVAL * restrict * restrict * restrict p3;} data_pointer;} * *restrict * restrict) out_buff;
/*
	    FOREACH cell IN gridCELL2D | cell{0..0}
	    {
		int cb,ce;
		get_indices_c(g,cell.block,&cb,&ce);
		memcpy(&var[cell],&CellBuffer1D[cell.block*g->blkSize], ce);
	    }
	    */
            {
                size_t min_block = g->mpi_rank == (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
                size_t max_block = g->mpi_rank < (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
                for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
                    for (size_t cell_index = (0); cell_index < (g->blkSize); cell_index++) {
                        var->data_pointer.p2[(block_index)][(cell_index)] = CellBuffer1D[(block_index * g->blkSize) + cell_index];
                    }
                }
            }
        } else {
            Error = nc_get_var_float(ncid_r, varid, &EdgeBuffer1D[0]);
            struct {
                char *name;
                int loc;
                int dim;
                union {
                    GVAL *restrict * restrict p2;
                    GVAL *restrict * restrict * restrict p3;
                } data_pointer;
            } *var = (struct {
                      char *name; int loc; int dim; union {
                      GVAL * restrict * restrict p2; GVAL * restrict * restrict * restrict p3;} data_pointer;} * *restrict * restrict) out_buff;
/*
	    FOREACH edge IN gridEDGE2D | edge{0..0}
	    {
		int eb,ee;
		get_indices_e(g,edge.block,&eb,&ee);
		memcpy(&var[edge],&EdgeBuffer1D[edge.block*g->blkSize], ee);
	    }
	    */
            {
                size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
                size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
                for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
                    for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                        var->data_pointer.p2[(block_index)][(edge_index)] = EdgeBuffer1D[(block_index * g->blkSize) + edge_index];
                    }
                }
            }
        }
    }
    PrintErrorAndExit(Error);
}

void io_read_start()
{
}

// does nothing in this implementation
void io_write_define(GRID * g, const char *variable_name, GVAL * out_buff, datatype_t file_type, GRID_VAR_POSITION var_pos, GRID_VAR_DIMENSION var_dim, io_var_t * out_varid)
{
    int ret, Error;
    int *Dims_2D = Dims_2D_Cell;
    int *Dims_1D = Dims_1D_Cell;
    if (var_pos == GRID_POS_EDGE) {
        Dims_2D = Dims_2D_Edge;
        Dims_1D = Dims_1D_Edge;
    }
//TODO: implement other datatypes!
    assert(file_type == FLOAT32);
    if (var_dim == GRID_DIM_3D) {
        Error = nc_def_var(ncid_w, variable_name, NC_FLOAT, 2, Dims_2D, &out_varid->id);
        PrintErrorAndExit(Error);
    } else {
        Error = nc_def_var(ncid_w, variable_name, NC_FLOAT, 1, Dims_1D, &out_varid->id);
        PrintErrorAndExit(Error);
    }
    out_varid->memory = out_buff;
    out_varid->Dim = var_dim;
    out_varid->Pos = var_pos;
    out_varid->Name = variable_name;
// Write out attributes
    Error = nc_put_att_int(ncid_w, out_varid->id, "position", NC_INT, 1, (int *) &(out_varid->Pos));
    PrintErrorAndExit(Error);
    Error = nc_put_att_int(ncid_w, out_varid->id, "dimension", NC_INT, 1, (int *) &(out_varid->Dim));
    PrintErrorAndExit(Error);
}

void io_write_registration_complete(GRID * g)
{
    int Error = nc_enddef(ncid_w);
    PrintErrorAndExit(Error);
}

void io_write_announce(GRID * g, io_var_t * var)
{
// NOTE: check if variables is already in the queue
    int Cur = WriteQueue.Tail;
    int isLooped = WriteQueue.isLooped;
    int isEnqueued = 0;
    while (!(isLooped == 0 && Cur == WriteQueue.Head))  // NOTE: while !isEmpty
    {
        int Id = WriteQueue.Elements[Cur].id;
        if (Id == var->id)
            isEnqueued = 1;
        if (Cur == ArrayCount(WriteQueue.Elements) - 1)
            isLooped = 0;
        Cur = ((Cur + 1) % (ArrayCount(WriteQueue.Elements)));
    }
// NOTE: if not, enqueue it
    if (!isEnqueued) {
        WriteQueue.Elements[WriteQueue.Head] = *var;
        WriteQueue.Head = (WriteQueue.Head + 1);
        if (WriteQueue.Head == ArrayCount(WriteQueue.Elements)) {
            WriteQueue.Head = 0;
            WriteQueue.isLooped = 1;
        }
        assert(!(WriteQueue.isLooped && (WriteQueue.Head == WriteQueue.Tail)));
    }
}

// NOTE: if this fires,
// the queue is to small.
// either increase the size
// or move to better data structure
void io_write_start(GRID * g)
{
    while (!(WriteQueue.isLooped == 0 && WriteQueue.Head == WriteQueue.Tail))   // NOTE: while !isEmpty
    {
        double StartTime = stime();
        io_var_t IoVar = WriteQueue.Elements[WriteQueue.Tail];
// is not used anywhere. remove?
        int ret, Error;
        int BlkCnt = g->cBlkCnt;
        if (IoVar.Pos == GRID_POS_EDGE) {
            BlkCnt = g->eBlkCnt;
        }
// Write out data
        if (IoVar.Dim == GRID_DIM_3D) {
            size_t Start[2];
            size_t Count[2];
            if (IoVar.Pos == GRID_POS_CELL) {
                struct {
                    char *name;
                    int loc;
                    int dim;
                    union {
                        GVAL *restrict * restrict p2;
                        GVAL *restrict * restrict * restrict p3;
                    } data_pointer;
                } *p = IoVar.memory;
/* TODO: why can't I use memcpy
		FOREACH cell IN grid | cell{0..0}
		{
		    int cb,ce;
		    get_indices_c(g,cell.block,&cb,&ce);
		                                  
		    memcpy(&Buffer[cell.height][cell.block*g->blkSize], &p[cell], ce);
		}
		*/
// Buffer = CellBuffer2D;                  
                {
                    size_t min_block = g->mpi_rank == (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
                    size_t max_block = g->mpi_rank < (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
                    for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
                        for (size_t height_index = (0); height_index < (g->height); height_index++) {
                            for (size_t cell_index = (0); cell_index < (g->blkSize); cell_index++) {
                                CellBuffer2D[height_index][(block_index * g->blkSize) + cell_index] = p->data_pointer.p3[(block_index)][(height_index)][(cell_index)];
                            }
                        }
                    }
                }
                for (int h = 0; h < g->height; h++) {
                    size_t Start[] = { h, 0 };
                    size_t Count[] = { 1, g->cellCount };
                    Error = nc_put_vara_float(ncid_w, IoVar.id, Start, Count, &CellBuffer2D[h][0]);
                }
            } else              //Error = nc_put_var_float(ncid_w, IoVar.id, &CellBuffer2D[0][0]);
            {
                struct {
                    char *name;
                    int loc;
                    int dim;
                    union {
                        GVAL *restrict * restrict p2;
                        GVAL *restrict * restrict * restrict p3;
                    } data_pointer;
                } *p = IoVar.memory;
/*
		FOREACH edge IN grid | edge{0..0}
		{
		    int eb,ee;
		    get_indices_e(g,edge.block,&eb,&ee);
		    memcpy(&EdgeBuffer2D[edge.height][edge.block*g->blkSize], &p[edge], ee);
		}
		*/
                {
                    size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
                    size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
                    for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
                        for (size_t height_index = (0); height_index < (g->height); height_index++) {
                            for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                                EdgeBuffer2D[height_index][(block_index * g->blkSize) + edge_index] = p->data_pointer.p3[(block_index)][(height_index)][(edge_index)];
                            }
                        }
                    }
                }
                for (int h = 0; h < g->height; h++) {
                    size_t Start[] = { h, 0 };
                    size_t Count[] = { 1, g->edgeCount };
                    Error = nc_put_vara_float(ncid_w, IoVar.id, Start, Count, &EdgeBuffer2D[h][0]);
                }
            }
        } else                  //Error = nc_put_var_float(ncid_w, IoVar.id, &Buffer[0][0]);
        {
            if (IoVar.Pos == GRID_POS_CELL) {
                struct {
                    char *name;
                    int loc;
                    int dim;
                    union {
                        GVAL *restrict * restrict p2;
                        GVAL *restrict * restrict * restrict p3;
                    } data_pointer;
                } *p = IoVar.memory;
/*
		FOREACH cell IN gridCELL2D | cell{0..0}
		{
		    int cb,ce;
		    get_indices_c(g,cell.block,&cb,&ce);
		    memcpy(&Buffer[cell.block*g->blkSize], &p[cell], ce);
		}
		*/
                {
                    size_t min_block = g->mpi_rank == (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
                    size_t max_block = g->mpi_rank < (0) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->cBlkCnt - 1) / (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->cBlkCnt % (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
                    for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
                        for (size_t cell_index = (0); cell_index < (g->blkSize); cell_index++) {
                            CellBuffer1D[(block_index * g->blkSize) + cell_index] = p->data_pointer.p2[(block_index)][(cell_index)];
                        }
                    }
                }
                Error = nc_put_var_float(ncid_w, IoVar.id, &CellBuffer1D[0]);
            } else {
                struct {
                    char *name;
                    int loc;
                    int dim;
                    union {
                        GVAL *restrict * restrict p2;
                        GVAL *restrict * restrict * restrict p3;
                    } data_pointer;
                } *p = IoVar.memory;
/*
		FOREACH edge IN gridEDGE2D | edge{0..0}
		{
		    int eb,ee;
		    get_indices_e(g,edge.block,&eb,&ee);
		    memcpy(&Buffer[edge.block*g->blkSize], &p[edge], ee);
		}
		*/
                {
                    size_t min_block = g->mpi_rank == (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
                    size_t max_block = g->mpi_rank < (0) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? 0 : g->mpi_rank == (g->eBlkCnt - 1) / (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ? g->eBlkCnt % (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : (((g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
#pragma omp parallel for
                    for (size_t block_index = (min_block); block_index < (max_block); block_index++) {
                        for (size_t edge_index = (0); edge_index < (g->blkSize); edge_index++) {
                            EdgeBuffer1D[(block_index * g->blkSize) + edge_index] = p->data_pointer.p2[(block_index)][(edge_index)];
                        }
                    }
                }
                Error = nc_put_var_float(ncid_w, IoVar.id, &EdgeBuffer1D[0]);
            }
        }
        PrintErrorAndExit(Error);
        if (WriteQueue.Tail == ArrayCount(WriteQueue.Elements) - 1)
            WriteQueue.isLooped = 0;
        WriteQueue.Tail = ((WriteQueue.Tail + 1) % (ArrayCount(WriteQueue.Elements)));
        double ElapsedTime = stime() - StartTime;
        printf("ElapsedTime: %f\n", ElapsedTime);
    }
}

void io_write_finalize(GRID * g)
{
    FreeIOBuffers(g);
    nc_close(ncid_w);
}
