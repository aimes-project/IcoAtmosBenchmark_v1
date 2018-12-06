#include <mpi.h>

#include <sys/time.h>
#include <stdint.h>
#include "option.h"
#include "grid.h"
#include "component.h"
#include "io.h"
#include "io_test.h"
double mtime()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double r = 1000000.0 * tv.tv_sec + tv.tv_usec;
    return r / 1000000.0;
}

static int timesteps = 100;
static double maxtime = 5.0;
static double mintime = 0.0;
static int gridsize = 4;
static int gridheight = 64;
static char *io_init_file;
static char *io_output_file;
static int iostep = 3;
static int *test_netcdf_io = 0;
static option_help options[] = { {'T', "timesteps", "The number of timesteps.", OPTION_OPTIONAL_ARGUMENT, 'd', &timesteps}, {'G', "gridsize", "The dimension of the grid in G*G.H", OPTION_OPTIONAL_ARGUMENT, 'd', &gridsize}, {'H', "gridheight", "The dimension of the grid in G*G*H.", OPTION_OPTIONAL_ARGUMENT, 'd', &gridheight}, {'M', "maxtime", "The maximum runtime in seconds.", OPTION_OPTIONAL_ARGUMENT, 'F', &maxtime}, {'m', "mintime", "The minimum time in seconds.", OPTION_OPTIONAL_ARGUMENT, 'F', &mintime}, {'R', "io_init_file", "The I/O initialization file (to read initial values).", OPTION_REQUIRED_ARGUMENT, 's', &io_init_file}, {'W', "io_output_file", "The output file.", OPTION_REQUIRED_ARGUMENT, 's', &io_output_file}, {'t', "test_netcdf_io", "Tests the writing and reading of netcdf files with random values.", OPTION_OPTIONAL_ARGUMENT, 'd', &test_netcdf_io}, LAST_OPTION };

extern MODEL_COMPONENT com1;
extern MODEL_COMPONENT com2;
extern MODEL_COMPONENT com3;
extern MODEL_COMPONENT com4;
extern MODEL_COMPONENT com5;
int main(int argc, char * *argv)
{
    int printhelp = 0;
    parseOptions(argc, argv, options, &printhelp);
    if (printhelp != 0) {
        print_help(options, 0);
        exit(0);
    }
    GRID *g = malloc(sizeof(GRID));
    {
        MPI_Init(NULL, NULL);
        MPI_Comm_size(MPI_COMM_WORLD, &g->mpi_world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &g->mpi_rank);
    }
    init_grid(g, gridsize, gridheight);
//currently grid size for one dimension
// Netcdf IO test
    if (test_netcdf_io) {
        TestIO(g);
        {
            MPI_Finalize();
        }
        return 0;
    }
//const char *io_init_file = "test_gv_temp2";
//const char *io_output_file = "output.cdf";
    io_read_init(g, io_init_file);
    io_write_init(g, io_output_file);
    if (!com1.loaded)
        com1.init(g);
/*if(!com2.loaded)com2.init(g);
    if(!com3.loaded)com3.init(g);
    if(!com4.loaded)com4.init(g);*/
    if (!com5.loaded)
        com5.init(g);
    io_write_registration_complete(g);
    io_read_start();
    uint64_t ts;
    double dt = 0;
    double st = mtime();
    {
        for (ts = 0; ((ts < timesteps) || (dt < mintime)) && dt < maxtime; ts++) {
            com5.compute(g);
            com1.io(g);
            com5.io(g);
            io_write_start(g);
            dt = mtime() - st;
        }
    }
/*com1.compute(g);
        com2.compute(g);
        com3.compute(g);
        com4.compute(g);*/
// TODO: io step
//com2.io(g);
//com3.io(g);
//com4.io(g);
/*uint64_t cs=com1.checksum(g)+com2.checksum(g)+com3.checksum(g)+com4.checksum(g);
    uint64_t gcs=0;
    MPI_Reduce(&cs, &gcs, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0,MPI_COMM_WORLD);
    if(g->mpi_rank==0)printf("checksum=%X\n",gcs);*/
    if (g->mpi_rank == 0) {
        double gflops, mem;
        gflops = (              /*com1.flops(g) + com2.flops(g) + com3.flops(g) + com4.flops(g) + */
                     com5.flops(g)) * ts / (1024 * 1024 * 1024) / dt;
        mem = (com1.memory(g));
/* + com2.memory(g) + com3.memory(g) + com4.memory(g) + com5.memory(g)*/
        printf("%ld,%f,%f,%f\n", ts, dt, gflops, mem);
        uint64_t cs = com1.checksum(g) + com5.checksum(g);
        uint64_t gcs = 0;
        MPI_Reduce(&cs, &gcs, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        if (g->mpi_rank == 0)
            printf("checksum=%X\n", gcs);
    }
    com1.cleanup(g);
/*com2.cleanup(g);
    com3.cleanup(g);
    com4.cleanup(g);*/
    com5.cleanup(g);
    io_write_finalize(g);
    {
        MPI_Finalize();
    }
}
