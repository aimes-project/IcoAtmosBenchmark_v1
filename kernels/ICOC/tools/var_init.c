#include <mpi.h>

extern int local_cell_blocks;
extern int local_edge_blocks;
#include <stdint.h>
#include "../io.h"
#include "../option.h"
#include "../grid.h"
static char  * Filename;
static int gridsize = 4;
static int gridheight = 64;
static option_help options[] = {{'G' , "gridsize" , "The dimension of the grid in G*G.H" , OPTION_OPTIONAL_ARGUMENT , 'd' , &gridsize} , {'H' , "gridheight" , "The dimension of the grid in G*G*H." , OPTION_OPTIONAL_ARGUMENT , 'd' , &gridheight} , {'f' , "filename" , "Path to the output file." , OPTION_REQUIRED_ARGUMENT , 's' , &Filename} , LAST_OPTION};
void Init_gv_temp(GRID  * g)
{
struct {
char  * name;
int loc;
int dim;
union {
GVAL  * restrict * restrict p2;
GVAL  * restrict * restrict * restrict p3;
} data_pointer;
}  * gv_temp;
{
int num_blocks = local_cell_blocks ?  local_cell_blocks : ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
 gv_temp = malloc( 24) ;
 gv_temp->name = "gv_temp" ;
 gv_temp->loc = 0 ;
 gv_temp->dim = 3 ;
 gv_temp->data_pointer.p3 = malloc( ( num_blocks * g->height * g->blkSize) * sizeof(GVAL) + ( num_blocks * g->height) * sizeof(char  * ) + ( num_blocks) * sizeof(char  * )) ;
char  * pos = (char  * ) gv_temp->data_pointer.p3 + num_blocks * sizeof(char  * );
char  * pos2 = (char  * ) gv_temp->data_pointer.p3 + num_blocks * sizeof(char  * ) + num_blocks * g->height * sizeof(char  * );
for(int b = 0; b < num_blocks ; b++) {
 gv_temp->data_pointer.p3[ b] = (GVAL  * * ) pos ;
 pos += g->height * sizeof(char  * ) ;
for(int k = 0; k < g->height ; k++) {
 gv_temp->data_pointer.p3[ b][ k] = (GVAL  * ) pos2 ;
 pos2 += g->blkSize * sizeof(GVAL) ;
for(int c = 0; c < g->blkSize ; c++) {
 gv_temp->data_pointer.p3[ b][ k][ c] = (GVAL) 0 ;
}
}
}
}
io_var_t io_gv_temp;
{
size_t min_block = g->mpi_rank == ( 0) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  0 % ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
size_t max_block = g->mpi_rank < ( 0) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > ( g->cBlkCnt - 1) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  0 : g->mpi_rank == ( g->cBlkCnt - 1) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  g->cBlkCnt % ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  g->cBlkCnt % ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
 #pragma omp parallel for
for(size_t block_index = ( min_block); block_index < ( max_block) ; block_index++) {
for(size_t height_index = ( 0); height_index < ( g->height) ; height_index++) {
for(size_t cell_index = ( 0); cell_index < ( g->blkSize) ; cell_index++) {
if (  block_index == 0 && cell_index == 0 && g->mpi_rank == 0 )  gv_temp->data_pointer.p3[ ( block_index)][ ( height_index)][ ( cell_index)] = 100.0f ; else  gv_temp->data_pointer.p3[ ( block_index)][ ( height_index)][ ( cell_index)] = 0.0f ;
}
}
}
}
 io_write_define( g , "gv_temp" , (GVAL  * ) gv_temp , FLOAT32 , GRID_POS_CELL , GRID_DIM_3D , &io_gv_temp) ;
 io_write_announce( g , &io_gv_temp) ;
}
void Init_gv_ind2Dparam(GRID  * g)
{
struct {
char  * name;
int loc;
int dim;
union {
GVAL  * restrict * restrict p2;
GVAL  * restrict * restrict * restrict p3;
} data_pointer;
}  * gv_ind2Dparam;
{
int num_blocks = local_edge_blocks ?  local_edge_blocks : ( ( ( g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
 gv_ind2Dparam = malloc( 24) ;
 gv_ind2Dparam->name = "gv_ind2Dparam" ;
 gv_ind2Dparam->loc = 1 ;
 gv_ind2Dparam->dim = 2 ;
 gv_ind2Dparam->data_pointer.p2 = malloc( ( num_blocks * g->blkSize) * sizeof(GVAL) + ( num_blocks) * sizeof(char  * )) ;
char  * pos = (char  * ) gv_ind2Dparam->data_pointer.p2 + num_blocks * sizeof(char  * );
for(int b = 0; b < num_blocks ; b++) {
 gv_ind2Dparam->data_pointer.p2[ b] = (GVAL  * ) pos ;
 pos += g->blkSize * sizeof(GVAL) ;
for(int e = 0; e < g->blkSize ; e++) {
 gv_ind2Dparam->data_pointer.p2[ b][ e] = (GVAL) 0 ;
}
}
}
io_var_t io_gv_ind2Dparam;
{
size_t min_block = g->mpi_rank == ( 0) / ( ( ( g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  0 % ( ( ( g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
size_t max_block = g->mpi_rank < ( 0) / ( ( ( g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > ( g->eBlkCnt - 1) / ( ( ( g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  0 : g->mpi_rank == ( g->eBlkCnt - 1) / ( ( ( g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  g->eBlkCnt % ( ( ( g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  g->eBlkCnt % ( ( ( g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : ( ( ( g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : ( ( ( g->eBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
 #pragma omp parallel for
for(size_t block_index = ( min_block); block_index < ( max_block) ; block_index++) {
for(size_t edge_index = ( 0); edge_index < ( g->blkSize) ; edge_index++) {
 gv_ind2Dparam->data_pointer.p2[ ( block_index)][ ( edge_index)] = 1.0f ;
}
}
}
 io_write_define( g , "gv_ind2Dparam" , (GVAL  * ) gv_ind2Dparam , FLOAT32 , GRID_POS_EDGE , GRID_DIM_2D , &io_gv_ind2Dparam) ;
 io_write_announce( g , &io_gv_ind2Dparam) ;
}
void Init_gv_o8param(GRID  * g)
{
struct {
char  * name;
int loc;
int dim;
union {
GVAL  * restrict * restrict p2;
GVAL  * restrict * restrict * restrict p3;
} data_pointer;
}  * gv_o8param0;
{
int num_blocks = local_cell_blocks ?  local_cell_blocks : ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
 gv_o8param0 = malloc( 24) ;
 gv_o8param0->name = "gv_o8param0" ;
 gv_o8param0->loc = 0 ;
 gv_o8param0->dim = 2 ;
 gv_o8param0->data_pointer.p2 = malloc( ( num_blocks * g->blkSize) * sizeof(GVAL) + ( num_blocks) * sizeof(char  * )) ;
char  * pos = (char  * ) gv_o8param0->data_pointer.p2 + num_blocks * sizeof(char  * );
for(int b = 0; b < num_blocks ; b++) {
 gv_o8param0->data_pointer.p2[ b] = (GVAL  * ) pos ;
 pos += g->blkSize * sizeof(GVAL) ;
for(int c = 0; c < g->blkSize ; c++) {
 gv_o8param0->data_pointer.p2[ b][ c] = (GVAL) 0 ;
}
}
}
io_var_t io_gv_o8param0;
{
size_t min_block = g->mpi_rank == ( 0) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  0 % ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
size_t max_block = g->mpi_rank < ( 0) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > ( g->cBlkCnt - 1) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  0 : g->mpi_rank == ( g->cBlkCnt - 1) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  g->cBlkCnt % ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  g->cBlkCnt % ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
 #pragma omp parallel for
for(size_t block_index = ( min_block); block_index < ( max_block) ; block_index++) {
for(size_t cell_index = ( 0); cell_index < ( g->blkSize) ; cell_index++) {
 gv_o8param0->data_pointer.p2[ ( block_index)][ ( cell_index)] = 1.0f / 3.0f ;
}
}
}
 io_write_define( g , "gv_o8param0" , (GVAL  * ) gv_o8param0 , FLOAT32 , GRID_POS_CELL , GRID_DIM_2D , &io_gv_o8param0) ;
 io_write_announce( g , &io_gv_o8param0) ;
struct {
char  * name;
int loc;
int dim;
union {
GVAL  * restrict * restrict p2;
GVAL  * restrict * restrict * restrict p3;
} data_pointer;
}  * gv_o8param1;
{
int num_blocks = local_cell_blocks ?  local_cell_blocks : ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
 gv_o8param1 = malloc( 24) ;
 gv_o8param1->name = "gv_o8param1" ;
 gv_o8param1->loc = 0 ;
 gv_o8param1->dim = 2 ;
 gv_o8param1->data_pointer.p2 = malloc( ( num_blocks * g->blkSize) * sizeof(GVAL) + ( num_blocks) * sizeof(char  * )) ;
char  * pos = (char  * ) gv_o8param1->data_pointer.p2 + num_blocks * sizeof(char  * );
for(int b = 0; b < num_blocks ; b++) {
 gv_o8param1->data_pointer.p2[ b] = (GVAL  * ) pos ;
 pos += g->blkSize * sizeof(GVAL) ;
for(int c = 0; c < g->blkSize ; c++) {
 gv_o8param1->data_pointer.p2[ b][ c] = (GVAL) 0 ;
}
}
}
io_var_t io_gv_o8param1;
{
size_t min_block = g->mpi_rank == ( 0) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  0 % ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
size_t max_block = g->mpi_rank < ( 0) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > ( g->cBlkCnt - 1) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  0 : g->mpi_rank == ( g->cBlkCnt - 1) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  g->cBlkCnt % ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  g->cBlkCnt % ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
 #pragma omp parallel for
for(size_t block_index = ( min_block); block_index < ( max_block) ; block_index++) {
for(size_t cell_index = ( 0); cell_index < ( g->blkSize) ; cell_index++) {
 gv_o8param1->data_pointer.p2[ ( block_index)][ ( cell_index)] = 1.0f / 3.0f ;
}
}
}
 io_write_define( g , "gv_o8param1" , (GVAL  * ) gv_o8param1 , FLOAT32 , GRID_POS_CELL , GRID_DIM_2D , &io_gv_o8param1) ;
 io_write_announce( g , &io_gv_o8param1) ;
struct {
char  * name;
int loc;
int dim;
union {
GVAL  * restrict * restrict p2;
GVAL  * restrict * restrict * restrict p3;
} data_pointer;
}  * gv_o8param2;
{
int num_blocks = local_cell_blocks ?  local_cell_blocks : ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
 gv_o8param2 = malloc( 24) ;
 gv_o8param2->name = "gv_o8param2" ;
 gv_o8param2->loc = 0 ;
 gv_o8param2->dim = 2 ;
 gv_o8param2->data_pointer.p2 = malloc( ( num_blocks * g->blkSize) * sizeof(GVAL) + ( num_blocks) * sizeof(char  * )) ;
char  * pos = (char  * ) gv_o8param2->data_pointer.p2 + num_blocks * sizeof(char  * );
for(int b = 0; b < num_blocks ; b++) {
 gv_o8param2->data_pointer.p2[ b] = (GVAL  * ) pos ;
 pos += g->blkSize * sizeof(GVAL) ;
for(int c = 0; c < g->blkSize ; c++) {
 gv_o8param2->data_pointer.p2[ b][ c] = (GVAL) 0 ;
}
}
}
io_var_t io_gv_o8param2;
{
size_t min_block = g->mpi_rank == ( 0) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  0 % ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
size_t max_block = g->mpi_rank < ( 0) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > ( g->cBlkCnt - 1) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  0 : g->mpi_rank == ( g->cBlkCnt - 1) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  g->cBlkCnt % ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  g->cBlkCnt % ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
 #pragma omp parallel for
for(size_t block_index = ( min_block); block_index < ( max_block) ; block_index++) {
for(size_t cell_index = ( 0); cell_index < ( g->blkSize) ; cell_index++) {
 gv_o8param2->data_pointer.p2[ ( block_index)][ ( cell_index)] = 1.0f / 3.0f ;
}
}
}
 io_write_define( g , "gv_o8param2" , (GVAL  * ) gv_o8param2 , FLOAT32 , GRID_POS_CELL , GRID_DIM_2D , &io_gv_o8param2) ;
 io_write_announce( g , &io_gv_o8param2) ;
}
void Init_gv_o8par2(GRID  * g)
{
struct {
char  * name;
int loc;
int dim;
union {
GVAL  * restrict * restrict p2;
GVAL  * restrict * restrict * restrict p3;
} data_pointer;
}  * gv_o8par2;
{
int num_blocks = local_cell_blocks ?  local_cell_blocks : ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
 gv_o8par2 = malloc( 24) ;
 gv_o8par2->name = "gv_o8par2" ;
 gv_o8par2->loc = 0 ;
 gv_o8par2->dim = 3 ;
 gv_o8par2->data_pointer.p3 = malloc( ( num_blocks * g->height * g->blkSize) * sizeof(GVAL) + ( num_blocks * g->height) * sizeof(char  * ) + ( num_blocks) * sizeof(char  * )) ;
char  * pos = (char  * ) gv_o8par2->data_pointer.p3 + num_blocks * sizeof(char  * );
char  * pos2 = (char  * ) gv_o8par2->data_pointer.p3 + num_blocks * sizeof(char  * ) + num_blocks * g->height * sizeof(char  * );
for(int b = 0; b < num_blocks ; b++) {
 gv_o8par2->data_pointer.p3[ b] = (GVAL  * * ) pos ;
 pos += g->height * sizeof(char  * ) ;
for(int k = 0; k < g->height ; k++) {
 gv_o8par2->data_pointer.p3[ b][ k] = (GVAL  * ) pos2 ;
 pos2 += g->blkSize * sizeof(GVAL) ;
for(int c = 0; c < g->blkSize ; c++) {
 gv_o8par2->data_pointer.p3[ b][ k][ c] = (GVAL) 0 ;
}
}
}
}
io_var_t io_gv_o8par2;
{
size_t min_block = g->mpi_rank == ( 0) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  0 % ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : 0;
size_t max_block = g->mpi_rank < ( 0) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) || g->mpi_rank > ( g->cBlkCnt - 1) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  0 : g->mpi_rank == ( g->cBlkCnt - 1) / ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  g->cBlkCnt % ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) ?  g->cBlkCnt % ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size) : ( ( ( g->cBlkCnt) + g->mpi_world_size - 1) / g->mpi_world_size);
 #pragma omp parallel for
for(size_t block_index = ( min_block); block_index < ( max_block) ; block_index++) {
for(size_t height_index = ( 0); height_index < ( g->height) ; height_index++) {
for(size_t cell_index = ( 0); cell_index < ( g->blkSize) ; cell_index++) {
 gv_o8par2->data_pointer.p3[ ( block_index)][ ( height_index)][ ( cell_index)] = 1.0f / 2.0f ;
}
}
}
}
 io_write_define( g , "gv_o8par2" , (GVAL  * ) gv_o8par2 , FLOAT32 , GRID_POS_CELL , GRID_DIM_3D , &io_gv_o8par2) ;
 io_write_announce( g , &io_gv_o8par2) ;
}
int main(int argc , char  * * argv)
{
int PrintHelp = 0;
 parseOptions( argc , argv , options , &PrintHelp) ;
if (  PrintHelp ) {
 print_help( options , 0) ;
 exit( 0) ;
}
GRID  * g = malloc( sizeof(GRID));
{
 MPI_Init( NULL , NULL) ;
 MPI_Comm_size( MPI_COMM_WORLD , &g->mpi_world_size) ;
 MPI_Comm_rank( MPI_COMM_WORLD , &g->mpi_rank) ;
}
 init_grid( g , gridsize , gridheight) ;
 io_write_init( g , Filename) ;
 Init_gv_temp( g) ;
 Init_gv_ind2Dparam( g) ;
 Init_gv_o8param( g) ;
 Init_gv_o8par2( g) ;
 io_write_registration_complete( g) ;
 io_write_start( g) ;
 io_write_finalize( g) ;
{
 MPI_Finalize() ;
}
}