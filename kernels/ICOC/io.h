#ifndef TEST_IO_H

#define TEST_IO_H
#include "memory.h"
typedef enum { FLOAT16, FLOAT32, FLOAT64 } datatype_t;
/*
 * Contains all the information about the variable
 */
typedef struct {
    int id;
// an module specific ID, this can later be a struct or sth...
    void *memory;
    datatype_t memtype;
// Do we need these?
    datatype_t storagetype;
// TODO: can the DSL do this book keeping? 
    GRID_VAR_POSITION Pos;
    GRID_VAR_DIMENSION Dim;
    const char *Name;
} io_var_t;
#define ArrayCount(array) sizeof((array))/sizeof((array[0]))
// TODO: maybe use a dynamic data-structure
typedef struct {
    int Head;
    int Tail;
    int isLooped;
    io_var_t Elements[32];
} write_queue;
/*
 * The model will call read_init exactly once before all init() methods of the modules are executed
 */
void io_read_init(GRID * g, const char *filename);
/**
 * Reads the variable from the open file...
 * TODO: in practise the register_function will actually block and read the variable
 */
void io_read_register(GRID * g, const char *variable_name, GVAL * out_buff, datatype_t file_type, datatype_t memory_type, GRID_VAR_POSITION var_pos, GRID_VAR_DIMENSION var_dim);
/*
 * Shall start all asynchronous I/O... TODO in practice this is a NULL operation, because the read_register has done the reading already
 * The model will call it exactly once after all init() functions have been executed.
 */
void io_read_start();
void io_write_init(GRID * g, const char *filename);
/**
 * Register the variable to write to, this needs to be done *once* for the whole model
 */
void io_write_define(GRID * g, const char *variable_name, GVAL * out_buff, datatype_t file_type, GRID_VAR_POSITION var_pos, GRID_VAR_DIMENSION var_dim, io_var_t * out_varid);
void io_write_registration_complete(GRID * g);
/**
 * prepares and enqueues variables for the io phase
 */
void io_write_announce(GRID * g, io_var_t * var);
/*
 * Writes out all the enqueued variables
 * TODO: maybe a more fitting name? 
 */
void io_write_start(GRID * g);
void io_write_finalize(GRID * g);
/*
 * May be need a finalize function to clean and finish I/O related activities
 */
#endif
