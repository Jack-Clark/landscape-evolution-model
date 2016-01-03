
#include "lem.h"
#include "Data.h"
#include "watershed.h"
#include "flowroutines.h"

#ifndef floodingdriverH
#define floodingdriverH

void outputEdges(const char* filename, int table_size, int* table_counter, int* hash_table, int list_size);

void outputEdges(const char* filename, int table_size, int* table_counter, double* hash_table, int list_size);

void floodingDriver(dim3 dimGrid, dim3 dimBlock, Data* data, Data* device, Catchment* catchments, int ncell_x, int ncell_y, int cell_size, int iter);



#endif
