#ifndef flowroutinesH
#define flowroutinesH

#include "lem.h"

typedef int TYPE;

__constant__ const TYPE code[] = {64 , 128, 1, 2, 4, 8, 16, 32, 0};


__global__ void MultipleFlowDir(int *mask, double *zs, double *slopes, int *aspects, int cell_size, int ncell_x, int ncell_y, int *dx, int *dy, int last);

__global__ void singleFlowDir(int *mask, double *zs, double *slopes, int *aspects, int cell_size,	int ncell_x, int ncell_y, int *dx, int *dy, int last);

__global__ void route_plateaus(int *mask, double *zs, int *aspects, float* shortest_paths,	int ncell_x, int ncell_y, int *dx, int *dy, int *change_flag, double *lowHeight);

__global__ void slope_plateaus(int *mask, double *zs, float* shortest_paths, int ncell_x, int ncell_y, double *lowHeight, double *slopes, int cell_size);

__global__ void remap(int *mask, int *aspects, int ncell_x, int ncell_y);

__global__ void comp_shortest_paths_sink(int *mask, double *zs, int *aspects, float* shortest_paths, int ncell_x, int ncell_y, int *dx, int *dy, int *change_flag);

__global__ void shortest_paths_plateaus_init(int *mask, double *zs, double *slopes, int *aspects, float* shortest_paths, int cell_size, int ncell_x, int ncell_y, double *lowHeight);

__global__ void flow_boundary(int *mask, double *zs, double *slopes, int *aspects, int ncell_x, int ncell_y, int *dx, int *dy);



#endif
