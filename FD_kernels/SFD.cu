#include "SFD.h"

//******************************************************
//**
//** Main Code for SFD calculation
//** 
//****
//** This code sets up and calls all kernels
//**
//******************************************************
void cuFlowDirection(Data* data, Data* device, Catchment* catchments, int iter)

{
  cudaEvent_t start, stop;
  float time;
  size_t freenow, total;

 // start the timer for SFD;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  int xmove [9] = {0,1,1,1,0,-1,-1,-1,0};
  int ymove [9] = {-1,-1,0,1,1,1,0,-1,0};
  int *dx;
  int *dy;
  dx = &xmove[0];
  dy = &ymove[0];

  //Data device;
  data->dx = dx;
  data->dy = dy;

  int block_ncell_y = 16;
  int block_ncell_x = 16;
  int fullsize;
  int ncell_x = data->mapInfo.width;
  int ncell_y = data->mapInfo.height;
  int csize = (int) data->mapInfo.cellsize;

  //dim3 dimGrid(ncell_x/block_ncell_x + 1, ncell_y/block_ncell_y + 1);

  dim3 dimGrid(ncell_x/block_ncell_x, ncell_y/block_ncell_y);
  dim3 dimBlock(block_ncell_x, block_ncell_y);
  fullsize= ncell_x * ncell_y;

  // declare transient memory to be freed on exit
  data->shortest_paths = (float*) malloc(fullsize * sizeof(float));
  data->fdmod = (int *) malloc(sizeof(int) * data->mapInfo.height * data->mapInfo.width);
  // *** Set up memory space to hold watershed id's
  data->watershed_id = (int*) malloc(ncell_x * ncell_y * sizeof(int));
  fprintf(data->outlog, "FD: dev mem allocate shortest_path, fdmod, watershed_id :%s\n", cudaGetErrorString(cudaGetLastError()));

  cudaMemcpy((void *) device->dx, data->dx, 9 * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy((void *) device->dy, data->dy, 9 * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy((void *) device->dem, data->dem, fullsize * sizeof(double), cudaMemcpyHostToDevice);
  fprintf(data->outlog, "FD: memcopy dx, dy, dem :%s\n", cudaGetErrorString(cudaGetLastError()));

  //********************************************************************
  //** Work out initial flow directions
  //********************************************************************
  // Each thread takes one cell and asigns a flow direction to it based
  // on its lowest neighbour
  //
  // NW  N NE    7  0  1
  //  W  *  E    6  8  2
  // SW  S SE    5  4  3

  singleFlowDir<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->SlopePtr, device->fd, csize, ncell_x, ncell_y, device->dx, device->dy, 0);
  fprintf(data->outlog,"FD: first single flow :%s\n", cudaGetErrorString(cudaGetLastError()));
  //********************************************************************
  //** Set flow direction for edge cells
  //********************************************************************
  // Mark all edge cells as flowing out from the DEM : THIS IS NOT NEEDED WHEN THE MASK IS INSIDE THE GRID BOUNDARIES


 // flow_boundary<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->SlopePtr, device->fd, data->mapInfo.width, data->mapInfo.height, device->dx, device->dy);

  // non-flat cells (aspect != 8) have shortest path = 0 all others have 
  // shortest path = large number
  // All values of lowHeight are set to 0

  shortest_paths_plateaus_init<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->SlopePtr, device->fd, device->shortest_paths, csize /*data->mapInfo.cellsize*/, data->mapInfo.width, data->mapInfo.height, device->lowHeight);
  fprintf(data->outlog,"FD: first shortest_path_plateaus_init :%s\n", cudaGetErrorString(cudaGetLastError()));
// Set up a flag to indicate if we have had any change since the last iteration?

  int *change_flag_h, *change_flag_d;
  change_flag_h = (int *)malloc(sizeof(int));
  cudaMalloc((void **)&change_flag_d, sizeof(int));
  int count = 0;

  do
  {
    // Clear flag
    *change_flag_h = 0;
    cudaMemcpy(change_flag_d, change_flag_h, sizeof(int), cudaMemcpyHostToDevice);

    // Each cell looks at the cells around itself finds the one with the same 
    // height as itself and shortest path to exit and flows towards that cell
    // When we adopt the shortest route we adopt its lowHeight value too
    route_plateaus<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->fd, device->shortest_paths,  data->mapInfo.width, data->mapInfo.height, device->dx, device->dy, change_flag_d, device->lowHeight);
    if (count <3) fprintf(data->outlog,"FD: first route_plateaus_init :%s\n", cudaGetErrorString(cudaGetLastError()));
	cudaMemcpy(change_flag_h, change_flag_d, sizeof(int), cudaMemcpyDeviceToHost);  // copy back flag

    //printf("Any change in plateaus? %d\n", *change_flag_h);
    count ++;

    //if(*change_flag_h == 1)
  } while(*change_flag_h == 1); // while at least some cells have been modified, repeat

  fprintf(data->outlog,"FD: route plateaus :%s\n", cudaGetErrorString(cudaGetLastError()));
  fprintf(data->outlog,"FD: Routed plateaus in %d iterations. Now starting flooding driver... \n", count);


  	floodingDriver(dimGrid, dimBlock, data, device, catchments, data->mapInfo.width, data->mapInfo.height, csize /*data->mapInfo.cellsize*/, iter);
  	fprintf(data->outlog,"FD: errors after return from flooding driver :%s\n", cudaGetErrorString(cudaGetLastError()));


  	cudaMemcpy(data->SlopePtr, device->SlopePtr, data->mapInfo.height * data->mapInfo.width * sizeof(double), cudaMemcpyDeviceToHost);

	// now remap the aspects for the Accumulation routines when using SFD only
	remap<<<dimGrid, dimBlock>>>(device->mask, device->fd,data->mapInfo.width, data->mapInfo.height )  ;
	fprintf(data->outlog,"FD: errors after remap :%s\n", cudaGetErrorString(cudaGetLastError()));
	cudaMemcpy(data->fd, device->fd, data->mapInfo.height * data->mapInfo.width * sizeof(int), cudaMemcpyDeviceToHost);

	//		MultipleFlowDir<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->SlopePtr, device->fd, csize, ncell_x, ncell_y, device->dx, device->dy, 0);
	//		cudaMemcpy(data->fd, device->fd, data->mapInfo.height * data->mapInfo.width * sizeof(int), cudaMemcpyDeviceToHost);

  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);
  cudaEventDestroy(start);	
  cudaEventDestroy(stop);

  cudaFree(change_flag_d);

  free(data->shortest_paths);
  free(data->fdmod);
  free(change_flag_h);


  printf("Time to complete flow routing routine algorithm %.6f s\n", time / 1000.0);
  fprintf(data->outlog,"FD: time to complete flow routing routine algorithm %.6f s\n", time / 1000.0);

  cudaMemGetInfo(&freenow, &total);
  fprintf(data->outlog,"FD: Memory on CUDA card free at end of routing: %d total: %d\n",freenow/1024,total/1024);
}
