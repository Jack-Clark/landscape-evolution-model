#include "contribA.h"

__global__ void getcontribarea(int *mask, int *fd, int *contribA, int *ok, int *secondOK, int *progressed, int rows, int cols)
{
  int self;
  int irow = blockIdx.y * blockDim.y + threadIdx.y;
  int icol = blockIdx.x * blockDim.x + threadIdx.x;

  if (irow >= rows || icol >= cols) // || irow < 0 || icol < 0)
    return;

  self = irow * cols + icol;

  if (mask[self] == 0) return; // don't calculate if not in catchment(s) of interest

  int nie, nise, nis, nisw, niw, ninw, nin, nine;
  
  if (ok[self] != 0) return;

  if (secondOK[self] == 1) {
    ok[self] = 1;
    *progressed = 1;
    return;
  }
  int accum = 1;
  
  nie  = self        + 1 ;
  nise = self + cols + 1 ;
  nis  = self + cols     ;
  nisw = self + cols - 1 ;
  niw  = self        - 1 ;
  ninw = self - cols - 1 ;
  nin  = self - cols     ;
  nine = self - cols + 1 ;

   
    if (icol < cols - 1 && fd[nie] & WEST) {
        if (!ok[nie]) return;
        
        accum += contribA[nie];
    }
    if (icol < cols - 1 && irow < rows - 1 && fd[nise] & NORTHWEST) {
        if (!ok[nise]) return;
        
        accum += contribA[nise];
    }
    if (irow < rows - 1 && fd[nis] & NORTH) {
        if (!ok[nis]) return;
        
        accum += contribA[nis];
    }
    if (icol > 0 && irow < rows - 1 && fd[nisw] & NORTHEAST) {
        if (!ok[nisw]) return;
        
        accum += contribA[nisw];
    }
    if (icol > 0 && fd[niw] & EAST) {
        if (!ok[niw]) return;
        
        accum += contribA[niw];
    }
    if (icol > 0 && irow > 0 && fd[ninw] & SOUTHEAST) {
        if (!ok[ninw]) return;
        
        accum += contribA[ninw];
    }
    if (irow > 0 && fd[nin] & SOUTH) {
        if (!ok[nin]) return;
        
        accum += contribA[nin];
    }
    if (irow > 0 && icol < cols - 1 && fd[nine] & SOUTHWEST) {
        if (!ok[nine]) return;
        
        accum += contribA[nine];
    }
    contribA[self] = accum;
    secondOK[self] = 1;
    *progressed = 1;
}

int calccontribA(Data* data, Data* device, Catchment* catchments)

{

  int block_ncell_y = 16;
  int block_ncell_x = 16;
  int ncell_x = data->mapInfo.width;
  int ncell_y = data->mapInfo.height;

  dim3 dimGrid(ncell_x/block_ncell_x + 1, ncell_y/block_ncell_y + 1);
  dim3 dimBlock(block_ncell_x, block_ncell_y);

  int	fullsize= data->mapInfo.width * data->mapInfo.height;
  //data->contribA = (int*) malloc(sizeof(int) * data->mapInfo.width* data->mapInfo.height);
		  for (int i = 0; i < data->mapInfo.width * data->mapInfo.height; i++) {
			  data->contribA[i] = 0;
		  }
    
  cudaMalloc((void **) &device->contribA, data->mapInfo.width* data->mapInfo.height * sizeof(int));
  /*printf("contribA:%s\n", cudaGetErrorString(cudaGetLastError()));

  int cpucontribA_max = 0;
	for (int i = 0; i < ncell_y; i++) {
		for (int j = 0; j < ncell_x; j++) {
			if (data->contribA[i * ncell_x + j] > cpucontribA_max)
				{
				cpucontribA_max = data->contribA[i * ncell_x + j];

				}
		}
	}
	printf("max contribA before: %d\n", cpucontribA_max);
*/
  cudaMemcpy(device->contribA, data->contribA, data->mapInfo.width * data->mapInfo.height * sizeof(int), cudaMemcpyHostToDevice);
  //printf("copy empty contribA:%s\n", cudaGetErrorString(cudaGetLastError()));



  int *ok;
  cudaMalloc((void **) &ok, data->mapInfo.width* data->mapInfo.height * sizeof(int));
  //printf("ok:%s\n", cudaGetErrorString(cudaGetLastError()));
  
  int *secondOK;
  cudaMalloc((void **) &secondOK, data->mapInfo.width * data->mapInfo.height * sizeof(int));
 // printf("secondOK:%s\n", cudaGetErrorString(cudaGetLastError()));
  
  int *okGrid_h = (int*) malloc(sizeof(int) * data->mapInfo.width * data->mapInfo.height);
  //printf(":%s\n", cudaGetErrorString(cudaGetLastError()));
  
  for (int i = 0; i < data->mapInfo.width * data->mapInfo.height; i++) {
    okGrid_h[i] = 0;
  }
  
  int *progressed;
  cudaMalloc((void **) &progressed, sizeof(int));
  //printf("progressed:%s\n", cudaGetErrorString(cudaGetLastError()));

  cudaMemcpy(ok, okGrid_h, data->mapInfo.width * data->mapInfo.height * sizeof(int),   cudaMemcpyHostToDevice); // set ok grid to zero
  cudaMemcpy(secondOK, okGrid_h, data->mapInfo.width * data->mapInfo.height * sizeof(int), cudaMemcpyHostToDevice); // set secondOK grid to zero

  int consecutiveZero = 0;
  int *progress = (int*) malloc(sizeof(int));
  *progress = 0;

  do {
	   *progress = 0;
  	   cudaMemcpy(progressed, progress, sizeof(int), cudaMemcpyHostToDevice);   
	
	   getcontribarea<<<dimGrid, dimBlock>>>(device->mask, device->fd, device->contribA, ok, secondOK, progressed, data->mapInfo.height, data->mapInfo.width); // the accumulation algorithm without partitioning or weights
	   
	   cudaMemcpy(progress, progressed, sizeof(int), cudaMemcpyDeviceToHost);

	   if (*progress == 0) {
			//printf("Progress was not made \n");
			consecutiveZero ++;
		}
	  else
			consecutiveZero = 0;

	} while (consecutiveZero < 5);

  cudaMemcpy(data->contribA, device->contribA, data->mapInfo.width * data->mapInfo.height * sizeof(int), cudaMemcpyDeviceToHost); // copy back contributing area grid in no of cells

  data->cpucontribA_max = 0;
	for (int i = 0; i < ncell_y; i++) {
		for (int j = 0; j < ncell_x; j++) {
			if (data->contribA[i * ncell_x + j] > data->cpucontribA_max)
				{
				data->cpucontribA_max = data->contribA[i * ncell_x + j];

				}
		}
	}
	printf("max contribA: %d\n", data->cpucontribA_max);



//  cudaFree(device->contribA);  // do not free as we will use this.
  cudaFree(ok);
  cudaFree(secondOK);
  cudaFree(progressed);

  free(okGrid_h);
  free(progress);
 // free(data->shortest_paths);

  return 0;
}
