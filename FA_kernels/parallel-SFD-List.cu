#include "FA_SFD.h"


__global__ void kernelfunction_SFD_NoPart_List(int *mask, int *fd, double *fa, int rows, int cols, double *weights, unsigned int* pkgprogressd, unsigned int* left)
{
	//int self = blockIdx.y*cols*BLOCKROWS + blockIdx.x*BLOCKCOLS + threadIdx.y*cols + threadIdx.x;
	int irow = blockIdx.y * blockDim.y + threadIdx.y;
	int icol = blockIdx.x * blockDim.x + threadIdx.x;
	int maxSize = rows * cols;

	if (irow >= rows || icol >= cols)
		return;

	int self = irow * cols + icol;
	if (mask[self] == 0) return; // don't calculate if not in catchment(s) of interest
	int nie, nise, nis, nisw, niw, ninw, nin, nine;

	double accum = 1.0 * weights[self];

	nie  = self        + 1 ;
	nise = self + cols + 1 ;
	nis  = self + cols     ;
	nisw = self + cols - 1 ;
	niw  = self        - 1 ;
	ninw = self - cols - 1 ;
	nin  = self - cols     ;
	nine = self - cols + 1 ;

	unsigned int myPos = 0;

	if (icol < cols - 1 && fd[nie] & WEST) {
		if (fa[nie] < 0) {
			myPos = atomicInc(pkgprogressd, maxSize);
		    left[myPos] = self;
			return;
		}
		accum += fa[nie];
	}
	if (icol < cols - 1 && irow < rows - 1 && fd[nise] & NORTHWEST) {
		if (fa[nise] < 0)  {
			myPos = atomicInc(pkgprogressd, maxSize);
		    left[myPos] = self;
			return;
		}

		accum += fa[nise];
	}
	if (irow < rows - 1 && fd[nis] & NORTH) {
		if (fa[nis] < 0)  {
			myPos = atomicInc(pkgprogressd, maxSize);
		    left[myPos] = self;
			return;
		}
		accum += fa[nis];
	}
	if (icol > 0 && irow < rows - 1 && fd[nisw] & NORTHEAST) {
		if (fa[nisw] < 0)  {
			myPos = atomicInc(pkgprogressd, maxSize);
		    left[myPos] = self;
			return;
		}

		accum += fa[nisw];
	}
	if (icol > 0 && fd[niw] & EAST) {
		if (fa[niw] < 0)  {
			myPos = atomicInc(pkgprogressd, maxSize);
		    left[myPos] = self;
			return;
		}

		accum += fa[niw];
	}
	if (icol > 0 && irow > 0 && fd[ninw] & SOUTHEAST) {
		if (fa[ninw] < 0)  {
			myPos = atomicInc(pkgprogressd, maxSize);
		    left[myPos] = self;
			return;
		}

		accum += fa[ninw];
	}
	if (irow > 0 && fd[nin] & SOUTH) {
		if (fa[nin] < 0)  {
			myPos = atomicInc(pkgprogressd, maxSize);
		    left[myPos] = self;
			return;
		}

		accum += fa[nin];
	}
	if (irow > 0 && icol < cols - 1 && fd[nine] & SOUTHWEST) {
		if (fa[nine] < 0)  {
			myPos = atomicInc(pkgprogressd, maxSize);
		    left[myPos] = self;
			return;
		}

		accum += fa[nine];
	}
	fa[self] = accum;
}

__global__ void kernelfunction_SFD_NoPart_ListProgress(int *mask, int *fd, double *fa, int rows, int cols, double *weights, unsigned int* pkgprogressd, unsigned int* left, unsigned int* had, int hadsize)
{
	int pos = blockIdx.y*cols*blockDim.y + blockIdx.x*blockDim.x + threadIdx.y*cols + threadIdx.x;
	int maxSize = rows * cols;

	int nie, nise, nis, nisw, niw, ninw, nin, nine;

	if (pos >= hadsize)
		return;

	int self = had[pos];
//	if (mask[self] == 0) return; // don't calculate if not in catchment(s) of interest

	double accum = 1.0 * weights[self];

	nie  = self        + 1 ;
	nise = self + cols + 1 ;
	nis  = self + cols     ;
	nisw = self + cols - 1 ;
	niw  = self        - 1 ;
	ninw = self - cols - 1 ;
	nin  = self - cols     ;
	nine = self - cols + 1 ;

	int myPos = 0;

	int icol = self % cols;
	int irow = self / cols;

	if (icol < cols - 1 && fd[nie] & WEST) {
		if (fa[nie] < 0) {
			myPos = atomicInc(pkgprogressd, maxSize);
		    left[myPos] = self;
			return;
		}
		accum += fa[nie];
	}
	if (icol < cols - 1 && irow < rows - 1 && fd[nise] & NORTHWEST) {
		if (fa[nise] < 0)  {
			myPos = atomicInc(pkgprogressd, maxSize);
		    left[myPos] = self;
			return;
		}

		accum += fa[nise];
	}
	if (irow < rows - 1 && fd[nis] & NORTH) {
		if (fa[nis] < 0)  {
			myPos = atomicInc(pkgprogressd, maxSize);
		    left[myPos] = self;
			return;
		}
		accum += fa[nis];
	}
	if (icol > 0 && irow < rows - 1 && fd[nisw] & NORTHEAST) {
		if (fa[nisw] < 0)  {
			myPos = atomicInc(pkgprogressd, maxSize);
		    left[myPos] = self;
			return;
		}

		accum += fa[nisw];
	}
	if (icol > 0 && fd[niw] & EAST) {
		if (fa[niw] < 0)  {
			myPos = atomicInc(pkgprogressd, maxSize);
		    left[myPos] = self;
			return;
		}

		accum += fa[niw];
	}
	if (icol > 0 && irow > 0 && fd[ninw] & SOUTHEAST) {
		if (fa[ninw] < 0)  {
			myPos = atomicInc(pkgprogressd, maxSize);
		    left[myPos] = self;
			return;
		}

		accum += fa[ninw];
	}
	if (irow > 0 && fd[nin] & SOUTH) {
		if (fa[nin] < 0)  {
			myPos = atomicInc(pkgprogressd, maxSize);
		    left[myPos] = self;
			return;
		}

		accum += fa[nin];
	}
	if (irow > 0 && icol < cols - 1 && fd[nine] & SOUTHWEST) {
		if (fa[nine] < 0)  {
			myPos = atomicInc(pkgprogressd, maxSize);
		    left[myPos] = self;
			return;
		}

		accum += fa[nine];
	}
	fa[self] = accum;

}


int process_SFD_NoPart_List(Data* data, Data* device, int iter) {
	printf("In process\n");

    int gridRows = data->mapInfo.height;
	int gridColumns = data->mapInfo.width;
	int grid1 = gridColumns / (BLOCKCOLS );
	int grid2 = gridRows / (BLOCKROWS );
	int fullsize = gridRows * gridColumns;

	unsigned int *progress_d;
	checkCudaErrors(cudaMalloc((void **) &progress_d, sizeof(unsigned int)) );
	unsigned int *progress_h = (unsigned int*) malloc(sizeof(unsigned int));
	*progress_h = 0;

	checkCudaErrors(cudaMemcpy(progress_d, progress_h, sizeof(unsigned int), cudaMemcpyHostToDevice));

	// Pull straight back
//	checkCudaErrors(cudaMemcpy(progress_h, progress_d. sizeof(unsigned int), cudaMemcpyDeviceToHost));


	unsigned int* list1;
	unsigned int* list2;
	unsigned int* listT;
	checkCudaErrors(cudaMalloc((void **) &list1, fullsize * sizeof(unsigned int)) ); //FIXME - need to figure out a way to set this value
	checkCudaErrors(cudaMalloc((void **) &list2, fullsize * sizeof(unsigned int)) ); //FIXME - need to figure out a way to set this value
	//checkCudaErrors(cudaMalloc((void **) &listT, sizeof(unsigned int)) );

	//printf("ListT = %x\n", listT);

	dim3 dimGrid(grid1, grid2);
	dim3 dimBlock(BLOCKCOLS, BLOCKROWS);

	//printf(":%s\n", cudaGetErrorString(cudaGetLastError()));


	// first run a kernel to solve those cells which are on the edges and produce the first list

	//__global__ void kernelfunction_SFD_NoPart_List(int *mask, int *fd, double *fa, int rows, int cols, double *weights, unsigned int* pkgprogressd, int* left)

	kernelfunction_SFD_NoPart_List<<<dimGrid, dimBlock>>>(device->mask, device->fd, device->fa, gridRows, gridColumns, device->runoffweight, progress_d, list1);

	//printf(":First kernel function %s\n", cudaGetErrorString(cudaGetLastError()));


	// get the size of the array
	checkCudaErrors(cudaMemcpy(progress_h, progress_d, sizeof(unsigned int), cudaMemcpyDeviceToHost));

	unsigned int lastTot = gridRows * gridColumns;

	//printf("starting loop left=%d\n", *progress_h);
	//printf("max left = %d\n", lastTot);

	unsigned int* temp = (unsigned int*) malloc(sizeof(unsigned int));

	int nextBlocks;

	while (*progress_h > 0) { // while the array still has elements
		// Swap list 1 and list 2
		listT = list1;
		list1 = list2;
		list2 = listT;
		//printf("Cells left to process = %d\n", *progress_h);
		if (*progress_h > lastTot) {
			printf("The number of incorrect cells should be coming down!\n");
			scanf("%d", &lastTot);
		}
		lastTot = *progress_h;
		// call kernel to process the list
		nextBlocks = *progress_h / 128 + 1;

		// reset the value of progressed before restarting

		*temp = 0;
		checkCudaErrors(cudaMemcpy(progress_d, temp, sizeof(unsigned int), cudaMemcpyHostToDevice) );
		//__global__ void kernelfunction_SFD_NoPart_ListProgress(int *fd, double *fa, int rows, int cols, double *weights, unsigned int* pkgprogressd, int* left, int* had, int hadsize)
		kernelfunction_SFD_NoPart_ListProgress<<<nextBlocks, 128>>>(device->mask, device->fd, device->fa, gridRows, gridColumns, device->runoffweight, progress_d, list1, list2, *progress_h);

		//cudaDeviceSynchronize();
		//printf(":loop kernel function %s\n", cudaGetErrorString(cudaGetLastError()));

		// get the new cell count
		checkCudaErrors(cudaMemcpy(progress_h, progress_d, sizeof(unsigned int), cudaMemcpyDeviceToHost) );
		//printf("Left= %d\n", *progress_h);
	}

	free(temp);

	// Copy flow accumulation back
	int count = 0;

	checkCudaErrors(cudaMemcpy(data->fa, device->fa, fullsize * sizeof(double),   cudaMemcpyDeviceToHost));
	fprintf(data->outlog, "FA: FA memcopy back operation 3:%s\n", cudaGetErrorString(cudaGetLastError()));

	double FA_max;
	int FAindex = 0;
	double cpuFA_max = 0.0;

	if (iter == 1) // cpu calculation otherwise we cannot locate the outletcell index
	{
		for (int i = 0; i < gridRows; i++) {
			for (int j = 0; j < gridColumns; j++) {
				if (data->fa[i * gridColumns + j] > cpuFA_max)
					{
					cpuFA_max = data->fa[i * gridColumns + j];
					FAindex = i * gridColumns + j;
					}
			}
		}
		data->FA_max = cpuFA_max;
		data->outletcellidx = FAindex; // this is the outlet cell which will be maintained throughout the simulation
	} else // do it faster using GPU in all subsequent iterations
		{
			thrust::device_ptr<double> max_FA = thrust::device_pointer_cast(device->fa);
			FA_max = thrust::reduce(max_FA, max_FA + fullsize, (double) 0, thrust::maximum<double>());
			data->FA_max = FA_max;
		}

	fprintf(data->outlog, "FA: Maximum FA is  %.6f s\n\n", data->FA_max);
	fprintf(data->outlog, "FA: Outletcell index is  %d s\n\n", data->outletcellidx);

	printf("Maximum FA is  %.6f s\n\n", data->FA_max);

	for (int i = 0; i < gridRows; i++) {
		for (int j = 0; j < gridColumns; j++) {
				if (data->fa[i * gridColumns + j] < 0) {
				//if (count < 20)
					//printf("[%d,%d] = %8.7f\n", i, j, data->fa[i * ncell_x + j]);
				count++;
				}
		}
	}
	fprintf(data->outlog, "FA: Bad value count (i.e. not in catchment(s) = %d\n", count);

	cudaFree(list1);
	cudaFree(list2);
	//cudaFree(listT);
	cudaFree(progress_d);

	free(progress_h);


	return 1;
}

void correctflow_SFD_NoPart_List(Data* data, Data* device, int iter) {
	if (cudaSuccess != cudaSetDevice(CUDA_DEVICE)) {
		printf("Unable to access CUDA card\n");
		return ;
	}

	int x;

	//printf("nrows = %d ncols = %d\n", nrows, ncols);
	cudaEvent_t start, stop;
	float time;
    int gridRows = data->mapInfo.height;
	int gridColumns = data->mapInfo.width;
	int fullsize = gridRows * gridColumns;

	// Set all values to 0.0
	for (x = 0; x < gridRows * gridColumns; ++x) {
		data->fa[x] = -1.0;
	}

	fprintf(data->outlog, "FA: set fagrid values to -1\n");
	cudaMemcpy(device->fa, data->fa, fullsize * sizeof(double), cudaMemcpyHostToDevice);
	fprintf(data->outlog, "FA: FA memcopy operation :%s\n", cudaGetErrorString(cudaGetLastError()));

	fprintf(data->outlog, "FA: Calling process\n");

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	process_SFD_NoPart_List(data, device, iter);

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	printf("Time to complete FA_SFD_list : %.6f s\n", time / 1000.0);
}
