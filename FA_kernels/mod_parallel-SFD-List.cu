#include "FA_SFD.h"
#include <assert.h>

__global__ void kernelfunction_SFD_Initital_Compute_Deps_And_Calculate_Zero_Dep_Cells(int *mask, int *fd, double *fa, int rows, int cols, double *weights, unsigned int* pkgprogressd, int *dep) {

	// The number of neighbour cells flowing into this cell
	int depCount = 0;

	// number of neighbours that have had their FA calculated
 	int numNeighboursReady = 0;

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

	if (icol < cols - 1 && fd[nie] & WEST) {
		depCount++;
		if (fa[nie] >= 0) {
			numNeighboursReady++;
			accum += fa[nie];
		}
	}
	if (icol < cols - 1 && fd[nise] & NORTHWEST) {
		depCount++;
		if (fa[nise] >= 0) {
			numNeighboursReady++;
			accum += fa[nise];
		} 
	}
	if (icol < cols - 1 && fd[nis] & NORTH) {
		depCount++;
		if (fa[nis] >= 0) {
			numNeighboursReady++;
			accum += fa[nis];
		}
	}
	if (icol < cols - 1 && fd[nisw] & NORTHEAST) {
		depCount++;
		if (fa[nisw] >= 0) {
			numNeighboursReady++;
			accum += fa[nisw];
		}
	}
	if (icol < cols - 1 && fd[niw] & EAST) {
		depCount++;
		if (fa[niw] >= 0) {
			numNeighboursReady++;
			accum += fa[niw];
		}
	}
	if (icol < cols - 1 && fd[ninw] & SOUTHEAST) {
		depCount++;
		if (fa[ninw] >= 0) {
			numNeighboursReady++;
			accum += fa[ninw];
		} 
	}
	if (icol < cols - 1 && fd[nin] & SOUTH) {
		depCount++;
		if (fa[nis] >= 0) {
			numNeighboursReady++;
			accum += fa[nis];
		}
	}
	if (icol < cols - 1 && fd[nine] & SOUTHWEST) {
		depCount++;
		if (fa[nine] >= 0) {
			numNeighboursReady++;
			accum += fa[nine];
		}
	}

	if((depCount - numNeighboursReady) == 0) {
		fa[self] = accum;
	} else {
		atomicInc(pkgprogressd, maxSize);
	}
	dep[self] = depCount;
}


__global__ void kernelfunction_SFD_Calculate_Zero_Dependency_Cells(int *mask, int *fd, double *fa, int rows, int cols, double *weights, unsigned int* pkgprogressd) {

	int irow = blockIdx.y * blockDim.y + threadIdx.y;
	int icol = blockIdx.x * blockDim.x + threadIdx.x;
	int maxSize = rows * cols;

	if (irow >= rows || icol >= cols)
		return;

	int self = irow * cols + icol;
	if (mask[self] == 0) return; // don't calculate if not in catchment(s) of interest

	if(fa[self] >= 0) return;

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

	if (icol < cols - 1 && fd[nie] & WEST) {
		if (fa[nie] < 0) {
			atomicInc(pkgprogressd, maxSize);
			return;
		}
		accum += fa[nie];
	}
	if (icol < cols - 1 && fd[nise] & NORTHWEST) {
		if (fa[nise] < 0) {
			atomicInc(pkgprogressd, maxSize);
			return;
		}
		accum += fa[nise];
	}
	if (icol < cols - 1 && fd[nis] & NORTH) {
		if (fa[nis] < 0) {
			atomicInc(pkgprogressd, maxSize);
			return;
		}
		accum += fa[nis];
	}
	if (icol < cols - 1 && fd[nisw] & NORTHEAST) {
		if (fa[nisw] < 0) {
			atomicInc(pkgprogressd, maxSize);
			return;
		}
		accum += fa[nisw];
	}
	if (icol < cols - 1 && fd[niw] & EAST) {
		if (fa[niw] < 0) {
			atomicInc(pkgprogressd, maxSize);
			return;
		}
		accum += fa[niw];
	}
	if (icol < cols - 1 && fd[ninw] & SOUTHEAST) {
		if (fa[ninw] < 0) {
			atomicInc(pkgprogressd, maxSize);
			return;
		}
		accum += fa[ninw];
	}
	if (icol < cols - 1 && fd[nin] & SOUTH) {
		if (fa[nis] < 0) {
			atomicInc(pkgprogressd, maxSize);
			return;
		}
		accum += fa[nis];
	}
	if (icol < cols - 1 && fd[nine] & SOUTHWEST) {
		if (fa[nine] < 0) {
			atomicInc(pkgprogressd, maxSize);
			return;
		}
		accum += fa[nine];
	}

	fa[self] = accum;
}


__global__ void kernelfunction_SFD_Calculate_Single_Dependency_Chains(int *mask, int *fd, double *fa, int rows, int cols, double *weights, unsigned int* pkgprogressd, int *dep, int *neighbourOffset) {

	int irow = blockIdx.y * blockDim.y + threadIdx.y;
	int icol = blockIdx.x * blockDim.x + threadIdx.x;
	int maxIndex = (rows * cols) - 1;

	if (irow >= rows || icol >= cols)
		return;

	int self = irow * cols + icol;
	
	if (mask[self] == 0) return; // don't calculate if not in catchment(s) of interest

	if(fa[self] > 0) return;

	// 1D index of the cell that the current cell flows into
	assert(fd[self] >= 0 && "\nThis kernel expects all cells to have a positive flow direction assigned. Behaviour is undefined for other values.\n");
	int nextCellInFlow = self + neighbourOffset[__ffs(fd[self]) - 1];
	int currentCell = self;

	while(nextCellInFlow < maxIndex && nextCellInFlow >= 0 && dep[nextCellInFlow] == 1 && fa[nextCellInFlow] < 0) {
		fa[nextCellInFlow] = 1.0 * weights[nextCellInFlow];
		fa[nextCellInFlow] += fa[currentCell];
		currentCell = nextCellInFlow;
		nextCellInFlow = currentCell + neighbourOffset[__ffs(fd[currentCell]) - 1];
		atomicDec(pkgprogressd, 0);
	}
}



int mod_process_SFD_NoPart_List(Data* data, Data* device, int iter) {
	printf("In process\n");

    int gridRows = data->mapInfo.height;
	int gridColumns = data->mapInfo.width;
	int grid1 = gridColumns / (BLOCKCOLS );
	int grid2 = gridRows / (BLOCKROWS );
	int fullsize = gridRows * gridColumns;

	// The number of cells remaining with no FA
	unsigned int progress_h = 0;
	unsigned int *progress_d;
	checkCudaErrors(cudaMalloc((void **) &progress_d, sizeof(progress_h)));
	checkCudaErrors(cudaMemcpy(progress_d, &progress_h, sizeof(progress_h), cudaMemcpyHostToDevice));

	/* An array that will store the number of dependencies: dependencyMap[x], for a given cell x, 
	   as calculated in kernelfunction_SFD_Compute_Deps_And_Resolve_Zero_Deps                  */
	int *dependencyMap;
	checkCudaErrors(cudaMalloc((void **) &dependencyMap, fullsize * sizeof(int)));

	/* NeighbourOffset_h stores the precomputed neighbour offsets that can be added to a cell's 
	   index to find each one of it's neighbours. The index of the offset corresponds to the neighbour's
	   location relative to the current cell index. To get the offset for a neighbour in a given direction,
	   use log2(direction)-1 to compute the index for example, to get the WEST neighbour, simply calculate
	   log2(WEST)-1 = log2(32)-1 = 4. TODO: Refactor this into a device function */
	int neighbourOffset_h[] = {1, gridColumns+1, gridColumns, gridColumns-1, -1, -gridColumns-1, -gridColumns, -gridColumns+1};
	int *neighbourOffset_d;
	checkCudaErrors(cudaMalloc((void **) &neighbourOffset_d, sizeof(neighbourOffset_h)/sizeof(int)));
	checkCudaErrors(cudaMemcpy(neighbourOffset_d, neighbourOffset_h, sizeof(neighbourOffset_h)/sizeof(int), cudaMemcpyHostToDevice));

	dim3 dimGrid(grid1, grid2);
	dim3 dimBlock(BLOCKCOLS, BLOCKROWS);
	unsigned int iterationCount = 0;

	kernelfunction_SFD_Initital_Compute_Deps_And_Calculate_Zero_Dep_Cells<<<dimGrid, dimBlock>>>(device->mask, device->fd, device->fa, gridRows, gridColumns, device->runoffweight, progress_d, dependencyMap);
	assert(cudaGetLastError() == cudaSuccess):

	kernelfunction_SFD_Calculate_Single_Dependency_Chains<<<dimGrid, dimBlock>>>(device->mask, device->fd, device->fa, gridRows, gridColumns, device->runoffweight, progress_d, dependencyMap, neighbourOffset_d);
	assert(cudaGetLastError() == cudaSuccess):

	checkCudaErrors(cudaMemcpy(progress_h, progress_d, sizeof(progress_h), cudaMemcpyDeviceToHost));

	iterationCount++;
	unsigned int lastTot = gridRows * gridColumns;
	unsigned int temp;

	/* All subsequent kernels are called inside the while loop*/
	while (progress_h > 0) { // while the array still has elements
		iterationCount++;
		assert(progress_h > lastTot && "\nERROR: The number of incomplete cells has not fallen. This is impossible.\n");

		lastTot = progress_h;

		temp = 0;
		checkCudaErrors(cudaMemcpy(progress_d, &temp, sizeof(temp), cudaMemcpyHostToDevice));


		kernelfunction_SFD_Calculate_Zero_Dependency_Cells<<<dimGrid, dimBlock>>>(device->mask, device->fd, device->fa, gridRows, gridColumns, device->runoffweight, progress_d);
		assert(cudaGetLastError() == cudaSuccess):

		kernelfunction_SFD_Calculate_Single_Dependency_Chains<<<dimGrid, dimBlock>>>(device->mask, device->fd, device->fa, gridRows, gridColumns, device->runoffweight, progress_d, dependencyMap, neighbourOffset_d);
		assert(cudaGetLastError() == cudaSuccess):

		checkCudaErrors(cudaMemcpy(&progress_h, progress_d, sizeof(progress_h), cudaMemcpyDeviceToHost));
	}
	printf("\nThe total number of FA iterations was: %d\n", iterationCount);
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
				count++;
				}
		}
	}
	fprintf(data->outlog, "FA: Bad value count (i.e. not in catchment(s) = %d\n", count);

	cudaFree(progress_d);
	cudaFree(dependencyMap);
	cudaFree(neighbourOffset_d);

	return 1;
}

void mod_correctflow_SFD_NoPart_List(Data* data, Data* device, int iter) {
	if (cudaSuccess != cudaSetDevice(CUDA_DEVICE)) {
		printf("Unable to access CUDA card\n");
		return ;
	}

	int x;
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

	mod_process_SFD_NoPart_List(data, device, iter);

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	printf("Time to complete FA_SFD_mod_kernels : %.6f s\n", time / 1000.0);
}
