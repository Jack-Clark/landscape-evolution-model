/**
 * @file process_SFD_block_level_single_chain.cu
 * @author Jack Clark (jack.clark@durham.ac.uk)
 * @date 27/9/2016
 * @brief An algorithm for computing SFD flow accumulation on a GPU. Very very slow algorithm ~10-20x slower than NoPart_List algorithm.
 **/

#include "FA_SFD.h"
#include <assert.h>

__global__ void find_deps_and_calculate_zero_dep_cells(int *mask, int *fd, double *fa, int rows, int cols, double *weights, unsigned int* pkgprogressd, int *dep) {

	int irow = blockIdx.y * blockDim.y + threadIdx.y;
	int icol = blockIdx.x * blockDim.x + threadIdx.x;
	int maxSize = rows * cols;

	if (irow >= rows || icol >= cols)
		return;

	int self = irow * cols + icol;

	dep[self] = -1;

	int depCount = 0;

 	int numNeighboursReady = 0;

	if (mask[self] == 0) { // don't calculate if not in catchment(s) of interest
		return;
	}

	if (fa[self] >= 0) return;

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
	if (icol < cols - 1 && irow < rows - 1 && fd[nise] & NORTHWEST) {
		depCount++;
		if (fa[nise] >= 0) {
			numNeighboursReady++;
			accum += fa[nise];
		} 
	}
	if (irow < rows - 1 && fd[nis] & NORTH) {
		depCount++;
		if (fa[nis] >= 0) {
			numNeighboursReady++;
			accum += fa[nis];
		}
	}
	if (icol > 0 && irow < rows - 1 && fd[nisw] & NORTHEAST) {
		depCount++;
		if (fa[nisw] >= 0) {
			numNeighboursReady++;
			accum += fa[nisw];
		}
	}
	if (icol > 0 && fd[niw] & EAST) {
		depCount++;
		if (fa[niw] >= 0) {
			numNeighboursReady++;
			accum += fa[niw];
		}
	}
	if (icol > 0 && irow > 0 && fd[ninw] & SOUTHEAST) {
		depCount++;
		if (fa[ninw] >= 0) {
			numNeighboursReady++;
			accum += fa[ninw];
		} 
	}
	if (irow > 0 && fd[nin] & SOUTH) {
		depCount++;
		if (fa[nin] >= 0) {
			numNeighboursReady++;
			accum += fa[nin];
		}
	}
	if (irow > 0 && icol < cols - 1 && fd[nine] & SOUTHWEST) {
		depCount++;
		if (fa[nine] >= 0) {
			numNeighboursReady++;
			accum += fa[nine];
		}
	}

	if(depCount == numNeighboursReady) {
		fa[self] = accum;
	} else {
		atomicInc(pkgprogressd, maxSize);
	}
	dep[self] = depCount;
}


__global__ void calculate_zero_dependency_cells(int *mask, int *fd, double *fa, int rows, int cols, double *weights, unsigned int* pkgprogressd) {

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
	if (icol < cols - 1 && irow < rows - 1 && fd[nise] & NORTHWEST) {
		if (fa[nise] < 0) {
			atomicInc(pkgprogressd, maxSize);
			return;
		}
		accum += fa[nise];
	}
	if (irow < rows - 1 && fd[nis] & NORTH) {
		if (fa[nis] < 0) {
			atomicInc(pkgprogressd, maxSize);
			return;
		}
		accum += fa[nis];
	}
	if (icol > 0 && irow < rows - 1 && fd[nisw] & NORTHEAST) {
		if (fa[nisw] < 0) {
			atomicInc(pkgprogressd, maxSize);
			return;
		}
		accum += fa[nisw];
	}
	if (icol > 0 && fd[niw] & EAST) {
		if (fa[niw] < 0) {
			atomicInc(pkgprogressd, maxSize);
			return;
		}
		accum += fa[niw];
	}
	if (icol > 0 && irow > 0 && fd[ninw] & SOUTHEAST) {
		if (fa[ninw] < 0) {
			atomicInc(pkgprogressd, maxSize);
			return;
		}
		accum += fa[ninw];
	}
	if (irow > 0 && fd[nin] & SOUTH) {
		if (fa[nin] < 0) {
			atomicInc(pkgprogressd, maxSize);
			return;
		}
		accum += fa[nin];
	}
	if (irow > 0 && icol < cols - 1 && fd[nine] & SOUTHWEST) {
		if (fa[nine] < 0) {
			atomicInc(pkgprogressd, maxSize);
			return;
		}
		accum += fa[nine];
	}

	fa[self] = accum;
}


__global__ void calculate_single_dependency_chains(int *mask, int *fd, double *fa, int rows, int cols, double *weights, unsigned int* pkgprogressd, int *dep, int *neighbourOffset) {


	int irow = blockIdx.y * blockDim.y + threadIdx.y;
	int icol = blockIdx.x * blockDim.x + threadIdx.x;
	int maxIndex = (rows * cols) - 1;

	if (irow >= rows || icol >= cols)
		return;

	int self = irow * cols + icol;

	if (mask[self] == 0) return; // don't calculate if not in catchment(s) of interest

	if(fa[self] < 0) return;

	int currentCell = self;
	int nextCellInFlow = currentCell + neighbourOffset[__ffs(fd[currentCell]) - 1];

	int chainLength = 1;

	while(nextCellInFlow < maxIndex && nextCellInFlow >= 0 && mask[nextCellInFlow] != 0 && dep[nextCellInFlow] == 1 && fa[nextCellInFlow] < 0) {
		fa[nextCellInFlow] = 1.0 * weights[nextCellInFlow];
		fa[nextCellInFlow] += fa[currentCell];
		currentCell = nextCellInFlow;
		nextCellInFlow = currentCell + neighbourOffset[__ffs(fd[currentCell]) - 1];
		atomicSub(pkgprogressd, 1);
		chainLength++;
	}
}


int process_SFD_global_level_single_chain(Data* data, Data* device, int iter) {
	printf("In process\n");

    int rows = data->mapInfo.height;
	int cols = data->mapInfo.width;
	int gridCols = cols / BLOCKCOLS;
	int gridRows = rows / BLOCKROWS;
	int totalCells = rows * cols;

	// The number of cells remaining with no FA
	unsigned int numCellsRemaining_h = 0;
	unsigned int *numCellsRemaining_d;
	checkCudaErrors(cudaMalloc((void **) &numCellsRemaining_d, sizeof(numCellsRemaining_h)));
	checkCudaErrors(cudaMemcpy(numCellsRemaining_d, &numCellsRemaining_h, sizeof(numCellsRemaining_h), cudaMemcpyHostToDevice));

	/* An array that will store the number of dependencies: dependencyMap[x], for a given cell x, 
	   as calculated in kernelfunction_SFD_Compute_Deps_And_Resolve_Zero_Deps                  */
	int *dependencyMap;
	checkCudaErrors(cudaMalloc((void **) &dependencyMap, totalCells * sizeof(int)));

	/* NeighbourOffset_h stores the precomputed neighbour offsets that can be added to a cell's 
	   index to find each one of it's neighbours. The index of the offset corresponds to the neighbour's
	   location relative to the current cell index. To get the offset for a neighbour in a given direction,
	   use log2(direction)-1 to compute the index for example, to get the WEST neighbour, simply calculate
	   log2(WEST)-1 = log2(32)-1 = 4. TODO: Refactor this into a device function */
	int neighbourOffset_h[] = {1, cols+1, cols, cols-1, -1, -cols-1, -cols, -cols+1};
	int *neighbourOffset_d;
	checkCudaErrors(cudaMalloc((void **) &neighbourOffset_d, sizeof(neighbourOffset_h)/sizeof(int)));
	checkCudaErrors(cudaMemcpy(neighbourOffset_d, neighbourOffset_h, sizeof(neighbourOffset_h)/sizeof(int), cudaMemcpyHostToDevice));

	dim3 dimGrid(gridCols, gridRows);
	dim3 dimBlock(BLOCKCOLS, BLOCKROWS);

	find_deps_and_calculate_zero_dep_cells<<<dimGrid, dimBlock>>>(device->mask, device->fd, device->fa, rows, cols, device->runoffweight, numCellsRemaining_d, dependencyMap);

	calculate_single_dependency_chains<<<dimGrid, dimBlock>>>(device->mask, device->fd, device->fa, rows, cols, device->runoffweight, numCellsRemaining_d, dependencyMap, neighbourOffset_d);

	checkCudaErrors(cudaMemcpy(&numCellsRemaining_h, numCellsRemaining_d, sizeof(numCellsRemaining_h), cudaMemcpyDeviceToHost));

	unsigned int temp;

	unsigned int lastTot = rows * cols;

	while (numCellsRemaining_h > 0) {

		assert(numCellsRemaining_h < lastTot);

		lastTot = numCellsRemaining_h;

		temp = 0;

		checkCudaErrors(cudaMemcpy(numCellsRemaining_d, &temp, sizeof(temp), cudaMemcpyHostToDevice));

		calculate_zero_dependency_cells<<<dimGrid, dimBlock>>>(device->mask, device->fd, device->fa, rows, cols, device->runoffweight, numCellsRemaining_d);

		calculate_single_dependency_chains<<<dimGrid, dimBlock>>>(device->mask, device->fd, device->fa, rows, cols, device->runoffweight, numCellsRemaining_d, dependencyMap, neighbourOffset_d);

		checkCudaErrors(cudaMemcpy(&numCellsRemaining_h, numCellsRemaining_d, sizeof(numCellsRemaining_h), cudaMemcpyDeviceToHost));
	}


	checkCudaErrors(cudaMemcpy(data->fa, device->fa, totalCells * sizeof(double),   cudaMemcpyDeviceToHost));
	fprintf(data->outlog, "FA: FA memcopy back operation 3:%s\n", cudaGetErrorString(cudaGetLastError()));

	double FA_max;
	int FAindex = 0;
	double cpuFA_max = 0.0;

	if (iter == 1) {
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (data->fa[i * cols + j] > cpuFA_max) {
					cpuFA_max = data->fa[i * cols + j];
					FAindex = i * cols + j;
				}
			}
		}
		data->FA_max = cpuFA_max;
		data->outletcellidx = FAindex; // this is the outlet cell which will be maintained throughout the simulation
	} else {
			thrust::device_ptr<double> max_FA = thrust::device_pointer_cast(device->fa);
			FA_max = thrust::reduce(max_FA, max_FA + totalCells, (double) 0, thrust::maximum<double>());
			data->FA_max = FA_max;
	}

	fprintf(data->outlog, "FA: Maximum FA is  %.6f s\n\n", data->FA_max);
	fprintf(data->outlog, "FA: Outletcell index is  %d s\n\n", data->outletcellidx);

	printf("Maximum FA is  %.6f s\n\n", data->FA_max);

	int count = 0;

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (data->fa[i * cols + j] < 0) {
				count++;
			}
		}
	}
	fprintf(data->outlog, "FA: Bad value count (i.e. not in catchment(s) = %d\n", count);

	cudaFree(numCellsRemaining_d);
	cudaFree(dependencyMap);
	cudaFree(neighbourOffset_d);

	return 1;
}
