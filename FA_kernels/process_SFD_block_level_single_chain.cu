/**
 * @file process_SFD_block_level_single_chain.cu
 * @author Jack Clark (jack.clark@durham.ac.uk)
 * @date 1/10/2016
 * @brief An algorithm for computing SFD flow accumulation on a GPU. ~2x slower than NoPart_List algorithm.
 **/

#include "FA_SFD.h"

#define BLOCK_SIZE 512

__global__ void find_deps_and_compute_zero_dep_cells(int *mask, int *fd, double *fa, int rows, int cols, double *weights, unsigned int* pkgprogressd, unsigned int* left, int * dep)
{
	int irow = blockIdx.y * blockDim.y + threadIdx.y;
	int icol = blockIdx.x * blockDim.x + threadIdx.x;
	int maxSize = rows * cols;

	if (irow < rows && icol < cols) {

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

		int dependencyCount = 0;

		int readyFlag = 1;

		if (icol < cols - 1 && fd[nie] & WEST) {
			dependencyCount++;
			if (fa[nie] < 0 && readyFlag) {
				myPos = atomicInc(pkgprogressd, maxSize);
			    left[myPos] = self;
				readyFlag = 0;
			}
			accum += fa[nie];
		}
		if (icol < cols - 1 && irow < rows - 1 && fd[nise] & NORTHWEST) {
			dependencyCount++;
			if (fa[nise] < 0 && readyFlag)  {
				myPos = atomicInc(pkgprogressd, maxSize);
			    left[myPos] = self;
				readyFlag = 0;
			}
			accum += fa[nise];
		}
		if (irow < rows - 1 && fd[nis] & NORTH) {
			dependencyCount++;
			if (fa[nis] < 0 && readyFlag)  {
				myPos = atomicInc(pkgprogressd, maxSize);
			    left[myPos] = self;
				readyFlag = 0;
			}
			accum += fa[nis];
		}
		if (icol > 0 && irow < rows - 1 && fd[nisw] & NORTHEAST) {
			dependencyCount++;
			if (fa[nisw] < 0 && readyFlag)  {
				myPos = atomicInc(pkgprogressd, maxSize);
			    left[myPos] = self;
				readyFlag = 0;
			}
			accum += fa[nisw];
		}
		if (icol > 0 && fd[niw] & EAST) {
			dependencyCount++;
			if (fa[niw] < 0 && readyFlag)  {
				myPos = atomicInc(pkgprogressd, maxSize);
			    left[myPos] = self;
				readyFlag = 0;
			}
			accum += fa[niw];
		}
		if (icol > 0 && irow > 0 && fd[ninw] & SOUTHEAST) {
			dependencyCount++;
			if (fa[ninw] < 0 && readyFlag)  {
				myPos = atomicInc(pkgprogressd, maxSize);
			    left[myPos] = self;
				readyFlag = 0;
			}
			accum += fa[ninw];
		}
		if (irow > 0 && fd[nin] & SOUTH) {
			dependencyCount++;
			if (fa[nin] < 0 && readyFlag)  {
				myPos = atomicInc(pkgprogressd, maxSize);
			    left[myPos] = self;
				readyFlag = 0;
			}
			accum += fa[nin];
		}
		if (irow > 0 && icol < cols - 1 && fd[nine] & SOUTHWEST) {
			dependencyCount++;
			if (fa[nine] < 0 && readyFlag)  {
				myPos = atomicInc(pkgprogressd, maxSize);
			    left[myPos] = self;
				readyFlag = 0;
			}
			accum += fa[nine];
		}

		assert(dependencyCount >= 0);
		assert(dependencyCount <= 8);

		if(readyFlag) fa[self] = accum;
		
		dep[self] = dependencyCount;
		assert(dep[self] >= 0);
		assert(dep[self] <= 8);
	}	
}


__global__ void Assign_BlockIds_To_Threads(int cols, unsigned int* had, int hadsize, int *blockIds) {
	int pos = blockIdx.y*cols*blockDim.y + blockIdx.x*blockDim.x + threadIdx.y*cols + threadIdx.x;

	if(pos < hadsize) {
		int self = had[pos];
		blockIds[self] = pos/blockDim.x;
		assert(blockIds[self] >= 0);
		assert(blockIds[self] < hadsize/blockDim.x + 1);
	}

}

__global__ void calculate_block_single_chains(int *mask, int *fd, double *fa, int rows, int cols, double *weights, unsigned int* pkgprogressd, unsigned int* left, unsigned int* had, int hadsize, int * dep, int * neighbourOffset, int *blockIds)
{
	int pos = blockIdx.y*cols*blockDim.y + blockIdx.x*blockDim.x + threadIdx.y*cols + threadIdx.x;
	
	__shared__ int shmem_offset[9];

	if(threadIdx.x <= 8) {
		shmem_offset[threadIdx.x] = neighbourOffset[threadIdx.x];
	}

	__syncthreads();

	int maxSize = rows * cols;

	int nie, nise, nis, nisw, niw, ninw, nin, nine;
	int self;
	int readyFlag;
	int myPos;

	if (pos < hadsize) {

		self = had[pos];

		double accum = 1.0 * weights[self];

		nie  = self        + 1 ;
		nise = self + cols + 1 ;
		nis  = self + cols     ;
		nisw = self + cols - 1 ;
		niw  = self        - 1 ;
		ninw = self - cols - 1 ;
		nin  = self - cols     ;
		nine = self - cols + 1 ;

		int icol = self % cols;
		int irow = self / cols;

		readyFlag = 1;

		if (icol < cols - 1 && fd[nie] & WEST && readyFlag) {
			if (fa[nie] < 0) {
				readyFlag = 0;
			}
			accum += fa[nie];
		}
		if (icol < cols - 1 && irow < rows - 1 && fd[nise] & NORTHWEST && readyFlag) {
			if (fa[nise] < 0)  {
				readyFlag = 0;
			}
			accum += fa[nise];
		}
		if (irow < rows - 1 && fd[nis] & NORTH && readyFlag) {
			if (fa[nis] < 0)  {
				readyFlag = 0;
			}
			accum += fa[nis];
		}
		if (icol > 0 && irow < rows - 1 && fd[nisw] & NORTHEAST && readyFlag) {
			if (fa[nisw] < 0)  {
				readyFlag = 0;
			}
			accum += fa[nisw];
		}
		if (icol > 0 && fd[niw] & EAST && readyFlag) {
			if (fa[niw] < 0)  {
				readyFlag = 0;
			}
			accum += fa[niw];
		}
		if (icol > 0 && irow > 0 && fd[ninw] & SOUTHEAST && readyFlag) {
			if (fa[ninw] < 0)  {
				readyFlag = 0;
			}
			accum += fa[ninw];
		}
		if (irow > 0 && fd[nin] & SOUTH && readyFlag) {
			if (fa[nin] < 0)  {
				readyFlag = 0;
			}
			accum += fa[nin];
		}
		if (irow > 0 && icol < cols - 1 && fd[nine] & SOUTHWEST && readyFlag) {
			if (fa[nine] < 0)  {
				readyFlag = 0;
			}
			accum += fa[nine];
		}

		if(readyFlag) fa[self] = accum;
	}

	__syncthreads();

	if(pos < hadsize) {

		int currentCell = self;
		assert(dep[currentCell] >= 0);
		assert(dep[currentCell] <= 8);
		assert(fd[currentCell] >= 0);
		assert(fd[currentCell] <= 128);
		int nextCellInFlow = currentCell + shmem_offset[__ffs(fd[currentCell])];

		while(nextCellInFlow < maxSize && nextCellInFlow >= 0 && blockIdx.x == blockIds[nextCellInFlow] && mask[nextCellInFlow] != 0 && dep[nextCellInFlow] == 1 && fa[nextCellInFlow] < 0 && readyFlag) {
			fa[nextCellInFlow] = 1.0 * weights[nextCellInFlow];
			fa[nextCellInFlow] += fa[currentCell];
			currentCell = nextCellInFlow;
			assert(fd[currentCell] >= 0);
			assert(fd[currentCell] <= 128);
			nextCellInFlow = currentCell + shmem_offset[__ffs(fd[currentCell])];
		}
	}

	__syncthreads();

	if(pos < hadsize) {
		myPos = 0;
		if(fa[self] < 0) {
			myPos = atomicInc(pkgprogressd, maxSize);
			left[myPos] = self;
		}
	}
}	



int process_SFD_block_level_single_chain(Data* data, Data* device, int iter) {
	printf("In process\n");

    int rows = data->mapInfo.height;
	int cols = data->mapInfo.width;
	int gridCols = cols / BLOCKCOLS;
	int gridRows = rows / BLOCKROWS;
	int totalCells = rows * cols;

	unsigned int *numCellsRemaining_d;
	checkCudaErrors(cudaMalloc((void **) &numCellsRemaining_d, sizeof(unsigned int)) );
	unsigned int *numCellsRemaining_h = (unsigned int*) malloc(sizeof(unsigned int));
	*numCellsRemaining_h = 0;

	checkCudaErrors(cudaMemcpy(numCellsRemaining_d, numCellsRemaining_h, sizeof(unsigned int), cudaMemcpyHostToDevice));

	unsigned int* list1;
	unsigned int* list2;
	unsigned int* listT;
	checkCudaErrors(cudaMalloc((void **) &list1, totalCells * sizeof(unsigned int)) ); //FIXME - need to figure out a way to set this value
	checkCudaErrors(cudaMalloc((void **) &list2, totalCells * sizeof(unsigned int)) ); //FIXME - need to figure out a way to set this value

	/* An array that will store the number of dependencies: dependencyMap[x], for a given cell x, 
	   as calculated in kernelfunction_SFD_Compute_Deps_And_Resolve_Zero_Deps                  */
	int *dependencyMap;
	checkCudaErrors(cudaMalloc((void **) &dependencyMap, totalCells * sizeof(int)));

	/* NeighbourOffset_h stores the precomputed neighbour offsets that can be added to a cell's 
	   index to find each one of it's neighbours. The index of the offset corresponds to the neighbour's
	   location relative to the current cell index. To get the offset for a neighbour in a given direction,
	   use log2(direction)-1 to compute the index, for example, to get the WEST neighbour, simply calculate
	   log2(WEST)-1 = log2(32)-1 = 4. */
	int neighbourOffset_h[] = {-(cols * rows) - 1, 1, cols+1, cols, cols-1, -1, -cols-1, -cols, -cols+1};
	int *neighbourOffset_d;
	checkCudaErrors(cudaMalloc((void **) &neighbourOffset_d, sizeof(neighbourOffset_h)/sizeof(int)));
	checkCudaErrors(cudaMemcpy(neighbourOffset_d, neighbourOffset_h, sizeof(neighbourOffset_h)/sizeof(int), cudaMemcpyHostToDevice));

	int *blockIds;
	checkCudaErrors(cudaMalloc((void **) &blockIds, sizeof(int) * totalCells));

	dim3 dimGrid(gridCols, gridRows);
	dim3 dimBlock(BLOCKCOLS, BLOCKROWS);

	find_deps_and_compute_zero_dep_cells<<<dimGrid, dimBlock>>>(device->mask, device->fd, device->fa, rows, cols, device->runoffweight, numCellsRemaining_d, list1, dependencyMap);

	checkCudaErrors(cudaMemcpy(numCellsRemaining_h, numCellsRemaining_d, sizeof(unsigned int), cudaMemcpyDeviceToHost));

	unsigned int lastTot = rows * cols;

	unsigned int* temp = (unsigned int*) malloc(sizeof(unsigned int));

	int nextBlocks;

	while (*numCellsRemaining_h > 0) {

		listT = list1;
		list1 = list2;
		list2 = listT;

		assert(*numCellsRemaining_h < lastTot);
		
		lastTot = *numCellsRemaining_h;

		nextBlocks = *numCellsRemaining_h / BLOCK_SIZE + 1;

		*temp = 0;

		checkCudaErrors(cudaMemcpy(numCellsRemaining_d, temp, sizeof(unsigned int), cudaMemcpyHostToDevice) );

		Assign_BlockIds_To_Threads<<<nextBlocks, BLOCK_SIZE>>>(cols, list2, *numCellsRemaining_h, blockIds);
			
		calculate_block_single_chains<<<nextBlocks, BLOCK_SIZE>>>(device->mask, device->fd, device->fa, rows, cols, device->runoffweight, numCellsRemaining_d, list1, list2, *numCellsRemaining_h, dependencyMap, neighbourOffset_d, blockIds);

		checkCudaErrors(cudaMemcpy(numCellsRemaining_h, numCellsRemaining_d, sizeof(unsigned int), cudaMemcpyDeviceToHost) );
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
	printf("Count: %d\n", count);
	fprintf(data->outlog, "FA: Bad value count (i.e. not in catchment(s) = %d\n", count);

	cudaFree(list1);
	cudaFree(list2);
	cudaFree(neighbourOffset_d);
	cudaFree(blockIds);
	cudaFree(dependencyMap);
	cudaFree(numCellsRemaining_d);
	free(temp);
	free(numCellsRemaining_h);

	return 1;
}