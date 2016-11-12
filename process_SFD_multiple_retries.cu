/**
 * @file process_SFD_multiple_retries.cu
 * @author Jack Clark (jack.clark@durham.ac.uk)
 * @date 8/10/2016
 * @brief multiple retries algorithm for computing SFD flow accumulation on a GPU. ~2x faster than NoPart_List algorithm.
 **/

#include "FA_SFD.h"

#define MULTI_RETRY_BLOCK_SIZE 128
#define FA_CELL_RETRIES 25

__global__ void SFD_Multiple_Retries_Initial(int *mask, int *fd, double *fa, int rows, int cols, double *weights, unsigned int* pkgprogressd, unsigned int* left)
{
	int irow = blockIdx.y * blockDim.y + threadIdx.y;
	int icol = blockIdx.x * blockDim.x + threadIdx.x;
	int maxSize = rows * cols;

	if (irow >= rows || icol >= cols)
		return;

	int self = irow * cols + icol;
	if (mask[self] == 0) return;
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

__global__ void SFD_Multiple_Retries_Progress_All_Shmem(int *mask, int *fd, double *fa, int rows, int cols, double *weights, unsigned int* pkgprogressd, unsigned int* left, unsigned int* had, int hadsize)
{
	int pos = blockIdx.y*cols*blockDim.y + blockIdx.x*blockDim.x + threadIdx.y*cols + threadIdx.x;

	int nie, nise, nis, nisw, niw, ninw, nin, nine, self, myPos, irow, icol, readyFlag, maxSize, faHasBeenCalculated, i;
	double accum;

	__shared__ int shmem_fd[MULTI_RETRY_BLOCK_SIZE * 8];

	__shared__ double shmem_fa[MULTI_RETRY_BLOCK_SIZE * 9];

	if (pos < hadsize) {

		self = had[pos];

		nie  = self        + 1 ;
		nise = self + cols + 1 ;
		nis  = self + cols     ;
		nisw = self + cols - 1 ;
		niw  = self        - 1 ;
		ninw = self - cols - 1 ;
		nin  = self - cols     ;
		nine = self - cols + 1 ;

		icol = self % cols;
		irow = self / cols;

		faHasBeenCalculated = 0;

		shmem_fd[threadIdx.x * 8]     = fd[nie];
		shmem_fd[threadIdx.x * 8 + 1] = fd[nise];
		shmem_fd[threadIdx.x * 8 + 2] = fd[nis];
		shmem_fd[threadIdx.x * 8 + 3] = fd[nisw];
		shmem_fd[threadIdx.x * 8 + 4] = fd[niw];
		shmem_fd[threadIdx.x * 8 + 5] = fd[ninw];
		shmem_fd[threadIdx.x * 8 + 6] = fd[nin];
		shmem_fd[threadIdx.x * 8 + 7] = fd[nine];

	}

	for(i = 0; i < FA_CELL_RETRIES; i++) {

		if(pos < hadsize && !faHasBeenCalculated) {

			shmem_fa[threadIdx.x * 9]     = fa[self];
			shmem_fa[threadIdx.x * 9 + 1] = fa[nie];
			shmem_fa[threadIdx.x * 9 + 2] = fa[nise];
			shmem_fa[threadIdx.x * 9 + 3] = fa[nis];
			shmem_fa[threadIdx.x * 9 + 4] = fa[nisw];
			shmem_fa[threadIdx.x * 9 + 5] = fa[niw];
			shmem_fa[threadIdx.x * 9 + 6] = fa[ninw];
			shmem_fa[threadIdx.x * 9 + 7] = fa[nin];
			shmem_fa[threadIdx.x * 9 + 8] = fa[nine];

			accum = 1.0 * weights[self];
			readyFlag = 1;
			
			if (icol < cols - 1 && shmem_fd[threadIdx.x * 8] & WEST && readyFlag) {
				if (shmem_fa[threadIdx.x * 9 + 1] < 0) {
					readyFlag = 0;
				}
				accum += shmem_fa[threadIdx.x * 9 + 1];
			}
			
			if (icol < cols - 1 && irow < rows - 1 && shmem_fd[threadIdx.x * 8 + 1] & NORTHWEST && readyFlag) {
				if (shmem_fa[threadIdx.x * 9 + 2] < 0)  {
					readyFlag = 0;
				}
				accum += shmem_fa[threadIdx.x * 9 + 2];
			}
			
			if (irow < rows - 1 && shmem_fd[threadIdx.x * 8 + 2] & NORTH && readyFlag) {
				if (shmem_fa[threadIdx.x * 9 + 3] < 0)  {
					readyFlag = 0;
				}
				accum += shmem_fa[threadIdx.x * 9 + 3];
			}
			
			if (icol > 0 && irow < rows - 1 && shmem_fd[threadIdx.x * 8 + 3] & NORTHEAST && readyFlag) {
				if (shmem_fa[threadIdx.x * 9 + 4] < 0)  {
					readyFlag = 0;
				}
				accum += shmem_fa[threadIdx.x * 9 + 4];
			}
			
			if (icol > 0 && shmem_fd[threadIdx.x * 8 + 4] & EAST && readyFlag) {
				if (shmem_fa[threadIdx.x * 9 + 5] < 0)  {
					readyFlag = 0;
				}
				accum += shmem_fa[threadIdx.x * 9 + 5];
			}
			
			if (icol > 0 && irow > 0 && shmem_fd[threadIdx.x * 8 + 5] & SOUTHEAST && readyFlag) {
				if (shmem_fa[threadIdx.x * 9 + 6] < 0)  {
					readyFlag = 0;
				}
				accum += shmem_fa[threadIdx.x * 9 + 6];
			}
			
			if (irow > 0 && shmem_fd[threadIdx.x * 8 + 6] & SOUTH && readyFlag) {
				if (shmem_fa[threadIdx.x * 9 + 7] < 0)  {
					readyFlag = 0;
				}
				accum += shmem_fa[threadIdx.x * 9 + 7];
			}
			
			if (irow > 0 && icol < cols - 1 && shmem_fd[threadIdx.x * 8 + 7] & SOUTHWEST && readyFlag) {
				if (shmem_fa[threadIdx.x * 9 + 8] < 0)  {
					readyFlag = 0;
				}
				accum += shmem_fa[threadIdx.x * 9 + 8];
			}
			
			if(readyFlag) {
				shmem_fa[threadIdx.x * 9] = accum;
				fa[self] = shmem_fa[threadIdx.x * 9];
				faHasBeenCalculated = 1;
			}

		}

		__syncthreads();
	}

	if(pos < hadsize) {
		if(!faHasBeenCalculated) {
			myPos = 0;
			maxSize = rows * cols;
			myPos = atomicInc(pkgprogressd, maxSize);
			left[myPos] = self;
		}
			
	}

}

__global__ void SFD_Multiple_Retries_Progress_FD_Shmem(int *mask, int *fd, double *fa, int rows, int cols, double *weights, unsigned int* pkgprogressd, unsigned int* left, unsigned int* had, int hadsize)
{
	int pos = blockIdx.y*cols*blockDim.y + blockIdx.x*blockDim.x + threadIdx.y*cols + threadIdx.x;

	int nie, nise, nis, nisw, niw, ninw, nin, nine, self, myPos, irow, icol, readyFlag, maxSize, faHasBeenCalculated, i;
	double accum;

	__shared__ int shmem_fd[MULTI_RETRY_BLOCK_SIZE * 8];

	if (pos < hadsize) {

		self = had[pos];

		nie  = self        + 1 ;
		nise = self + cols + 1 ;
		nis  = self + cols     ;
		nisw = self + cols - 1 ;
		niw  = self        - 1 ;
		ninw = self - cols - 1 ;
		nin  = self - cols     ;
		nine = self - cols + 1 ;

		icol = self % cols;
		irow = self / cols;

		faHasBeenCalculated = 0;

		shmem_fd[threadIdx.x * 8]     = fd[nie];
		shmem_fd[threadIdx.x * 8 + 1] = fd[nise];
		shmem_fd[threadIdx.x * 8 + 2] = fd[nis];
		shmem_fd[threadIdx.x * 8 + 3] = fd[nisw];
		shmem_fd[threadIdx.x * 8 + 4] = fd[niw];
		shmem_fd[threadIdx.x * 8 + 5] = fd[ninw];
		shmem_fd[threadIdx.x * 8 + 6] = fd[nin];
		shmem_fd[threadIdx.x * 8 + 7] = fd[nine];

	}

	for(i = 0; i < FA_CELL_RETRIES; i++) {

		if(pos < hadsize && !faHasBeenCalculated) {

			accum = 1.0 * weights[self];
			readyFlag = 1;
			
			if (icol < cols - 1 && shmem_fd[threadIdx.x * 8] & WEST && readyFlag) {
				if (fa[nie] < 0) {
					readyFlag = 0;
				}
				accum += fa[nie];
			}
			
			if (icol < cols - 1 && irow < rows - 1 && shmem_fd[threadIdx.x * 8 + 1] & NORTHWEST && readyFlag) {
				if (fa[nise] < 0)  {
					readyFlag = 0;
				}
				accum += fa[nise];
			}
			
			if (irow < rows - 1 && shmem_fd[threadIdx.x * 8 + 2] & NORTH && readyFlag) {
				if (fa[nis] < 0)  {
					readyFlag = 0;
				}
				accum += fa[nis];
			}
			
			if (icol > 0 && irow < rows - 1 && shmem_fd[threadIdx.x * 8 + 3] & NORTHEAST && readyFlag) {
				if (fa[nisw] < 0)  {
					readyFlag = 0;
				}
				accum += fa[nisw];
			}
			
			if (icol > 0 && shmem_fd[threadIdx.x * 8 + 4] & EAST && readyFlag) {
				if (fa[niw] < 0)  {
					readyFlag = 0;
				}
				accum += fa[niw];
			}
			
			if (icol > 0 && irow > 0 && shmem_fd[threadIdx.x * 8 + 5] & SOUTHEAST && readyFlag) {
				if (fa[ninw] < 0)  {
					readyFlag = 0;
				}
				accum += fa[ninw];
			}
			
			if (irow > 0 && shmem_fd[threadIdx.x * 8 + 6] & SOUTH && readyFlag) {
				if (fa[nin] < 0)  {
					readyFlag = 0;
				}
				accum += fa[nin];
			}
			
			if (irow > 0 && icol < cols - 1 && shmem_fd[threadIdx.x * 8 + 7] & SOUTHWEST && readyFlag) {
				if (fa[nine] < 0)  {
					readyFlag = 0;
				}
				accum += fa[nine];
			}
			
			if(readyFlag) {
				fa[self] = accum;
				faHasBeenCalculated = 1;
			}

		}

		__syncthreads();
	}

	if(pos < hadsize) {
		if(!faHasBeenCalculated) {
			myPos = 0;
			maxSize = rows * cols;
			myPos = atomicInc(pkgprogressd, maxSize);
			left[myPos] = self;
		}
			
	}

}


int process_SFD_Multiple_Retries(Data* data, Data* device, int iter) {
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
	checkCudaErrors(cudaMalloc((void **) &list1, totalCells * sizeof(unsigned int)) );
	checkCudaErrors(cudaMalloc((void **) &list2, totalCells * sizeof(unsigned int)) );

	dim3 dimGrid(gridCols, gridRows);
	dim3 dimBlock(BLOCKCOLS, BLOCKROWS);

	SFD_Multiple_Retries_Initial<<<dimGrid, dimBlock>>>(device->mask, device->fd, device->fa, rows, cols, device->runoffweight, numCellsRemaining_d, list1);

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

		nextBlocks = *numCellsRemaining_h / MULTI_RETRY_BLOCK_SIZE + 1;

		*temp = 0;

		checkCudaErrors(cudaMemcpy(numCellsRemaining_d, temp, sizeof(unsigned int), cudaMemcpyHostToDevice) );

		SFD_Multiple_Retries_Progress_FD_Shmem<<<nextBlocks, MULTI_RETRY_BLOCK_SIZE>>>(device->mask, device->fd, device->fa, rows, cols, device->runoffweight, numCellsRemaining_d, list1, list2, *numCellsRemaining_h);

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
		data->outletcellidx = FAindex; // This is the outlet cell which will be maintained throughout the simulation
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

	cudaFree(list1);
	cudaFree(list2);
	cudaFree(numCellsRemaining_d);
	free(temp);
	free(numCellsRemaining_h);

	return 0;
}