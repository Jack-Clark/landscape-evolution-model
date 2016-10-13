
#include "FA_SFD.h"
__global__ void kernelfunction_scaleFA(int *mask, double *fa, int *contribA, int rows, int cols, int totalA)
{
	int irow = blockIdx.y * blockDim.y + threadIdx.y;
	int icol = blockIdx.x * blockDim.x + threadIdx.x;

	if (irow >= rows || icol >= cols)
		return;

	int self = irow * cols + icol;

	if (mask[self] == 0) return; // don't calculate if not in catchment(s) of interest

	fa[self] = fa[self] * (contribA[self] / totalA);

	}

__global__ void kernelfunction_FA(int *mask, int *fd, double *fa, int *progressed, int rows, int cols, double *weights)
{
	int irow = blockIdx.y * blockDim.y + threadIdx.y;
	int icol = blockIdx.x * blockDim.x + threadIdx.x;

	if (irow >= rows || icol >= cols)
		return;

	int self = irow * cols + icol;

	if (mask[self] == 0) return; // don't calculate if not in catchment(s) of interest

	int nie, nise, nis, nisw, niw, ninw, nin, nine;

	if (fa[self] >= 0) return;

	double accum = weights[self];

	nie  = self        + 1 ;
	nise = self + cols + 1 ;
	nis  = self + cols     ;
	nisw = self + cols - 1 ;
	niw  = self        - 1 ;
	ninw = self - cols - 1 ;
	nin  = self - cols     ;
	nine = self - cols + 1 ;

	if (icol < cols - 1 && (fd[nie] & WEST)) {
		if (fa[nie] < 0) return;

		accum += fa[nie];
	}
	if (icol < cols - 1 && irow < rows - 1 && (fd[nise] & NORTHWEST)) {
		if (fa[nise] < 0) return;

		accum += fa[nise];
	}
	if (irow < rows - 1 && (fd[nis] & NORTH)) {
		if (fa[nis] < 0) return;

		accum += fa[nis];
	}
	if (icol > 0 && irow < rows - 1 && (fd[nisw] & NORTHEAST)) {
		if (fa[nisw] < 0) return;

		accum += fa[nisw];
	}
	if (icol > 0 && (fd[niw] & EAST)) {
		if (fa[niw] < 0) return;

		accum += fa[niw];
	}
	if (icol > 0 && irow > 0 && (fd[ninw] & SOUTHEAST)) {
		if (fa[ninw] < 0) return;

		accum += fa[ninw];
	}
	if (irow > 0 && (fd[nin] & SOUTH)) {
		if (fa[nin] < 0) return;

		accum += fa[nin];
	}
	if (irow > 0 && icol < cols - 1 && (fd[nine] & SOUTHWEST)) {
		if (fa[nine] < 0) return;

		accum += fa[nine];
	}
	fa[self] = accum;
	*progressed = 1;
}


int scale_FA(int* mask, int* contribA, double* fa, int totalA, int rows, int cols, int grid1, int grid2, int grid3, int blockRows, int blockColumns, int dimBlock3)
{
	dim3 dimGrid(grid1, grid2);
	dim3 dimBlock(blockColumns, blockRows);

	kernelfunction_scaleFA<<<dimGrid, dimBlock>>>(mask, fa, contribA, rows, cols, totalA);

	return 0;
}




int process_FA(Data* data, int* mask, int *fd, double *fa, double *runoffweight, int gridRows, int gridColumns, int dimGrid1, int dimGrid2, int dimGrid3, int blockRows, int blockColumns, int dimBlock3)
{
	// Create space for data on card
	int *progressed;
	cudaMalloc((void **) &progressed, sizeof(int));

	int *progress = (int*) malloc(sizeof(int));
	*progress = 0;

	dim3 dimGrid(dimGrid1, dimGrid2);
	dim3 dimBlock(blockColumns, blockRows);
	int consecutiveZero = 0;


	int count = 0;
	do {
		*progress = 0;
		cudaMemcpy(progressed, progress, sizeof(int), cudaMemcpyHostToDevice);

		kernelfunction_FA<<<dimGrid, dimBlock>>>(mask, fd, fa, progressed, gridRows, gridColumns, runoffweight); //send the wieghts package to the kernel
		if (count <3) fprintf(data->outlog, "FA: kernelfunctionFA :%s\n", cudaGetErrorString(cudaGetLastError()));
		//printf("kernel fuction called:%s\n", cudaGetErrorString(cudaGetLastError()));
		count ++;
		cudaMemcpy(progress, progressed, sizeof(int), cudaMemcpyDeviceToHost);
		//printf("copy progressed:%s\n", cudaGetErrorString(cudaGetLastError()));
		if (*progress == 0) {
			//printf("Progress was not made %d \n", count);
			consecutiveZero ++;
		}
		else
			consecutiveZero = 0;
	} while (consecutiveZero < 2);

	cudaFree(progressed);

	return *progress;
}

void accumulateflowSFD(Data* data, Data* device, int iter) {
    int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	int grid1 = data->mapInfo.width  / (BLOCKCOLS );
	int grid2 = data->mapInfo.height / (BLOCKROWS );
	int fullsize = ncell_x * ncell_y;

	process_FA(data, device->mask, device->fd, device->fa, device->runoffweight, ncell_y, ncell_x, grid1, grid2, 1, BLOCKROWS, BLOCKCOLS, 1);

	// Copy flow accumulation back
	int count = 0;

	cudaMemcpy(data->fa, device->fa, fullsize * sizeof(double),   cudaMemcpyDeviceToHost);
	fprintf(data->outlog, "FA: FA memcopy back operation 3:%s\n", cudaGetErrorString(cudaGetLastError()));

	printf("FA: FA memcopy back operation 3:%s\n", cudaGetErrorString(cudaGetLastError()));


	double FA_max;
	int FAindex = 0;
	double cpuFA_max = 0.0;

	if (iter == 1) // cpu calculation otherwise we cannot locate the outletcell index
	{
		for (int i = 0; i < ncell_y; i++) {
			for (int j = 0; j < ncell_x; j++) {
				if (data->fa[i * ncell_x + j] > cpuFA_max)
					{
					cpuFA_max = data->fa[i * ncell_x + j];
					FAindex = i * ncell_x + j;
					}
			}
		}
		data->FA_max = cpuFA_max;
		data->outletcellidx = FAindex; // this is the outlet cell which will be maintained throughout the simulation
	} else // do it faster using GPU in all subsequent iterations
		{
			int fullsize = data->mapInfo.width * data->mapInfo.height;
			thrust::device_ptr<double> max_FA = thrust::device_pointer_cast(device->fa);
			FA_max = thrust::reduce(max_FA, max_FA + fullsize, (double) 0, thrust::maximum<double>());
			data->FA_max = FA_max;
		}

	fprintf(data->outlog, "FA: Maximum FA is  %.6f s\n\n", data->FA_max);
	fprintf(data->outlog, "FA: Outletcell index is  %d \n\n", data->outletcellidx);

	printf("Maximum FA is  %.6f s\n\n", data->FA_max);

	for (int i = 0; i < ncell_y; i++) {
		for (int j = 0; j < ncell_x; j++) {
				if (data->fa[i * ncell_x + j] < 0) {
				//if (count < 20)
					//printf("[%d,%d] = %8.7f\n", i, j, data->fa[i * ncell_x + j]);
				count++;
				}
		}
	}
	fprintf(data->outlog, "FA: Bad value count (i.e. not in catchment(s) = %d\n", count);

}
