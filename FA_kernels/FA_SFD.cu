
#include "FA_SFD.h"

void correctflow_SFD(Data* data, Data* device, int iter) {
	if (cudaSuccess != cudaSetDevice(CUDA_DEVICE)) {
		printf("Unable to access CUDA card\n");
		return ;
	}

	int x;

	cudaEvent_t start, stop;
	float time;
    int rows = data->mapInfo.height;
	int cols = data->mapInfo.width;
	int totalCells = rows * cols;

	for (x = 0; x < rows * cols; ++x) {
		data->fa[x] = -1.0;
	}

	fprintf(data->outlog, "FA: set fagrid values to -1\n");
	cudaMemcpy(device->fa, data->fa, totalCells * sizeof(double), cudaMemcpyHostToDevice);
	fprintf(data->outlog, "FA: FA memcopy operation :%s\n", cudaGetErrorString(cudaGetLastError()));

	fprintf(data->outlog, "FA: Calling process\n");

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	// Change this function to whichever version of the algorithm you want to run.
	process_SFD_Multiple_Retries(data, device, iter);

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	printf("Time to complete FA_SFD_list : %.6f s\n", time / 1000.0);
}
