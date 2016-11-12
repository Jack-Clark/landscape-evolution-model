
#include "FA_SFD.h"

int correctflow_SFD(Data* data, Data* device, int iter, int algorithmID) {
	if (cudaSuccess != cudaSetDevice(CUDA_DEVICE)) {
		printf("Unable to access CUDA card\n");
		return 1;
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

	switch(algorithmID) {

		case 1:
			process_SFD_NoPart_List(data, device, iter);
			break;
	
		case 2:
			process_SFD_block_level_single_chains(data, device, iter);
			break;

		case 3:
			process_SFD_global_level_single_chains(data, device, iter);
			break;

		case 4:
			process_SFD_Multiple_Retries(data, device, iter);
			break;

		default:
			fprintf(data->outlog, "\nInvalid algorithmID. Exiting...\n");
			printf("\nInvalid algorithmID. Exiting...\n");
			return 1;
	}

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	printf("Time to complete FA_SFD_list : %.6f s\n", time / 1000.0);

	return 0;
}
