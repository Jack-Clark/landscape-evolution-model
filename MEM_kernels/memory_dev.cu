
#include "memory_dev.h"

void setdevicespace_FD(Data* data, Data* device)
{
	 size_t freenow, total;
	 int fullsize;
	 int ncell_x = data->mapInfo.width;
	 int ncell_y = data->mapInfo.height;
	 fullsize= ncell_x * ncell_y;

	  cudaMalloc((void**) &device->fd, fullsize * sizeof(int));
	  cudaMalloc((void **)&(device->dx), 9 * sizeof(int));
	  cudaMalloc((void **)&(device->dy), 9 * sizeof(int));
	  cudaMalloc((void **)&(device->shortest_paths), fullsize * sizeof(float));
	  cudaMalloc((void**)&(device->lowHeight), fullsize * sizeof(double));
	  cudaMalloc((void **) &device->watershed_id, ncell_x * ncell_y * sizeof(int));

	  cudaMemGetInfo(&freenow, &total);
	  fprintf(data->outlog, "Memory on CUDA card free after FD space allocated: %d total: %d \n",freenow/1024,total/1024);
	  fprintf(data->outlog, "FD: setdevicespace:%s\n", cudaGetErrorString(cudaGetLastError()));

}


void cleardevicespace_FD(Data* data, Data* device)
{
	size_t freenow, total;

		cudaFree(device->fd);
		cudaFree(device->dx);
		cudaFree(device->dy);
		cudaFree(device->shortest_paths);
		cudaFree(device->lowHeight);
		cudaFree(device->watershed_id);
		fprintf(data->outlog, "FD: error after FD clear :%s\n", cudaGetErrorString(cudaGetLastError()));

		cudaMemGetInfo(&freenow, &total);
		fprintf(data->outlog, "FD: Memory on CUDA card free after FD space freed: %d total: %d \n\n",freenow/1024,total/1024);

}

void setdevicespace_FA(Data* data, Data* device)
{
	int full_size;
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	full_size= ncell_x * ncell_y;

	cudaMalloc( (void**) &device->runoffweight, full_size * sizeof(double));
	cudaMalloc( (void**) &device->fa, full_size * sizeof(double));
	cudaMalloc((void**) &device->fd, full_size * sizeof(int));

	cudaMalloc( (void**) &device->stonePtr, full_size * sizeof(double));
	cudaMalloc( (void**) &device->TotBPtr, full_size * sizeof(double));
	cudaMalloc( (void**) &device->soilMPtr, full_size * sizeof(double));

	// now copy the necessary data - these will not overlap becasue they are all on the same stream

	//checkCudaErrors(cudaSetDevice(0));
	//checkCudaErrors( cudaMemcpy( device->fa, data->fa, full_size * sizeof(double), cudaMemcpyHostToDevice)) ;
	checkCudaErrors( cudaMemcpy( device->fd, data->fd, full_size * sizeof(int), cudaMemcpyHostToDevice)) ;
	checkCudaErrors( cudaMemcpy( device->runoffweight, data->runoffweight, full_size * sizeof(double), cudaMemcpyHostToDevice)) ;
	checkCudaErrors( cudaMemcpy( device->stonePtr, data->stonePtr, full_size * sizeof(double), cudaMemcpyHostToDevice)) ;
	checkCudaErrors( cudaMemcpy( device->TotBPtr, data->TotBPtr, full_size * sizeof(double), cudaMemcpyHostToDevice)) ;
	checkCudaErrors( cudaMemcpy( device->soilMPtr, data->soilMPtr, full_size * sizeof(double), cudaMemcpyHostToDevice)) ;

	fprintf(data->outlog, "FA: setdevicespace_FA:%s\n", cudaGetErrorString(cudaGetLastError()));
}

void cleardevicespace_FA(Data* data, Data* device)
{
	size_t freenow, total;

	cudaFree(device->fd);
	cudaFree(device->runoffweight);
	cudaFree(device->fa);

	cudaFree(device->contribA); // free it here as it is no longer needed)

	cudaMemGetInfo(&freenow, &total);
	fprintf(data->outlog, "FA: Memory on CUDA card free after FA space freed: %d total: %d \n\n",freenow/1024,total/1024);
}


void setdevicespace_Process(Data* data, Data* device)
{
	size_t freenow, total;
	int full_size;
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	full_size= ncell_x * ncell_y;



		cudaMalloc( (void**) &device->fa,       full_size * sizeof(double));
		cudaMalloc( (void**) &device->fd,       full_size * sizeof(int));
		cudaMalloc( (void**) &device->dz,       full_size * sizeof(double)); // create room for product dz
		cudaMalloc( (void**) &device->finesPtr, full_size * sizeof(double));
		cudaMalloc( (void**) &device->soilTPtr, full_size * sizeof(double));
		cudaMalloc( (void**) &device->nutPtr,   full_size * sizeof(double));
		cudaMalloc( (void**) &device->soilBPtr, full_size * sizeof(double));
		cudaMalloc( (void**) &device->eroPtr,   full_size * sizeof(double));
		cudaMalloc( (void**) &device->geliPtr,  full_size * sizeof(double));
		cudaMalloc( (void**) &device->inciPtr,  full_size * sizeof(double));
		cudaMalloc( (void**) &device->depoPtr,  full_size * sizeof(double));
		cudaMalloc( (void**) &device->weatherC, full_size * sizeof(double));
		cudaMalloc( (void**) &device->weatherP, full_size * sizeof(double));

		fprintf(data->outlog, "MOD: setdevicespace_Process :%s\n", cudaGetErrorString(cudaGetLastError()));


		// stones, TotBio, soilM plus dem, slope and mask still on device
		checkCudaErrors( cudaMemcpy ( device->fa,       data->fa,         full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->fd,       data->fd,         full_size * sizeof(int),    cudaMemcpyHostToDevice) );

		checkCudaErrors( cudaMemcpy ( device->SlopePtr,  data->SlopePtr,  full_size * sizeof(double), cudaMemcpyHostToDevice) );

		checkCudaErrors( cudaMemcpy ( device->finesPtr, data->finesPtr,   full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->soilTPtr, data->soilTPtr,   full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->nutPtr,   data->nutPtr,     full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->soilBPtr, data->soilBPtr,   full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->eroPtr,   data->eroPtr,     full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->geliPtr,  data->geliPtr,    full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->inciPtr,  data->inciPtr,    full_size * sizeof(double), cudaMemcpyHostToDevice) );
		checkCudaErrors( cudaMemcpy ( device->depoPtr,  data->depoPtr,    full_size * sizeof(double), cudaMemcpyHostToDevice) );

		checkCudaErrors( cudaMemcpy ( device->SlopePtr,  data->SlopePtr,  full_size * sizeof(double), cudaMemcpyHostToDevice) );

		fprintf(data->outlog, "MOD: Matrix memcopy operations :%s\n", cudaGetErrorString(cudaGetLastError()));

		cudaMemGetInfo(&freenow, &total);
		fprintf(data->outlog, "MOD: Memory on CUDA card free after model matrix space allocated: %d total: %d \n",freenow/1024,total/1024);
}

void cleardevicespace_Process(Data* data, Data* device)
{
	size_t freenow, total;

	cudaFree(device->fa);
	cudaFree(device->fd);

	cudaFree(device->dz);

	cudaFree(device->finesPtr);
	cudaFree(device->soilTPtr);
	cudaFree(device->nutPtr);
	cudaFree(device->soilBPtr);
	cudaFree(device->eroPtr);
	cudaFree(device->geliPtr);
	cudaFree(device->inciPtr);
	cudaFree(device->depoPtr);
	cudaFree(device->weatherC);
	cudaFree(device->weatherP);

	// free after being left at end of FA routines.
	cudaFree(device->stonePtr);
	cudaFree(device->TotBPtr);
	cudaFree(device->soilMPtr);

	cudaMemGetInfo(&freenow, &total);
	fprintf(data->outlog, "MOD: Memory on CUDA card free after model space freed: %d total: %d \n",freenow/1024,total/1024);
}

int copyMask(Data* data, Data* device)
{

	 int fullsize;
	 int ncell_x = data->mapInfo.width;
	 int ncell_y = data->mapInfo.height;
	 fullsize= ncell_x * ncell_y;

	 cudaMalloc( (void**) &device->mask, fullsize * sizeof(int)); // create space for the mask

	 cudaMemcpy(device->mask, data->mask, fullsize * sizeof(int), cudaMemcpyHostToDevice);  // copy back flag
	 fprintf(data->outlog, "Mask data sent to device \n");

	thrust::device_ptr<int> activecells = thrust::device_pointer_cast(device->mask);
	data->activecells  = thrust::count(activecells, activecells + fullsize, 1);
	printf("No of active cells = %d \n", data->activecells);

	return 0;
}



int createDeviceSpace(Data* data, Data* device)
{
	size_t freenow, total;

	 int fullsize;
	 int ncell_x = data->mapInfo.width;
	 int ncell_y = data->mapInfo.height;
	 fullsize= ncell_x * ncell_y;

	  cudaMalloc((void **)&(device->dem), fullsize* sizeof(double));
	  cudaMalloc((void **)&(device->SlopePtr), fullsize * sizeof(double));

	  cudaMalloc((void **)&(device->summary), fullsize * sizeof(double));

	fprintf(data->outlog,"Allocated DEM and slope matrices on device :%s\n", cudaGetErrorString(cudaGetLastError()));

	cudaMemGetInfo(&freenow, &total);
	fprintf(data->outlog,"Memory on CUDA card free after device DEM and slope grids allocated: %d total: %d \n",freenow/1024,total/1024);

	printf("Device space created \n");
	return 0;
}


int clearDeviceSpace(Data* data, Data* device)
{
	size_t freenow, total;

		cudaFree(device->dem);
		cudaFree(device->SlopePtr);
		cudaFree(device->summary);


	cudaMemGetInfo(&freenow, &total);
	printf("Memory on CUDA card free after DEM and slope device grids space freed: %d total: %d \n",freenow/1024,total/1024);
	fprintf(data->outlog,"Memory on CUDA card free after DEM and slope device grids space freed: %d total: %d \n",freenow/1024,total/1024);


	free(data->watershed_id);


	return 0;
}



int zerogrids(Data* data)
{

	memset(data->eroPtr, 0.0, sizeof(data->eroPtr));
	memset(data->geliPtr, 0.0, sizeof(data->eroPtr));
	memset(data->inciPtr, 0.0, sizeof(data->inciPtr));
	memset(data->depoPtr, 0.0, sizeof(data->depoPtr));

	return 0;
}
