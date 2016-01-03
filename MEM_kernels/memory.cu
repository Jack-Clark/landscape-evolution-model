
#include "memory.h"


int createcontribAspace(Data* data)
{
	int fullsize;
	int dataSize;
	fullsize =  data->mapInfo.width * data->mapInfo.height;
	dataSize = fullsize * sizeof(int);
	data->contribA = (int *) malloc(dataSize);
	fprintf(data->outlog,"Host memory allocation for contribA  \n");
	return 0;
}


int clearcontribAspace(Data* data)
{
	free(data->contribA);
	//free(data->watershed_id); // need to clear this?

	return 0;
}

int createProcessMatrices(Data* data)
{
  int fullsize;
  int dataSize;
  int dataSizeInt;
  fullsize =  data->mapInfo.width * data->mapInfo.height;
  dataSize = fullsize * sizeof(double);

// these are the static grids in which data is stored from one iteration to the next ie. these are ONLY freed at the end of the simulation

  checkCudaErrors(cudaMallocHost((void **)&data->fd, sizeof(int) * data->mapInfo.height * data->mapInfo.width));
  fprintf(data->outlog, "Flow direction space on host allocated \n");

  checkCudaErrors(cudaMallocHost((void **)&data->fa,           dataSize));
  fprintf(data->outlog, "Flow accumulation space on host allocated \n");

  checkCudaErrors(cudaMallocHost((void **)&data->SlopePtr,     dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->runoffweight, dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->stonePtr,     dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->finesPtr,     dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->soilMPtr,     dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->soilBPtr,     dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->soilTPtr,     dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->nutPtr,       dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->TotBPtr,      dataSize));

  checkCudaErrors(cudaMallocHost((void **)&data->eroPtr,       dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->geliPtr,      dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->inciPtr,      dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->depoPtr,      dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->dz,           dataSize));

  checkCudaErrors(cudaMallocHost((void **)&data->weatherC,     dataSize));
  checkCudaErrors(cudaMallocHost((void **)&data->weatherP,     dataSize));

  fprintf(data->outlog, "All hosts matrices memory allocated \n");

  return 0;

}


int deleteProcessMatrices(Data* data)
{

	  checkCudaErrors(cudaFreeHost(data->dem));
	  checkCudaErrors(cudaFreeHost(data->fd));
	  checkCudaErrors(cudaFreeHost(data->fa));
	  checkCudaErrors(cudaFreeHost(data->SlopePtr));

	  checkCudaErrors(cudaFreeHost(data->runoffweight));
	  checkCudaErrors(cudaFreeHost(data->stonePtr));
	  checkCudaErrors(cudaFreeHost(data->finesPtr));
	  checkCudaErrors(cudaFreeHost(data->soilMPtr));
	  checkCudaErrors(cudaFreeHost(data->soilBPtr));
	  checkCudaErrors(cudaFreeHost(data->soilTPtr));
	  checkCudaErrors(cudaFreeHost(data->nutPtr));
	  checkCudaErrors(cudaFreeHost(data->TotBPtr));

	  checkCudaErrors(cudaFreeHost(data->eroPtr));
	  checkCudaErrors(cudaFreeHost(data->geliPtr));
	  checkCudaErrors(cudaFreeHost(data->inciPtr));
	  checkCudaErrors(cudaFreeHost(data->depoPtr));
	  checkCudaErrors(cudaFreeHost(data->dz));
	  checkCudaErrors(cudaFreeHost(data->weatherC));
	  checkCudaErrors(cudaFreeHost(data->weatherP));

	  fprintf(data->outlog, "All hosts matrices memory freed \n");

	return 0;
}


//////////////////////////////////////////////////////////////////////////////
// Setup store for catchment data ( needed for summary outputs etc)
//////////////////////////////////////////////////////////////////////////////

int createCatchmentSpace(Data* data, Catchment* Catchments) {
	//allocate space for catchment data and selective list and set values to zero
	Catchments->watershed_id = (int *) calloc(sizeof(int) , data->mapInfo.height * data->mapInfo.width);
	Catchments->mask = (int *) calloc(sizeof(int),  data->mapInfo.height * data->mapInfo.width); // all mask values set to zero

	fprintf(data->outlog, "Catchment space allocated \n");
	return 0;
}


