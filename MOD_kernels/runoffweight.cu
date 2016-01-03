
#include "runoffweight.h"
/*!

  Calculates the actual runoff from individual cell in m.

  OLD VERSION:: runoff_coeff = (1.0 - 0.0085 * (stonePtr[self])) * exp (-0.0025 * (TotBPtr[self]) / 76.5)	* tanh (0.6 + (soilMPtr[self]) / 0.2);

  runoff_coeff = (1.0 - 0.0085 * (stonePtr[self])) * exp (-0.0025 * (TotBPtr[self]) / 76.5)	; // newer version is simplified?

  runoff = ppn * runoff_coeff;  ppn is depth of rain in m, runoff is therefore a depth in m

*/



__global__ void calcrunoffweight(int ncell_x, int ncell_y, double ppn, int* mask, double* stonePtr, double* TotBPtr, double* soilMPtr, double* runoffweight)
{
  double runoff = 0.0;
  double runoff_coeff;

  int irow, icol;
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;

  // If outside of DEM return (nothing to do)
  if(icol >= ncell_x || irow >= ncell_y)
    return;

  int self = irow * ncell_x + icol;
	  if (mask[self] == 0)
		  {
		  runoffweight[self] = 0.;
		  return; // don't calculate if not in catchment(s) of interest
		  }

  runoff_coeff = (1.0 - 0.0085 * (stonePtr[self])) * exp (-0.0025 * (TotBPtr[self]) / 76.5)	* tanh (0.6 + (soilMPtr[self]) / 0.2);

 // runoff_coeff = (1.0 - 0.0085 * (stonePtr[self])) * exp (-0.0025 * (TotBPtr[self]) / 76.5)	; // newer version is simplified?


  runoff = ppn * runoff_coeff; // ppn is depth of rain in m, runoff is therefore a depth in m

  runoffweight[self] = runoff * 0.04; //actual runoff depth per cell in m


  soilMPtr[self] += ppn * (1. - runoff_coeff);

	if ((soilMPtr[self]) > 0.5)
				{
					soilMPtr[self] = 0.5;
				}
	if (runoffweight[self] == 0) printf("weight = 0 IN CATCHMENT");

}



int computeRunOffWeights(Data* data, Data* device)
{
	double currentRain;
	double ppn;

	data->cellarea = data->mapInfo.cellsize * data->mapInfo.cellsize ;


	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	int fullsize = ncell_x * ncell_y;

	dim3  threads( BLOCKCOLS, BLOCKROWS);
	dim3  grid(ncell_x/BLOCKCOLS, ncell_y/ BLOCKROWS);

	currentRain = data->rain; // rain is in mm
	ppn =  currentRain/1000  ;  // depth per cell scaled now in m

	fprintf(data->outlog, "FA: ppn (mm) = %12.8lf \n", ppn);
	

	calcrunoffweight<<<grid, threads>>>(ncell_x, ncell_y, ppn, device->mask, device->stonePtr, device->TotBPtr, device->soilMPtr, device->runoffweight);
	fprintf(data->outlog, "FA: calcrunoffweight:%s\n", cudaGetErrorString(cudaGetLastError()));

	thrust::device_ptr<double> runoff_d = thrust::device_pointer_cast(device->runoffweight);
	thrust::device_ptr<double> summary_d = thrust::device_pointer_cast(device->summary);

	//predicate function is_not_zero defined in header
	thrust::copy_if(runoff_d, runoff_d + fullsize, summary_d, is_not_zero());

	double maxRO;
	double minRO;

	maxRO = thrust::reduce(summary_d, summary_d + data->activecells, (double) 0, thrust::maximum<double>());
	minRO = thrust::reduce(summary_d, summary_d + data->activecells, (double) 10, thrust::minimum<double>());

	fprintf(data->outlog, "FA: Max RunOff: %f,  Min RunOff: %f \n", maxRO, minRO);

	cudaMemcpy(data->runoffweight, device->runoffweight, sizeof(double)* ncell_x* ncell_y, cudaMemcpyDeviceToHost); // here just for checking
	fprintf(data->outlog, "FA: runoff weights calculated :%s\n", cudaGetErrorString(cudaGetLastError()));

	thrust::fill(summary_d, summary_d + data->activecells, 0.0); // reset the summary grid

	return 1;
}
