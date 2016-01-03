// Simple Process Kernels for calculating diffuse and concentrated erosion, returning them in erosion and incision matrices

#include "eroincidep.h"



__global__ void diff_erosion(int *mask, double *ePtr_d, double *FAPtr_d, double *SlopePtr_d, double *finesPtr_d, double *nutPtr_d,  double *soilTPtr_d, double *stonePtr_d, double *TotBPtr_d,
				int soil_struct, int profile_permeability, int ncell_x, int ncell_y)
{

	/*!
	  Diffuse erosion by unconcentrated overland flow
	  erodibility is the k factor in the USLE.  Need to check whether a zero value is likely
	  The cell flow is the FA corrected for cell size (100. / cell_size) * (100. / cell_size). It is not clear whether this is really appropriate
	  This is based on a modified Musgrave approach. Calculation is based upon Mitchell and Bubenzer, 1980, Zhang et al., 2002
	  See Wainwright 2005 p.106
	  OF is runoff depth in mm and Slope is surface slope in m per m. according to Wainwright
	*/



	int icol = blockIdx.x * blockDim.x + threadIdx.x; //column
	int irow = blockIdx.y * blockDim.y + threadIdx.y; //row

	// simple out of bounds check
	if( icol >= ncell_x  || irow >= ncell_y )
		return;

	int self = irow * ncell_x + icol;
	if (mask[self] == 0)
		{
		ePtr_d[self] = 0.0;
		return; // don't calculate if not in catchment(s) of interest
		}


	double erodibility = (2.77e-6 * pow (finesPtr_d[self], 1.14) * (12. - nutPtr_d[self])) + (0.043 * (soil_struct - 2.)) + (0.033 * (profile_permeability - 3.));
	// both soil_strict and profile_permeability are unused i.e. both are zero at all times

	if (erodibility < 0.0) erodibility = 0.0;
	double OF = (FAPtr_d[self]) * 1000 ; // now in mm;
	//double threshold = 600. ;//* (100./ncell_x) * (100./ncell_x);
	double threshold = 6. ;//* (100./ncell_x) * (100./ncell_x);

	if (OF > threshold)            // limit to when concentrated erosion sets in
		OF = threshold;

	double slope = (SlopePtr_d[self]);
	//if (slope > 2.0) slope = 2.0;

	double answer = ((  erodibility * OF * OF * pow (slope, 1.67 ) *
						(pow ((110. - (stonePtr_d[self])), 0.877) / 61.7023))) *
								exp (-0.07 * (TotBPtr_d[self]) / 76.5) ;


	if (answer > (soilTPtr_d[self])) answer = soilTPtr_d[self];
	ePtr_d[self] = answer;

}

__global__ void conc_erosion( int *mask, double *inciPtr, double *FAPtr, double *SlopePtr, double *TotBPtr, double *soilTPtr, int ncell_x, int ncell_y)
{
	/*!
	 * Concentrated erosion by channelized flow [mm]. Simuluated using a modified stream-power approach
	 * see Eq [3] of Wainwright 1995 p. 106
	 * the concentrated erosion parameter k2 set at 0.00005 seems low - need to check
	 */
	int icol = blockIdx.x * blockDim.x + threadIdx.x; //column
	int irow = blockIdx.y * blockDim.y + threadIdx.y; //row

	// simple out of bounds check
	if( icol >= ncell_x  || irow >= ncell_y )
		return;

	int self = irow * ncell_x + icol;
	if (mask[self] == 0)
		{
		inciPtr[self] = 0.0;
		return; // don't calculate if not in catchment(s) of interest
		}

	//double flow_thresh = 6. ;// * (100./ncell_x) * (100./ncell_x);
	double flow_thresh = 6. ;// * (100./ncell_x) * (100./ncell_x);
	double flow_A = 0.00035;  // value for k2 concentrated erosion erodibility - from what is value derived? changed from 0.0005 7/12/15
	double flow_B = 2.;
	//double max_incision = 0.005;
	double max_incision = SlopePtr[self] * 20 ; //note 20 here is cell size

	//if (max_incision > 0.01) max_incision = 0.01;

	if (max_incision > 0.0004) max_incision = 0.0004;// from 0.002

	//double f_temp = (( FAPtr[self] ) * 1000 / 25)  - flow_thresh;
	double f_temp =  (FAPtr[self] * 1000) - flow_thresh;
	
	if (f_temp < 0.0)
	{
		f_temp = 0.0;
		inciPtr[self] = 0.0;
	}
	else
	{
		//    but still account for vegetation effects
		if (SlopePtr[self] == 0.0)
			inciPtr[self]= flow_A * pow (f_temp, flow_B) * 0.01 * exp (-0.07 * (TotBPtr[self]) / 76.5);
		else
			inciPtr[self] = flow_A * pow (f_temp, flow_B) * (SlopePtr[self]) * exp (-0.07 * (TotBPtr[self]) / 76.5);

		//if ((inciPtr[self]) > (max_incision + ( soilTPtr[self])) )	inciPtr[self] = (max_incision + (soilTPtr[self]));

		if ( inciPtr[self] > max_incision)	inciPtr[self] = max_incision;
	}
}

__global__ void gelifluction( int *mask, double* geliPtr, double temperature, double *eroPtr, double *soilTPtr, double *SlopePtr, int ncell_x, int ncell_y, double *FA, double FAMax)
{
	/*! JW solifluxion component using version of Anderson (Culling) model Q=-k dz/dx where k=0.5*beta*f*zeta
	JW beta base value = 0.01-0.03 --> 0.02 (Anderson pers comm),
	JW f = 0.9236-0.0601*MAT-0.005MAT^2 (Matsuoka)
	JW zeta = 17.011-0.1745*MAT-0.0487MAT^2 (Grossi et al.)
	JW multiplying together gives fourth-order polynomial between -20�C<=MAT<=8.75�C, k set to 0.001 outside this range
    */

	int icol = blockIdx.x * blockDim.x + threadIdx.x; //column
	int irow = blockIdx.y * blockDim.y + threadIdx.y; //row

	// simple out of bounds check
	if( icol >= ncell_x  || irow >= ncell_y )
		return;

	int self = irow * ncell_x + icol;
	if (mask[self] == 0)
		{
		geliPtr[self] = 0.0;
		return; // give nodatavalue if not in catchment(s) of interest
		}

	double kappa = 0.00001;

	if (temperature >= -20. && temperature <= 8.87)
	{
		kappa = 0.15711 - 0.011835 * temperature - 0.001195 * temperature * temperature + 0.000038 * temperature * temperature * temperature  + 0.000002 * temperature * temperature * temperature * temperature;
	}

	double slopethreshold;
	double scaler = 1 - (FA[self] / 600);
	if (FA[self] > 600) scaler = 0.0; // turn gelifluction off in channel cells.

	slopethreshold = SlopePtr[self];
	if (slopethreshold > 1.) slopethreshold = 1.0; // placed this constraint back in

	geliPtr[self] = kappa * slopethreshold * scaler; // * sfx_d;  //JW sfx is the scaling parameter passed via the parameter file or defaulting to 1. if not present

	if ((eroPtr[self] + geliPtr[self]) > soilTPtr[self])
	{
		geliPtr[self] = soilTPtr[self] - eroPtr[self];
	}

	if (geliPtr[self] < 0.0) geliPtr[self] = 0.0;
}


__global__ void accumulate(int *mask, int* AspectPtr, double* eroPtr, double* geliPtr, double* inciPtr, double* depoPtr, int *ok, int *secondOK, int *progressed, int rows, int cols, double ddk, double dck, double dgk,
	   double* finesPtr, double* stonePtr, double* nutPtr, double* soilTPtr )
{
	/*!
	 * The deposition routine uses a version of correct flow but pulls sediment from the feeding cells
	 * Each cells waits for the surrounding cells to be calculated and then pulls the appropriate amount of
	 * diffuse and concentrated erosion dependent upon the ddk and dck values (because they are pulling this is 1-ddk and 1-dck
	 * By recording which cell has supplied [from] the calibre mix of stones and fines can be updated proportionally
	 *
	 */
  int from;
  int irow = blockIdx.y * blockDim.y + threadIdx.y;
  int icol = blockIdx.x * blockDim.x + threadIdx.x;

  if (irow >= rows || icol >= cols) // || irow < 0 || icol < 0)
    return;

	int self = irow * cols + icol;
	if (mask[self] == 0)
		{
		depoPtr[self] = 0.0;
		return; // don't calculate if not in catchment(s) of interest
		} // don't calculate if not in catchment(s) of interest

  from = self; // for cells which receive no additions

  int nie, nise, nis, nisw, niw, ninw, nin, nine;

  //self = blockIdx.y*cols*BLOCKROWS + blockIdx.x*BLOCKCOLS + threadIdx.y*cols + threadIdx.x;
  //ok[self] = 1;
  //secondOK[self] = 2;
  //return;
  if (ok[self] != 0) return;

  if (secondOK[self] == 1) {
    ok[self] = 1;
    *progressed = 1;
    return;
  }

  double diffuse = 0. ;
  double conc =  0. ;
  double geli = 0. ;

  nie  = self        + 1 ;
  nise = self + cols + 1 ;
  nis  = self + cols     ;
  nisw = self + cols - 1 ;
  niw  = self        - 1 ;
  ninw = self - cols - 1 ;
  nin  = self - cols     ;
  nine = self - cols + 1 ;


    if (icol < cols - 1 && (AspectPtr[nie] & WEST)) {
        if (!ok[nie]) return; // check to see if the cell from which you receive sediment is OK i.e. has a finished computation if not return

        diffuse += eroPtr[nie] * (1-ddk) ;
		conc += inciPtr[nie] * (1-dck) ;
		geli += geliPtr[nie] * (1-dgk);
		from = nie;
	}
    if (icol < cols - 1 && irow < rows - 1 && (AspectPtr[nise] & NORTHWEST)) {
        if (!ok[nise]) return;

        diffuse += eroPtr[nise] * (1-ddk);
		conc += inciPtr[nise] * (1-dck);
		geli += geliPtr[nie] * (1-dgk);
		from = nise;
    }
    if (irow < rows - 1 && (AspectPtr[nis] & NORTH)) {
        if (!ok[nis]) return;

        diffuse += eroPtr[nis] * (1-ddk);
		conc += inciPtr[nis] * (1-dck);
		geli += geliPtr[nie] * (1-dgk);
		from = nis;
    }
    if (icol > 0 && irow < rows - 1 && (AspectPtr[nisw] & NORTHEAST)) {
        if (!ok[nisw]) return;

        diffuse += (eroPtr[nisw] * (1-ddk) );
		conc += (inciPtr[nisw] * (1-dck));
		geli += geliPtr[nie] * (1-dgk);
		from = nisw;
    }
    if (icol > 0 && (AspectPtr[niw] & EAST)) {
        if (!ok[niw]) return;

        diffuse += (eroPtr[niw] * (1-ddk));
		conc += (inciPtr[niw] * (1-dck));
		geli += geliPtr[nie] * (1-dgk);
		from = niw;
    }
    if (icol > 0 && irow > 0 && (AspectPtr[ninw] & SOUTHEAST)) {
        if (!ok[ninw]) return;

        diffuse += (eroPtr[ninw] * (1-ddk));
		conc += (inciPtr[ninw] * (1-dck) );
		geli += geliPtr[nie] * (1-dgk);
		from = ninw;
    }
    if (irow > 0 && (AspectPtr[nin] & SOUTH)) {
        if (!ok[nin]) return;

        diffuse += (eroPtr[nin] * (1-ddk) );
		conc += (inciPtr[nin] * (1-dck) );
		geli += geliPtr[nie] * (1-dgk);
		from = nin;
    }
    if (irow > 0 && icol < cols - 1 && (AspectPtr[nine] & SOUTHWEST) ) {
        if (!ok[nine]) return;

        diffuse += (eroPtr[nine] * (1-ddk));
		conc += (inciPtr[nine] * (1-dck));
		geli += geliPtr[nie] * (1-dgk);
		from = nine;
    }

	//depoPtr[self] = ((eroPtr[self] * ddk) + (inciPtr[self] * dck)) + diffuse + conc;

	// I think this should be :
	depoPtr[self] = ((eroPtr[self]+diffuse) * ddk) + ((inciPtr[self]+conc) * dck) + ((geliPtr[self]+geli) * dgk);

    secondOK[self] = 1;
    *progressed = 1;


    // replace original deposit routine to calculate new proportions of fines, stones and nutrients
	double weight = (soilTPtr[self]) + (depoPtr[self]) ;
	if (weight != 0.0)
	{
		double old = (soilTPtr[self]) / weight ;
		double input = (depoPtr[from]) /weight;
		stonePtr[self] =   old * (stonePtr[self])   + input * (stonePtr[from]) ;
		if (stonePtr[self] < 0.01) stonePtr[self] = 0.01;
		else if (stonePtr[self] > 99.9) stonePtr[self] = 99.9;

		finesPtr[self] = 100 - stonePtr[self];  // now same as calculation in weathering
		//finesPtr[self] = old * (finesPtr[self]) + input * (finesPtr[from]);

		nutPtr[self] =   old * (nutPtr[self])   + input * (nutPtr[from]) ;
	}

}

__global__ void surface_change (int *mask, double* dz_d, double* eroPtr_d, double* inciPtr_d, double* depoPtr_d, double* geliPtr_d, int rows, int cols)
{

  int irow = blockIdx.y * blockDim.y + threadIdx.y;
  int icol = blockIdx.x * blockDim.x + threadIdx.x;

  if (irow >= rows || icol >= cols) // || irow < 0 || icol < 0)
    return;

  int self = irow * cols + icol;
	if (mask[self] == 0)
		{
		dz_d[self] = 0.0;
		return;
		} // give nodatavalue if not in catchment(s) of interest

  dz_d[self] = (depoPtr_d[self] - eroPtr_d[self] - inciPtr_d[self] - geliPtr_d[self]) ;
}

__global__ void weathering (int *mask, double temperature, double rain, double *FAPtr_d, int type, double *weatherC_d, double *weatherP_d, double *soilTPtr_d,
	                                 double* finesPtr_d, double* stonePtr_d, double* soilMPtr_d, double* eroPtr_d, double* dz_d, double cell_size, int rows, int cols)
{

  double surface_flow;
  double evap_trans;
  int irow = blockIdx.y * blockDim.y + threadIdx.y;
  int icol = blockIdx.x * blockDim.x + threadIdx.x;

  if (irow >= rows || icol >= cols) // || irow < 0 || icol < 0)
    return;

  int self = irow * cols + icol;


  if (mask[self] == 0) return; // don't calculate if not in catchment(s) of interest


  surface_flow = FAPtr_d[self] * 1000; // / 5;  // converted back to mm here
  if (rain + surface_flow <= PE_d)
			evap_trans = rain + surface_flow;
  else
			evap_trans = rain + surface_flow - PE_d;


	switch (type)
	{
        case 0:  // White regression method (see Dreybrodt)
	   	    weatherC_d[self] = 6.3 + 0.049 * (rain - evap_trans);
	   	    break;

        case 1:  // Dreybrodt as used by Kaufmann and Braun 2001 Terra Nova 13, 313-20.
			double PCO2;
			if (soilTPtr_d[self] > 0.1)
				PCO2 = pow (10., -2.5 + ((temperature - 15.) / 10.));
			else
				PCO2 = pow (10., -3.5 + ((temperature - 15.) / 10.));

			double Ca2eq = pow (PCO2 * ((K1_d * KC_d * KH_d) / (4. * K2_d * gammCa_d *  gammHCO3_d * gammHCO3_d)), (1. / 3.));

			//weatherC_d[self] = 0.0148148148 * Ca2eq * surface_flow / (cell_size * cell_size);
			weatherC_d[self] = 0.0148148148 * Ca2eq * surface_flow / 250000; //denominator now has same scalar as the original 100*100, changed from 10k on 30/11/15

			break;
	}
	if (weatherC_d[self]>0.001) weatherC_d[self]=0.001; // new constraint added 13/10/15

	weatherP_d[self] = P0_d * exp (-(P1_d) * soilTPtr_d[self]);

	// following makes the assumption of 2% solid residue from chemical
	//    weathering
	//double dz = depoPtr_d[self] - eroPtr_d[self] - inciPtr_d[self];
	soilTPtr_d[self] = soilTPtr_d[self] + dz_d[self] + weatherC_d[self] * 0.02 + weatherP_d[self];

	if (soilTPtr_d[self] <= 0.0)
	{
		soilTPtr_d[self] = 0.0;
        stonePtr_d[self] = 0.0;
        finesPtr_d[self] = 0.0;
        soilMPtr_d[self] = 0.02;
	}

	else
	{
		if (weatherP_d[self] > 0)
			stonePtr_d[self] += stonePtr_d[self] * (tanh (((eroPtr_d[self] / soilTPtr_d[self]) - 0.5) / 0.25) * -tanh (((weatherC_d[self] / weatherP_d[self]) - 1.) / 0.25));
		else
			stonePtr_d[self] += stonePtr_d[self] * tanh (((eroPtr_d[self] / soilTPtr_d[self]) - 0.5) / 0.25);

		if (stonePtr_d[self] < 0.01) stonePtr_d[self] = 0.01;
		else if (stonePtr_d[self] > 99.9) stonePtr_d[self] = 99.9;

		finesPtr_d[self] = 100. - stonePtr_d[self];





		soilMPtr_d[self] -= (1.e-3 * evap_trans / soilTPtr_d[self]); // 1.e-3 converts ET to m


		if (soilMPtr_d[self] < 0.02) soilMPtr_d[self] = 0.02;
        else if (soilMPtr_d[self] > 0.5) // assume constant porosity of 50% for the moment
		                   soilMPtr_d[self]  = 0.5;  // but don't feedback water sitting on surface at end of year yet

	}

	//weatherC_d[self] = 0.0;



}



void calc_diff_erosion(Data* data, Data* device)
{
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	int full_size = ncell_x * ncell_y;

	dim3  threads( 16, 16);
	dim3  grid(ncell_x/BLOCKCOLS, ncell_y/ BLOCKROWS);

	diff_erosion<<< grid, threads >>>( device->mask, device->eroPtr, device->fa, device->SlopePtr, device->finesPtr, device->nutPtr, device->soilTPtr, device->stonePtr, device->TotBPtr, data->soil_struct, data->profile_permeability, ncell_x, ncell_y);
	fprintf(data->outlog, "MOD: diff_erosion :%s\n", cudaGetErrorString(cudaGetLastError()));
	checkCudaErrors( cudaMemcpy ( data->eroPtr,   device->eroPtr,   full_size * sizeof(double), cudaMemcpyDeviceToHost) );

	thrust::device_ptr<double> difftot_d = thrust::device_pointer_cast(device->eroPtr);
	cudaSetDevice(0);
	data->totE = thrust::reduce(difftot_d, difftot_d + full_size, (double) 0);

	if (data->totE == NAN)
	{
		printf("data from totE is NaN \n");
		exit(0);
	}

}

void calc_conc_erosion(Data* data, Data* device)
{
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	int full_size = ncell_x * ncell_y;

	dim3  threads( 16, 16);
	dim3  grid(ncell_x/BLOCKCOLS, ncell_y/ BLOCKROWS);

	conc_erosion<<< grid, threads >>>( device->mask, device->inciPtr, device->fa, device->SlopePtr, device->TotBPtr, device->soilTPtr, ncell_x, ncell_y);
	fprintf(data->outlog, "MOD: conc_erosion :%s\n", cudaGetErrorString(cudaGetLastError()));
	checkCudaErrors( cudaMemcpy ( data->eroPtr,   device->eroPtr,   full_size * sizeof(double), cudaMemcpyDeviceToHost) );

	thrust::device_ptr<double> incitot_d = thrust::device_pointer_cast(device->inciPtr);
	cudaSetDevice(0);
	data->totI = thrust::reduce(incitot_d, incitot_d + full_size, (double) 0);
	//printf("total Incision from thrust is %10.8lf \n", data->totI);


	if (data->totI == NAN) {
		printf("data from totI is NaN \n");
		exit(0);
	}
}

void calc_gelifluction(Data* data, Data* device)
{
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	int full_size = ncell_x * ncell_y;

	const double sfxh = 1.0;
	const double* hsfx = &sfxh;

	if( !(cudaMemcpyToSymbol(sfx_d, hsfx,              sizeof(double), 0 , cudaMemcpyHostToDevice) == cudaSuccess)) printf ("failed \n");

	dim3  threads( 16, 16);
	dim3  grid(ncell_x/BLOCKCOLS, ncell_y/ BLOCKROWS);

	printf("FA_max before gelifluction %f \n", data->FA_max);

	//__global__ void gelifluction( int *mask, double* geliPtr, double temperature, double *eroPtr, double *soilTPtr, double *SlopePtr, int ncell_x, int ncell_y)
	gelifluction<<< grid, threads >>>( device->mask, device->geliPtr, data->temp, device->eroPtr, device->soilTPtr, device->SlopePtr,  ncell_x, ncell_y, device->fa, data->FA_max);

	fprintf(data->outlog, "MOD: gelifluction :%s\n", cudaGetErrorString(cudaGetLastError()));
	checkCudaErrors( cudaMemcpy ( data->geliPtr,   device->geliPtr,   full_size * sizeof(double), cudaMemcpyDeviceToHost) );

	thrust::device_ptr<double> gelitot_d = thrust::device_pointer_cast(device->geliPtr);
	cudaSetDevice(0);
	data->totG = thrust::reduce(gelitot_d, gelitot_d + full_size, (double) 0);
	//printf("total Incision from thrust is %10.8lf \n", data->totI);


	if (data->totG == NAN) {
		printf("data from totG is NaN \n");
		exit(0);
	}
}

void calc_sedflux(Data* data, Data* device)
{

  cudaEvent_t start, stop;
  float time;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  //  former process routine now merged
  int gridRows = data->mapInfo.height;
  int gridColumns = data->mapInfo.width;
  int full_size = gridColumns * gridRows;

  int *ok;
  int *secondOK;
  int *okGrid_h = (int*) malloc(sizeof(int) * gridRows * gridColumns);
  int *progressed;

  cudaMalloc((void **) &ok, gridRows * gridColumns * sizeof(int));
  cudaMalloc((void **) &secondOK, gridRows * gridColumns * sizeof(int));
  cudaMalloc((void **) &progressed, sizeof(int));

  for (int i = 0; i < gridRows * gridColumns; i++) {
    okGrid_h[i] = 0;
  }

// Copy data to card
  cudaMemcpy(ok, okGrid_h, gridRows * gridColumns * sizeof(int),   cudaMemcpyHostToDevice);
  cudaMemcpy(secondOK, okGrid_h, gridRows * gridColumns * sizeof(int), cudaMemcpyHostToDevice);

  int *progress = (int*) malloc(sizeof(int));
  *progress = 0;

  dim3  dimBlock( BLOCKCOLS, BLOCKROWS);
  dim3  dimGrid(gridColumns/BLOCKCOLS, gridRows/ BLOCKROWS);

  int consecutiveZero = 0;

 // transport distance exponents

 	//		double ddk = 0.95 * (23. / cell_size2); // travel distance parameter for diffuse flow -- normalized relative to 23 m plot
	//      double dck = 0.05 * (100. / cell_size2); // travel distance parameter for concentrated flow -- normalized relative to 100 m

  // is the aim is to export ~5% of diffuse and 95% of concentrated ?

  			double ddk = 0.95; // proportion deposited
    		double dck = 0.5; // proportion deposited
    		double dgk = 0.98; // solifluxion does not move with this parameter set to 1.

			if (ddk > 0.95)	ddk = 0.95;
			if (dck > 0.9)	dck = 0.9;

  do {

    *progress = 0;

	cudaMemcpy(progressed, progress, sizeof(int), cudaMemcpyHostToDevice);


	accumulate<<<dimGrid, dimBlock>>>(device->mask, device->fd, device->eroPtr, device->geliPtr, device->inciPtr, device->depoPtr,
										ok, secondOK, progressed, gridRows, gridColumns, ddk, dck, dgk, device->finesPtr, device->stonePtr, device->nutPtr, device->soilTPtr);

    cudaMemcpy(progress, progressed, sizeof(int), cudaMemcpyDeviceToHost);

		if (*progress == 0)
			{
   // printf("Progress was not made \n");
				consecutiveZero ++;
			}
		else
				consecutiveZero = 0;

	} while (consecutiveZero < 5);


  cudaMemcpy(data->depoPtr, device->depoPtr, gridRows * gridColumns * sizeof(double),   cudaMemcpyDeviceToHost);

// Free the CUDA copies of the package data
  cudaFree(ok);
  cudaFree(secondOK);
  cudaFree(progressed);

/* Free host copies of the package data */
  free(okGrid_h);
  free(progress);

  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  printf("time in sed accumulation %.6f s\n", time / 1000.0);

	thrust::device_ptr<double> deptot_d = thrust::device_pointer_cast(device->depoPtr);
	cudaSetDevice(0);
	data->totD = thrust::reduce(deptot_d, deptot_d + full_size, (double) 0);
	//printf("total Dep from thrust is %10.8lf \n", data->totD);
	if (data->totD == NAN) {
		printf("data from totD is NaN \n");
		exit(0);
	}

}

void calc_dz(Data* data, Data* device)
{
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;
	int full_size = ncell_x * ncell_y;

	dim3  threads( 16, 16);
	dim3  grid(ncell_x/BLOCKCOLS, ncell_y/ BLOCKROWS);

	// Calculate dz due to erosion and depositional processes
	surface_change<<<grid, threads>>>(device->mask, device->dz, device->eroPtr, device->inciPtr, device->depoPtr, device->geliPtr, data->mapInfo.height, data->mapInfo.width);
	fprintf(data->outlog, "MOD: surface_change :%s\n", cudaGetErrorString(cudaGetLastError()));

	//thrust::device_ptr<double> max_dz = thrust::device_pointer_cast(device->dz);
	//cudaSetDevice(0);


	//double* maxdz = thrust::max_element(thrust::device, device->dz, device->dz + full_size );
	//printf("\nMaximum dz = %10.5lf \n", *maxdz);



}

void calc_weathering(Data* data, Data* device)
{
	int ncell_x = data->mapInfo.width;
	int ncell_y = data->mapInfo.height;

	double temperature = data->temp ; // current iteration temperature
	double rain = data->rain ; // current iteration rainfall

	const double PEh = 1.44 * pow (temperature + 6., 2);;
	const double P0h = 7.7e-5;
	const double P1h = 2.3 ;

	const double* hPE = &PEh;
	const double* hP0 = &P0h;
	const double* hP1 = &P1h;

	// Based upon Dreybrodt as used by Kaufmann and Braun 2001 Terra Nova 13, 313-20.
	// Do these calculations once before entering kernel and transfer values into shared memory
    double t_abs = temperature + 273.16;
    double t2 = t_abs * t_abs;
    double ltemp = log10 (t_abs);
    double DHA = -0.4883 + 8.074e-4 * temperature;
    double DHB = -0.3241 + 1.600e-4 * temperature;
    double sqI = sqrt (0.1);
    const double gammCa = pow (10, (-4. * DHA * sqI) / (1. + 5.e-8 * DHB * sqI));
    const double gammHCO3 = pow (10, (-DHA * sqI) /  (1. + 5.4e-8 * DHB * sqI));
    const double K1 = pow (10, -356.3094 - 0.06091964 * t_abs + 21834.37 / t_abs + 126.8339 * ltemp -  1684915 / t2);
    const double K2 = pow (10, -107.8871 - 0.03252849 * t_abs +  5151.79 / t_abs + 38.92561 * ltemp - 563713.9 / t2);
    const double KC = pow (10, -171.9065 - 0.077993 * t_abs +  2839.319 / t_abs + 71.595 * ltemp);
    const double KH = pow (10, 108.3865 + 0.01985076 * t_abs - 6919.53 / t_abs - 40.45154 * ltemp +  669365.0 / t2);

    const double* gammCa_h = &gammCa ;
    const double* gammHCO3_h = &gammHCO3;
    const double* K1_h = &K1;
    const double* K2_h = &K2;
    const double* KC_h = &KC;
    const double* KH_h = &KH;

	if( !(cudaMemcpyToSymbol(PE_d, hPE,              sizeof(double), 0 , cudaMemcpyHostToDevice) == cudaSuccess)) printf ("failed \n");
	if( !(cudaMemcpyToSymbol(P0_d, hP0,              sizeof(double), 0 , cudaMemcpyHostToDevice) == cudaSuccess)) printf ("failed \n");
	if( !(cudaMemcpyToSymbol(P1_d, hP1,              sizeof(double), 0 , cudaMemcpyHostToDevice) == cudaSuccess)) printf ("failed \n");
	if( !(cudaMemcpyToSymbol(gammCa_d, gammCa_h ,    sizeof(double), 0 , cudaMemcpyHostToDevice) == cudaSuccess)) printf ("failed \n");
	if( !(cudaMemcpyToSymbol(gammHCO3_d, gammHCO3_h, sizeof(double), 0 , cudaMemcpyHostToDevice) == cudaSuccess)) printf ("failed \n");
	if( !(cudaMemcpyToSymbol(K1_d, K1_h,             sizeof(double), 0 , cudaMemcpyHostToDevice) == cudaSuccess)) printf ("failed \n");
	if( !(cudaMemcpyToSymbol(K2_d, K2_h,             sizeof(double), 0 , cudaMemcpyHostToDevice) == cudaSuccess)) printf ("failed \n");
	if( !(cudaMemcpyToSymbol(KC_d, KC_h,             sizeof(double), 0 , cudaMemcpyHostToDevice) == cudaSuccess)) printf ("failed \n");
	if( !(cudaMemcpyToSymbol(KH_d, KH_h,             sizeof(double), 0 , cudaMemcpyHostToDevice) == cudaSuccess)) printf ("failed \n");

	int type = 1;

	dim3  threads( BLOCKCOLS, BLOCKROWS);
	dim3  grid(ncell_x/BLOCKCOLS, ncell_y/ BLOCKROWS);

	  //  former process routine now merged
	  int gridRows = data->mapInfo.height;
	  int gridColumns = data->mapInfo.width;
	  int full_size = gridColumns * gridRows;

	//__global__ void weathering2 (int *mask, double temperature, double rain, double *FAPtr_d, int type, double *weatherC_d, double *weatherP_d, double *soilTPtr_d,
		//                                 double* finesPtr_d, double* stonePtr_d, double* soilMPtr_d, double* eroPtr_d, double* dz_d, double cell_size, int rows, int cols)
	weathering<<<grid, threads>>>(device->mask, temperature, rain, device->fa, type, device->weatherC, device->weatherP, device->soilTPtr, device->finesPtr, device->stonePtr, device->soilMPtr,
		                              device->eroPtr, device->dz, data->mapInfo.cellsize, data->mapInfo.height, data->mapInfo.width);


	thrust::device_ptr<double> weatherC_d = thrust::device_pointer_cast(device->weatherC);
	cudaSetDevice(0);
	data->totweatherC = thrust::reduce(weatherC_d, weatherC_d + full_size, (double) 0);
	//printf("total Incision from thrust is %10.8lf \n", data->totI);

	thrust::device_ptr<double> weatherP_d = thrust::device_pointer_cast(device->weatherP);
	cudaSetDevice(0);
	data->totweatherP = thrust::reduce(weatherP_d, weatherP_d + full_size, (double) 0);
	//printf("total Incision from thrust is %10.8lf \n", data->totI);





	//printf("MOD: weathering :%s\n", cudaGetErrorString(cudaGetLastError()));
	fprintf(data->outlog, "MOD: weathering :%s\n", cudaGetErrorString(cudaGetLastError()));

}

