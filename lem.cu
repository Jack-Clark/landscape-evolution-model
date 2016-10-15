
#include "lem.h" // this contains headers for all standard and add-in libraries

#include "config.h"
#include "MapInfo.h"
#include "Data.h"
#include "memory_dev.h"
#include "memory.h"

#include "SFD.h" // gives access to flow direction routine
#include "runoffweight.h" // gives access to runoffweight calculations (losses)
#include "contribA.h"
#include "FA_SFD.h" //gives access to flow accumulation algorithm
#include "erosion.h" // gives access to sediment erosion, transport, deposition model

int main(int argc, char* argv[]) {

	const char* clim_file;
	const char* file;
	const char* maskfile;

	char* heightfile = (char*) malloc(sizeof(char) *100);
	char* SFDfile = (char*) malloc(sizeof(char) *100);
    char* FAfile = (char*) malloc(sizeof(char) *100);
	char* erofile = (char*) malloc(sizeof(char) *100);
	char* incifile = (char*) malloc(sizeof(char) *100);
	char* gelifile = (char*) malloc(sizeof(char) *100);
	char* depofile = (char*) malloc(sizeof(char) *100);
	char* slopefile = (char*) malloc(sizeof(char) *100);

	char* finesfile = (char*) malloc(sizeof(char) *100);
	char* stonesfile = (char*) malloc(sizeof(char) *100);
	char* totbiofile = (char*) malloc(sizeof(char) *100);
	char* soilTfile = (char*) malloc(sizeof(char) *100);
	char* nutfile = (char*) malloc(sizeof(char) *100);

	char* wCfile = (char*) malloc(sizeof(char) *100);
	char* wPfile = (char*) malloc(sizeof(char) *100);

	char* catchmap = (char*) malloc(sizeof(char) *100);
	char* catchmask = (char*) malloc(sizeof(char) *100);
	char* contrib = (char*) malloc(sizeof(char) *100);
	char* rivermask = (char*) malloc(sizeof(char) *100);

	char* logfile = (char*) malloc(sizeof(char) *100);

	clock_t startrun;
	clock_t endrun;
	double Ptime;
	
	size_t freenow, total;

	Data data;
	Data device;
	Catchment catchments;

	data.usethisgpu = 1; // this sets which GPU to use - remember to alter the target exe name to match

	data.flowdirectiontype = 0; // 0 for SFD and 1 for MFD
	data.WhichGrid = 8 ;// currently only works up to 25 million cell grid
	data.max_iterations = 5;
	data.outfiles = 1; // only option at present!

	// new parameters to adjust vegetation response
	data.temp_inflex = 6. ;
	data.temp_sens = 6.;

	sprintf(logfile, "Output%d/logfile_%dm_%d.txt", data.usethisgpu, data.WhichGrid, data.max_iterations);

	data.logfileOut = logfile;
	data.outlog = fopen(data.logfileOut, "w");
	if(!data.outlog) {
		printf("\nThe file: %s failed to open.", data.logfileOut);
		exit(0);
	}
	fprintf(data.outlog, "logfile:  grid:%dm_   Max_iterations:%d. \n",data.WhichGrid, data.max_iterations);
	cudaMemGetInfo(&freenow, &total);
	fprintf(data.outlog,"Memory on CUDA card free at start of iteration: %d total: %d\n",freenow/1024,total/1024);
	
	startrun = clock();
	switch (data.WhichGrid)
	{

		case 8:
			file = "inputgrids/thames20m.asc"; // this is the modified 20m DEM i.e. given an outlet cell
			maskfile = "inputgrids/20mMask.tif"; // this is the mask file
			//totgridsize = 8666880  2928 cols, 2960 rows
			break;
	}
	clim_file = "inputdata/SSTbase.txt";
		

	fprintf(data.outlog,"Landscape file   = %s\n", file);
	//fprintf(data.outlog,"Output directory = %s\n", OUTPUT_PATH);
	fprintf(data.outlog,"Climate file     = %s\n", clim_file);

	ReadClimateDataFromFile(&data, clim_file); // read in appropriate number of data pairs, temp and rain and allocate sufficient memory


	// Set up the dem and other data sets
	readgdalDEMfromFile(&data, file);  // the host malloc for dem is in here
	readgdalmaskfromFile(&data, maskfile) ; // read in the datamask
	copyMask(&data, &device); // copy mask to device only needs to be done once

	createCatchmentSpace(&data, &catchments);

	createProcessMatrices(&data); // allocate sufficient memory on host
	setProcessMatrices(&data);  // initialise values for attributes

	// set runin conditions
	  data.rain = 300; //550; //300;
	  data.temp = 0.0; // 2.0; //0. ;
	  printf("current rain %f \n", data.rain);

	createcontribAspace(&data); // create space for contribA data

	// Perform the iterations
    printf("Performing the simulation\n");

    fprintf(data.outlog, "Performing the simulation\n");

    fflush(data.outlog);

    data.max_iterations = 135000;

    int start_pos_sim; // start position in climate file
    int end_pos_sim;
    int runiniter; // how many on default rain and temp set above
    int totalnumiter;

    runiniter = 5;
    start_pos_sim = 0; //  140k dataset

    data.start_year = (int *) malloc(sizeof(int));
    *data.start_year = start_pos_sim;

    end_pos_sim = 5;

    totalnumiter = runiniter + (end_pos_sim-start_pos_sim);


	for (int i = 1; i < totalnumiter + 1; i++) {

	printf("\n***********************\n");\
	printf("Starting Iteration %d\n", i);
	printf("***********************\n");

	fprintf(data.outlog,"\n***********************\n");
	fprintf(data.outlog,"Starting Iteration %d\n", i);
	fprintf(data.outlog,"***********************\n");


		createDeviceSpace(&data, &device);


		if (i >= runiniter ){
			setClimate(&data, start_pos_sim );
			start_pos_sim ++;
		}


		// undertake flow routing and compute catchments and contributing areas
		setdevicespace_FD(&data, &device);
			cuFlowDirection(&data, &device, &catchments, i);
			calccontribA(&data, &device, &catchments); // calculate contributing area for each cell
		cleardevicespace_FD(&data, &device);



		setdevicespace_FA(&data, &device);  // load matrices for runoffweight calculation
		computeRunOffWeights(&data, &device);


		//if (data.flowdirectiontype == 0) {
		// ****************************************  SFD   ****************************************************
		printf("Entering FA: \n");
		//accumulateflowSFD(&data, &device, i);
		correctflow_SFD(&data, &device, i);
		cleardevicespace_FA(&data, &device);

			//correctflow_SFD_NoPart_List(data.fa, data.fd, data.mapInfo.width, data.mapInfo.height, data.runoffweight, data.WhichGrid);
				// To use more than one GPU for SFD use...
			//multicorrectflow(data.fa, data.fd, data.mapInfo.width, data.mapInfo.height, data.runoffweight, data.WhichGrid);


		//}
		//else
		// ****************************************  MFD   *****************************************************
		//{
		//	correctmfdflow(&data);
			//getLastCudaError("Anything here in lem.cu .\n");
		//}

		setdevicespace_Process(&data, &device);
			erosionGPU(&data, &device, &catchments, i);
	    cleardevicespace_Process(&data, &device);

		if ( i==10 || i%1000==0 )
		{
			// Write the recalculated DEM
			sprintf(heightfile, "Output%d/%d_height.tif",data.usethisgpu,  i);
			writegdalGRIDtoFile(&data, &catchments, heightfile, 0, 0);
			
			// Write the recalculated SFD
			if (data.flowdirectiontype == 0) sprintf(SFDfile, "Output%d/%d_SFD.tif",data.usethisgpu,  i);
			if (data.flowdirectiontype == 1) sprintf(SFDfile, "Output%d/MFD_%dm_%d.tif", data.usethisgpu, data.WhichGrid, i);
			writegdalGRIDtoFile(&data, &catchments, SFDfile, 1, 0);

			//sprintf(contrib, "output1/contrib_%dm_%d.tif", data.WhichGrid, i);
			//writegdalGRIDtoFile(&data, &catchments, contrib, 1, 3);

			// Write the recalculated FA
			sprintf(FAfile, "Output%d/%d_FA.tif", data.usethisgpu,  i);
			writegdalGRIDtoFile(&data, &catchments, FAfile, 0, 1);
					
			// Write the recalculated erosion
			sprintf(erofile, "Output%d/%d_ero.tif", data.usethisgpu,  i);
			writegdalGRIDtoFile(&data, &catchments, erofile, 0, 2);
					
			// Write the recalculated incision
			sprintf(incifile, "Output%d/%d_inci.tif", data.usethisgpu, i);
			writegdalGRIDtoFile(&data, &catchments, incifile, 0, 3);
					
			// Write the recalculated deposition
			sprintf(depofile, "Output%d/%d_depo.tif", data.usethisgpu, i);
			writegdalGRIDtoFile(&data, &catchments, depofile, 0, 4);
					
			// Write the recalculated slope
			sprintf(slopefile, "Output%d/%d_slope.tif", data.usethisgpu,  i);
			writegdalGRIDtoFile(&data, &catchments, slopefile, 0, 5);

			// Write the recalculated fines
			sprintf(finesfile, "Output%d/%d_fines.tif", data.usethisgpu,  i);
			writegdalGRIDtoFile(&data, &catchments, finesfile, 0, 6);

			// Write the recalculated stone
			sprintf(stonesfile, "Output%d/%d_stone.tif", data.usethisgpu, i);
			writegdalGRIDtoFile(&data, &catchments, stonesfile, 0, 7);

			sprintf(totbiofile, "Output%d/%d_totbio.tif", data.usethisgpu,  i);
			//writegdalGRIDtoFile(&data, &catchments, totbiofile, 0, 8);

			sprintf(soilTfile, "Output%d/%d_soilT.tif", data.usethisgpu,  i);
			writegdalGRIDtoFile(&data, &catchments, soilTfile, 0, 9);

			sprintf(nutfile, "Output%d/%d_nut.tif", data.usethisgpu,  i);
			writegdalGRIDtoFile(&data, &catchments, nutfile, 0, 10);

			sprintf(wCfile, "Output%d/%d_wC.tif", data.usethisgpu,  i);
			writegdalGRIDtoFile(&data, &catchments, wCfile, 0, 11);

			sprintf(wPfile, "Output%d/%d_wP.tif",data.usethisgpu,  i);
			writegdalGRIDtoFile(&data, &catchments, wPfile, 0, 12);

			sprintf(gelifile, "Output%d/%d_geli.tif",data.usethisgpu,  i);
			writegdalGRIDtoFile(&data, &catchments, gelifile, 0, 13);
		}	

		printf("Finished Iteration %d\n", i);
		fprintf(data.outlog,"Finished Iteration %d\n", i);
		//printf("Here\n");
		printf("\nIteration logged\n");
		writeSummaryDataToFile(&data, data.outfiles, i); // note summary output is every iteration with argument of 1, 2, 5 or 10 only
		printf("\nSummary logged\n");
		zerogrids(&data);
		//printf("herer2\n");
		clearDeviceSpace(&data, &device);


		cudaMemGetInfo(&freenow, &total);
		fprintf(data.outlog,"Memory on CUDA card free at end of iteration: %d total: %d\n",freenow/1024,total/1024);

	} // end of iterations

	printf("%s\n", "\nEnded main iteration loop\n");
	endrun = clock();
	Ptime = (double) (endrun-startrun)/ CLOCKS_PER_SEC ;
	fprintf(data.outlog, "Total simulation time of : %20.17lf \n\n", Ptime);

	cudaFree(device.mask);
	deleteProcessMatrices(&data);  //delete memory on host
	clearcontribAspace(&data);
	free(data.start_year);
	fclose(data.outlog);

}
