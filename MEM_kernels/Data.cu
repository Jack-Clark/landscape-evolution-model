
#include "Data.h"

int writelogfile(Data* data)
{
	data->outlog = fopen(data->logfileOut, "a");
	fprintf(data->outlog, "logfile:  grid:%dm_   Max_iterations:%d. \n");
	fclose(data->outlog);
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
// GRID functions
//////////////////////////////////////////////////////////////////////////////

int ReadClimateDataFromFile(Data* data, const char* clim_file)
{
	char a[100];
	int year;
	float changetemp;
	float changerain;

	// create sufficient space for climate data
	data->rainchange = (double *) calloc (data->max_iterations, sizeof(double));
	data->tempchange = (double *) calloc (data->max_iterations, sizeof(double));
	data->year =       (int *)    calloc (data->max_iterations, sizeof(int));

	// now read in the data
	FILE* in = fopen(clim_file, "r");
    if (!clim_file)
        {
			printf("Climate file does not exist \n ");
	        exit (0);
        }

	// read the header
	fscanf(in, "%s", a); // year
	fscanf(in, "%s", a); // ppt
	fscanf(in, "%s", a); // temp

	for(int i = 0; i < data->max_iterations; i++)
	{

		fscanf(in, "%d", &year);
		data->year[i] = (int) year;
		fscanf(in, "%f", &changerain);
		data->rainchange[i] = (double) changerain;
		fscanf(in, "%f", &changetemp);
		data->tempchange[i] = (double) changetemp;
	}

	fclose(in);
	printf("climate file read \n");
	fprintf(data->outlog,"climate file read \n");
	return 0;
}

int setClimate(Data* data, int iteration)
{	

	//double mean_rain = 850. ;
	//double mean_temp = 10.0 ;

	data->last_rain = data->rain; // set up defaults
	data->last_temp = data->temp; 
	//data->rain = mean_rain + mean_rain * data->rainchange[iteration-1];
	//data->temp = mean_temp + data->tempchange[iteration-1];

	data->rain = data->rainchange[iteration-1];
	data->temp = data->tempchange[iteration-1];



	fprintf(data->outlog,"rain: %f. temp: %f, last_rain: %f, last temp: %f %f \n", data->rain, data->temp, data->last_rain, data->last_temp);

	return 0;
}



int readgdalDEMfromFile(Data* data, const char* fileName) {
	
	// Required by gdal n.b. change to c notation
    GDALDatasetH  poDataset;
    GDALRasterBandH  poBand;
    GDALDriverH  poDriver;
    double  adfGeoTransform[6]; //contains projection data
	int  nodatavalue;
	int nXSize, nYSize;
	

	GDALAllRegister();
		
    poDataset = GDALOpen( fileName, GA_ReadOnly );
      if( poDataset == NULL )
      {
        printf("File does not exist");
        exit (0);
      }

	poDriver = GDALGetDatasetDriver ( poDataset );
    
	fprintf(data->outlog, "Driver: %s/%s\n",
			GDALGetDriverShortName( poDriver ),
			GDALGetDriverLongName (poDriver ) );

	poBand = GDALGetRasterBand ( poDataset, 1);
	data->mapInfo.width = GDALGetRasterBandXSize( poBand );
	data->mapInfo.height = GDALGetRasterBandYSize( poBand );

	fprintf(data->outlog, "Size is %d %d\n",
		data->mapInfo.width, 
		data->mapInfo.height);


    if( GDALGetGeoTransform( poDataset, adfGeoTransform ) == CE_None )
    {
		data->mapInfo.cellsize = adfGeoTransform[1] ;
		data->mapInfo.xllcorner = adfGeoTransform[0] ;
		data->mapInfo.yllcorner = adfGeoTransform[3] ;
    }

    fprintf(data->outlog,"Cell Size = (%.6f)\n",data->mapInfo.cellsize );

	poBand = GDALGetRasterBand ( poDataset, 1);
	data->mapInfo.nodata = GDALGetRasterNoDataValue(poBand, &nodatavalue);
	fprintf(data->outlog, "No Data Value = (%.6f)\n",data->mapInfo.nodata );

	// Read the values into a buffer zs
	nXSize = GDALGetRasterBandXSize( poBand );
	nYSize = GDALGetRasterBandYSize( poBand );
	checkCudaErrors(cudaMallocHost((void **)&data->dem, data->mapInfo.height * data->mapInfo.width * sizeof(double)));

	printf("nXSize : %d, nYSize : %d \n", nXSize, nYSize);

    // ************ Read in height data from source ************************
    GDALRasterIO( poBand, GF_Read, 0, 0, nXSize, nYSize, data->dem, nXSize, nYSize, GDT_Float64, 0, 0 );
	printf( "Data Read In Complete \n");
	fprintf(data->outlog,"DEM Data Read In Complete \n");
	return 1;
}


int readgdalmaskfromFile(Data* data, const char* fileName) {

	// Required by gdal n.b. change to c notation
    GDALDatasetH  poDataset;
    GDALRasterBandH  poBand;
    GDALDriverH  poDriver;

    int  nodatavalue;
	int nXSize, nYSize;


	GDALAllRegister();

    poDataset = GDALOpen( fileName, GA_ReadOnly );
      if( poDataset == NULL )
      {
        printf("File does not exist");
        exit(0);
      }

	poDriver = GDALGetDatasetDriver ( poDataset );

	fprintf(data->outlog,"Driver: %s/%s\n",
			GDALGetDriverShortName( poDriver ),
			GDALGetDriverLongName (poDriver ) );

	poBand = GDALGetRasterBand ( poDataset, 1);
	nodatavalue = GDALGetRasterNoDataValue(poBand, &nodatavalue);
	printf( "No Data Value = (%d)\n",nodatavalue );

	// Read the values into a buffer zs
	nXSize = GDALGetRasterBandXSize( poBand );
	nYSize = GDALGetRasterBandYSize( poBand );
	data->mask = (int *) CPLMalloc(sizeof(int)*nXSize*nYSize);

	fprintf(data->outlog,"nXSize : %d, nYSize : %d \n", nXSize, nYSize);

    // ************ Read in height data from source ************************
    GDALRasterIO( poBand, GF_Read, 0, 0, nXSize, nYSize, data->mask, nXSize, nYSize, GDT_Int32  , 0, 0 );
	printf( "Mask Data Read In Complete \n");
	fprintf(data->outlog,"Mask Data Read In Complete \n");
	return 1;
}


//int writegdalGRIDtoFile(Data* data, const char* fileName, int whichtype, int what) {
int writegdalGRIDtoFile(Data* data, Catchment* catchments, char* fileName, int whichtype, int what) {

	double transformcell = - (data->mapInfo.cellsize) ;
	double adfGeoTransform[6] = {data->mapInfo.xllcorner, data->mapInfo.cellsize, 0, data->mapInfo.yllcorner, 0, transformcell };


	//printf(" %f, %f, %f, %f, %f, %f ", adfGeoTransform[0],adfGeoTransform[1],adfGeoTransform[2],adfGeoTransform[3],adfGeoTransform[4],adfGeoTransform[5]);

	GDALDatasetH DDstDS; // for doubles
	GDALDatasetH IntDstDS; // for integers
	
    OGRSpatialReferenceH hSRS;
    char *pszSRS_WKT = NULL;

    char **Options = NULL;
	const char *pszFormat = "GTiff";
    GDALDriverH hDriver = GDALGetDriverByName( pszFormat );
	GDALRasterBandH  DstBand, DstBandFD;
	
	if( hDriver == NULL )  exit( 1 );

	Options = CSLSetNameValue( Options, "COMPRESS","LZW");

    hSRS = OSRNewSpatialReference( NULL );
    //OSRSetWellKnownGeogCS( hSRS, "EPSG:27700" );
    OSRImportFromEPSG(hSRS, 27700);
    OSRExportToWkt( hSRS, &pszSRS_WKT );
    OSRDestroySpatialReference( hSRS );

	switch (whichtype) {
		case 0: // write doubles
		DDstDS = GDALCreate( hDriver, fileName, data->mapInfo.width, data->mapInfo.height, 1, GDT_Float64, Options );
		GDALSetGeoTransform( DDstDS, adfGeoTransform );


	    GDALSetProjection( DDstDS, pszSRS_WKT );
	    CPLFree( pszSRS_WKT );

		DstBand = GDALGetRasterBand( DDstDS, 1 );
		GDALSetRasterNoDataValue(DstBand, -9999);


			switch (what) {
				case 0: // write the DEM
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->dem, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 ); 
				break;
				
				case 1: // write FA
					GDALSetRasterNoDataValue(DstBand, -1);
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->fa, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 ); 
				break;
				
				case 2: // write erosion (diffuse)
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->eroPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;
				
				case 3: // write incision (con concentrated)
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->inciPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 ); 
				break;
				
				case 4: // write deposition
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->depoPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 ); 
				break;

				case 5: // write slope
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->SlopePtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 ); 
				break;

				case 6: // write fines
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->finesPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 7: // write stone
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->stonePtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 8: // write total Bio
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->TotBPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 9: // write soil thickness
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->soilTPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 10: // write nutrients
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->nutPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 11: // write chemical weathering
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->weatherC, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 12: // write Physical weathering
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->weatherP, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;

				case 13: // write solifluction
					GDALRasterIO( DstBand, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->geliPtr, data->mapInfo.width, data->mapInfo.height, GDT_Float64, 0, 0 );
				break;
			}	
		

		GDALClose (DDstDS);
		
		break;
		
		case 1: // write the integers
		IntDstDS = GDALCreate( hDriver, fileName, data->mapInfo.width, data->mapInfo.height, 1, GDT_Int32, Options );
		GDALSetGeoTransform( IntDstDS, adfGeoTransform );

	    GDALSetProjection( IntDstDS, pszSRS_WKT );
	    CPLFree( pszSRS_WKT );


		DstBandFD = GDALGetRasterBand( IntDstDS, 1 );
			switch (what) {
			case 0: // write the SFD file
					GDALRasterIO( DstBandFD, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->fd, data->mapInfo.width, data->mapInfo.height, GDT_Int32, 0, 0 );
			break;
			
			case 1: // write the catchment id map
				GDALRasterIO( DstBandFD, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, catchments->watershed_id, data->mapInfo.width, data->mapInfo.height, GDT_Int32, 0, 0 ); 
			break;

			case 2: // write the catchment mask
				GDALRasterIO( DstBandFD, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, catchments->mask, data->mapInfo.width, data->mapInfo.height, GDT_Int32, 0, 0 ); 
			break;

			case 3: // write the contributing area
				GDALRasterIO( DstBandFD, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->contribA, data->mapInfo.width, data->mapInfo.height, GDT_Int32, 0, 0 ); 
			break;

			case 4: // write the rivermask
			//	GDALRasterIO( DstBandFD, GF_Write, 0, 0, data->mapInfo.width, data->mapInfo.height, data->mask, data->mapInfo.width, data->mapInfo.height, GDT_Int32, 0, 0 );
			break;
			}
		GDALClose (IntDstDS);
		break;
	}

	return 1;
}


int writeSummaryDataToFile(Data* data, int howmany, int iteration)
{
	
	// outfiles given csv extension for wasy read into excel etc for plotting
	char* outfilename = (char*) malloc(sizeof(char) *100); 
	sprintf(outfilename, "summarydata/modelsum%d/model%d.csv", data->usethisgpu, data->usethisgpu);
	const char* outfile1 = outfilename;	

	FILE* out1; 

			if (iteration == 1)
			{
				out1 = fopen(outfile1, "w");
				if(!out1) {
					printf("\nThe file %s couldn't be opened for writing", outfile1);
					exit(0);
				}
				fprintf(out1, "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s \n",
						      "iteration","year", "temp", "rain", "diffuse", "concentrated", "gelifluction", "chemicalW", "physicalW", "deposition", "totbio", "FA_max", "sedoutflow", "waterout", "sinkcount");
				fclose(out1);
			}

	int thisyear;
	if (iteration <500 ) thisyear = -500 + iteration;
	if (iteration >= 500) thisyear = data->year[(*data->start_year) + iteration -500];


	out1=fopen(outfile1,"a");
	double sedoutflow = ( (data->totE + data->totI + data->totG - data->totD) / 5069059 )  ; // outputs average per cell in mm, erosion already expressed in mm?
	printf("Ave erosion rate mm/a : %10.8lf \n", sedoutflow);
	double totalrunoff = ((data->FA_max) * 1000) / 5069059 ; // water depth is in m [see runoffweights] therefore * 1000 for mm div no of cells to get rainfall equivalent
	fprintf(out1, "%d, %d, %f, %f, %10.8lf, %10.81lf, %10.8lf, %10.8lf, %10.8lf, %10.8lf, %10.8lf, %10.3lf, %10.8lf,  %10.5lf, %d \n",
			        iteration, thisyear, data->temp, data->rain, data->totE, data->totI, data->totG, data->totweatherC, data->totweatherP, data->totD, data->totbio, data->FA_max,  sedoutflow, totalrunoff, data->sinkcount);
	fclose(out1);

	free(outfilename);
	return 0;
}



//////////////////////////////////////////////////////////////////////////////
// Process Data Matrices
//////////////////////////////////////////////////////////////////////////////

int setProcessMatrices(Data* data)
{
  int fullsize;
  int x;

  double calinitBio;
  
  fullsize =  data->mapInfo.width * data->mapInfo.height;
  
  // first reading for actual model is -6.8, ~53% ppt


  calinitBio = ( 1000.0 + (5.0 * (data->rain - 100.0) ) ) / 25;  // scaled for 20m cells?

  for (x = 0; x < fullsize; x++)
  {
	  if (data->mask[x] == 1) // inside the catchment(s)
	  {
		  data->runoffweight[x] = 	  1.0 ;
		  data->soilTPtr[x] = 		  1.0 ;
		  data->stonePtr[x] = 		 20.0 ;
		  data->finesPtr[x] = 		 80.0 ;
		  data->soilMPtr[x] = 		  0.5 ;
		  data->nutPtr[x] = 		 10.0 ;
		  data->soilBPtr[x] = 		  0.0 ;
		  data->TotBPtr[x] = 		calinitBio;

	  } else // if not in catchment(s) set values to zero
	  {
		  data->runoffweight[x] = 	  0.0 ;
		  data->soilTPtr[x] = 		  0.0 ;
		  data->stonePtr[x] = 		  0.0 ;
		  data->finesPtr[x] = 		  0.0 ;
		  data->soilMPtr[x] = 		  0.0 ;
		  data->nutPtr[x] =   		  0.0 ;
		  data->soilBPtr[x] = 		  0.0 ;
		  data->TotBPtr[x] =  		  0.0 ;
	  }
		  data->eroPtr[x] = 		  0.0 ;
		  data->geliPtr[x] = 		  0.0 ;
		  data->inciPtr[x] = 		  0.0 ;
		  data->depoPtr[x] = 		  0.0 ;
		  data->weatherC[x] = 		  0.0 ;
		  data->weatherP[x] = 		  0.0 ;
  }

		  data->Qthreshold = 		 600.0; // threshold depth for diffuse/concentrated erosion ?mm - this needs confirmation
		  data->soil_struct = 		   2.0;
		  data->profile_permeability = 3.0;

  printf("All process grids now allocated \n\n");
  fprintf(data->outlog, "All process grids now allocated and default values set \n\n");
  fprintf(data->outlog, "runoffweight: %f \nsoilThickness: %f \nstone%: %f, fines%: %f \nsoilMoisture: %f, nutrients: %f, soilBio: %f, TotalBio: %f \n",
		               data->runoffweight[3000000], data->soilTPtr[3000000], data->stonePtr[3000000], data->finesPtr[3000000], data->soilMPtr[3000000], data->nutPtr[3000000], data->soilBPtr[3000000], data->TotBPtr[3000000]);

  return 0;
}
