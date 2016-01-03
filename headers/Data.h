
#ifndef DATA
#define DATA

#include "lem.h"
#include "MapInfo.h"


/*! This is the main data structure used to pass values across the sub-models */



typedef struct DataS {
    double* dem;
	int* fd;
	int* fdmod;
	double* fa;
	int* catchmentarea;
	int usethisgpu;

	MapInfo mapInfo;

	double cellarea;
	double* runoffweight ;
	double* SlopePtr; // slope grid
	int cpucontribA_max;

	double* finesPtr /*! fines_content% */ ;
	double* stonePtr /*! stone content% */ ;
	double* soilTPtr /*! soil thickness [mm] */ ;
	double* soilMPtr /*! soil moisture % */ ;
	double* soilBPtr  /*! soil biomass % */ ;
	double* oldsoilBPtr /*! old soil biomass % */ ;
	double* nutPtr   /*! nutrient content % */ ;
	double* TotBPtr  /*! total biomass calculated load in kg per ??*/ ;

	double* eroPtr /*! diffuse erosion grid */ ;
	double* geliPtr /*! gelifluction grid */ ;
	double* inciPtr /*! concentrated erosion grid */ ;
	double* depoPtr /*! deposition grid */ ;

	double totE /*! total annual diffuse erosion */ ;
	double totI /*! total annual concentrated erosion */ ;
	double totD /*! total annual deposition */ ;
	double totG /*! total annual gelifluction */ ;
	double totbio /*! total annual biomass */ ;
	double totweatherP;
	double totweatherC;

	double* summary /*! used to perform reduction stats using thurst */;

	int soil_struct  /*! soil structure unused at present */ ;
	int profile_permeability  /*! soil permeability unused at present */ ;

	double* weatherP /*! physical weathering */ ;
	double* weatherC /*! chemical weathering */ ;
	//double* evapPtr;

	int max_iterations;
	int WhichGrid;
	int flowdirectiontype;
	int outfiles;
	double  adfGeoTransform[6]; // needed for gdal write operations to contain projection parameters

	int* dx;
	int* dy;
	double* dz;
	float* shortest_paths;
    int* watershed_id;
    double* lowHeight;

	int* contribA;
	double Qthreshold;

	int outletcellidx;

	double rain;
	double last_rain;
	double* rainchange;
	double temp;
	double last_temp;
	double* tempchange;
	int* year;
	int* start_year;

	int* mask;
	int* rivermask;

	int sinkcount;

	const char* logfileOut;
	FILE* outlog;

	double FA_max;


	double temp_inflex;
	double temp_sens;

	int activecells;
} Data, DataCopy;


typedef struct CatchmentS {

	int* watershed_id;
	int* catchmentarea;
	int* toptencatchmentArea;
	int* toptencatchmentID;
	int* toptenrootcell;
	int* mask;

} Catchment;



struct is_not_zero
{
	__host__ __device__ bool operator() (double x)
	{
	return (x != 0.0) ;
	}
};

struct is_not_negative
{
	__host__ __device__ bool operator() (double x)
	{
	return (x >= 0.0) ;
	}
};

int writelogfile(Data* data);

int ReadClimateDataFromFile(Data* data, const char* clim_file);

int setClimate(Data* data, int iterations);

int zerogrids(Data* data);

int readgdalDEMfromFile(Data* data, const char* fileName);

int readgdalmaskfromFile(Data* data, const char* fileName) ;

int writegdalGRIDtoFile(Data* data, Catchment* catchments, char* fileName, int whichtype, int what) ;

int writeSummaryDataToFile(Data* data, int howmany, int iteration);

int setProcessMatrices(Data* data);

#endif
