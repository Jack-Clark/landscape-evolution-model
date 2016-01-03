#ifndef DEMTYPE
#define DEMTYPE

typedef struct MapInfoS {
	int width;
	int height;
	double xllcorner;
	double yllcorner;
	double cellsize;
	int data_is_double;
	double nodata;
} MapInfo;

int calculateRealXY(double* x, double* y, MapInfo* map, int r, int c);

#endif
