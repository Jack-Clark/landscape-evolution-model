#include "MapInfo.h"

// Compute the actual x and y coordinates of a cell
// xllcorner, yllcorner is the top left hand corner of the DEM
// To compute the actual coordinates of a cell we perform
// x = xllcorner + cellsize * r
// y = yllcorner + cellsize * c
int calculateRealXY(double* x, double* y, MapInfo* map, int r, int c) {
  *x = map->cellsize * c + map->xllcorner;
  *y =  ((map->height* map->cellsize)+map->yllcorner) - (map->cellsize * r);
  return 1;
}