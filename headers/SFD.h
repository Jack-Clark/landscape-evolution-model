#ifndef SFDH
#define SFDH

#include "lem.h"
#include "Data.h"
#include "Directions.h"
#include "watershed.h"
#include "floodingdriver.h"

void cuFlowDirection(Data* data, Data* device, Catchment* catchments, int iter);

#endif
