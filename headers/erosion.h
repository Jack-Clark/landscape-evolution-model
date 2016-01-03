
#ifndef EROSION_H_
#define EROSION_H_

#include "lem.h"
#include "Data.h"
#include "eroincidep.h"
#include "config.h"
#include "Directions.h"


void erosionGPU( Data* data, Data* device, Catchment* catchments, int iter);


#endif /* EROSION_H_ */
