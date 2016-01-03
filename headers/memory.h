/*
 * memory
 *
 *  Created on: Jan 31, 2014
 *      Author: ndm12
 */

#ifndef MEMORY_
#define MEMORY_

#include "lem.h"
#include "Data.h"

int createcontribAspace(Data* data);
int clearcontribAspace(Data* data);

int createProcessMatrices(Data* data);
int deleteProcessMatrices(Data* data);

int createCatchmentSpace(Data* data, Catchment* Catchments);

#endif /* MEMORY_ */
