/*
 * updates.h
 *
 *  Created on: Feb 25, 2014
 *      Author: ndm12
 */
#include "lem.h"
#include "Data.h"

#ifndef UPDATES_H_
#define UPDATES_H_

void update_newSurface(Data* data, Data* device, int iter);
void update_nutrients(Data* data, Data* device);
void update_vegetation(Data* data, Data* device);

#endif /* UPDATES_H_ */
