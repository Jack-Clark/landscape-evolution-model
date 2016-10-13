/*
 * FA_SFD.h
 *
 *  Created on: Feb 11, 2014
 *      Author: ndm12
 */

#ifndef FA_SFD_H_
#define FA_SFD_H_

#include "lem.h"
#include "Data.h"
#include "Directions.h"
#include "config.h"

void correctflow_SFD(Data* data, Data* device, int iter);

int process_SFD_Multiple_Retries(Data* data, Data* device, int iter);

int process_SFD_NoPart_List(Data* data, Data* device, int iter);

int process_SFD_block_level_single_chain(Data* data, Data* device, int iter);

int process_SFD_global_level_single_chain(Data* data, Data* device, int iter);

void accumulateflowSFD(Data* data, Data* device, int iter);

#endif /* FA_SFD_H_ */
