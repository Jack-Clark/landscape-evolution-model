/*
 * FA_SFD.h
 *
 *  Created on: Feb 11, 2014
 *      Author: ndm12
 */

#ifndef FA_SFD_H_
#define FA_SFD_H_

/* TODO: Refactor header inclusion. These should be included in the source files as needed, not bundled into one header. 
   Only Data.h is required in this header. */
#include "lem.h"
#include "Data.h"
#include "Directions.h"
#include "config.h"

int correctflow_SFD(Data* data, Data* device, int iter, int algorithmID);

int process_SFD_Multiple_Retries(Data* data, Data* device, int iter);

int process_SFD_NoPart_List(Data* data, Data* device, int iter);

int process_SFD_block_level_single_chains(Data* data, Data* device, int iter);

int process_SFD_global_level_single_chains(Data* data, Data* device, int iter);

#endif /* FA_SFD_H_ */
