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

//#include "pad.h"
//#include "partition.h"
//#include "MGPU_correct.h"
//#include "mfd_correct.h"
//#include "parallel_sfd_noPart.h"
//#include "parallel_sfd_noPart_noOK.h"
//#include "parallel_sfd_noPart_List.h"


void accumulateflowSFD(Data* data, Data* device, int iter);

void correctflow_SFD_NoPart_List(Data* data, Data* device, int iter);


#endif /* FA_SFD_H_ */
