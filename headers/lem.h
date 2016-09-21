/*
 * lem.h
 *
 *  Created on: Feb 18, 2014
 *      Author: ndm12
 */

#ifndef LEM_H_
#define LEM_H_

// Skeleton outline of CUDA enabled Landscape Evolution Model
#include <string>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <time.h>
#include <float.h>
#include <math.h>

// The below files are part of the GDAL library, which needs to be downloaded and installed
#include "gdal.h" // note this is a different header file to the c++ version
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_core.h"
#include "ogr_srs_api.h"

// found in ${CUDA_HOME}/inc
#include <cuda.h>
#include <device_functions.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// found in ${CUDA_HOME}/samples/common/inc. These should be <> not "", and then included during compilation.
#include <helper_functions.h>
#include <helper_timer.h>
#include <helper_cuda.h>

#include <thrust/system_error.h>
#include <thrust/device_ptr.h>
#include <thrust/functional.h>
#include <thrust/reduce.h>
#include <thrust/copy.h>
#include <thrust/count.h>
#include <thrust/fill.h>
#include <thrust/extrema.h>
#include <thrust/execution_policy.h>



#endif /* LEM_H_ */
