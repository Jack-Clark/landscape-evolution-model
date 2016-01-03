#include "flowroutines.h"


//  Before MFD can be called sinks and flats must have been routed. Should they be carved for a number of iterations using SFD?
//  Slopes in all directions need to be stored for use in accumulation partitioning
//  In MFD we must abandon the 0-7 directions in favour of the byte notation

// NW  N NE    32  64   128
//  W  *  E    16   0   1
// SW  S SE     8   4   2



__global__ void MultipleFlowDir(int *mask, double *zs, double *slopes, int *aspects, int cell_size,
				int ncell_x, int ncell_y, int *dx, int *dy, int last)
{
  int cellx, celly; // oldx, oldy;
  int icol, irow, dcell;  // dtemp, xx, yy
  double smax, stemp; //runoff;
  float dc;
  int aspect;
  int oldaspect;
  int morethanSFD;
  double slope;
  int tempaspect;
  
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;
  if(icol >= ncell_x || irow >= ncell_y) // if outside of DEM - nothing to do
    return;
  
 // int self = irow * ncell_x + icol;
 // if (mask[self] != 1) {
//	  aspects[self] = 0;
//	  return; // don't calculate if not in catchment(s) of interest
 // }

  aspect = 0;
  morethanSFD = 0;

  //smax = -1.e38;
  for (dcell = 0; dcell < 8; dcell ++) {
	
    cellx = icol + dx [dcell];
    celly = irow + dy [dcell];
    if (cellx >= 0 && cellx < ncell_x && celly >= 0 && celly < ncell_y) { // for each of my neighbours
      if (dcell % 2 == 0)
        dc = 1.0;   //even directions are cardinal
      else
        dc = 1.41;  //odd directions are diagonals
      
	  slope = (zs[irow * ncell_x + icol] - zs[celly * ncell_x + cellx]) / (cell_size * dc);

	  // code for all downslope neighbours
	  if (slope > 0.0 ) {
		  aspect += code[dcell];
		  morethanSFD++;
	  }
      
	  // slopesMFD[irow * ncell_x + icol]  - need to think about how to store these
    }
  }

  if (morethanSFD < 2) {
	  oldaspect = code[(aspects[irow * ncell_x + icol] )] ; // recode existing values as we go
	  aspects[irow * ncell_x + icol] = oldaspect;
  }

  // apply MFD ONLY to cells which have a slope to more than one downslope member as SFD routing has already been done.
  
  if (morethanSFD > 1) aspects[irow * ncell_x + icol] = aspect; // only set for non-plateau cells
  
 }

//**************************************************
// Work out which direction I should flow in:

// NW  N NE    7  0  1
//  W  *  E    6  8  2
// SW  S SE    5  4  3


__global__ void singleFlowDir(int *mask, double *zs, double *slopes, int *aspects, int cell_size,
				int ncell_x, int ncell_y, int *dx, int *dy, int last)
{
  int cellx, celly; // oldx, oldy;
  int icol, irow, dcell;  // dtemp, xx, yy
  double smax, stemp; //runoff;
  float dc;
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;

  if(icol >= ncell_x || irow >= ncell_y) // if outside of DEM - nothing to do
    return;

  int self = irow * ncell_x + icol;

  if (mask[self] == 0) {
	  aspects[self] = 8;
	  return; // we can force the non catchment cells to go through plateau routing which will make them flow to the edge when edge routing has been done.
  }

  smax = -1.e38;
  for (dcell = 0; dcell < 8; dcell ++) {
    cellx = icol + dx [dcell];
    celly = irow + dy [dcell];
    if (cellx >= 0 && cellx < ncell_x && celly >= 0 && celly < ncell_y) { // for each of my neighbours
      if (dcell % 2 == 0)
        dc = 1.0;   //even directions are cardinal
      else
        dc = 1.41;  //odd directions are diagonals
      stemp = (zs[irow * ncell_x + icol] - zs[celly * ncell_x + cellx]) / (cell_size * dc);

      if ((mask[celly * ncell_x + cellx]) != 1) stemp = -1.0; // if the target cell is outside the catchment(s) then set temp slope to negative

      if (stemp > smax && stemp > 0) { // if this slope is bigger than any we've found so far...
        smax = stemp;
	    slopes[irow * ncell_x + icol] = stemp;
        aspects[irow * ncell_x + icol] = dcell;
      }
    }
  }
  if (smax == -1.e38) { // no outlet cells 
    aspects[irow * ncell_x + icol] = 8;
    slopes[irow * ncell_x + icol] = 0.0;
  }
}




//**************************************************************
// Route water over plateaus
//**************************************************************
// find the cell around self which has same height and a lower distance to the exit 
// route to that cell
// Q: Is this thread safe?
// Could one cell update happening at the same time as another cause erronious computation?

__global__ void route_plateaus(int *mask, double *zs, int *aspects, float* shortest_paths,
				int ncell_x, int ncell_y, int *dx, int *dy, int *change_flag, double *lowHeight)
{
	int irow, icol, dcell, cellx, celly;
	float dc, dis;
	irow = blockIdx.y * blockDim.y + threadIdx.y;
	icol = blockIdx.x * blockDim.x + threadIdx.x;

  // If outside of DEM return (nothing to do)
	if(icol >= ncell_x || irow >= ncell_y)
		return;

	int self = irow * ncell_x + icol;
	if (mask[self] !=1) return; // don't calculate if not in catchment(s) of interest

  // Get the current shortest path for this cell
	float min_distance = shortest_paths[irow * ncell_x + icol];
	
	// If it's zero then we're not part of a plateau so finish
	if(min_distance == 0)
		return;

	int flow_nei_idx = -1;
	double this_ele = zs[irow * ncell_x + icol];
	for (dcell = 0; dcell < 8; dcell ++)
	{
		cellx = icol + dx [dcell];
		celly = irow + dy [dcell];
		if (cellx >= 0 && cellx < ncell_x && celly >= 0 && celly < ncell_y)
		{
			dc = 1.0;
			if (dcell % 2 == 0)
			{
				dc = 1.0;   //even directions
										//are cardinal
			}
			else
			{
			  dc = 1.41;  //odd directions are diagonals
			}
			dis = shortest_paths[celly * ncell_x + cellx] + dc;
			
			if (this_ele == zs[celly * ncell_x + cellx] && min_distance > dis) 
		    {
				flow_nei_idx = dcell;
				min_distance = dis;
				*change_flag = 1;  // don't care how many cells have changed - just that cells have changed
				aspects[irow * ncell_x + icol] = flow_nei_idx;  // record aspect to this lower cell and the new distance
				shortest_paths[irow * ncell_x + icol] = min_distance;
				lowHeight[irow * ncell_x + icol] = lowHeight[celly * ncell_x + cellx];
			}
		}
	}
}


//**************************************************************
// Compute the slope for cells in the plateau
//**************************************************************
__global__ void slope_plateaus(int *mask, double *zs, float* shortest_paths, int ncell_x, int ncell_y, double *lowHeight, double *slopes, int cell_size)
{
  int irow, icol; // dcell, cellx, celly
  //float dis; // dc
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;

  // If outside of DEM return (nothing to do)
  if(icol >= ncell_x || irow >= ncell_y)
    return;

  int self = irow * ncell_x + icol;
  if (mask[self] != 1) return; // don't calculate if not in catchment(s) of interest

  // Get the current shortest path for this cell
  float min_distance = shortest_paths[irow * ncell_x + icol];

  // if we're not part of the plateau then nothing to do
  if (min_distance == 0)
    return;

  // Get the height at which this cell flows to (when it leaves the plateau)
  double low = lowHeight[irow * ncell_x + icol];
  
  // Get the height of this cell
  double height = zs[irow * ncell_x + icol];
  
  // The slope is (height - low) / (min_distance * cell_size)
  slopes[irow * ncell_x + icol] = (height - low) / ((min_distance * cell_size)*10000000);
  //slopes[irow * ncell_x + icol] = (height - low) / ((min_distance * cell_size));
  //slopes[irow * ncell_x + icol] = DBL_MIN;
}


__global__ void remap(int *mask, int *aspects, int ncell_x, int ncell_y)
{
  int irow, icol; // dcell, cellx, celly
  //float dis; // dc
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;

    // If outside of DEM return (nothing to do)
  if(icol >= ncell_x || irow >= ncell_y)
    return;

  //int self = irow * ncell_x + icol;
  //if (mask[self] != 1)
//	  {
//	  aspects[self] = 0;
//	  return; // don't calculate if not in catchment(s) of interest
//	  }


  aspects[irow * ncell_x + icol] = code [aspects[irow * ncell_x + icol] ];

}




//**************************************************************
//** Compute shortest paths in sinks
//**************************************************************
// Same as route water over plateaus except we don't change the 
// aspect
// Q: Do we need to do this again? The contentse of shortest_paths
// at the end of this will be the same as from the previous 
// run...?
// Q: Depth first works on plateaus as at least one cell in the 
// plateau will have a 0 value for shortest_path (ie this is the 
// exit cell from the plateau. However, for a sink there is no cell
// in the sink with shortest_path value of 0 as all cells around the 
// sink are higher than the cells in the sink. Will this do anything?

__global__ void comp_shortest_paths_sink(int *mask, double *zs, int *aspects, float* shortest_paths,
				int ncell_x, int ncell_y, int *dx, int *dy, int *change_flag)
{
	int irow, icol, dcell, cellx, celly;
	float dc, dis;
	irow = blockIdx.y * blockDim.y + threadIdx.y;
	icol = blockIdx.x * blockDim.x + threadIdx.x;

	if(icol >= ncell_x || irow >= ncell_y) // if we're outside of the DEM - nothing to do so return
		return;

	//int self = irow * ncell_x + icol;
	//if (mask[self] != 1) return; // don't calculate if not in catchment(s) of interest



	float min_distance = shortest_paths[irow * ncell_x + icol];
	if(min_distance == 0)  // if we're not part of a flat
	  return;

	int flow_nei_idx = -1;
	int this_ele = zs[irow * ncell_x + icol];
	for (dcell = 0; dcell < 8; dcell ++)
	{
		cellx = icol + dx [dcell];
		celly = irow + dy [dcell];
		if (cellx >= 0 && cellx < ncell_x && celly >= 0 && celly < ncell_y)
		{
			dc = 1.0;
			if (dcell % 2 == 0)
			{
			  dc = 1.0;   //even directions
							      //are cardinal
			}
			else
			{
				dc = 1.41;  //odd directions
                    //are diagonals
			}
			dis = shortest_paths[celly * ncell_x + cellx] + dc;
			if((((this_ele == zs[celly * ncell_x + cellx]) ||
			  this_ele - zs[celly * ncell_x + cellx] < 0.000001 && this_ele > zs[celly * ncell_x + cellx]) || (this_ele - zs[celly * ncell_x + cellx] > -0.000001 &&
				this_ele < zs[celly * ncell_x + cellx])) && min_distance > dis)
				// Q: again - can this be replaced with:
				// if (abs(this_ele - zs[celly * ncell_x + cellx) < 0.000001 && min_distance > dis)
			{
				flow_nei_idx = dcell;
				min_distance = dis;
				*change_flag = 1;  // don't care which one has changed - just that there's been change
				shortest_paths[irow * ncell_x + icol] = min_distance;
			}
		}
	}
}


//***************************************************
// Set up for working out shortest paths on plateaus
//*****
// If cell is not a plateau (aspect = 8) then set shortest path to 0
// else set shortest path to large number
// Q: is +ve infinity a valid number in CUDA - if so better to use that rather than DBL_MAX

__global__ void shortest_paths_plateaus_init(int *mask, double *zs, double *slopes, int *aspects, float* shortest_paths, int cell_size,
				int ncell_x, int ncell_y, double *lowHeight)
{
	int irow, icol;
	irow = blockIdx.y * blockDim.y + threadIdx.y;
	icol = blockIdx.x * blockDim.x + threadIdx.x;
	if(icol >= ncell_x || irow >= ncell_y)
		return;
	
	//int self = irow * ncell_x + icol;
	//if (mask[self] != 1) return; // don't calculate if not in catchment(s) of interest


	// Set the lowHeight to 0 for each cell
	lowHeight[irow * ncell_x + icol] = 0.0;
	if(aspects[irow * ncell_x + icol] != 8)  
	{
		shortest_paths[irow * ncell_x + icol] = 0;
		return;
	}
	else {
//	shortest_paths[irow * ncell_x + icol] = DBL_MAX - 1.0;
		shortest_paths[irow * ncell_x + icol] = 20000000;}
}


//********************************************
// Set what happens on the boundary of the DEM
//********************************************

__global__ void flow_boundary(int *mask, double *zs, double *slopes, int *aspects, int ncell_x, int ncell_y, int *dx, int *dy) //directions altered to correspond better with slim
{
    int irow, icol, dcell, cellx, celly;
    irow = blockIdx.y * blockDim.y + threadIdx.y;
        icol = blockIdx.x * blockDim.x + threadIdx.x;
        
        if(icol >= ncell_x || irow >= ncell_y)
        return;
        double  smax ;
        int dtemp;
        
        if(icol == 0) //western boundary
        {
                smax = 0.00000000001; 
                dtemp = 6 ; // make it flow out if no  neighbour in the grid is lower i.e. west
                for (dcell = 0; dcell < 5; dcell ++)
                        {
                                cellx = icol + dx [dcell];
                                celly = irow + dy [dcell];
                        
                                if (slopes[irow * ncell_x + icol] > smax)
                                {
                                        smax = slopes[irow * ncell_x + icol];
                                        dtemp = aspects[irow * ncell_x + icol];
                                }
                        }
                
                aspects[irow * ncell_x + icol] = dtemp;
                slopes[irow * ncell_x + icol] = smax;
        }

        if(icol == ncell_x - 1) //Eastern boundary
    {
                smax = 0.00000000001; 
                dtemp = 2 ; // if no lower neighbours send it east

                for (dcell = 4; dcell < 8; dcell ++) // was 7 i.e. could not go north!
                        {
                                cellx = icol + dx [dcell];
                                celly = irow + dy [dcell];
                        
                                if (slopes[irow * ncell_x + icol] > smax)
                                {
                                        smax = slopes[irow * ncell_x + icol];
                                        dtemp = aspects[irow * ncell_x + icol];
                                }
                        }
                
                aspects[irow * ncell_x + icol] = dtemp;
                slopes[irow * ncell_x + icol] = smax;
        }

    if(irow == 0) //Northern Boundary
        {
                smax = 0.00000000001; 
                dtemp = 0 ; // if no lower neighbours send it north

                for (dcell = 3; dcell < 6; dcell ++) // should be 2 - cannot go east!
                        {
                                cellx = icol + dx [dcell];
                                celly = irow + dy [dcell];
                        
                                if (slopes[irow * ncell_x + icol] > smax)
                                {
                                        smax = slopes[irow * ncell_x + icol];
                                        dtemp = aspects[irow * ncell_x + icol];
                                }
                        }
                aspects[irow * ncell_x + icol] = dtemp;
                slopes[irow * ncell_x + icol] = smax;
        }

        if (irow == ncell_y-1) // Southern Boundary 
        {
                smax = 0.00000000001; //make sure there is a way out
                dtemp = 4; // same as original! if no lower neighbours send it south

                for (dcell = 0; dcell < 3; dcell ++) // was <2
                        {
                                cellx = icol + dx [dcell];
                                celly = irow + dy [dcell];
                        
                                if (slopes[irow * ncell_x + icol] > smax)
                                {
                                        smax = slopes[irow * ncell_x + icol];
                                        dtemp = aspects[irow * ncell_x + icol];
                                }
                        }
                for (dcell = 6; dcell < 8; dcell ++) // was dcell = 7;
                        {
                                cellx = icol + dx [dcell];
                                celly = irow + dy [dcell];
                                if (slopes[irow * ncell_x + icol] > smax)
                                {
                                        smax = slopes[irow * ncell_x + icol];
                                        dtemp = aspects[irow * ncell_x + icol];
                                }
                        }

                aspects[irow * ncell_x + icol] = dtemp;
                slopes[irow * ncell_x + icol] = smax;
        }
        
}
