#include "watershed.h"

using namespace std;

watershed_edge::watershed_edge()
{
	node1 = 0;
	node2 = 0;
	weight = 0.0;
}

watershed_edge::~watershed_edge()
{
}

watershed_edge::watershed_edge(int n1, int n2, HASH_TYPE w)
{
    node1 = n1;
    node2 = n2;
    weight = w;
}

int watershed_edge::get_node1()
{
    return node1;
}

int watershed_edge::get_node2()
{
    return node2;
}

HASH_TYPE watershed_edge::get_weight() const
{
    return weight;
}

void watershed_edge::set_node1(int n)
{
    node1 = n;
}

void watershed_edge::set_node2(int n)
{
    node2 = n;
}

void watershed_edge::set_weight(HASH_TYPE w)
{
    weight = w;
}

bool watershed_edge::operator<(const watershed_edge &other)
{
    return ((*this).weight < other.weight);
}

union_set::union_set(int s)
{
    size = s;
    parent = (int*) calloc(size, sizeof(int));
    rank = (int*) calloc(size, sizeof(int));
}

union_set::~union_set()
{
    free(parent);
    free(rank);
}

bool union_set::in_set(int x)
{
    if(x == 0)
        return true;
    if(x > 0 && x < size)
        if(parent[x] > 0)
            return true;
    return false;
}

int union_set::find_root(int x)
{
    if(parent[x] != x)
        return find_root(parent[x]);
    else return x;
}

void union_set::make_set(int x)
{
    if(x > 0 && x < size)
    {
        parent[x] = x;
        rank[x] = 0;
    }
}

void union_set::merge(int a, int b)
{
    int aroot = find_root(a);
    int broot = find_root(b);
    if(aroot == broot)
        return;
    if(rank[aroot] > rank[broot])
        parent[broot] = aroot;
    else
    {
        parent[aroot] = broot;
        if(rank[aroot] == rank[broot])
            rank[broot]++;
    }
}


__device__ double atomicMin(double* address, double val) {
	// HACK!!
	// As there is currently no atomic min for doubles use the following...
	double min = *address;
	double assumed;
	double old = *address;
    //double val = *ptr;
    do {       // If we have a value lower than that stored at the moment
		assumed = old;
		min = fmin(old, val);
		if (min == old)
			break;
        old = __longlong_as_double(atomicCAS((unsigned long long int*)address, __double_as_longlong(assumed), __double_as_longlong(min)));
    //    // When finished if it worked val will contain the old value and ptr will hold the new low value
    //	// if it fails - something else got in there first! - we get the new value held in prt back
    } while (assumed != old); // keep going until val <= min
    return old;
}
//*******************************************************
//** Initialise for Watersheds
//*******************************************************
// If we're a flat cell then we're part of a sink
// Use own cell id as watershed value
// If cell is not part of a watershed point at the cell 
// this cell flows into

__global__ void init_watershed(int *mask, int *aspects, int *watershed_id, double *zs, float *shortest_paths, int *dx, int *dy, int *counter, int ncell_x, int ncell_y)
{
  int irow, icol;
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;
  if(icol >= ncell_x || irow >= ncell_y) // if outside of DEM nothing to do
    return;

	int self = irow * ncell_x + icol;
	if (mask[self] != 1)
	{
	    watershed_id[self] = 0;
		return;
	}
		
  int this_idx = irow * ncell_x + icol;
  //watershed_id[this_idx] = -1;
  if((icol == 0) || (icol == ncell_x - 1) || (irow == 0) || (irow == ncell_y - 1)) // if we're on the outside of the DEM use a special watershed (ncell_x*ncell_y)
  {
    watershed_id[this_idx] = 0;
	return;
  }

  if(aspects[this_idx] == 8) // this is a sink as all plateaus are now routed 
    watershed_id[this_idx] = atomicAdd(counter,1) + 1;
  else {
    // find cell we flow into
    int dir = aspects[this_idx];
    int cellx = icol + dx[dir];
    int celly = irow + dy[dir];
    // check we are still in the DEM
    if (cellx < 0 || cellx >= ncell_x || celly < 0 || celly >= ncell_y)
      watershed_id[this_idx] = ncell_x * ncell_y; // use special watershed
    else 
      // save this as a negative number so we don't confuse it with a watershed
      // each value is shifted down by 1 so that index 0 isn't used (it could be a 
      // watershed in the DEM)
      watershed_id[this_idx] = -(celly * ncell_x + cellx + 1);
  }
}


// This routine is as above but does not treat the outer cells as a 'special case' 
__global__ void final_watershed(int *mask, int *aspects, int *watershed_id, double *zs, float *shortest_paths, int *dx, int *dy, int *counter, int ncell_x, int ncell_y)
{
  int irow, icol;
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;
  if(icol >= ncell_x || irow >= ncell_y) // if outside of DEM nothing to do
    return;

	int self = irow * ncell_x + icol;
	if (mask[self] != 1) return; // don't calculate if not in catchment(s) of interest

		
  int this_idx = irow * ncell_x + icol;
  if(aspects[this_idx] == 8) // this is a sink as all plateaus are now routed 
    watershed_id[this_idx] = atomicAdd(counter,1) + 1;
  else {
    int dir = aspects[this_idx];
    int cellx = icol + dx[dir];
    int celly = irow + dy[dir];
    if (cellx < 0 || cellx >= ncell_x || celly < 0 || celly >= ncell_y)
      watershed_id[this_idx] = ncell_x * ncell_y; // use special watershed
    else 
      watershed_id[this_idx] = -(celly * ncell_x + cellx + 1);
  }
}


//**********************************************************
// Process Watersheds
//**********************************************************
// If any cell around me part of the same plateau has a higher
// watershed index use that one instead of mine.
// Q: Is this thread safe? What if two neigbouring cells update
// at the same time?

__global__ void init_watershed_sink(int *mask, int *aspects, int *watershed_id, double *zs, int *dx, int *dy, int ncell_x, int ncell_y, int *changed) {
  int irow, icol;
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;
  if(icol >= ncell_x || irow >= ncell_y) // Outside DEM - return as no work to do...
    return;


	int self = irow * ncell_x + icol;
	if (mask[self] != 1) return; // don't calculate if not in catchment(s) of interest


  int this_idx = irow * ncell_x + icol;
  if(aspects[this_idx] == 8) { // if we're part of a flat 
    int w_id = -1;
    for (int dcell = 0; dcell < 8; dcell ++) { // for each of the cells around me...
      int cellx = icol + dx [dcell];
      int celly = irow + dy [dcell];
      if (cellx >= 0 && cellx < ncell_x && celly >= 0 && celly < ncell_y) { // check we're still inside of the DEM
        int idx = celly * ncell_x + cellx;
        if(zs[idx] == zs[this_idx] && watershed_id[idx] > watershed_id[this_idx]) {
        // if there's a watershed next to me (within the same plateau) with a higher index then use the index of that watershed
          w_id = watershed_id[celly * ncell_x + cellx];
          *changed = 1; // something has chaneged...
        }
      }
    }
    if(w_id != -1)
      watershed_id[this_idx] = w_id; // I've got a new watershed value - so update...
  }
}

//*****************************************************************
//** Identify Watersheds
//*****************************************************************
// for each cell which doesn't have a watershed yet...
// assign watershed value for cell I flow into (if assigned)
// if we don't have a watershed value from this point at the 
// cell that the cell I point at points to...
__global__ void identify_watershed(int *mask, int *aspects, int *watershed_id, int *dx, int *dy, int ncell_x, int ncell_y, int *changed)
{
  int irow, icol;
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;
  if(icol >= ncell_x || irow >= ncell_y)   // If outside of DEM nothing to do...
    return;

//	int self = irow * ncell_x + icol;
//	if (mask[self] != 1) return; // don't calculate if not in catchment(s) of interest

  //int this_idx = irow * ncell_x + icol;
  if(watershed_id[irow * ncell_x + icol] < 0) // if I don't have a watershed value yet...
  { // point to the cell that the cell I'm pointing at points to
    // each watershed is shifted down by 1 to deal with 0 - which could be a watershed
    watershed_id[irow * ncell_x + icol] = watershed_id[-watershed_id[irow * ncell_x + icol]-1];
    // this will either be a watershed value or a pointer to another cell (closer to the watershed)
    *changed = 1; // something updated so need to go around again
  }
}




//****************************************************************
//** Initialise the hash table
//****************************************************************
//
__global__ void init_hash_table(HASH_TYPE *hash_table, int *table_counter, int table_size, int list_size)
{
  int irow, icol;
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;
  if(icol > list_size || irow > table_size) // If I'm not in the space for the hash table return...
    return;

  hash_table[irow * list_size * 3 + icol] = -1.0;     // set all values to -1
  hash_table[irow * list_size * 3 + icol + 1] = -1.0;
  hash_table[irow * list_size * 3 + icol + 2] = -1.0;
  table_counter[irow] = -1;
}

//****************************************************************
//** Identify watershed boundaries
//****************************************************************
// If we're part of a boundary between two watersheds:
// Try to find boundary indexes in the hash table - if present ensure height for pair is minimum of old value and this value
// If not in hash table, add to hash table

__global__ void identify_watershed_boundary_write_back(int *mask, int *watershed_id, double *zs, HASH_TYPE *hash_table, int *table_counter, int table_size, int list_size, int *dx, int *dy, int ncell_x, int ncell_y, int *is_succeeded, int* counter)
{
  int irow, icol;
  irow = blockIdx.y * blockDim.y + threadIdx.y;
  icol = blockIdx.x * blockDim.x + threadIdx.x;
  if(icol >= ncell_x || irow >= ncell_y) // If I'm outside of the DEM nothing to do so return...
    return;

	int self = irow * ncell_x + icol;
	if (mask[self] != 1) return;

  int dcell, cellx, celly;
  int this_idx = irow * ncell_x + icol;
  double height = -1.0;
  //*is_succeeded = 0;

  for (dcell = 0; dcell < 8; dcell ++) // for each cell around me...
  {
    cellx = icol + dx [dcell];
    celly = irow + dy [dcell];
    if (cellx >= 0 && cellx < ncell_x && celly >= 0 && celly < ncell_y)
    {
      //atomicAdd(counter, 1);
      int idx = celly * ncell_x + cellx;
      int this_id = watershed_id[this_idx]; // my watershed index
      int nei_id = watershed_id[idx];       // watershed index of neighbour
      if(this_id != nei_id) // if I'm part of a boundary between watersheds...
      {
        if(zs[this_idx] >= zs[idx]) // if I'm the higher cell of the boundary
          height = zs[this_idx];// store my height
        else
	  continue; // Ignore here - when neighbour is central cell it will be processed there...
	//height = zs[idx]; // store height of other cell
        atomicAdd(counter, 1);
        int hash_val = (this_id + nei_id) % table_size;  // compute hash value for this boundary
        //int *list = hash_table[hash_val];
        int found = 0;
        // try to find this hash value in those stored already...
        for(int i = 0; i < table_counter[hash_val] * 3; i += 3)
        {
          int label1 = (int) hash_table[hash_val * (list_size * 3) + i];
          int label2 = (int) hash_table[hash_val * (list_size * 3) + i + 1];
          //HASH_TYPE edge_height = (HASH_TYPE) hash_table[hash_val * (list_size * 3) + i + 2];
          if((label1 == this_id && label2 == nei_id) || (label1 == nei_id && label2 == this_id))
          {
            // if this boundary pair exists already...
            // make sure we have the lower height of the previous one and this one.
            atomicMin(hash_table + (hash_val * (list_size * 3) + i + 2), height);
            found = 1; // we've found this pair
          }
        }
        if(found == 0) // did we find this?
        {  // no
          // ask for a unique storage space...
          int list_idx = (atomicAdd(table_counter + hash_val, 1) + 1) * 3;
          if(list_idx < list_size * 3) // is there enough space to store another item in the hash?
          { // yes
            //list[list_idx] = this_id;
            //list[list_idx + 1] = nei_id;
            //list[list_idx + 2] = height;
            // Store the three peices of information about this watershed boundary
            hash_table[hash_val * (list_size * 3) + list_idx] = (HASH_TYPE) this_id;
            hash_table[hash_val * (list_size * 3) + list_idx + 1] = (HASH_TYPE) nei_id;
            hash_table[hash_val * (list_size * 3) + list_idx + 2] = (HASH_TYPE) height;
          }
          else
           *is_succeeded = 0; // mark that we couldn't store this data...
        } 
      }
    }
  }
  //*is_succeeded = 0;
}

__global__ void raise_watershed(int *mask, int *watershed_id, HASH_TYPE *raise, double *zs, int ncell_x, int ncell_y)
{
    int irow, icol;
    irow = blockIdx.y * blockDim.y + threadIdx.y;
	icol = blockIdx.x * blockDim.x + threadIdx.x;
	if(icol >= ncell_x || irow >= ncell_y)
        return;

	int self = irow * ncell_x + icol;
	if (mask[self] != 1) return; // don't calculate if not in catchment(s) of interest

    int this_idx = irow * ncell_x + icol;
    int label = watershed_id[this_idx];
    if(zs[this_idx] < raise[label])
        zs[this_idx] = raise[label];
}

bool operator<(const watershed_edge& a, const watershed_edge &b)
{
    return (a.get_weight()) < (b.get_weight());
}

//*******************************************************************
//** Initialise and sort watershed boundary labels
//*******************************************************************
//
int init_and_sort_watershed_bound_labels(Data* data, HASH_TYPE *hash_table, int *table_counter, int table_size, int list_size, watershed_edge *sorted_list, int sorted_list_size, int *max_watershed_p)
{
  //printf("Allocate memory for edge list.\n");
  int max_watershed = 0;
  int max_ws_index = -1;
  int current_idx = 0;
  //watershed_edge* memory_index_of_last;
  // store the labels of all watersheds that have been found
  
  //int* watershed_labels = (int*) malloc(sizeof(int) * sorted_list_size);
    int* watershed_labels = (int*) calloc (sorted_list_size, sizeof(int));


  //printf("Memory allocated for watersheds = %d\n", sorted_list_size);
  for(int i = 0; i < table_size; i++)
  {
    int last_start = current_idx;
    int edge_num = table_counter[i] + 1;
    //printf("i: %d\n", i);
    //printf("i = %d, Boundary to be exceeded? %d\n", i, current_idx + counter > sorted_list_size);
    for(int j = 0; j < edge_num * 3; j += 3)
    {
      //printf("i = %d j = %d\n", i, j);
      
	  int n1 =  (int) hash_table[i * list_size * 3 + j];  // added casting 
      int n2 =  (int) hash_table[i * list_size * 3 + j + 1];
	  
      HASH_TYPE w = hash_table[i * list_size * 3 + j + 2];

	  	  if (n1 < 0 ) 
	  {
		  printf("n1 < 0");
//		  _getch();
	  }



      // perform a merge sort to see if this value already exists in the list...
      int placed = 0;
      for (int pos = last_start; pos < current_idx; pos++) {
        if (n1 == sorted_list[pos].get_node1() && n2 == sorted_list[pos].get_node2() ||
            n2 == sorted_list[pos].get_node1() && n1 == sorted_list[pos].get_node2() ) 
		{
          // This pair already exists in the list so just check to see if we have a 
          // lower value for weight
	      //printf("found (%d,%d,%d) at %d\n", n1, n2, w, pos);
          if (w < sorted_list[pos].get_weight()) {
            sorted_list[pos].set_weight(w);
          } // otherwise the value is larger. In eather case we're finished looking
          placed = 1;
          break;
		}
      }

      // if we didn't find the watershed pair in there already...
      if (placed == 0) 
	  {
        sorted_list[current_idx].set_node1(n1);
        sorted_list[current_idx].set_node2(n2);
        sorted_list[current_idx].set_weight(w);
        //memory_index_of_last = &sorted_list[current_idx];

        current_idx++;

        // if placed == 1 then we've already stored these
        int found_n1 = 0;
        int found_n2 = 0;
        for (int pp = 0; pp < max_watershed; pp++) 
		{
          if (watershed_labels[pp] == n1)
            found_n1 = 1;
          if (watershed_labels[pp] == n2)
            found_n2 = 1;
          if (found_n1 && found_n2)
            break;
		}
        if (found_n1 == 0) 
		{
			watershed_labels[max_watershed] = n1;
			max_watershed++;
				if (max_ws_index < n1)
				max_ws_index = n1;
		}
		if (found_n2 == 0) 
		{
			watershed_labels[max_watershed] = n2;
			max_watershed++;
			if (max_ws_index < n2)
			 max_ws_index = n2;
		}
      }
    }
  }
  
  *max_watershed_p = max_ws_index;
  fprintf(data->outlog,"FD: Sort edge list.\n");
  fprintf(data->outlog,"FD: Sorted_list_size: %d, Current_idx: %d\n", sorted_list_size, current_idx);
  //qsort(sorted_list, sorted_list_size, sizeof(watershed_edge*), compare_watershed_edge);
 // printf("Memory_index_of_last = %u\n", memory_index_of_last);
  fprintf(data->outlog,"FD: sorted_list + sorted_list_size = %u\n", sorted_list + sorted_list_size);
  fprintf(data->outlog,"FD: current idx = : %d \n", current_idx);

  sort(sorted_list, sorted_list + current_idx);
  
  // check the sorted list
   for(int i = 0; i < sorted_list_size; i++) 
   {
	  double height1 = sorted_list[i].get_weight();
	  int    label1a = sorted_list[i].get_node1();
      int    label2a = sorted_list[i].get_node2();

	  if (label1a < 0 ) 
	  {
		  fprintf(data->outlog,"FD: label1 < 0");
//		  _getch();
	  }
   }

   fprintf(data->outlog,"FD: List sorted.\n");
//  for (int i = 0; i < current_idx; i++) {
  //  printf("[%d] W1: %d W2: %d H: %f\n", i, sorted_list[i].get_node1(), sorted_list[i].get_node2(), sorted_list[i].get_weight());
// }
   fprintf(data->outlog,"FD: Number of watersheds = %d\n", max_watershed);
   data->sinkcount = max_watershed;
  
  //_getch();
  free(watershed_labels);
  return current_idx;
}

