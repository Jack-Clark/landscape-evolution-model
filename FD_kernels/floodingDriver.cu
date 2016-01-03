#include "floodingdriver.h"


using namespace std;

void outputEdges(const char* filename, int table_size, int* table_counter, int* hash_table, int list_size) {
  FILE* out = fopen(filename, "w");
  for(int i = 0; i < table_size; i++)
  {
    if(table_counter[i] > 0)
      // go through each entry within the table
      for(int j = 0; j < (table_counter[i] + 1) * 3; j+=3)
      {
        fprintf(out, "Watershed 1: %d, Watershed 2: %d, Height: %d\n", hash_table[i * list_size * 3 + j], hash_table[i * list_size * 3 + j + 1], hash_table[i * list_size * 3 + j + 2]);
      }
  }
  fclose(out);
}

void outputEdges(const char* filename, int table_size, int* table_counter, double* hash_table, int list_size) {
  FILE* out = fopen(filename, "w");
  for(int i = 0; i < table_size; i++)
  {
    if(table_counter[i] > 0)
      // go through each entry within the table
      for(int j = 0; j < (table_counter[i] + 1) * 3; j+=3)
      {
        fprintf(out, "Watershed 1: %d, Watershed 2: %d, Height: %f\n", (int) hash_table[i * list_size * 3 + j], (int) hash_table[i * list_size * 3 + j + 1], hash_table[i * list_size * 3 + j + 2]);
      }
  }
  fclose(out);
}


void floodingDriver(dim3 dimGrid, dim3 dimBlock, Data* data, Data* device, Catchment* catchments, int ncell_x, int ncell_y, int cell_size, int iter) {
#ifdef OUTPUT
  printf("********************\n");
  printf("Computing watershed.\n");
  printf("********************\n");
#endif
  int block_ncell_y = 16;
  int block_ncell_x = 16;
  int *counter;
  int *counter_d;
  int count = 0;
  int *change_flag_d;
  int *change_flag_h;
//  size_t freenow, total;

  // ** Flats (0) can still exist - but these are now parts of sinks...
  // ** so deal with these  
  counter = (int*) malloc(sizeof(int)); //cleared
  *counter = 0;
  cudaMalloc((void**)&counter_d, sizeof(int));
  cudaMemcpy(counter_d, counter, sizeof(int), cudaMemcpyHostToDevice);

  // Identify each sink cell and assign it a watershed value
  init_watershed<<<dimGrid, dimBlock>>>(device->mask, device->fd, device->watershed_id, device->dem, device->shortest_paths, device->dx, device->dy, counter_d, ncell_x, ncell_y);
  fprintf(data->outlog,"FD: init_watershed :%s\n", cudaGetErrorString(cudaGetLastError()));
  change_flag_h = (int*) malloc(sizeof(int));//cleared
  cudaMalloc((void**) &change_flag_d, sizeof(int));

  do {
    *change_flag_h = 0; // clear flag
    cudaMemcpy(change_flag_d, change_flag_h, sizeof(int), cudaMemcpyHostToDevice);

    //** Process Watersheds
    //**********************************************************************
    // Reduce number of watersheds - neighbouring cells in the same plateau
    // adopt the higher watershed values...

    init_watershed_sink<<<dimGrid, dimBlock>>>(device->mask, device->fd, device->watershed_id, device->dem, device->dx, device->dy, ncell_x, ncell_y, change_flag_d);
    if(count <3) fprintf(data->outlog,"FD: init_watershed_sink :%s\n", cudaGetErrorString(cudaGetLastError()));
	cudaMemcpy(change_flag_h, change_flag_d, sizeof(int), cudaMemcpyDeviceToHost);

	count ++;

  } while(*change_flag_h == 1);  // while we're still doing things...

  fprintf(data->outlog,"FD: Processed watersheds in %d iterations\n", count);

  count = 0;
  do {
    // Reset change flag
    *change_flag_h = 0;
    cudaMemcpy(change_flag_d, change_flag_h, sizeof(int), cudaMemcpyHostToDevice);
    // *********************************************************************
    // ** identify watershed
    // *********************************************************************
    // Assign watershed values to cells that don't currently have values
    // assume same watershed value as cell I flow into...
    // Q: Does this allocate a watershed value to all cells in the DEM?

    identify_watershed<<<dimGrid, dimBlock>>>(device->mask, device->fd, device->watershed_id, device->dx, device->dy, ncell_x, ncell_y, change_flag_d);
    if(count <3) fprintf(data->outlog,"FD: identify_watersheds :%s\n", cudaGetErrorString(cudaGetLastError()));
    cudaMemcpy(change_flag_h, change_flag_d, sizeof(int), cudaMemcpyDeviceToHost);

    count ++;
  } while(*change_flag_h == 1);

  fprintf(data->outlog,"FD: Identified watersheds in %d iterations\n", count);
  // From here on in we're doing a hash table solution...
  // Is this a good idea...?

  // The values here need to be adjusted for each DEM size
  int table_size = 8000; //32767;
  int list_size = 800;
  HASH_TYPE *hash_table;
  HASH_TYPE *hash_table_d;
  int *is_succeeded, *is_succeeded_d, *table_counter, *table_counter_d;
  is_succeeded = (int*) malloc(sizeof(int)); //cleared
  *is_succeeded = 1;

  hash_table = (HASH_TYPE*) calloc(table_size*list_size*3, sizeof(HASH_TYPE)); //cleared
  table_counter = (int*) malloc(sizeof(int) * table_size); //cleared

//  int getmem = sizeof(HASH_TYPE) * table_size * list_size * 3;

  cudaMalloc((void **)&hash_table_d, sizeof(HASH_TYPE) * table_size * list_size * 3);
  cudaMalloc((void **)&table_counter_d, sizeof(int) * table_size);
  cudaMalloc((void **)&is_succeeded_d, sizeof(int));

  dim3 dimTableGrid(list_size/block_ncell_x + 1, table_size/block_ncell_y + 1);

  cudaMemcpy((void*) is_succeeded_d, (void*) is_succeeded, sizeof(int), cudaMemcpyHostToDevice);

  // Set all values in hash table to -1
  //** Initialise hash table...
  // Set all values in hash table to -1
  init_hash_table<<<dimTableGrid, dimBlock>>>(hash_table_d, table_counter_d, table_size, list_size);
  fprintf(data->outlog,"FD: init_hash_table :%s\n", cudaGetErrorString(cudaGetLastError()));

  // Place boundary watershed indexes into hash table along with heights
  identify_watershed_boundary_write_back<<<dimGrid, dimBlock>>>(device->mask, device->watershed_id, device->dem, hash_table_d, table_counter_d, table_size, list_size, device->dx, device->dy, ncell_x, ncell_y, is_succeeded_d, counter_d);
  fprintf(data->outlog,"FD: identify_watershed_boundary_write_back  :%s\n", cudaGetErrorString(cudaGetLastError()));
  // Copy hash table back to main memory...
  cudaMemcpy((void*) is_succeeded, (void*) is_succeeded_d, sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy((void*) counter, (void*) counter_d, sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy((void*) table_counter, (void*) table_counter_d, sizeof(int) * table_size, cudaMemcpyDeviceToHost);
  cudaMemcpy((void *) hash_table, hash_table_d, sizeof(HASH_TYPE) * table_size * list_size * 3, cudaMemcpyDeviceToHost);
  fprintf(data->outlog,"FD: idevice hash memcopies :%s\n", cudaGetErrorString(cudaGetLastError()));

  int edgeNum = 0;

  // go through each table entry... to count up how many pairs we have.
  for(int i = 0; i < table_size; i++) {
    if(table_counter[i] > 0)
      // go through each entry within the table
      for(int j = 0; j < (table_counter[i] + 1) * 3; j+=3) {
        edgeNum++;
      }
  }

  int sorted_list_size = 0;
  int max_watershed = -1;
  for(int i = 0; i < table_size; i++) {
    sorted_list_size += table_counter[i] + 1;
  }
  // create new array just large enough to store the edges
  fprintf(data->outlog,"FD: sorted list size1 = %d \n", sorted_list_size);

  watershed_edge *sorted_list = new watershed_edge[sorted_list_size];  // cleared

  sorted_list_size = init_and_sort_watershed_bound_labels(data, hash_table, table_counter, table_size, list_size, sorted_list, sorted_list_size, &max_watershed);
 
  int cccount = 0;
  int nncount = 0;
  for(int i = 0; i < sorted_list_size; i++) {
	  double height1 = sorted_list[i].get_weight();
	  int    label1a = sorted_list[i].get_node1();
      int    label2a = sorted_list[i].get_node2();
	  double raise_d1 =  height1;
	  if (raise_d1 == sorted_list[i].get_weight()) {  // we were an integer to start with
		  cccount ++;
	  }
	  else
		nncount ++;
  }

  union_set *set = new union_set(max_watershed + 1); //cleared

  fprintf(data->outlog,"FD: max_watershed %d \n", max_watershed);

  int *done;
  done = (int*) calloc(max_watershed + 1, sizeof(int)); //cleared

  HASH_TYPE *raise, *raise_d;
  raise = (HASH_TYPE*) calloc(max_watershed + 1, sizeof(HASH_TYPE)); //cleared
  
  cudaMalloc((void**)&raise_d, (max_watershed + 1) * sizeof(HASH_TYPE));
  cudaDeviceSynchronize();
  
  done[0] = 1;

  for(int i = 0; i < sorted_list_size; i++) {
    int label1, label2, root1, root2;
	HASH_TYPE height;
	
    label1 = sorted_list[i].get_node1();
    label2 = sorted_list[i].get_node2();
    height = (HASH_TYPE) sorted_list[i].get_weight();
    if(!(set->in_set(label1)))
      set->make_set(label1);
    if(!(set->in_set(label2)))
      set->make_set(label2);
    root1 = set->find_root(label1);
    root2 = set->find_root(label2);
    // If both labels are done, ignore the edge
    if(done[root1] && done[root2]) {
      continue;
    }
    // If only one label is done, assigne the other one to be done, and the other watershed is to be flooded
    if(done[root1] || done[root2]) {
      if(done[root1]) {
	done[root2] = 1;
	raise[root2] = height;
      }
      else {
	done[root1] = 1;
	raise[root1] = height;
      }

      continue;
    }
    //* If both labels are not done, merge them
    if(root1 != root2) {
      set->merge(root1, root2);
      raise[root2] = raise[root1] = height;
    }
  }
  for(int i = 0; i < max_watershed + 1; i++) {
    if(raise[i] != 0)
      {
	raise[i] = raise[set->find_root(i)];
      }
  }

  int ccount = 0;
  int ncount = 0;
  for(int ii = 0; ii < max_watershed; ii++) {
	  int raise_i = (int) raise[ii];
	  double raise_d = (double) raise_i;
	  if (raise_d == raise[ii]) {  // we were an integer to start with
		  if (raise_i > 0)
		  ccount ++;
	  }
	  else
		ncount ++;
  }

  cudaMemcpy((void*)raise_d, (void*) raise, (max_watershed + 1) * sizeof(HASH_TYPE), cudaMemcpyHostToDevice);

  raise_watershed<<<dimGrid, dimBlock>>>(device->mask, device->watershed_id, raise_d, device->dem, ncell_x, ncell_y);
  fprintf(data->outlog,"FD: raise_watershed  :%s\n", cudaGetErrorString(cudaGetLastError()));

  cudaDeviceSynchronize();

  singleFlowDir<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->SlopePtr, device->fd, cell_size, ncell_x, ncell_y, device->dx, device->dy, 1);
  fprintf(data->outlog,"FD: second single flow direction  :%s\n", cudaGetErrorString(cudaGetLastError()));

  flow_boundary<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->SlopePtr, device->fd, ncell_x, ncell_y, device->dx, device->dy);
  fprintf(data->outlog,"FD: second flow boundary  :%s\n", cudaGetErrorString(cudaGetLastError()));

  //printf("after flow boundary:%s\n", cudaGetErrorString(cudaGetLastError()));
 
  shortest_paths_plateaus_init<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->SlopePtr, device->fd, device->shortest_paths, cell_size, ncell_x, ncell_y, device->lowHeight);
  fprintf(data->outlog,"FD: second shortest_paths_plateaus  :%s\n", cudaGetErrorString(cudaGetLastError()));

  
  do {

    *change_flag_h = 0;
    cudaMemcpy(change_flag_d, change_flag_h, sizeof(int), cudaMemcpyHostToDevice);
    route_plateaus<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->fd, device->shortest_paths,  ncell_x, ncell_y, device->dx, device->dy, change_flag_d, device->lowHeight);
    cudaMemcpy(change_flag_h, change_flag_d, sizeof(int), cudaMemcpyDeviceToHost);
  } while(*change_flag_h == 1);
  fprintf(data->outlog,"FD: after route_plateaus:%s\n", cudaGetErrorString(cudaGetLastError()));

  // Now work out the slopes on cells in a plateaux
  slope_plateaus<<<dimGrid, dimBlock>>>(device->mask, device->dem, device->shortest_paths, ncell_x, ncell_y, device->lowHeight, device->SlopePtr, cell_size);
  fprintf(data->outlog,"FD: slope pleateaus  :%s\n", cudaGetErrorString(cudaGetLastError()));


//  Clean up memory usage

  free(counter);
  free(change_flag_h);
  free(is_succeeded);
  free(hash_table);
  free(table_counter);
  free(raise);
  free(done);

  delete [] sorted_list;
  delete  set;

  cudaFree(counter_d);
  cudaFree(change_flag_d);
  cudaFree(raise_d);
  cudaFree(hash_table_d);
  cudaFree(table_counter_d);
  cudaFree(is_succeeded_d);

  
}
