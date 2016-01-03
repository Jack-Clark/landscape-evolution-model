#ifndef watershedH
#define watershedH

#include "lem.h"
#include "config.h"
#include "flowroutines.h"
#include "SFD.h"

class watershed_edge
{
    private:
        int node1;
        int node2;
        HASH_TYPE weight;
    
	public:
        watershed_edge();
        watershed_edge(int n1, int n2, HASH_TYPE w);
        ~watershed_edge();
        int get_node1();
        int get_node2();
        void set_node1(int n);
        void set_node2(int n);
		void set_weight(HASH_TYPE w);
        HASH_TYPE get_weight() const;
        bool operator<(const watershed_edge &other);
};

class union_set
{
    private:
        int *parent;
        int *rank;
        int size;
    public:
        union_set(int s);
        ~union_set();
        bool in_set(int x);
        void make_set(int x);
        int find_root(int x);
        void merge(int a, int b);
};

__device__ double atomicMin(double* address, double val);

__global__ void init_watershed(int *mask, int *aspects, int *watershed_id, double *zs, float *shortest_paths, int *dx, int *dy, int *counter, int ncell_x, int ncell_y);

__global__ void final_watershed(int *mask, int *aspects, int *watershed_id, double *zs, float *shortest_paths, int *dx, int *dy, int *counter, int ncell_x, int ncell_y);

__global__ void init_watershed_sink(int *mask, int *aspects, int *watershed_id, double *zs, int *dx, int *dy, int ncell_x, int ncell_y, int *changed);

__global__ void identify_watershed(int *mask, int *aspects, int *watershed_id, int *dx, int *dy, int ncell_x, int ncell_y, int *changed);

__global__ void init_hash_table(HASH_TYPE *hash_table, int *table_counter, int table_size, int list_size);

__global__ void identify_watershed_boundary_write_back(int *mask, int *watershed_id, double *zs, HASH_TYPE *hash_table, int *table_counter, int table_size, int list_size, int *dx, int *dy, int ncell_x, int ncell_y, int *is_succeeded, int* counter);

__global__ void raise_watershed(int *mask, int *watershed_id, HASH_TYPE *raise, double *zs, int ncell_x, int ncell_y);

int init_and_sort_watershed_bound_labels(Data* data, HASH_TYPE *hash_table, int *table_counter, int table_size, int list_size, watershed_edge *sorted_list, int sorted_list_size, int *max_watershed_p);


#endif
