#ifndef __GASAL_H__
#define __GASAL_H__


#include <cstdlib>
#include <cstdint>
#include "/usr/local/cuda-9.2/targets/x86_64-linux/include/cuda_runtime.h"

#ifndef HOST_MALLOC_SAFETY_FACTOR
#define HOST_MALLOC_SAFETY_FACTOR 5
#endif

#define CHECKCUDAERROR(error, cuda_api_call) \
		do{\
			cudaError_t err = error;\
			if (cudaSuccess != err ) { \
				fprintf(stderr, "[GPUseed CUDA ERROR:] Error in %s: %s(CUDA error no.=%d). Line no. %d in file %s\n", cuda_api_call, cudaGetErrorString(err), err,  __LINE__, __FILE__); \
				exit(EXIT_FAILURE);\
			}\
		}while(0)\

#ifdef __CUDA_ARCH__
#define CHECKCUBERROR(error, cub_call) \
		do{\
			err = error;\
			if (cudaSuccess != err ) { \
				printf("[GPUseed CUB ERROR:] Error in CUB library function %s: ", cub_call); \
				CuDebugExit(err); \
			}\

		}while(0)\

#else
#define CHECKCUBERROR(error, cub_call) \
		do{\
			cudaError_t err = error;\
					if (cudaSuccess != err ) { \
						printf("[GPUseed CUB ERROR:] Error in CUB library function %s: ", cub_call); \
						CuDebugExit(err); \
					}\
		}while(0)\

#endif

inline int CudaCheckKernelLaunch()
{
	cudaError err = cudaGetLastError();
	if ( cudaSuccess != err )
	{
		return -1;
	}

	return 0;
}


enum comp_start{
	WITHOUT_START,
	WITH_START
};

// Generic enum for ture/false. Using this instead of bool to generalize templates out of Int values for secondBest.
// Can be usd more generically, for example for WITH_/WITHOUT_START.
enum Bool{
	FALSE,
	TRUE
};

enum data_source{
	NONE,
	QUERY,
	TARGET,
	BOTH
};

enum algo_type{
	UNKNOWN,
	GLOBAL,
	SEMI_GLOBAL,
	LOCAL,
	BANDED,
	MICROLOCAL,
	FIXEDBAND
};

enum operation_on_seq{
	FORWARD_NATURAL,
	REVERSE_NATURAL,
	FORWARD_COMPLEMENT,
	REVERSE_COMPLEMENT,
};

// data structure of linked list to allow extension of memory on host side.
struct host_batch{
	uint8_t *data;
	uint32_t offset;
	struct host_batch* next;
};
typedef struct host_batch host_batch_t;

// Data structure to hold results. Can be instantiated for host or device memory (see res.cpp)
struct gasal_res{
	int32_t *aln_score;
	int32_t *query_batch_end;
	int32_t *target_batch_end;
	int32_t *query_batch_start;
	int32_t *target_batch_start;
};
typedef struct gasal_res gasal_res_t;

//stream data
typedef struct {
	uint8_t *unpacked_query_batch;
	uint32_t *packed_fow_query_batch;
	uint32_t *packed_rev_query_batch;
	uint32_t *query_batch_offsets;
	uint32_t *query_batch_lens;
	uint32_t *gpu_thread_query_num;
	uint32_t *gpu_thread_query_idx;
	uint32_t *n_smems_fow;
	uint32_t *n_smems_rev;
	uint32_t *n_smems_fow_rev;
	uint32_t *n_seed_occ;
	uint32_t *n_smems_sum;
	uint64_t *n_seed_occ_sum;
	uint32_t *n_queries_processed;

	uint2 *seed_intervals;
	int2 *seed_query_pos;
	uint2 *seed_intervals_compact;
	int2 *seed_query_pos_compact;
	uint32_t *n_query_seeds;
	uint32_t *host_n_query_seeds;
	uint32_t *n_query_seeds_scan;
	uint32_t *seed_ref_pos;
	uint32_t *host_seed_ref_pos;
	uint32_t *seed_sa_idx;

	uint32_t *n_smems_max;

	uint32_t *n_smems_scan;
	uint32_t *n_seed_occ_scan;
	uint32_t *is_smem_flag;

	void *cub_scan_temp_storage;
	size_t cub_scan_storage_bytes;
	void *cub_select_temp_storage;
	size_t cub_select_storage_bytes;
	void *cub_sort_temp_storage;
	size_t cub_sort_storage_bytes;

	
	host_batch_t *extensible_host_unpacked_query_batch;


	uint8_t *host_query_op;
	uint8_t *query_op;

	uint32_t *host_query_batch_offsets;

	uint32_t *host_query_batch_lens;


	gasal_res_t *host_res; // the results that can be read on host - THE STRUCT IS ON HOST SIDE, ITS CONTENT IS ON HOST SIDE.
	gasal_res_t *device_cpy; // a struct that contains the pointers to the device side - THE STRUCT IS ON HOST SIDE, but the CONTENT is malloc'd on and points to the DEVICE SIDE
	gasal_res_t *device_res; // the results that are written on device - THE STRUCT IS ON DEVICE SIDE, ITS CONTENT POINTS TO THE DEVICE SIDE.

	gasal_res_t *host_res_second; 
	gasal_res_t *device_res_second; 
	gasal_res_t *device_cpy_second;

	uint32_t max_query_batch_bytes;
	uint32_t max_n_queries;
	uint32_t max_n_seeds;
	uint32_t max_query_size;
	
	uint32_t actual_n_queries;

	cudaStream_t str;
	int is_free;

} gasal_gpu_storage_t;

//vector of streams
typedef struct {
	int n;
	gasal_gpu_storage_t *a;
}gasal_gpu_storage_v;


//match/mismatch and gap penalties
typedef struct{
	int32_t match;
	int32_t mismatch;
	int32_t gap_open;
	int32_t gap_extend;
} gasal_subst_scores;


#endif
