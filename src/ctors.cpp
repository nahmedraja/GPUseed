
#include "gasal.h"

#include "args_parser.h"

#include "host_batch.h"

#include "res.h"

#include "ctors.h"

#include "cub/cub/cub.cuh"


gasal_gpu_storage_v GPUseed_init_gpu_storage_v(int n_streams) {
	gasal_gpu_storage_v v;
	v.a = (gasal_gpu_storage_t*)calloc(n_streams, sizeof(gasal_gpu_storage_t));
	v.n = n_streams;
	return v;

}


void GPUseed_init_streams(gasal_gpu_storage_v *gpu_storage_vec, int max_query_batch_bytes, int max_n_queries, uint32_t max_n_seeds, int max_query_size, Parameters *params) {

	cudaError_t err;
	int i;
	max_query_batch_bytes = max_query_batch_bytes%8 == 0 ? max_query_batch_bytes: max_query_batch_bytes + (8 - max_query_batch_bytes%8);
	if (max_n_seeds > 4e9) {
		fprintf(stderr, "[GASAL ERROR:] max_n_seeds should <= 4e9\n");
		exit(EXIT_FAILURE);
	}
	max_n_seeds = max_n_seeds >= (2*(params->doRev ? 2 : 1) * max_query_batch_bytes) ? max_n_seeds : (2*(params->doRev ? 2 : 1) * max_query_batch_bytes);
	for (i = 0; i < gpu_storage_vec->n; i++) {

		gpu_storage_vec->a[i].extensible_host_unpacked_query_batch = gasal_host_batch_new(max_query_batch_bytes, 0);

		CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].unpacked_query_batch), max_query_batch_bytes * sizeof(uint8_t)), "cudaMalloc");



		int max_packed_query_batch_size = max_query_batch_bytes%8 == 0 ? max_query_batch_bytes/8 : max_query_batch_bytes/8 + (8 - max_query_batch_bytes/8);

		CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].packed_fow_query_batch), max_query_batch_bytes/8 * sizeof(uint32_t)), "cudaMalloc");
		if (params->doRev) CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].packed_rev_query_batch), max_query_batch_bytes/8 * sizeof(uint32_t)), "cudaMalloc");


		CHECKCUDAERROR(cudaHostAlloc(&(gpu_storage_vec->a[i].host_query_batch_lens), max_n_queries * sizeof(uint32_t), cudaHostAllocDefault), "cudaHostAlloc");
		CHECKCUDAERROR(cudaHostAlloc(&(gpu_storage_vec->a[i].host_query_batch_offsets), max_n_queries * sizeof(uint32_t), cudaHostAllocDefault), "cudaHostAlloc");

		CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].query_batch_lens), max_n_queries * sizeof(uint32_t)), "cudaMalloc");

		CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].query_batch_offsets), max_n_queries * sizeof(uint32_t)), "cudaMalloc");

		CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].n_smems_fow), max_n_queries*sizeof(uint32_t)), "cudaMalloc");
		if (params->doRev) CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].n_smems_rev), max_n_queries*sizeof(uint32_t)), "cudaMalloc");

		CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].n_smems_scan),(max_n_queries+1)*sizeof(uint32_t)), "cudaMalloc");
		CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].n_seed_occ_scan), ((params->doRev ? 2 : 1)*(max_query_batch_bytes) + 1)*sizeof(uint32_t)), "cudaMalloc");

		CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].is_smem_flag), max_n_seeds*sizeof(uint32_t)), "cudaMalloc");

		CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].n_seed_occ_sum), sizeof(uint64_t)), "cudaMalloc");
		CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].n_smems_max), (params->doRev ? 2 : 1)*sizeof(uint32_t)), "cudaMalloc");
		CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].n_queries_processed), sizeof(uint32_t)), "cudaMalloc");

		CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].seed_query_pos), max_n_seeds*sizeof(int2)), "cudaMalloc");

		CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].seed_query_pos_compact), max_n_seeds*sizeof(int2)), "cudaMalloc");

		if(params->doMem){
			CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].seed_ref_pos), max_n_seeds*sizeof(uint32_t)), "cudaMalloc");
			CHECKCUDAERROR(cudaHostAlloc(&(gpu_storage_vec->a[i].host_seed_ref_pos), max_n_seeds*sizeof(uint32_t), cudaHostAllocDefault), "cudaHostAlloc");
		}



		gpu_storage_vec->a[i]. cub_select_temp_storage = NULL;
		gpu_storage_vec->a[i]. cub_select_storage_bytes = 0;
		CHECKCUBERROR(cub::DeviceSelect::Flagged(gpu_storage_vec->a[i].cub_select_temp_storage, gpu_storage_vec->a[i].cub_select_storage_bytes, (uint32_t*)(gpu_storage_vec->a[i].seed_intervals), (uint16_t*)(gpu_storage_vec->a[i].is_smem_flag), (uint32_t*)(gpu_storage_vec->a[i].seed_intervals_compact), gpu_storage_vec->a[i].n_smems_sum, 2*(params->doRev ? 2 : 1)*(max_query_batch_bytes)), "DeviceSelect::Flagged");
		CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].cub_select_temp_storage),gpu_storage_vec->a[i].cub_select_storage_bytes), "cudaMalloc");


		gpu_storage_vec->a[i]. cub_sort_temp_storage = NULL;
		gpu_storage_vec->a[i]. cub_sort_storage_bytes = 0;
		int n_bits_max_query_size = (int)ceil(log2((double)max_query_size));
		if (params->doMem) {
			CHECKCUBERROR(cub::DeviceSegmentedRadixSort::SortPairs(gpu_storage_vec->a[i].cub_sort_temp_storage, gpu_storage_vec->a[i].cub_sort_storage_bytes, (uint64_t*)(gpu_storage_vec->a[i].seed_query_pos_compact), (uint64_t*)(gpu_storage_vec->a[i].seed_query_pos), gpu_storage_vec->a[i].seed_sa_idx, gpu_storage_vec->a[i].seed_ref_pos,  max_n_seeds, max_n_queries, gpu_storage_vec->a[i].n_query_seeds_scan, gpu_storage_vec->a[i].n_query_seeds_scan + 1, 0, n_bits_max_query_size), "DeviceSegmentedRadixSort::SortPairs");

		} else {
			CHECKCUBERROR(cub::DeviceSegmentedRadixSort::SortPairs(gpu_storage_vec->a[i].cub_sort_temp_storage, gpu_storage_vec->a[i].cub_sort_storage_bytes, (uint64_t*)(gpu_storage_vec->a[i].seed_query_pos_compact), (uint64_t*)(gpu_storage_vec->a[i].seed_query_pos), (uint64_t*)(gpu_storage_vec->a[i].seed_intervals_compact), (uint64_t*)(gpu_storage_vec->a[i].seed_intervals),   (params->doRev ? 2 : 1)*(max_query_batch_bytes), max_n_queries, gpu_storage_vec->a[i].n_smems_scan, gpu_storage_vec->a[i].n_smems_scan + 1, 0, n_bits_max_query_size), "DeviceSegmentedRadixSort::SortPairs");
		}
		CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].cub_sort_temp_storage),gpu_storage_vec->a[i].cub_sort_storage_bytes), "cudaMalloc");

		gpu_storage_vec->a[i]. cub_scan_temp_storage = NULL;
		gpu_storage_vec->a[i]. cub_scan_storage_bytes = 0;
		CHECKCUBERROR(cub::DeviceScan::ExclusiveSum(gpu_storage_vec->a[i]. cub_scan_temp_storage, gpu_storage_vec->a[i]. cub_scan_storage_bytes, gpu_storage_vec->a[i].n_seed_occ, gpu_storage_vec->a[i].n_seed_occ_scan, (params->doRev ? 2 : 1)*(max_query_batch_bytes)), "DeviceScan::ExclusiveSum");
		CHECKCUDAERROR(cudaMalloc(&(gpu_storage_vec->a[i].cub_scan_temp_storage),gpu_storage_vec->a[i].cub_scan_storage_bytes), "cudaMalloc");
		




		CHECKCUDAERROR(cudaStreamCreate(&(gpu_storage_vec->a[i].str)), "cudaStreamCreate");
		gpu_storage_vec->a[i].is_free = 1;
		gpu_storage_vec->a[i].max_query_batch_bytes = max_query_batch_bytes;
		gpu_storage_vec->a[i].max_n_queries = max_n_queries;
		gpu_storage_vec->a[i].max_query_batch_bytes = max_query_batch_bytes;
		gpu_storage_vec->a[i].max_n_queries = max_n_queries;
		gpu_storage_vec->a[i].max_n_seeds = max_n_seeds;
		gpu_storage_vec->a[i].max_query_size = max_n_seeds;
	}
}




void GPUseed_destroy_streams(gasal_gpu_storage_v *gpu_storage_vec, Parameters *params) {

	cudaError_t err;

	int i;
	for (i = 0; i < gpu_storage_vec->n; i ++) {
		
		gasal_host_batch_destroy(gpu_storage_vec->a[i].extensible_host_unpacked_query_batch);




		if (gpu_storage_vec->a[i].host_query_batch_offsets != NULL) CHECKCUDAERROR(cudaFreeHost(gpu_storage_vec->a[i].host_query_batch_offsets), "cudaFreeHost");
		if (gpu_storage_vec->a[i].host_query_batch_lens != NULL) CHECKCUDAERROR(cudaFreeHost(gpu_storage_vec->a[i].host_query_batch_lens), "cudaFreeHost");





		if (gpu_storage_vec->a[i].unpacked_query_batch != NULL) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].unpacked_query_batch), "cudaFree");

		if (gpu_storage_vec->a[i].packed_fow_query_batch != NULL) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].packed_fow_query_batch), "cudaFree");
		if (gpu_storage_vec->a[i].packed_rev_query_batch != NULL && params->doRev) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].packed_rev_query_batch), "cudaFree");

		if (gpu_storage_vec->a[i].n_smems_fow != NULL) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].n_smems_fow), "cudaFree");
		if (gpu_storage_vec->a[i].n_smems_rev  != NULL && params->doRev) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].n_smems_rev), "cudaFree");
		if (gpu_storage_vec->a[i].n_smems_scan != NULL) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].n_smems_scan), "cudaFree");
		if (gpu_storage_vec->a[i].n_seed_occ_scan != NULL) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].n_seed_occ_scan), "cudaFree");
		if (gpu_storage_vec->a[i].is_smem_flag !=NULL)	CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].is_smem_flag), "cudaFree");
		if (gpu_storage_vec->a[i].n_seed_occ_sum != NULL) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].n_seed_occ_sum), "cudaFree");
		if (gpu_storage_vec->a[i].n_smems_max != NULL) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].n_smems_max), "cudaFree");
		if (gpu_storage_vec->a[i].n_queries_processed != NULL) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].n_queries_processed), "cudaFree");
		if (gpu_storage_vec->a[i].seed_query_pos != NULL) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].seed_query_pos), "cudaFree");
		if (gpu_storage_vec->a[i].seed_query_pos_compact != NULL) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].seed_query_pos_compact), "cudaFree");
		if (gpu_storage_vec->a[i].seed_ref_pos != NULL && params->doMem) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].seed_ref_pos), "cudaFree");
		if (gpu_storage_vec->a[i].host_seed_ref_pos != NULL && params->doMem) CHECKCUDAERROR(cudaFreeHost(gpu_storage_vec->a[i].host_seed_ref_pos), "cudaFreeHost");
		if (gpu_storage_vec->a[i].cub_select_temp_storage != NULL) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].cub_select_temp_storage), "cudaFree");
		if (gpu_storage_vec->a[i].cub_sort_temp_storage != NULL) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].cub_sort_temp_storage), "cudaFree");
		if (gpu_storage_vec->a[i].cub_scan_temp_storage != NULL) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].cub_scan_temp_storage), "cudaFree");


		if (gpu_storage_vec->a[i].query_batch_offsets != NULL) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].query_batch_offsets), "cudaFree");

		if (gpu_storage_vec->a[i].query_batch_lens != NULL) CHECKCUDAERROR(cudaFree(gpu_storage_vec->a[i].query_batch_lens), "cudaFree");


		if (gpu_storage_vec->a[i].str != NULL)CHECKCUDAERROR(cudaStreamDestroy(gpu_storage_vec->a[i].str), "cudaStreamDestroy");
	}



}


void GPUseed_destroy_gpu_storage_v(gasal_gpu_storage_v *gpu_storage_vec) {

	if(gpu_storage_vec->a != NULL) free(gpu_storage_vec->a);
}

