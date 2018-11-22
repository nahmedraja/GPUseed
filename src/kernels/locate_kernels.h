#ifndef __LOCAL_KERNEL_TEMPLATE__
#define __LOCAL_KERNEL_TEMPLATE__


__global__ void find_n_query_to_process (uint32_t *n_query_seeds, uint32_t *n_smems_scan, uint32_t *n_seed_occ, uint32_t *n_seed_occ_sum, uint32_t *n_smems_sum, uint32_t *n_queries_processed, void *cub_scan_temp_storage, size_t cub_scan_storage_bytes, uint32_t actual_n_queries, uint32_t max_n_seeds) {


	CHECKCUBERROR(cub::DeviceReduce::Sum(cub_scan_temp_storage, cub_scan_storage_bytes, n_seed_occ , n_seed_occ_sum,  n_smems_sum[0]), "DeviceReduce::Sum");
	if (n_seed_occ_sum[0] <= max_n_seeds) {
		n_queries_processed[0] = actual_n_queries;
		return;
	}
	else {

	}
}

#endif
