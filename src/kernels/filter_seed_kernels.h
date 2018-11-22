#ifndef KERNEL_BANDED
#define KERNEL_BANDED

__global__ void sum_arrays(uint32_t *in1, uint32_t *in2, uint32_t *out, int num_items) {
	int tid = (blockIdx.x*blockDim.x) + threadIdx.x;
	if (tid >= num_items) return;

	out[tid] = in1[tid] + in2[tid];

}


__global__ void filter_smem_intervals_fow_rev_kernel(uint2 *seed_intervals_fow_rev, int2 *seed_read_pos_fow_rev, uint32_t *n_smems_fow, uint32_t *n_smems_rev, uint32_t *n_smems_fow_rev_scan, uint32_t n_smems_max, int n_tasks) {

	int tid = (blockIdx.x*blockDim.x) + threadIdx.x;
	if (tid >= 2*n_tasks) return;
	if (tid < n_tasks) {
		int thread_read_num = tid/n_smems_max;
		int offset_in_read = tid - (thread_read_num*n_smems_max);
		if(offset_in_read >= n_smems_fow[thread_read_num] || offset_in_read == 0) return;
		int intv_idx = n_smems_fow_rev_scan[thread_read_num];
		int seed_begin_pos = seed_read_pos_fow_rev[intv_idx + offset_in_read].x;
		int comp_seed_begin_pos = seed_read_pos_fow_rev[intv_idx + offset_in_read - 1].x;
		if(seed_begin_pos == comp_seed_begin_pos ) {
			seed_intervals_fow_rev[intv_idx + offset_in_read - 1]= make_uint2 (1, 0);
			//seed_read_pos_fow[intv_idx + offset_in_read].y =  -1;
		}
	} else {
		tid = tid - n_tasks;
		int thread_read_num = tid/n_smems_max;
		int offset_in_read = tid - (thread_read_num*n_smems_max);
		if(offset_in_read >= n_smems_rev[thread_read_num] || offset_in_read == 0) return;
		int intv_idx = n_smems_fow_rev_scan[thread_read_num] + n_smems_fow[thread_read_num];
		int seed_begin_pos = seed_read_pos_fow_rev[intv_idx + offset_in_read].y;
		int comp_seed_begin_pos = seed_read_pos_fow_rev[intv_idx + offset_in_read - 1].y;
		if(seed_begin_pos == comp_seed_begin_pos) {
			seed_intervals_fow_rev[intv_idx + offset_in_read] =  make_uint2 (1, 0);

		}
	}

//	int intv_idx = read_offsets[thread_read_num] - (thread_read_num*min_seed_size);//thread_read_num * (MAX_READ_LENGTH - min_seed_size);
//	uint32_t thread_read_idx = read_idx[tid%n_tasks];
//	if (tid < n_tasks) {
//		int smems_num = n_smems_fow[thread_read_num];
//		if (thread_read_idx >= smems_num) return;
//		int i;
//		int2 seed_pos_on_read = seed_read_pos_fow[intv_idx + thread_read_idx];
//		for (i = 0; i < smems_num; i++) {
//			int2 comp_seed_pos_on_read = seed_read_pos_fow[intv_idx + i];
//			/*if ((seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//					|| (seed_pos_on_read.x == comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//					|| (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y))*/
//			if ((seed_pos_on_read.x == comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)){
//				pass = 0;intv_idx
//				//break;
//
//			}
//		}
//		if (pass) {
//			int smems_num_rev = n_smems_rev[thread_read_num];
//			int i;
//			for (i = 0; i < smems_num_rev; i++) {
//				int2 comp_seed_pos_on_read = seed_read_pos_rev[intv_idx + i];
//				if ((seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//								|| (seed_pos_on_read.x == comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//								|| (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y)){
//					pass = 0;
//					//break;
//
//				}
//			}
//		}
//		if (pass == 0) {
//			seed_intervals_fow[intv_idx + thread_read_idx] =  make_uint2 (1, 0) ;
//			seed_read_pos_fow[intv_idx + thread_read_idx] =  make_int2 (-1, -1) ;
//		}
//
//	} else {
//		int smems_num = n_smems_rev[thread_read_num];
//		if (thread_read_idx >= smems_num) return;
//		int i;
//		int2 seed_pos_on_read = seed_read_pos_rev[intv_idx + thread_read_idx];
//		for (i = 0; i < smems_num; i++) {
//			int2 comp_seed_pos_on_read = seed_read_pos_rev[intv_idx + i];
//			/*if ((seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//					|| (seed_pos_on_read.x == comp_seed_pos_on_read.x && s$
//					eed_pos_on_read.y < comp_seed_pos_on_read.y)
//					|| (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y))*/
//			if (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y){
//				pass = 0;
//				//break;
//
//			}
//		}
//		if (pass) {
//			int smems_num_fow = n_smems_fow[thread_read_num];
//			int i;
//			for (i = 0; i < smems_num_fow; i++) {
//				int2 comp_seed_pos_on_read = seed_read_pos_fow[intv_idx + i];
//				if ((seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//						|| (seed_pos_on_read.x == comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//						|| (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y)){
//					pass = 0;
//					//break;
//_compact_gpu
//				}
//			}
//		}
//		if (pass == 0) {
//			seed_intervals_rev[intv_idx + thread_read_idx] =  make_uint2 (1, 0) ;
//			seed_read_pos_rev[intv_idx + thread_read_idx] =  make_int2 (-1, -1) ;
//		}
//
//	}





	return;

}

__global__ void filter_smem_intervals_fow_kernel(uint2 *seed_intervals_fow_rev, int2 *seed_read_pos_fow_rev, uint32_t *n_smems_fow, uint32_t *n_smems_fow_rev_scan, uint32_t n_smems_max, int n_tasks) {

	int tid = (blockIdx.x*blockDim.x) + threadIdx.x;
	if (tid >= n_tasks) return;
	int thread_read_num = tid/n_smems_max;
	int offset_in_read = tid - (thread_read_num*n_smems_max);
	if(offset_in_read >= n_smems_fow[thread_read_num] || offset_in_read == 0) return;
	int intv_idx = n_smems_fow_rev_scan[thread_read_num];
	int seed_begin_pos = seed_read_pos_fow_rev[intv_idx + offset_in_read].x;
	int comp_seed_begin_pos = seed_read_pos_fow_rev[intv_idx + offset_in_read - 1].x;
	if(seed_begin_pos == comp_seed_begin_pos ) {
		seed_intervals_fow_rev[intv_idx + offset_in_read - 1]= make_uint2 (1, 0);
		//seed_read_pos_fow[intv_idx + offset_in_read].y =  -1;
	}

//	int intv_idx = read_offsets[thread_read_num] - (thread_read_num*min_seed_size);//thread_read_num * (MAX_READ_LENGTH - min_seed_size);
//	uint32_t thread_read_idx = read_idx[tid%n_tasks];
//	if (tid < n_tasks) {
//		int smems_num = n_smems_fow[thread_read_num];
//		if (thread_read_idx >= smems_num) return;
//		int i;
//		int2 seed_pos_on_read = seed_read_pos_fow[intv_idx + thread_read_idx];
//		for (i = 0; i < smems_num; i++) {
//			int2 comp_seed_pos_on_read = seed_read_pos_fow[intv_idx + i];
//			/*if ((seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//					|| (seed_pos_on_read.x == comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//					|| (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y))*/
//			if ((seed_pos_on_read.x == comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)){
//				pass = 0;intv_idx
//				//break;
//
//			}
//		}
//		if (pass) {
//			int smems_num_rev = n_smems_rev[thread_read_num];
//			int i;
//			for (i = 0; i < smems_num_rev; i++) {
//				int2 comp_seed_pos_on_read = seed_read_pos_rev[intv_idx + i];
//				if ((seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//								|| (seed_pos_on_read.x == comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//								|| (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y)){
//					pass = 0;
//					//break;
//
//				}
//			}
//		}
//		if (pass == 0) {
//			seed_intervals_fow[intv_idx + thread_read_idx] =  make_uint2 (1, 0) ;
//			seed_read_pos_fow[intv_idx + thread_read_idx] =  make_int2 (-1, -1) ;
//		}
//
//	} else {
//		int smems_num = n_smems_rev[thread_read_num];
//		if (thread_read_idx >= smems_num) return;
//		int i;
//		int2 seed_pos_on_read = seed_read_pos_rev[intv_idx + thread_read_idx];
//		for (i = 0; i < smems_num; i++) {
//			int2 comp_seed_pos_on_read = seed_read_pos_rev[intv_idx + i];
//			/*if ((seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//					|| (seed_pos_on_read.x == comp_seed_pos_on_read.x && s$
//					eed_pos_on_read.y < comp_seed_pos_on_read.y)
//					|| (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y))*/
//			if (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y){
//				pass = 0;
//				//break;
//
//			}
//		}
//		if (pass) {
//			int smems_num_fow = n_smems_fow[thread_read_num];
//			int i;
//			for (i = 0; i < smems_num_fow; i++) {
//				int2 comp_seed_pos_on_read = seed_read_pos_fow[intv_idx + i];
//				if ((seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//						|| (seed_pos_on_read.x == comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//						|| (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y)){
//					pass = 0;
//					//break;
//_compact_gpu
//				}
//			}
//		}
//		if (pass == 0) {
//			seed_intervals_rev[intv_idx + thread_read_idx] =  make_uint2 (1, 0) ;
//			seed_read_pos_rev[intv_idx + thread_read_idx] =  make_int2 (-1, -1) ;
//		}
//
//	}





	return;

}

__global__ void count_mem_occ_kernel(uint2 *seed_intervals_fow_rev, uint32_t *n_smems_fow_rev, uint32_t *n_seeds_fow_rev, uint32_t *n_smems_fow_rev_scan, uint32_t* n_ref_pos_fow_rev,  uint32_t n_smems_max, int n_tasks) {

	int tid = (blockIdx.x*blockDim.x) + threadIdx.x;
	if (tid >= n_tasks) return;

	int thread_read_num = tid/n_smems_max;
	int offset_in_read = tid - (thread_read_num*n_smems_max);
	if(offset_in_read >= n_smems_fow_rev[thread_read_num]) return;
	int intv_idx = n_smems_fow_rev_scan[thread_read_num];
	int n_intervals = seed_intervals_fow_rev[intv_idx + offset_in_read].y - seed_intervals_fow_rev[intv_idx + offset_in_read].x + 1;
	n_seeds_fow_rev[intv_idx + offset_in_read] = n_intervals;
	if (n_intervals > 0)  atomicAdd(&n_ref_pos_fow_rev[thread_read_num], n_intervals);

	return;

}
__global__ void filter_smem_intervals_kernel_wrapper(uint2 *seed_intervals_fow_rev_compact, int2 *seed_read_pos_fow_rev_compact, uint2 *seed_intervals_fow_rev,  int2 *seed_read_pos_fow_rev, uint32_t *n_smems_fow, uint32_t *n_smems_rev, uint32_t *n_smems_fow_rev, uint32_t* n_seeds_fow_rev, uint32_t *n_smems_fow_rev_scan,  uint32_t* n_ref_pos_fow_rev, uint32_t *n_smems_max, uint32_t *n_smems_sum_fow_rev, void *cub_sort_temp_storage, size_t cub_sort_storage_bytes, int total_reads, int n_bits_max_read_size, bool doRev) {

	uint32_t n_smems_max_val;
	if(doRev) n_smems_max_val = n_smems_max[0] > n_smems_max[1] ? n_smems_max[0] : n_smems_max[1];
	else n_smems_max_val = n_smems_max[0] > n_smems_max[1] ? n_smems_max[0] : n_smems_max[1];
	int n_tasks = n_smems_max_val*total_reads;
	n_smems_sum_fow_rev[0] = n_smems_sum_fow_rev[0]/2;
	int BLOCKDIM = 128;
	int N_BLOCKS = ((doRev ? 2 : 1)*n_tasks + BLOCKDIM - 1)/BLOCKDIM;


	if(doRev) filter_smem_intervals_fow_rev_kernel<<<N_BLOCKS, BLOCKDIM>>>(seed_intervals_fow_rev_compact, seed_read_pos_fow_rev_compact, n_smems_fow, n_smems_rev, n_smems_fow_rev_scan, n_smems_max_val, n_tasks);
	else filter_smem_intervals_fow_kernel<<<N_BLOCKS, BLOCKDIM>>>(seed_intervals_fow_rev_compact, seed_read_pos_fow_rev_compact, n_smems_fow, n_smems_fow_rev_scan, n_smems_max_val, n_tasks);
	cudaError_t filter_kernel_err = cudaGetLastError();
	if ( cudaSuccess != filter_kernel_err )
	{
			printf("[GPUseed CUDA ERROR:] Error in compute_intervals kernel: %s(CUDA error no.=%d). Line no. %d in file %s\n", cudaGetErrorString(filter_kernel_err), filter_kernel_err,  __LINE__, __FILE__);
			exit(EXIT_FAILURE);
	}

	CHECKCUBERROR(cub::DeviceSegmentedRadixSort::SortPairs(cub_sort_temp_storage, cub_sort_storage_bytes, (uint64_t*)seed_read_pos_fow_rev_compact, (uint64_t*)seed_read_pos_fow_rev, (uint64_t*)seed_intervals_fow_rev_compact, (uint64_t*)seed_intervals_fow_rev,  n_smems_sum_fow_rev[0], total_reads, n_smems_fow_rev_scan, n_smems_fow_rev_scan + 1, 0, n_bits_max_read_size), "DeviceSegmentedRadixSort::SortPairs");

	count_mem_occ_kernel<<<N_BLOCKS, BLOCKDIM>>>(seed_intervals_fow_rev, n_smems_fow_rev, n_seeds_fow_rev, n_smems_fow_rev_scan, n_ref_pos_fow_rev, n_smems_max_val << 1, 2*n_tasks);
	cudaError_t count_kernel_err = cudaGetLastError();
	if ( cudaSuccess != count_kernel_err )
	{
		printf("[GPUseed CUDA ERROR:] Error in compute_intervals kernel: %s(CUDA error no.=%d). Line no. %d in file %s\n", cudaGetErrorString(count_kernel_err), count_kernel_err,  __LINE__, __FILE__);
		exit(EXIT_FAILURE);
	}
}

__global__ void filter_mem_intervals_fow_rev_kernel(uint2 *seed_intervals_fow_rev_compact, int2 *seed_read_pos_fow_rev_compact, uint2 *seed_intervals_fow_rev, uint32_t *n_smems_fow, uint32_t *n_smems_rev, uint32_t *n_smems_fow_rev_scan, uint32_t n_smems_max, int n_tasks) {

	int tid = (blockIdx.x*blockDim.x) + threadIdx.x;
	if (tid >= 2*n_tasks) return;
	if (tid < n_tasks) {
		int thread_read_num = tid/n_smems_max;
		int offset_in_read = tid - (thread_read_num*n_smems_max);
		if(offset_in_read >= n_smems_fow[thread_read_num] || offset_in_read == 0) return;
		int intv_idx = n_smems_fow_rev_scan[thread_read_num];
		if(offset_in_read == n_smems_fow[thread_read_num] - 1){
			seed_intervals_fow_rev[intv_idx + offset_in_read] = seed_intervals_fow_rev_compact[intv_idx + offset_in_read];
		}
		int seed_begin_pos = seed_read_pos_fow_rev_compact[intv_idx + offset_in_read].x;
		int comp_seed_begin_pos = seed_read_pos_fow_rev_compact[intv_idx + offset_in_read - 1].x;
		if(seed_begin_pos == comp_seed_begin_pos && ((seed_intervals_fow_rev_compact[intv_idx + offset_in_read].y - seed_intervals_fow_rev_compact[intv_idx + offset_in_read].x + 1) == (seed_intervals_fow_rev_compact[intv_idx + offset_in_read - 1].y - seed_intervals_fow_rev_compact[intv_idx + offset_in_read - 1].x + 1))) {
			seed_intervals_fow_rev[intv_idx + offset_in_read - 1] =  make_uint2 (1, 0);			//seed_read_pos_fow[intv_idx + offset_in_read].y =  -1;
		} else {
			seed_intervals_fow_rev[intv_idx + offset_in_read - 1] = seed_intervals_fow_rev_compact[intv_idx + offset_in_read - 1];
		}
	} else {
		tid = tid - n_tasks;
		int thread_read_num = tid/n_smems_max;
		int offset_in_read = tid - (thread_read_num*n_smems_max);
		if(offset_in_read >= n_smems_rev[thread_read_num] /*|| offset_in_read == 0*/) return;
		int intv_idx = n_smems_fow_rev_scan[thread_read_num] + n_smems_fow[thread_read_num];
		if(offset_in_read == 0){
			seed_intervals_fow_rev[intv_idx + offset_in_read] = seed_intervals_fow_rev_compact[intv_idx + offset_in_read];
			return;
		}
		int seed_begin_pos = seed_read_pos_fow_rev_compact[intv_idx + offset_in_read].y;
		int comp_seed_begin_pos = seed_read_pos_fow_rev_compact[intv_idx + offset_in_read - 1].y;
		if(seed_begin_pos == comp_seed_begin_pos && ((seed_intervals_fow_rev_compact[intv_idx + offset_in_read].y - seed_intervals_fow_rev_compact[intv_idx + offset_in_read].x + 1) == (seed_intervals_fow_rev_compact[intv_idx + offset_in_read - 1].y - seed_intervals_fow_rev_compact[intv_idx + offset_in_read - 1].x + 1))) {
			seed_intervals_fow_rev[intv_idx + offset_in_read] =  make_uint2 (1, 0);			//seed_read_pos_fow[intv_idx + offset_in_read].y =  -1;
		} else {
			seed_intervals_fow_rev[intv_idx + offset_in_read] = seed_intervals_fow_rev_compact[intv_idx + offset_in_read];
		}
	}

//	int intv_idx = read_offsets[thread_read_num] - (thread_read_num*min_seed_size);//thread_read_num * (MAX_READ_LENGTH - min_seed_size);
//	uint32_t thread_read_idx = read_idx[tid%n_tasks];
//	if (tid < n_tasks) {
//		int smems_num = n_smems_fow[thread_read_num];
//		if (thread_read_idx >= smems_num) return;
//		int i;
//		int2 seed_pos_on_read = seed_read_pos_fow[intv_idx + thread_read_idx];
//		for (i = 0; i < smems_num; i++) {
//			int2 comp_seed_pos_on_read = seed_read_pos_fow[intv_idx + i];
//			/*if ((seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//					|| (seed_pos_on_read.x == comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//					|| (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y))*/
//			if ((seed_pos_on_read.x == comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)){
//				pass = 0;intv_idx
//				//break;
//
//			}
//		}
//		if (pass) {
//			int smems_num_rev = n_smems_rev[thread_read_num];
//			int i;
//			for (i = 0; i < smems_num_rev; i++) {
//				int2 comp_seed_pos_on_read = seed_read_pos_rev[intv_idx + i];
//				if ((seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//								|| (seed_pos_on_read.x == comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//								|| (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y)){
//					pass = 0;
//					//break;
//
//				}
//			}
//		}
//		if (pass == 0) {
//			seed_intervals_fow[intv_idx + thread_read_idx] =  make_uint2 (1, 0) ;
//			seed_read_pos_fow[intv_idx + thread_read_idx] =  make_int2 (-1, -1) ;
//		}
//
//	} else {
//		int smems_num = n_smems_rev[thread_read_num];
//		if (thread_read_idx >= smems_num) return;
//		int i;
//		int2 seed_pos_on_read = seed_read_pos_rev[intv_idx + thread_read_idx];
//		for (i = 0; i < smems_num; i++) {
//			int2 comp_seed_pos_on_read = seed_read_pos_rev[intv_idx + i];
//			/*if ((seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//					|| (seed_pos_on_read.x == comp_seed_pos_on_read.x && s$
//					eed_pos_on_read.y < comp_seed_pos_on_read.y)
//					|| (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y))*/
//			if (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y){
//				pass = 0;
//				//break;
//
//			}
//		}
//		if (pass) {
//			int smems_num_fow = n_smems_fow[thread_read_num];
//			int i;
//			for (i = 0; i < smems_num_fow; i++) {
//				int2 comp_seed_pos_on_read = seed_read_pos_fow[intv_idx + i];
//				if ((seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//						|| (seed_pos_on_read.x == comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//						|| (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y)){
//					pass = 0;
//					//break;
//_compact_gpu
//				}
//			}
//		}
//		if (pass == 0) {
//			seed_intervals_rev[intv_idx + thread_read_idx] =  make_uint2 (1, 0) ;
//			seed_read_pos_rev[intv_idx + thread_read_idx] =  make_int2 (-1, -1) ;
//		}
//
//	}





	return;

}

__global__ void filter_mem_intervals_fow_kernel(uint2 *seed_intervals_fow_rev_compact, int2 *seed_read_pos_fow_rev_compact, uint2 *seed_intervals_fow_rev, uint32_t *n_smems_fow, uint32_t *n_smems_fow_rev_scan, uint32_t n_smems_max, int n_tasks) {

	int tid = (blockIdx.x*blockDim.x) + threadIdx.x;
	if (tid >= n_tasks) return;
	int thread_read_num = tid/n_smems_max;
	int offset_in_read = tid - (thread_read_num*n_smems_max);
	if(offset_in_read >= n_smems_fow[thread_read_num] || offset_in_read == 0) return;
	int intv_idx = n_smems_fow_rev_scan[thread_read_num];
	if(offset_in_read == n_smems_fow[thread_read_num] - 1){
		seed_intervals_fow_rev[intv_idx + offset_in_read] = seed_intervals_fow_rev_compact[intv_idx + offset_in_read];
	}
	int seed_begin_pos = seed_read_pos_fow_rev_compact[intv_idx + offset_in_read].x;
	int comp_seed_begin_pos = seed_read_pos_fow_rev_compact[intv_idx + offset_in_read - 1].x;
	if(seed_begin_pos == comp_seed_begin_pos && ((seed_intervals_fow_rev_compact[intv_idx + offset_in_read].y - seed_intervals_fow_rev_compact[intv_idx + offset_in_read].x + 1) == (seed_intervals_fow_rev_compact[intv_idx + offset_in_read - 1].y - seed_intervals_fow_rev_compact[intv_idx + offset_in_read - 1].x + 1))) {
		seed_intervals_fow_rev[intv_idx + offset_in_read - 1] =  make_uint2 (1, 0);			//seed_read_pos_fow[intv_idx + offset_in_read].y =  -1;
	} else {
		seed_intervals_fow_rev[intv_idx + offset_in_read - 1] = seed_intervals_fow_rev_compact[intv_idx + offset_in_read - 1];
	}

//	int intv_idx = read_offsets[thread_read_num] - (thread_read_num*min_seed_size);//thread_read_num * (MAX_READ_LENGTH - min_seed_size);
//	uint32_t thread_read_idx = read_idx[tid%n_tasks];
//	if (tid < n_tasks) {
//		int smems_num = n_smems_fow[thread_read_num];
//		if (thread_read_idx >= smems_num) return;
//		int i;
//		int2 seed_pos_on_read = seed_read_pos_fow[intv_idx + thread_read_idx];
//		for (i = 0; i < smems_num; i++) {
//			int2 comp_seed_pos_on_read = seed_read_pos_fow[intv_idx + i];
//			/*if ((seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//					|| (seed_pos_on_read.x == comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//					|| (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y))*/
//			if ((seed_pos_on_read.x == comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)){
//				pass = 0;intv_idx
//				//break;
//
//			}
//		}
//		if (pass) {
//			int smems_num_rev = n_smems_rev[thread_read_num];
//			int i;
//			for (i = 0; i < smems_num_rev; i++) {
//				int2 comp_seed_pos_on_read = seed_read_pos_rev[intv_idx + i];
//				if ((seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//								|| (seed_pos_on_read.x == comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//								|| (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y)){
//					pass = 0;
//					//break;
//
//				}
//			}
//		}
//		if (pass == 0) {
//			seed_intervals_fow[intv_idx + thread_read_idx] =  make_uint2 (1, 0) ;
//			seed_read_pos_fow[intv_idx + thread_read_idx] =  make_int2 (-1, -1) ;
//		}
//
//	} else {
//		int smems_num = n_smems_rev[thread_read_num];
//		if (thread_read_idx >= smems_num) return;
//		int i;
//		int2 seed_pos_on_read = seed_read_pos_rev[intv_idx + thread_read_idx];
//		for (i = 0; i < smems_num; i++) {
//			int2 comp_seed_pos_on_read = seed_read_pos_rev[intv_idx + i];
//			/*if ((seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//					|| (seed_pos_on_read.x == comp_seed_pos_on_read.x && s$
//					eed_pos_on_read.y < comp_seed_pos_on_read.y)
//					|| (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y))*/
//			if (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y){
//				pass = 0;
//				//break;
//
//			}
//		}
//		if (pass) {
//			int smems_num_fow = n_smems_fow[thread_read_num];
//			int i;
//			for (i = 0; i < smems_num_fow; i++) {
//				int2 comp_seed_pos_on_read = seed_read_pos_fow[intv_idx + i];
//				if ((seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//						|| (seed_pos_on_read.x == comp_seed_pos_on_read.x && seed_pos_on_read.y < comp_seed_pos_on_read.y)
//						|| (seed_pos_on_read.x > comp_seed_pos_on_read.x && seed_pos_on_read.y == comp_seed_pos_on_read.y)){
//					pass = 0;
//					//break;
//_compact_gpu
//				}
//			}
//		}
//		if (pass == 0) {
//			seed_intervals_rev[intv_idx + thread_read_idx] =  make_uint2 (1, 0) ;
//			seed_read_pos_rev[intv_idx + thread_read_idx] =  make_int2 (-1, -1) ;
//		}
//
//	}





	return;

}

__global__ void count_mem_occ_fow_rev_kernel(uint2 *seed_intervals_fow_rev, int2 *seed_read_pos_fow_rev, uint32_t *n_smems_fow,  uint32_t *n_smems_rev, uint32_t *n_seeds_fow_rev, uint32_t *n_smems_fow_rev_scan,uint32_t *n_ref_pos_fow_rev, uint32_t n_smems_max, int n_tasks) {

	 int tid = (blockIdx.x*blockDim.x) + threadIdx.x;
	 if (tid >= 2*n_tasks) return;

	 if (tid < n_tasks) {
		 int thread_read_num = tid/n_smems_max;
		 int offset_in_read = tid - (thread_read_num*n_smems_max);
		 if(offset_in_read >= n_smems_fow[thread_read_num]) return;
		 int intv_idx = n_smems_fow_rev_scan[thread_read_num];
		 int n_intervals = seed_intervals_fow_rev[intv_idx + offset_in_read].y - seed_intervals_fow_rev[intv_idx + offset_in_read].x + 1;
		 int n_intervals_to_add = n_intervals;
		 int next_n_intervals = 0;
		 if (n_intervals > 0) {
			 int seed_read_pos_x = seed_read_pos_fow_rev[intv_idx + offset_in_read].x;
			 int p = 1;
			 while (seed_read_pos_x == seed_read_pos_fow_rev[intv_idx + offset_in_read + p].x && offset_in_read + p < n_smems_fow[thread_read_num]) {
				 next_n_intervals = seed_intervals_fow_rev[intv_idx + offset_in_read + p].y - seed_intervals_fow_rev[intv_idx + offset_in_read + p].x + 1;
				 if (next_n_intervals > 0) {
					 n_intervals_to_add = n_intervals - next_n_intervals;
					 break;
				 }

				 p++;
			 }
			 atomicAdd(&n_ref_pos_fow_rev[thread_read_num], n_intervals_to_add);
		 }

		 n_seeds_fow_rev[intv_idx + offset_in_read] = n_intervals_to_add;



	 } else {
		 tid = tid - n_tasks;
		 int thread_read_num = tid/n_smems_max;
		 int offset_in_read = tid - (thread_read_num*n_smems_max);
		 if(offset_in_read >= n_smems_rev[thread_read_num]) return;
		 int intv_idx = n_smems_fow_rev_scan[thread_read_num] + n_smems_fow[thread_read_num];
		 int n_intervals = seed_intervals_fow_rev[intv_idx + offset_in_read].y - seed_intervals_fow_rev[intv_idx + offset_in_read].x + 1;
		 int n_intervals_to_add = n_intervals;
		 int next_n_intervals = 0;
		 if (n_intervals > 0) {
			 int seed_read_pos_y = seed_read_pos_fow_rev[intv_idx + offset_in_read].y;
			 int p = 1;
			 while (seed_read_pos_y == seed_read_pos_fow_rev[intv_idx + offset_in_read - p].y && offset_in_read - p >= 0) {
				 next_n_intervals = seed_intervals_fow_rev[intv_idx + offset_in_read - p].y - seed_intervals_fow_rev[intv_idx + offset_in_read - p].x + 1;
				 if (next_n_intervals > 0) {
					 n_intervals_to_add = n_intervals - next_n_intervals;
					 break;
				 }

				 p++;
			 }
			 atomicAdd(&n_ref_pos_fow_rev[thread_read_num], n_intervals_to_add);
		 }

		 n_seeds_fow_rev[intv_idx + offset_in_read] = n_intervals_to_add;

	 }

	return;

}


__global__ void count_mem_occ_fow_kernel(uint2 *seed_intervals_fow_rev, int2 *seed_read_pos_fow_rev, uint32_t *n_smems_fow, uint32_t *n_seeds_fow_rev, uint32_t *n_smems_fow_rev_scan,uint32_t *n_ref_pos_fow_rev, uint32_t n_smems_max, int n_tasks) {

	 int tid = (blockIdx.x*blockDim.x) + threadIdx.x;
	 if (tid >= n_tasks) return;

	 int thread_read_num = tid/n_smems_max;
	 int offset_in_read = tid - (thread_read_num*n_smems_max);
	 if(offset_in_read >= n_smems_fow[thread_read_num]) return;
	 int intv_idx = n_smems_fow_rev_scan[thread_read_num];
	 int n_intervals = seed_intervals_fow_rev[intv_idx + offset_in_read].y - seed_intervals_fow_rev[intv_idx + offset_in_read].x + 1;
	 int n_intervals_to_add = n_intervals;
	 int next_n_intervals = 0;
	 if (n_intervals > 0) {
		 int seed_read_pos_x = seed_read_pos_fow_rev[intv_idx + offset_in_read].x;
		 int p = 1;
		 while (seed_read_pos_x == seed_read_pos_fow_rev[intv_idx + offset_in_read + p].x && offset_in_read + p < n_smems_fow[thread_read_num]) {
			 next_n_intervals = seed_intervals_fow_rev[intv_idx + offset_in_read + p].y - seed_intervals_fow_rev[intv_idx + offset_in_read + p].x + 1;
			 if (next_n_intervals > 0) {
				 n_intervals_to_add = n_intervals - next_n_intervals;
				 break;
			 }

			 p++;
		 }
		 atomicAdd(&n_ref_pos_fow_rev[thread_read_num], n_intervals_to_add);
	 }

	 n_seeds_fow_rev[intv_idx + offset_in_read] = n_intervals_to_add;





	 return;

}
__global__ void filter_mem_intervals_kernel_wrapper(uint2 *seed_intervals_fow_rev_compact, int2 *seed_read_pos_fow_rev_compact, uint2 *seed_intervals_fow_rev, uint32_t *n_smems_fow, uint32_t *n_smems_rev, uint32_t *n_seeds_fow_rev, uint32_t *n_smems_fow_rev_scan, uint32_t *n_ref_pos_fow_rev, uint32_t *n_smems_max, uint32_t *n_smems_sum_fow_rev, int total_reads, bool doRev) {

	uint32_t n_smems_max_val;
	if(doRev) n_smems_max_val = n_smems_max[0] > n_smems_max[1] ? n_smems_max[0] : n_smems_max[1];
	else n_smems_max_val = n_smems_max[0] > n_smems_max[1] ? n_smems_max[0] : n_smems_max[1];
	int n_tasks = n_smems_max_val*total_reads;
	n_smems_sum_fow_rev[0] = n_smems_sum_fow_rev[0]/2;
	int BLOCKDIM = 128;
	int N_BLOCKS = ((doRev ? 2 : 1)*n_tasks + BLOCKDIM - 1)/BLOCKDIM;

	//filter_seed_intervals_gpu<<<N_BLOCKS, BLOCKDIM>>>(seed_intervals_fow_rev_compact, seed_read_pos_fow_rev_compact, n_smems_fow, n_smems_rev, n_smems_fow_rev_scan, n_smems_max_val, n_tasks);

	if (doRev) {
		filter_mem_intervals_fow_rev_kernel<<<N_BLOCKS, BLOCKDIM>>>(seed_intervals_fow_rev_compact, seed_read_pos_fow_rev_compact, seed_intervals_fow_rev, n_smems_fow, n_smems_rev, n_smems_fow_rev_scan, n_smems_max_val, n_tasks);
		count_mem_occ_fow_rev_kernel<<<N_BLOCKS, BLOCKDIM>>>(seed_intervals_fow_rev, seed_read_pos_fow_rev_compact, n_smems_fow, n_smems_rev, n_seeds_fow_rev, n_smems_fow_rev_scan, n_ref_pos_fow_rev, n_smems_max_val, n_tasks);
	} else{
		filter_mem_intervals_fow_kernel<<<N_BLOCKS, BLOCKDIM>>>(seed_intervals_fow_rev_compact, seed_read_pos_fow_rev_compact, seed_intervals_fow_rev, n_smems_fow, n_smems_fow_rev_scan, n_smems_max_val, n_tasks);
		count_mem_occ_fow_kernel<<<N_BLOCKS, BLOCKDIM>>>(seed_intervals_fow_rev, seed_read_pos_fow_rev_compact, n_smems_fow, n_seeds_fow_rev, n_smems_fow_rev_scan, n_ref_pos_fow_rev, n_smems_max_val, n_tasks);
	}

}
#endif
