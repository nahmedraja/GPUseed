#ifndef __KERNEL_GLOBAL__
#define __KERNEL_GLOBAL__


#define N_SHUFFLES 30
__global__ void compute_intervals_fow_rev_kernel(uint32_t *packed_read_batch_fow,  uint32_t *packed_read_batch_rev, uint32_t *read_sizes, uint32_t *read_offsets, uint2 *seed_intervals_fow_rev,
		int2 *seed_read_pos_fow_rev, uint32_t *read_num, uint32_t *read_idx, uint32_t *is_smem_fow_rev_flag, uint32_t *n_smems_fow,  uint32_t *n_smems_rev, int min_seed_size, GPUseed_ref_index ref_index, int n_tasks) {

	int tid = (blockIdx.x*blockDim.x) + threadIdx.x;
	if (tid >= 2*n_tasks) return;
	int thread_read_num = read_num[tid%n_tasks];
	int read_len = read_sizes[thread_read_num];
	int read_off = read_offsets[thread_read_num];
	//thread_read_num * (MAX_READ_LENGTH - min_seed_size);
	uint32_t thread_read_idx = read_idx[tid%n_tasks];
	int is_active = 0;
	int is_smem = 1;
	int is_shfl[N_SHUFFLES];
	int only_next_time = 0;
	uint32_t neighbour_active[N_SHUFFLES];
	uint32_t prev_intv_size[N_SHUFFLES];

	int m;
	for (m = 0; m < N_SHUFFLES; m++) {
		is_shfl[m] = tid%n_tasks ? 1 : 0;
		if (is_shfl[m]) is_shfl[m] = (thread_read_num == read_num[(tid%n_tasks) - (m+1)]) ? 1 : 0;
		if (is_shfl[m]) is_shfl[m] = ((tid%32) - m > 0) ? 1 : 0;
		prev_intv_size[m] = 0;
		neighbour_active[m] = 1;
	}

//	int is_shfl = tid%n_tasks ? 1 : 0;
//	if (is_shfl) is_shfl = (thread_read_num == read_num[(tid%n_tasks) - 1]) ? 1 : 0;
//	if (is_shfl) is_shfl = (tid%32) ? 1 : 0;
//	int is_shfl_2 = (thread_read_num == read_num[(tid%n_tasks) - 2]) ? 1 : 0;
//	if (is_shfl_2) is_shfl_2 = (tid%33) ? 1 : 0;
//	int is_shfl_3 = (thread_read_num == read_num[(tid%n_tasks) - 3]) ? 1 : 0;
//	if (is_shfl_3) is_shfl_3 = (tid%34) ? 1 : 0;
//	int is_shfl_4 = (thread_read_num == read_num[(tid%n_tasks) - 4]) ? 1 : 0;
//	if (is_shfl_4) is_shfl_4 = (tid%35) ? 1 : 0;
//	int is_shfl_5 = (thread_read_num == read_num[(tid%n_tasks) - 5]) ? 1 : 0;
//	if (is_shfl_5) is_shfl_5 = (tid%36) ? 1 : 0;
	int i, j;
	int base;
	bwtint_t l, u;
	if (tid < n_tasks) {
		int intv_idx = (2*(read_offsets[thread_read_num] - (thread_read_num*min_seed_size))) + read_len - min_seed_size - 1;
		int start = read_off&7;
		uint32_t *seq = &(packed_read_batch_fow[read_off >> 3]);
		uint32_t pre_calc_seed = 0;
		for (i = start + read_len - thread_read_idx - 1, j = 0; j < pre_calc_seed_len; i--, j++) {
			int reg_no = i >> 3;
			int reg_pos = i & 7;
			int reg = seq[reg_no];
			uint32_t base = (reg >> (28 - (reg_pos << 2)))&15;
			/*unknown bases*/
			if (base > 3) {
				break;
			}
			pre_calc_seed |= (base << (j<<1));
		}
		uint2 prev_seed_interval = j < pre_calc_seed_len ? make_uint2(1,0) : pre_calc_intervals[pre_calc_seed];
		int beg_i = i;
		if(prev_seed_interval.x <= prev_seed_interval.y) {
			is_active = 1;
			uint32_t curr_intv_size = prev_seed_interval.y - prev_seed_interval.x + 1;
//			uint32_t prev_intv_size = 0;
//			uint32_t prev_2_intv_size = 0;
//			uint32_t prev_3_intv_size = 0;
//			uint32_t prev_4_intv_size = 0;
//			uint32_t prev_5_intv_size = 0;
			l = prev_seed_interval.x, u = prev_seed_interval.y;
			for (; i >= start; i--) {
				/*get the base*/
				if (is_active) {
					prev_seed_interval = make_uint2(l,u);
					int reg_no = i >> 3;
					int reg_pos = i & 7;
					int reg = seq[reg_no];
					int base = (reg >> (28 - (reg_pos << 2)))&15;
					/*unknown bases*/
//					if (base > 3) {
//						is_active = 0;
//						break;
//					}

					uint2 intv = find_occ_gpu(bwt, prev_seed_interval.x - 1, prev_seed_interval.y, base);
					//calculate the range
					//l = L2_gpu[ch] + bwt_occ_gpu(bwt, prev_l - 1, ch) + 1;
					//u = L2_gpu[ch] + bwt_occ_gpu(bwt, prev_u, ch);
					l = L2_gpu[base] + intv.x + 1;
					u = L2_gpu[base] + intv.y;


					beg_i = (i == start) ? i - 1 : i;
				}
				//if (tid == 26 ||tid == 27 || tid == 28) printf("%d-->%d,%d, ",tid, u-l+1, itr);
//				uint32_t neighbour_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, 1);
//				uint32_t is_neighbour_active = __shfl_up_sync(0xFFFFFFFF, is_active, 1);
//				uint32_t neighbour_2_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, 2);
//				uint32_t is_neighbour_2_active = __shfl_up_sync(0xFFFFFFFF, is_active, 2);
//				uint32_t neighbour_3_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, 3);
//				uint32_t is_neighbour_3_active = __shfl_up_sync(0xFFFFFFFF, is_active, 3);
//				uint32_t neighbour_4_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, 4);
//				uint32_t is_neighbour_4_active = __shfl_up_sync(0xFFFFFFFF, is_active, 4);
//				uint32_t neighbour_5_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, 5);
//				uint32_t is_neighbour_5_active = __shfl_up_sync(0xFFFFFFFF, is_active, 5);

				for (m = 0; m <N_SHUFFLES; m++){
					uint32_t neighbour_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, m+1);
					uint32_t is_neighbour_active = __shfl_up_sync(0xFFFFFFFF, is_active, m+1);
					if(neighbour_active[m]) neighbour_active[m] = is_neighbour_active;
					if (is_shfl[m] && neighbour_active[m] && prev_intv_size[m] == neighbour_intv_size) {
						//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
						is_active = 0;
						is_smem = 0;
						break;
						//prev_seed_interval = make_uint2(m,m);
					}
					//if(is_shfl[m] == 0) break;
				}
				only_next_time = is_active ? only_next_time : only_next_time + 1;
				if(only_next_time == 2) break;


				for (m = N_SHUFFLES - 1; m >= 1; m--){
					prev_intv_size[m] = prev_intv_size[m-1];
				}
				prev_intv_size[0] = curr_intv_size;
//				if (is_shfl && is_neighbour_active && prev_intv_size == neighbour_intv_size) {
//					//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
//					is_active = 0;
//					//break;
//				}
//				if (is_shfl_2 && is_neighbour_2_active && prev_2_intv_size == neighbour_2_intv_size) {
//					//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
//					is_active = 0;
//					//break;
//				}
//
//				if (is_shfl_3 && is_neighbour_3_active && prev_3_intv_size == neighbour_3_intv_size) {
//					//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
//					is_active = 0;
//					//break;
//				}
//
//				if (is_shfl_4 && is_neighbour_4_active && prev_4_intv_size == neighbour_4_intv_size) {
//					//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
//					is_active = 0;
//					//break;
//				}
//				if (is_shfl_5 && is_neighbour_5_active && prev_5_intv_size == neighbour_5_intv_size) {
//					//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
//					is_active = 0;
//					//break;
//				}
//
//				prev_5_intv_size = prev_4_intv_size;
//				prev_4_intv_size = prev_3_intv_size;
//				prev_3_intv_size = prev_2_intv_size;
//				prev_2_intv_size = prev_intv_size;
//				prev_intv_size =  curr_intv_size;
				if (l > u || base > 3) {
					is_active = 0;
				}

				curr_intv_size =  l <= u ? u - l + 1 : curr_intv_size;

			}
		}
		if (read_len - thread_read_idx - beg_i + start - 1 >= min_seed_size && is_smem) {
			atomicAdd(&n_smems_fow[thread_read_num], 1);
			seed_intervals_fow_rev[intv_idx - thread_read_idx] = make_uint2(prev_seed_interval.x, prev_seed_interval.y);
			seed_read_pos_fow_rev[intv_idx - thread_read_idx] = make_int2 (beg_i - start + 1, read_len - thread_read_idx);
			is_smem_fow_rev_flag[intv_idx - thread_read_idx] = 0x00010001;
		}

	}
	else {
		int intv_idx = 2*(read_offsets[thread_read_num] - (thread_read_num*min_seed_size)) + read_len - min_seed_size;
		int start = read_off&7;
		uint32_t *seq = &(packed_read_batch_rev[read_off >> 3]);
		uint32_t pre_calc_seed = 0;
		for (i = start + thread_read_idx, j = 0; j < pre_calc_seed_len; i++, j++) {
			int reg_no = i >> 3;
			int reg_pos = i & 7;
			int reg = seq[reg_no];
			uint32_t base = (reg >> (28 - (reg_pos << 2)))&15;
			/*unknown bases*/
			if (base > 3) {

				break;
			}
			pre_calc_seed |= (base << (j<<1));
		}
		uint2 prev_seed_interval = j < pre_calc_seed_len ? make_uint2(1,0) : pre_calc_intervals[pre_calc_seed];
		int beg_i = i;
		if(prev_seed_interval.x <= prev_seed_interval.y) {
			is_active = 1;
			uint32_t curr_intv_size = prev_seed_interval.y - prev_seed_interval.x + 1;
//			uint32_t prev_intv_size = 0;
//			uint32_t prev_2_intv_size = 0;
//			uint32_t prev_3_intv_size = 0;
//			uint32_t prev_4_intv_size = 0;
//			uint32_t prev_5_intv_size = 0;
			l = prev_seed_interval.x, u = prev_seed_interval.y;
			for (; i < read_len + start; i++) {
				/*get the base*/
				if (is_active) {
					prev_seed_interval = make_uint2(l,u);
					int reg_no = i >> 3;
					int reg_pos = i & 7;
					int reg = seq[reg_no];
					int base = (reg >> (28 - (reg_pos << 2)))&15;
					/*unknown bases*/
//					if (base > 3) {
//						//break;
//						is_active = 0;
//					}

					uint2 intv = find_occ_gpu(bwt, prev_seed_interval.x - 1, prev_seed_interval.y, base);
					//calculate the range
					//l = L2_gpu[ch] + bwt_occ_gpu(bwt, prev_l - 1, ch) + 1;
					//u = L2_gpu[ch] + bwt_occ_gpu(bwt, prev_u, ch);
					l = L2_gpu[base] + intv.x + 1;
					u = L2_gpu[base] + intv.y;
//					if (l > u) {
//						//break;
//						is_active = 0;
//					}
					beg_i = i == (read_len + start - 1) ? read_len + start : i;
				}
//				uint32_t neighbour_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, 1);
//				uint32_t is_neighbour_active = __shfl_up_sync(0xFFFFFFFF, is_active, 1);
//				uint32_t neighbour_2_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, 2);
//				uint32_t is_neighbour_2_active = __shfl_up_sync(0xFFFFFFFF, is_active, 2);
//				uint32_t neighbour_3_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, 3);
//				uint32_t is_neighbour_3_active = __shfl_up_sync(0xFFFFFFFF, is_active, 3);
//				uint32_t neighbour_4_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, 4);
//				uint32_t is_neighbour_4_active = __shfl_up_sync(0xFFFFFFFF, is_active, 4);
//				uint32_t neighbour_5_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, 5);
//				uint32_t is_neighbour_5_active = __shfl_up_sync(0xFFFFFFFF, is_active, 5);
				for (m = 0; m < N_SHUFFLES; m++){
					uint32_t neighbour_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, m+1);
					uint32_t is_neighbour_active = __shfl_up_sync(0xFFFFFFFF, is_active, m+1);
					if(neighbour_active[m]) neighbour_active[m] = is_neighbour_active;
					if (is_shfl[m] && neighbour_active[m] && prev_intv_size[m] == neighbour_intv_size) {
						//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
						is_active = 0;
						is_smem = 0;
						break;
					}
					//if (is_shfl[m] == 0) break;
				}
				only_next_time = is_active ? only_next_time : only_next_time + 1;
				if(only_next_time == 2) break;

				for (m = N_SHUFFLES - 1 ; m >= 1; m--){
					prev_intv_size[m] = prev_intv_size[m-1];
				}
				prev_intv_size[0] = curr_intv_size;
//
//				if (is_shfl && is_neighbour_active && prev_intv_size == neighbour_intv_size) {
//					//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
//					is_active = 0;
//					//break;
//				}
//				if (is_shfl_2 && is_neighbour_2_active && prev_2_intv_size == neighbour_2_intv_size) {
//					//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
//					is_active = 0;
//					//break;
//				}

//
//				if (is_shfl_3 && is_neighbour_3_active && prev_3_intv_size == neighbour_3_intv_size) {
//					//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
//					is_active = 0;
//					//break;
//				}
//
//				if (is_shfl_4 && is_neighbour_4_active && prev_4_intv_size == neighbour_4_intv_size) {
//					//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
//					is_active = 0;
//					//break;
//				}
//				if (is_shfl_5 && is_neighbour_5_active && prev_5_intv_size == neighbour_5_intv_size) {
//					//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
//					is_active = 0;
//					//break;
//				}
//
//				prev_5_intv_size = prev_4_intv_size;
//				prev_4_intv_size = prev_3_intv_size;
//				prev_3_intv_size = prev_2_intv_size;
//				prev_2_intv_size = prev_intv_size;
//				prev_intv_size =  curr_intv_size;
				if (l > u || base > 3) {
					is_active = 0;
				}

				curr_intv_size =  l <= u ? u - l + 1 : curr_intv_size;
			}
			if (beg_i - start - thread_read_idx >= min_seed_size && is_smem) {
				atomicAdd(&n_smems_rev[thread_read_num], 1);
				seed_intervals_fow_rev[intv_idx + thread_read_idx] = make_uint2(prev_seed_interval.x, prev_seed_interval.y);
				//seed_read_pos_rev[intv_idx + thread_read_idx] =  make_int2 (read_len + start - beg_i, read_len - thread_read_idx) ;
				seed_read_pos_fow_rev[intv_idx + thread_read_idx] =  make_int2 (thread_read_idx|0x80000000, (beg_i - start)|0x80000000) ;
				is_smem_fow_rev_flag[intv_idx + thread_read_idx]=0x00010001;
			}
		}
//		seed_intervals_rev[intv_idx + thread_read_idx] = make_uint2(prev_seed_interval.x, prev_seed_interval.y);
//		//seed_read_pos_rev[intv_idx + thread_read_idx] =  make_int2 (read_len + start - beg_i, read_len - thread_read_idx) ;
//		seed_read_pos_rev[intv_idx + thread_read_idx] =  make_int2 (thread_read_idx, beg_i - start) ;
	}

	return;

}

__global__ void compute_intervals_fow_kernel(uint32_t *packed_read_batch, uint32_t *read_sizes, uint32_t *read_offsets, uint2 *seed_intervals,
		int2 *seed_read_pos, uint32_t *read_num, uint32_t *read_idx, uint32_t *is_smem_flag, uint32_t *n_smems, int min_seed_size, GPUseed_ref_index ref_index, int n_tasks) {

	int tid = (blockIdx.x*blockDim.x) + threadIdx.x;
	if (tid >= n_tasks) return;
	int thread_read_num = read_num[tid%n_tasks];
	int read_len = read_sizes[thread_read_num];
	int read_off = read_offsets[thread_read_num];
	//thread_read_num * (MAX_READ_LENGTH - min_seed_size);
	uint32_t thread_read_idx = read_idx[tid%n_tasks];
	int is_active = 0;
	int is_smem = 1;
	int is_shfl[N_SHUFFLES];
	int only_next_time = 0;
	uint32_t neighbour_active[N_SHUFFLES];
	uint32_t prev_intv_size[N_SHUFFLES];

	int m;
	for (m = 0; m < N_SHUFFLES; m++) {
		is_shfl[m] = tid%n_tasks ? 1 : 0;
		if (is_shfl[m]) is_shfl[m] = (thread_read_num == read_num[(tid%n_tasks) - (m+1)]) ? 1 : 0;
		if (is_shfl[m]) is_shfl[m] = ((tid%32) - m > 0) ? 1 : 0;
		prev_intv_size[m] = 0;
		neighbour_active[m] = 1;
	}

//	int is_shfl = tid%n_tasks ? 1 : 0;
//	if (is_shfl) is_shfl = (thread_read_num == read_num[(tid%n_tasks) - 1]) ? 1 : 0;
//	if (is_shfl) is_shfl = (tid%32) ? 1 : 0;
//	int is_shfl_2 = (thread_read_num == read_num[(tid%n_tasks) - 2]) ? 1 : 0;
//	if (is_shfl_2) is_shfl_2 = (tid%33) ? 1 : 0;
//	int is_shfl_3 = (thread_read_num == read_num[(tid%n_tasks) - 3]) ? 1 : 0;
//	if (is_shfl_3) is_shfl_3 = (tid%34) ? 1 : 0;
//	int is_shfl_4 = (thread_read_num == read_num[(tid%n_tasks) - 4]) ? 1 : 0;
//	if (is_shfl_4) is_shfl_4 = (tid%35) ? 1 : 0;
//	int is_shfl_5 = (thread_read_num == read_num[(tid%n_tasks) - 5]) ? 1 : 0;
//	if (is_shfl_5) is_shfl_5 = (tid%36) ? 1 : 0;
	int i, j;
	int base;
	bwtint_t l, u;
	int intv_idx = read_offsets[thread_read_num] - (thread_read_num*min_seed_size) + read_len - min_seed_size - 1;
	int start = read_off&7;
	uint32_t *seq = &(packed_read_batch[read_off >> 3]);
	uint32_t pre_calc_seed = 0;
	for (i = start + read_len - thread_read_idx - 1, j = 0; j < pre_calc_seed_len; i--, j++) {
		int reg_no = i >> 3;
		int reg_pos = i & 7;
		int reg = seq[reg_no];
		uint32_t base = (reg >> (28 - (reg_pos << 2)))&15;
		/*unknown bases*/
		if (base > 3) {
			break;
		}
		pre_calc_seed |= (base << (j<<1));
	}
	uint2 prev_seed_interval = j < pre_calc_seed_len ? make_uint2(1,0) : pre_calc_intervals[pre_calc_seed];
	int beg_i = i;
	if(prev_seed_interval.x <= prev_seed_interval.y) {
		is_active = 1;
		uint32_t curr_intv_size = prev_seed_interval.y - prev_seed_interval.x + 1;
		//			uint32_t prev_intv_size = 0;
		//			uint32_t prev_2_intv_size = 0;
		//			uint32_t prev_3_intv_size = 0;
		//			uint32_t prev_4_intv_size = 0;
		//			uint32_t prev_5_intv_size = 0;
		l = prev_seed_interval.x, u = prev_seed_interval.y;
		for (; i >= start; i--) {
			/*get the base*/
			if (is_active) {
				prev_seed_interval = make_uint2(l,u);
				int reg_no = i >> 3;
				int reg_pos = i & 7;
				int reg = seq[reg_no];
				int base = (reg >> (28 - (reg_pos << 2)))&15;
				/*unknown bases*/
				//					if (base > 3) {
				//						is_active = 0;
				//						break;
				//					}

				uint2 intv = find_occ_gpu(bwt, prev_seed_interval.x - 1, prev_seed_interval.y, base);
				//calculate the range
				//l = L2_gpu[ch] + bwt_occ_gpu(bwt, prev_l - 1, ch) + 1;
				//u = L2_gpu[ch] + bwt_occ_gpu(bwt, prev_u, ch);
				l = L2_gpu[base] + intv.x + 1;
				u = L2_gpu[base] + intv.y;


				beg_i = (i == start) ? i - 1 : i;
			}
			//if (tid == 26 ||tid == 27 || tid == 28) printf("%d-->%d,%d, ",tid, u-l+1, itr);
			//				uint32_t neighbour_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, 1);
			//				uint32_t is_neighbour_active = __shfl_up_sync(0xFFFFFFFF, is_active, 1);
			//				uint32_t neighbour_2_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, 2);
			//				uint32_t is_neighbour_2_active = __shfl_up_sync(0xFFFFFFFF, is_active, 2);
			//				uint32_t neighbour_3_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, 3);
			//				uint32_t is_neighbour_3_active = __shfl_up_sync(0xFFFFFFFF, is_active, 3);
			//				uint32_t neighbour_4_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, 4);
			//				uint32_t is_neighbour_4_active = __shfl_up_sync(0xFFFFFFFF, is_active, 4);
			//				uint32_t neighbour_5_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, 5);
			//				uint32_t is_neighbour_5_active = __shfl_up_sync(0xFFFFFFFF, is_active, 5);

			for (m = 0; m <N_SHUFFLES; m++){
				uint32_t neighbour_intv_size = __shfl_up_sync(0xFFFFFFFF, curr_intv_size, m+1);
				uint32_t is_neighbour_active = __shfl_up_sync(0xFFFFFFFF, is_active, m+1);
				if(neighbour_active[m]) neighbour_active[m] = is_neighbour_active;
				if (is_shfl[m] && neighbour_active[m] && prev_intv_size[m] == neighbour_intv_size) {
					//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
					is_active = 0;
					is_smem = 0;
					break;
					//prev_seed_interval = make_uint2(m,m);
				}
				//if(is_shfl[m] == 0) break;
			}
			only_next_time = is_active ? only_next_time : only_next_time + 1;
			if(only_next_time == 2) break;


			for (m = N_SHUFFLES - 1; m >= 1; m--){
				prev_intv_size[m] = prev_intv_size[m-1];
			}
			prev_intv_size[0] = curr_intv_size;
			//				if (is_shfl && is_neighbour_active && prev_intv_size == neighbour_intv_size) {
			//					//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
			//					is_active = 0;
			//					//break;
			//				}
			//				if (is_shfl_2 && is_neighbour_2_active && prev_2_intv_size == neighbour_2_intv_size) {
			//					//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
			//					is_active = 0;
			//					//break;
			//				}
			//
			//				if (is_shfl_3 && is_neighbour_3_active && prev_3_intv_size == neighbour_3_intv_size) {
			//					//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
			//					is_active = 0;
			//					//break;
			//				}
			//
			//				if (is_shfl_4 && is_neighbour_4_active && prev_4_intv_size == neighbour_4_intv_size) {
			//					//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
			//					is_active = 0;
			//					//break;
			//				}
			//				if (is_shfl_5 && is_neighbour_5_active && prev_5_intv_size == neighbour_5_intv_size) {
			//					//if (tid == 26 ||tid == 27 || tid == 28) printf("I am out thread_read_idx = %d, %d\n", thread_read_idx, i-start+1);
			//					is_active = 0;
			//					//break;
			//				}
			//
			//				prev_5_intv_size = prev_4_intv_size;
			//				prev_4_intv_size = prev_3_intv_size;
			//				prev_3_intv_size = prev_2_intv_size;
			//				prev_2_intv_size = prev_intv_size;
			//				prev_intv_size =  curr_intv_size;
			if (l > u || base > 3) {
				is_active = 0;
			}

			curr_intv_size =  l <= u ? u - l + 1 : curr_intv_size;

		}
	}
	if (read_len - thread_read_idx - beg_i + start - 1 >= min_seed_size && is_smem) {
		atomicAdd(&n_smems[thread_read_num], 1);
		seed_intervals[intv_idx - thread_read_idx] = make_uint2(prev_seed_interval.x, prev_seed_interval.y);
		seed_read_pos[intv_idx - thread_read_idx] = make_int2 (beg_i - start + 1, read_len - thread_read_idx);
		is_smem_flag[intv_idx - thread_read_idx] = 0x00010001;
	}



	return;

}
#endif
