#ifndef __KERNEL_SEQPAK__
#define __KERNEL_SEQPAK__

__global__ void pack_4bit_fow(uint32_t *read_batch, uint32_t* packed_read_batch_fow, int n_tasks) {

	int32_t i;
	const int32_t tid = (blockIdx.x * blockDim.x) + threadIdx.x; // thread index.
	uint32_t ascii_to_dna_reg_rev = 0x24031005;
	if (tid >= n_tasks) return;
	uint32_t *packed_read_batch_addr = &(read_batch[(tid << 1)]);
	uint32_t reg1 = packed_read_batch_addr[0]; //load 4 bases of the first sequence from global memory
	uint32_t reg2 = packed_read_batch_addr[1]; //load  another 4 bases of the S1 from global memory
	uint32_t pack_reg_4bit = 0;
	pack_reg_4bit |= ((ascii_to_dna_reg_rev >> ((reg1 & 7) << 2))&15)  << 28;        // ---
	pack_reg_4bit |=  ((ascii_to_dna_reg_rev >> (((reg1 >> 8) & 7) << 2))&15) << 24; //    |
	pack_reg_4bit |= ((ascii_to_dna_reg_rev >> (((reg1 >> 16) & 7) << 2))&15) << 20;//    |
	pack_reg_4bit |= ((ascii_to_dna_reg_rev >> (((reg1 >> 24) & 7) << 2))&15) << 16;//    |
	pack_reg_4bit |=  ((ascii_to_dna_reg_rev >> ((reg2 & 7) << 2))&15) << 12;        //     > pack data
	pack_reg_4bit |=  ((ascii_to_dna_reg_rev >> (((reg2 >> 8) & 7) << 2))&15) << 8;  //    |
	pack_reg_4bit |=  ((ascii_to_dna_reg_rev >> (((reg2 >> 16) & 7) << 2))&15) << 4; //    |
	pack_reg_4bit |= ((ascii_to_dna_reg_rev >> (((reg2 >> 24) & 7) << 2))&15);      //----
	packed_read_batch_fow[tid] = pack_reg_4bit; // write 8 bases of S1 packed into a unsigned 32 bit integer to global memory

}

__global__ void pack_4bit_rev(uint32_t *read_batch, uint32_t* packed_read_batch_rev, int n_tasks) {

	int32_t i;
	const int32_t tid = (blockIdx.x * blockDim.x) + threadIdx.x; // thread index.
	uint32_t ascii_to_dna_reg_rev = 0x14002035;
	if (tid >= n_tasks) return;
	uint32_t *packed_read_batch_addr = &(read_batch[(tid << 1)]);
	uint32_t reg1 = packed_read_batch_addr[0]; //load 4 bases of the first sequence from global memory
	uint32_t reg2 = packed_read_batch_addr[1]; //load  another 4 bases of the S1 from global memory
	uint32_t pack_reg_4bit = 0;
	pack_reg_4bit |= ((ascii_to_dna_reg_rev >> ((reg1 & 7) << 2))&15)  << 28;        // ---
	pack_reg_4bit |=  ((ascii_to_dna_reg_rev >> (((reg1 >> 8) & 7) << 2))&15) << 24; //    |
	pack_reg_4bit |= ((ascii_to_dna_reg_rev >> (((reg1 >> 16) & 7) << 2))&15) << 20;//    |
	pack_reg_4bit |= ((ascii_to_dna_reg_rev >> (((reg1 >> 24) & 7) << 2))&15) << 16;//    |
	pack_reg_4bit |=  ((ascii_to_dna_reg_rev >> ((reg2 & 7) << 2))&15) << 12;        //     > pack data
	pack_reg_4bit |=  ((ascii_to_dna_reg_rev >> (((reg2 >> 8) & 7) << 2))&15) << 8;  //    |
	pack_reg_4bit |=  ((ascii_to_dna_reg_rev >> (((reg2 >> 16) & 7) << 2))&15) << 4; //    |
	pack_reg_4bit |= ((ascii_to_dna_reg_rev >> (((reg2 >> 24) & 7) << 2))&15);      //----
	packed_read_batch_rev[tid] = pack_reg_4bit; // write 8 bases of S1 packed into a unsigned 32 bit integer to global memory

}

__global__ void assign_gpu_threads(uint32_t *thread_read_num, uint32_t *thread_read_idx, uint32_t *read_offsets, uint32_t* read_sizes, int min_seed_len, int n_tasks, int pack, int max_read_length) {

	int tid = (blockIdx.x*blockDim.x) + threadIdx.x;
	if (tid >= n_tasks) return;

	if (pack) {
		int read_len = read_sizes[tid]%8 ? read_sizes[tid] + (8 - read_sizes[tid]%8) : read_sizes[tid];
		uint32_t BLOCKDIM = (read_len >> 3) > 1024  ? 1024 : (read_len >> 3);
		uint32_t N_BLOCKS = ((read_len >> 3) + BLOCKDIM - 1) / BLOCKDIM;
//		cudaStream_t s;
//		cudaStreamCreateWithFlags(&s, cudaStreamNonBlocking);
//		int i;
//#pragma unroll
//		for (i = 0; i < read_len >> 3; i++){
//			thread_read_num[(read_offsets[tid] >> 3) + i] = tid;
//		}
//		assign_read_num<<<N_BLOCKS, BLOCKDIM, 0, s>>>(thread_read_num, read_offsets[tid] >> 3, tid, read_len >> 3);
//		cudaDeviceSynchronize();
	}
	else {
		int read_no = tid/max_read_length;
		int offset_in_read = tid - (read_no*max_read_length);
		if (offset_in_read >= read_sizes[read_no] - min_seed_len) return;
		thread_read_num[read_offsets[read_no] - (read_no*min_seed_len) + offset_in_read] = read_no;
		thread_read_idx[read_offsets[read_no] - (read_no*min_seed_len) + offset_in_read] = offset_in_read;
//		int i;
//		#pragma unroll
//		for (i = 0; i < read_sizes[tid] - min_seed_len; i++){
//			thread_read_num[read_offsets[tid] - (tid*min_seed_len) + i] = tid;
//			thread_read_idx[read_offsets[tid] - (tid*min_seed_len) + i] = i;
//		}
//		uint32_t BLOCKDIM = (read_len - min_seed_len) > 1024  ? 1024 : (read_len - min_seed_len);
//		uint32_t N_BLOCKS = ((read_len - min_seed_len) + BLOCKDIM - 1) / BLOCKDIM;
//		assign_threads_to_reads<<<N_BLOCKS, BLOCKDIM>>>(thread_read_num, thread_read_idx, read_offsets[tid] - (tid*min_seed_len), tid, (read_len - min_seed_len));
		//cudaDeviceSynchronize();
	}
	return;

}


#endif
