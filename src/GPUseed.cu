#include "gasal.h"
#include "args_parser.h"
#include "res.h"
#include "gasal_align.h"
#include "gasal_kernels.h"



inline void gasal_kernel_launcher(int32_t N_BLOCKS, int32_t BLOCKDIM, algo_type algo, comp_start start, gasal_gpu_storage_t *gpu_storage, int32_t actual_n_queries, int32_t k_band, data_source semiglobal_skipping_head, data_source semiglobal_skipping_tail, Bool secondBest)
{
	switch(algo)
	{
		
		KERNEL_SWITCH(LOCAL,		start, semiglobal_skipping_head, semiglobal_skipping_tail, secondBest);
		KERNEL_SWITCH(SEMI_GLOBAL,  start, semiglobal_skipping_head, semiglobal_skipping_tail, secondBest);		// MACRO that expands all 32 semi-global kernels
		KERNEL_SWITCH(GLOBAL,		start, semiglobal_skipping_head, semiglobal_skipping_tail, secondBest);
		KERNEL_SWITCH(BANDED,		start, semiglobal_skipping_head, semiglobal_skipping_tail, secondBest);

		default:
		break;

	}

}


//GASAL2 asynchronous (a.k.a non-blocking) alignment function
void gasal_aln_async(gasal_gpu_storage_t *gpu_storage, const uint32_t actual_query_batch_bytes, const uint32_t actual_n_queries, const uint32_t actual_max_query_size, GPUseed_ref_index ref_index,  int min_seed_size,  Parameters *params, uint32_t max_n_seeds = gpu_storage->max_n_seeds) {

	cudaError_t err;
	if (actual_n_queries <= 0) {
		fprintf(stderr, "[GASAL ERROR:] actual_n_queries <= 0\n");
		exit(EXIT_FAILURE);
	}
	if (actual_query_batch_bytes <= 0) {
		fprintf(stderr, "[GASAL ERROR:] actual_query_batch_bytes <= 0\n");
		exit(EXIT_FAILURE);
	}

	int actual_query_batch_bytes_8 = actual_query_batch_bytes%8 ? actual_query_batch_bytes + (8 - actual_query_batch_bytes%8) : actual_query_batch_bytes;


	if (actual_n_queries > gpu_storage->max_n_queries) {
			fprintf(stderr, "[GASAL ERROR:] actual_n_queries(%d) > max_n_queries(%d)\n", actual_n_queries, gpu_storage->max_n_queries);
			exit(EXIT_FAILURE);
	}

	//--------------if pre-allocated memory is less, allocate more--------------------------
	if (gpu_storage->max_query_batch_bytes < actual_query_batch_bytes_8) {

		int i = 2;
		while ( (gpu_storage->max_query_batch_bytes * i) < actual_query_batch_bytes) i++;
		gpu_storage->max_query_batch_bytes = gpu_storage->max_query_batch_bytes * i;

		fprintf(stderr, "[GASAL WARNING:] actual_query_batch_bytes(%d) > Allocated GPU memory (max_query_batch_bytes=%d). Therefore, allocating %d bytes on GPU (max_query_batch_bytes=%d). Performance may be lost if this is repeated many times.\n", actual_query_batch_bytes, gpu_storage->max_query_batch_bytes, gpu_storage->max_query_batch_bytes*i, gpu_storage->max_query_batch_bytes*i);


		if (gpu_storage->unpacked_query_batch != NULL) CHECKCUDAERROR(cudaFree(gpu_storage->unpacked_query_batch), "cudaFree");
		if (gpu_storage->packed_fow_query_batch != NULL) CHECKCUDAERROR(cudaFree(gpu_storage->packed_fow_query_batch), "cudaFree");
		if (gpu_storage->packed_rev_query_batch != NULL && params->doRev) CHECKCUDAERROR(cudaFree(gpu_storage->packed_fow_query_batch), "cudaFree");
		if (gpu_storage->n_seed_occ_scan != NULL) CHECKCUDAERROR(cudaFree(gpu_storage->n_seed_occ_scan), "cudaFree");
		if (gpu_storage->cub_select_temp_storage != NULL) CHECKCUDAERROR(cudaFree(gpu_storage->cub_select_temp_storage), "cudaFree");
		if (gpu_storage->cub_sort_temp_storage != NULL && !params->doMem) CHECKCUDAERROR(cudaFree(gpu_storage->cub_sort_temp_storage), "cudaFree");
		if (gpu_storage->cub_scan_temp_storage != NULL) CHECKCUDAERROR(cudaFree(gpu_storage->cub_scan_temp_storage), "cudaFree");


		CHECKCUDAERROR(cudaMalloc(&(gpu_storage->unpacked_query_batch), gpu_storage->max_query_batch_bytes * sizeof(uint8_t)), "cudaMalloc");
		CHECKCUDAERROR(cudaMalloc(&(gpu_storage->packed_fow_query_batch), (gpu_storage->max_query_batch_bytes/8) * sizeof(uint32_t)), "cudaMalloc");
		if(params->doRev)  CHECKCUDAERROR(cudaMalloc(&(gpu_storage->packed_rev_query_batch), (gpu_storage->max_query_batch_bytes/8) * sizeof(uint32_t)), "cudaMalloc");
		CHECKCUDAERROR(cudaMalloc(&(gpu_storage->n_seed_occ_scan), ((params->doRev ? 2 : 1)*(gpu_storage->max_query_batch_bytes) + 1)*sizeof(uint32_t)), "cudaMalloc");

		gpu_storage->cub_select_temp_storage = NULL;
		gpu_storage->cub_select_storage_bytes = 0;
		CHECKCUBERROR(cub::DeviceSelect::Flagged(gpu_storage->cub_select_temp_storage, gpu_storage->cub_select_storage_bytes, (uint32_t*)(gpu_storage->seed_intervals), (uint16_t*)(gpu_storage->is_smem_flag), (uint32_t*)(gpu_storage->seed_intervals_compact), gpu_storage->n_smems_sum, 2*(params->doRev ? 2 : 1)*(gpu_storage->max_query_batch_bytes), gpu_storage->str), "DeviceSelect::Flagged");
		CHECKCUDAERROR(cudaMalloc(&(gpu_storage->cub_select_temp_storage), gpu_storage->cub_select_storage_bytes), "cudaMalloc");

		if (!params->doMem) {
		gpu_storage->max_query_size = gpu_storage->max_query_size < actual_max_query_size ? actual_max_query_size : gpu_storage->max_query_size;
		gpu_storage->cub_sort_temp_storage = NULL;
		gpu_storage->cub_sort_storage_bytes = 0;
		int n_bits_max_query_size = (int)ceil(log2((double)gpu_storage->max_query_size));
		CHECKCUBERROR(cub::DeviceSegmentedRadixSort::SortPairs(gpu_storage->cub_sort_temp_storage, gpu_storage->cub_sort_storage_bytes, (uint64_t*)(gpu_storage->seed_query_pos_compact), (uint64_t*)(gpu_storage->seed_query_pos), (uint64_t*)(gpu_storage->seed_intervals_compact), (uint64_t*)(gpu_storage->seed_intervals),   (params->doRev ? 2 : 1)*(gpu_storage->max_query_batch_bytes), gpu_storage->max_n_queries, gpu_storage->n_smems_scan, gpu_storage->n_smems_scan + 1, 0, n_bits_max_query_size, gpu_storage->str), "DeviceSegmentedRadixSort::SortPairs");
		CHECKCUDAERROR(cudaMalloc(&(gpu_storage->cub_sort_temp_storage),gpu_storage->cub_sort_storage_bytes), "cudaMalloc");
		}

		gpu_storage->cub_scan_temp_storage = NULL;
		gpu_storage->cub_scan_storage_bytes = 0;
		CHECKCUBERROR(cub::DeviceScan::ExclusiveSum(gpu_storage->cub_scan_temp_storage, gpu_storage->cub_scan_storage_bytes, gpu_storage->n_seed_occ, gpu_storage->n_seed_occ_scan, (params->doRev ? 2 : 1)*(gpu_storage->max_query_batch_bytes), gpu_storage->str), "DeviceScan::ExclusiveSum");
		CHECKCUDAERROR(cudaMalloc(&(gpu_storage->cub_scan_temp_storage),gpu_storage->cub_scan_storage_bytes), "cudaMalloc");

	}


	if (gpu_storage->max_n_queries < actual_n_queries) {

		int i = 2;
		while ( (gpu_storage->max_n_queries * i) < actual_n_queries) i++;
		gpu_storage->max_n_queries = gpu_storage->max_n_queries * i;

		fprintf(stderr, "[GASAL WARNING:] actual_n_queries(%d) > max_n_queries(%d). Therefore, allocating memory for %d alignments on  GPU (max_n_queries=%d). Performance may be lost if this is repeated many times.\n", actual_n_queries, gpu_storage->max_n_queries, gpu_storage->max_n_queries*i, gpu_storage->max_n_queries*i);

		if (gpu_storage->query_batch_offsets != NULL) CHECKCUDAERROR(cudaFree(gpu_storage->query_batch_offsets), "cudaFree");
		if (gpu_storage->query_batch_lens != NULL) CHECKCUDAERROR(cudaFree(gpu_storage->query_batch_lens), "cudaFree");
		if (gpu_storage->n_smems_fow != NULL) CHECKCUDAERROR(cudaFree(gpu_storage->n_smems_fow), "cudaFree");
		if (gpu_storage->n_smems_rev != NULL && params->doRev) CHECKCUDAERROR(cudaFree(gpu_storage->n_smems_rev), "cudaFree");
		if (gpu_storage->n_smems_scan != NULL) CHECKCUDAERROR(cudaFree(gpu_storage->n_smems_scan), "cudaFree");
		if (gpu_storage->cub_sort_temp_storage != NULL ) CHECKCUDAERROR(cudaFree(gpu_storage->cub_sort_temp_storage), "cudaFree");

		CHECKCUDAERROR(cudaMalloc(&(gpu_storage->query_batch_lens), gpu_storage->max_n_queries * sizeof(uint32_t)), "cudaMalloc");
		CHECKCUDAERROR(cudaMalloc(&(gpu_storage->query_batch_offsets), gpu_storage->max_n_queries * sizeof(uint32_t)), "cudaMalloc");
		CHECKCUDAERROR(cudaMalloc(&(gpu_storage->n_smems_fow), gpu_storage->max_n_queries*sizeof(uint32_t)), "cudaMalloc");
		if (params->doRev) CHECKCUDAERROR(cudaMalloc(&(gpu_storage->n_smems_rev), gpu_storage->max_n_queries*sizeof(uint32_t)), "cudaMalloc");
		CHECKCUDAERROR(cudaMalloc(&(gpu_storage->n_smems_scan),(gpu_storage->max_n_queries+1)*sizeof(uint32_t)), "cudaMalloc");

		gpu_storage->max_query_size = gpu_storage->max_query_size < actual_max_query_size ? actual_max_query_size : gpu_storage->max_query_size;
		gpu_storage->cub_sort_temp_storage = NULL;
		gpu_storage->cub_sort_storage_bytes = 0;
		int n_bits_max_query_size = (int)ceil(log2((double)gpu_storage->max_query_size));
		if (params->doMem) {
			CHECKCUBERROR(cub::DeviceSegmentedRadixSort::SortPairs(gpu_storage->cub_sort_temp_storage, gpu_storage->cub_sort_storage_bytes, (uint64_t*)(gpu_storage->seed_query_pos_compact), (uint64_t*)(gpu_storage->seed_query_pos), gpu_storage->seed_sa_idx, gpu_storage->seed_ref_pos,  gpu_storage->max_n_seeds, gpu_storage->max_n_queries, gpu_storage->n_query_seeds_scan, gpu_storage->n_query_seeds_scan + 1, 0, n_bits_max_query_size, gpu_storage->str), "DeviceSegmentedRadixSort::SortPairs");

		} else {
			CHECKCUBERROR(cub::DeviceSegmentedRadixSort::SortPairs(gpu_storage->cub_sort_temp_storage, gpu_storage->cub_sort_storage_bytes, (uint64_t*)(gpu_storage->seed_query_pos_compact), (uint64_t*)(gpu_storage->seed_query_pos), (uint64_t*)(gpu_storage->seed_intervals_compact), (uint64_t*)(gpu_storage->seed_intervals),   (params->doRev ? 2 : 1)*(gpu_storage->max_query_batch_bytes), gpu_storage->max_n_queries, gpu_storage->n_smems_scan, gpu_storage->n_smems_scan + 1, 0, n_bits_max_query_size, gpu_storage->str), "DeviceSegmentedRadixSort::SortPairs");
		}
		CHECKCUDAERROR(cudaMalloc(&(gpu_storage->cub_sort_temp_storage),gpu_storage->cub_sort_storage_bytes), "cudaMalloc");



	}
	if (gpu_storage->max_query_size < actual_max_query_size) {

		gpu_storage->max_query_size = actual_max_query_size;

		if (gpu_storage->cub_sort_temp_storage != NULL ) CHECKCUDAERROR(cudaFree(gpu_storage->cub_sort_temp_storage), "cudaFree");
		gpu_storage->cub_sort_temp_storage = NULL;
		gpu_storage->cub_sort_storage_bytes = 0;
		int n_bits_max_query_size = (int)ceil(log2((double)gpu_storage->max_query_size));
		if (params->doMem) {
			CHECKCUBERROR(cub::DeviceSegmentedRadixSort::SortPairs(gpu_storage->cub_sort_temp_storage, gpu_storage->cub_sort_storage_bytes, (uint64_t*)(gpu_storage->seed_query_pos_compact), (uint64_t*)(gpu_storage->seed_query_pos), gpu_storage->seed_sa_idx, gpu_storage->seed_ref_pos,  gpu_storage->max_n_seeds, gpu_storage->max_n_queries, gpu_storage->n_query_seeds_scan, gpu_storage->n_query_seeds_scan + 1, 0, n_bits_max_query_size, gpu_storage->str), "DeviceSegmentedRadixSort::SortPairs");

		} else {
			CHECKCUBERROR(cub::DeviceSegmentedRadixSort::SortPairs(gpu_storage->cub_sort_temp_storage, gpu_storage->cub_sort_storage_bytes, (uint64_t*)(gpu_storage->seed_query_pos_compact), (uint64_t*)(gpu_storage->seed_query_pos), (uint64_t*)(gpu_storage->seed_intervals_compact), (uint64_t*)(gpu_storage->seed_intervals),   (params->doRev ? 2 : 1)*(gpu_storage->max_query_batch_bytes), gpu_storage->max_n_queries, gpu_storage->n_smems_scan, gpu_storage->n_smems_scan + 1, 0, n_bits_max_query_size, gpu_storage->str), "DeviceSegmentedRadixSort::SortPairs");
		}
		CHECKCUDAERROR(cudaMalloc(&(gpu_storage->cub_sort_temp_storage),gpu_storage->cub_sort_storage_bytes), "cudaMalloc");
	}

	//------------------------------------------

	//------------------------launch copying of sequence batches from CPU to GPU---------------------------

	// here you can track the evolution of your data structure processing with the printer: gasal_host_batch_printall(current);

	host_batch_t *current = gpu_storage->extensible_host_unpacked_query_batch;
	while (current != NULL)
	{
		if (current->next != NULL ) 
		{
			CHECKCUDAERROR(cudaMemcpyAsync( &(gpu_storage->unpacked_query_batch[current->offset]), 
											current->data, 
											current->next->offset - current->offset,
											cudaMemcpyHostToDevice, 
											gpu_storage->str ), "cudaMemcpyAsync" );
			
		} else {
			// it's the last page to copy
			CHECKCUDAERROR(cudaMemcpyAsync( &(gpu_storage->unpacked_query_batch[current->offset]), 
											current->data, 
											actual_query_batch_bytes - current->offset, 
											cudaMemcpyHostToDevice, 
											gpu_storage->str ), "cudaMemcpyAsync");
		}
		current = current->next;
	}

	// We could reverse-complement before packing, but we would get 2x more read-writes to memory.

	//----------------------launch copying of sequence offsets and lengths from CPU to GPU--------------------------------------
	CHECKCUDAERROR(cudaMemcpyAsync(gpu_storage->query_batch_lens, gpu_storage->host_query_batch_lens, actual_n_queries * sizeof(uint32_t), cudaMemcpyHostToDevice, gpu_storage->str), "cudaMemcpyAsync");
	CHECKCUDAERROR(cudaMemcpyAsync(gpu_storage->query_batch_offsets, gpu_storage->host_query_batch_offsets, actual_n_queries * sizeof(uint32_t), cudaMemcpyHostToDevice,  gpu_storage->str), "cudaMemcpyAsync");



	//-----------------------------------------------------------------------------------------------------------
	// TODO: Adjust the block size depending on the kernel execution.
	
	//-------------------------------------------launch packing kernel-------------------------------------------------
	int BLOCKDIM = 128;
	int N_BLOCKS = ((actual_query_batch_bytes_8>> 3)  + BLOCKDIM - 1)/BLOCKDIM;
	pack_4bit_fow<<<N_BLOCKS, BLOCKDIM, 0 , gpu_storage->str>>>((uint32_t*)gpu_storage->unpacked_query_batch, gpu_storage->packed_rev_query_batch, actual_query_batch_bytes_8 >> 3);
	cudaError_t pack_fow_kernel_err = cudaGetLastError();
	if ( cudaSuccess != pack_fow_kernel_err )
	{
		fprintf(stderr, "[GPUseed CUDA ERROR:]  Error in forward packing kernel: %s(CUDA error no.=%d). Line no. %d in file %s\n", cudaGetErrorString(pack_fow_kernel_err), pack_fow_kernel_err,  __LINE__, __FILE__);
		exit(EXIT_FAILURE);
	}

	if(params->doRev) {
		pack_4bit_rev<<<N_BLOCKS, BLOCKDIM, 0, gpu_storage->str>>>((uint32_t*)gpu_storage->unpacked_query_batch, gpu_storage->packed_rev_query_batch, actual_query_batch_bytes_8 >> 3);
		cudaError_t pack_rev_kernel_err = cudaGetLastError();
		if ( cudaSuccess != pack_rev_kernel_err )
		{
			fprintf(stderr, "[GPUseed CUDA ERROR:] Error in reverse packing kernel: %s(CUDA error no.=%d). Line no. %d in file %s\n", cudaGetErrorString(pack_rev_kernel_err), pack_rev_kernel_err,  __LINE__, __FILE__);
			exit(EXIT_FAILURE);
		}
	}

    
	//--------------------------------------------------------------------------------------------------------------------------

	N_BLOCKS = ((actual_n_queries*actual_max_query_size) + BLOCKDIM - 1)/BLOCKDIM;
	assign_gpu_threads<<<N_BLOCKS, BLOCKDIM, 0, gpu_storage->str>>>(gpu_storage->gpu_thread_query_num, gpu_storage->gpu_thread_query_idx, gpu_storage->query_batch_offsets, gpu_storage->query_batch_lens, min_seed_size, actual_n_queries*actual_max_query_size, 0, actual_max_query_size);

	cudaError_t assign_kernel_err = cudaGetLastError();
	if ( cudaSuccess != assign_kernel_err )
	{
		fprintf(stderr, "[GPUseed CUDA ERROR:] Error in assign_gpu_threads kernel: %s(CUDA error no.=%d). Line no. %d in file %s\n", cudaGetErrorString(assign_kernel_err), assign_kernel_err,  __LINE__, __FILE__);
		exit(EXIT_FAILURE);
	}

    //--------------------------------------launch alignment kernels--------------------------------------------------------------
	
	gpu_storage->seed_intervals = (uint2*)(&gpu_storage->seed_query_pos[(params->doRev ? 2 : 1)*(actual_query_batch_bytes - (actual_n_queries*min_seed_size))]);


	CHECKCUDAERROR(cudaMemsetAsync(gpu_storage->is_smem_flag, 0, (params->doRev ? 2 : 1)*(actual_query_batch_bytes - (actual_n_queries*min_seed_size))*sizeof(uint32_t), gpu_storage->str), "cudaMemsetAsync");
	CHECKCUDAERROR(cudaMemsetAsync(gpu_storage->n_smems_fow, 0, actual_n_queries*sizeof(uint32_t), gpu_storage->str), "cudaMemsetAsync");
	if(params->doRev) CHECKCUDAERROR(cudaMemsetAsync(gpu_storage->n_smems_rev, 0, actual_n_queries*sizeof(uint32_t), gpu_storage->str), "cudaMemsetAsync");

	int n_seed_cands = actual_query_batch_bytes - (actual_n_queries*params->min_seed_size);
	N_BLOCKS = (((params->doRev ? 2 : 1)*n_seed_cands) + BLOCKDIM - 1)/BLOCKDIM;
	if (params->doRev) {

		compute_intervals_fow_rev_kernel<<<N_BLOCKS, BLOCKDIM, 0, gpu_storage->str>>>(gpu_storage->packed_fow_query_batch, gpu_storage->packed_rev_query_batch,  gpu_storage->query_batch_lens, gpu_storage->query_batch_offsets, gpu_storage->seed_intervals,
				gpu_storage->seed_query_pos, gpu_storage->gpu_thread_query_num, gpu_storage->gpu_thread_query_idx, gpu_storage->is_smem_flag, gpu_storage->n_smems_fow, gpu_storage->n_smems_rev, min_seed_size , ref_index,  n_seed_cands);

	} else {
		N_BLOCKS = (n_seed_cands + BLOCKDIM - 1)/BLOCKDIM;
		compute_intervals_fow_kernel<<<N_BLOCKS, BLOCKDIM, 0, gpu_storage->str>>>(gpu_storage->packed_fow_query_batch, gpu_storage->query_batch_lens, gpu_storage->query_batch_offsets, gpu_storage->seed_intervals,
						gpu_storage->seed_query_pos, gpu_storage->gpu_thread_query_num, gpu_storage->gpu_thread_query_idx, gpu_storage->is_smem_flag, gpu_storage->n_smems_fow, min_seed_size , ref_index,  n_seed_cands);
	}

	cudaError_t comp_intv_kernel_err = cudaGetLastError();
	if ( cudaSuccess != comp_intv_kernel_err )
	{
		fprintf(stderr, "[GPUseed CUDA ERROR:] Error in compute_intervals kernel: %s(CUDA error no.=%d). Line no. %d in file %s\n", cudaGetErrorString(comp_intv_kernel_err), comp_intv_kernel_err,  __LINE__, __FILE__);
		exit(EXIT_FAILURE);
	}

	CHECKCUBERROR(cub::DeviceReduce::Max(gpu_storage->cub_scan_temp_storage, gpu_storage->cub_scan_storage_bytes, gpu_storage->n_smems_fow, gpu_storage->n_smems_max[0], actual_n_queries, gpu_storage->str), "DeviceReduce::Max");
	if (params->doRev) CHECKCUBERROR(cub::DeviceReduce::Max(gpu_storage->cub_scan_temp_storage, gpu_storage->cub_scan_storage_bytes, gpu_storage->n_smems_rev, gpu_storage->n_smems_max[1], actual_n_queries, gpu_storage->str), "DeviceReduce::Max");


	gpu_storage->n_smems_sum = &(gpu_storage->n_smems_scan[actual_n_queries]);

	CHECKCUBERROR(cub::DeviceSelect::Flagged(gpu_storage->cub_select_temp_storage, gpu_storage->cub_select_storage_bytes, (uint32_t*)gpu_storage->seed_intervals, (uint16_t*)gpu_storage->is_smem_flag, (uint32_t*)gpu_storage->seed_intervals_compact, gpu_storage->n_smems_sum, 2*2*(actual_query_batch_bytes - (actual_n_queries*min_seed_size)), gpu_storage->str), "DeviceSelect::Flagged");
	CHECKCUBERROR(cub::DeviceSelect::Flagged(gpu_storage->cub_select_temp_storage, gpu_storage->cub_select_storage_bytes, (uint32_t*)gpu_storage->seed_query_pos, (uint16_t*)gpu_storage->is_smem_flag, (uint32_t*)gpu_storage->seed_query_pos_compact, gpu_storage->n_smems_sum, 2*2*(actual_query_batch_bytes - (actual_n_queries*min_seed_size)), gpu_storage->str), "DeviceSelect::Flagged");


	if (params->doRev) {
		gpu_storage->n_smems_fow_rev = gpu_storage->query_batch_lens;
		N_BLOCKS = (actual_n_queries + BLOCKDIM - 1)/BLOCKDIM;
		sum_arrays<<<N_BLOCKS, BLOCKDIM, 0, gpu_storage->str>>>(gpu_storage->n_smems_fow, gpu_storage->n_smems_rev, gpu_storage->n_smems_fow_rev, actual_n_queries);
	}


	if (params->doRev) CHECKCUBERROR(cub::DeviceScan::ExclusiveSum(gpu_storage->cub_scan_temp_storage, gpu_storage->cub_scan_storage_bytes, gpu_storage->n_smems_fow_rev, gpu_storage->n_smems_scan, actual_n_queries, gpu_storage->str), "DeviceScan::ExclusiveSum");
	else CHECKCUBERROR(cub::DeviceScan::ExclusiveSum(gpu_storage->cub_scan_temp_storage, gpu_storage->cub_scan_storage_bytes, gpu_storage->n_smems_fow, gpu_storage->n_smems_scan, actual_n_queries, gpu_storage->str), "DeviceScan::ExclusiveSum");
     //-----------------------------------------------------------------------------------------------------------------------
	gpu_storage->n_query_seeds = gpu_storage->query_batch_offsets;

	CHECKCUDAERROR(cudaMemsetAsync(gpu_storage->n_query_seeds, 0, actual_n_queries*sizeof(uint32_t), gpu_storage->str), "cudaMemsetAsync");



	int n_bits_actual_max_query_size = (int)ceil(log2((double)actual_max_query_size));
	if (params->doSmem) {
		filter_smem_intervals_kernel_wrapper<<<1, 1, 0, gpu_storage->str>>>(gpu_storage->seed_intervals_compact, gpu_storage->seed_query_pos_compact, gpu_storage->seed_intervals, gpu_storage->seed_query_pos,   gpu_storage->n_smems_fow, gpu_storage->n_smems_rev, gpu_storage->n_smems_fow_rev,  gpu_storage->n_seed_occ, gpu_storage->n_smems_scan,  gpu_storage->n_query_seeds, gpu_storage->n_smems_max, gpu_storage->n_smems_sum, gpu_storage->cub_sort_temp_storage, gpu_storage->cub_sort_storage_bytes, actual_n_queries, n_bits_actual_max_query_size, params->doRev);
	}else if (params->doMem) {
		filter_mem_intervals_kernel_wrapper<<<1, 1, 0, gpu_storage->str>>>(gpu_storage->seed_intervals_compact, gpu_storage->seed_query_pos_compact, gpu_storage->seed_intervals, gpu_storage->n_smems_fow, gpu_storage->n_smems_rev, gpu_storage->n_seed_occ,  gpu_storage->n_smems_scan, gpu_storage->n_query_seeds, gpu_storage->n_smems_max, gpu_storage->n_smems_sum, actual_n_queries, params->doRev);

		int2 *swap =  gpu_storage->seed_query_pos_compact;
		gpu_storage->seed_query_pos_compact = gpu_storage->seed_query_pos;
		gpu_storage->seed_query_pos = swap;

	}




	//uint32_t host_n_smems_sum;
	//CHECKCUDAERROR(cudaMemcpyAsync(&host_n_smems_sum, gpu_storage->n_smems_sum, sizeof(uint32_t), cudaMemcpyDeviceToHost, gpu_storage->str), "cudaMemcpyAsync");
    //------------------------0launch the copying of alignment results from GPU to CPU--------------------------------------

	CHECKCUBERROR(cub::DeviceReduce::Sum(gpu_storage->cub_scan_temp_storage, gpu_storage->cub_scan_storage_bytes, gpu_storage->n_seed_occ , gpu_storage->n_seed_occ_sum,  host_n_smems_sum, gpu_storage->str), "DeviceReduce::Sum");



	find_n_query_to_process<<<1, 1, 0, gpu_storage->str>>> (gpu_storage->n_query_seeds, gpu_storage->n_smems_scan, gpu_storage->n_seed_occ, gpu_storage->n_seed_occ_sum, gpu_storage->n_smems_sum, gpu_storage->n_queries_processed, actual_n_queries, gpu_storage->max_n_seeds);


	CHECKCUBERROR(cub::DeviceScan::ExclusiveSum(gpu_storage->cub_scan_temp_storage, gpu_storage->cub_scan_storage_bytes, gpu_storage->n_seed_occ, gpu_storage->n_seed_occ_scan, host_n_smems_sum, gpu_storage->str), "DeviceScan::ExclusiveSum");
	//-----------------------------------------------------------------------------------------------------------------------
	
	uint64_t host_n_seed_occ_sum;
	CHECKCUDAERROR(cudaMemcpyAsync(&host_n_seed_occ_sum, gpu_storage->n_seed_occ_sum, sizeof(uint64_t), cudaMemcpyDeviceToHost, gpu_storage->str), "cudaMemcpyAsync");

	// not really needed to filter with params->secondBest, since all the pointers will be null and non-initialized.


    gpu_storage->is_free = 0; //set the availability of current stream to false

}


int gasal_is_aln_async_done(gasal_gpu_storage_t *gpu_storage) 
{
	cudaError_t err;
	if(gpu_storage->is_free == 1) return -2;//if no work is launched in this stream, return -2
	err = cudaStreamQuery(gpu_storage->str);//check to see if the stream is finished
	if (err != cudaSuccess ) {
		if (err == cudaErrorNotReady) return -1;
		else{
			fprintf(stderr, "[GASAL CUDA ERROR:] %s(CUDA error no.=%d). Line no. %d in file %s\n", cudaGetErrorString(err), err,  __LINE__, __FILE__);
			exit(EXIT_FAILURE);
		}
	}
	gpu_storage->is_free = 1;

	return 0;
}


void gasal_copy_subst_scores(gasal_subst_scores *subst){

	cudaError_t err;
	CHECKCUDAERROR(cudaMemcpyToSymbol(_cudaGapO, &(subst->gap_open), sizeof(int32_t), 0, cudaMemcpyHostToDevice));
	CHECKCUDAERROR(cudaMemcpyToSymbol(_cudaGapExtend, &(subst->gap_extend), sizeof(int32_t), 0, cudaMemcpyHostToDevice));
	int32_t gapoe = (subst->gap_open + subst->gap_extend);
	CHECKCUDAERROR(cudaMemcpyToSymbol(_cudaGapOE, &(gapoe), sizeof(int32_t), 0, cudaMemcpyHostToDevice));
	CHECKCUDAERROR(cudaMemcpyToSymbol(_cudaMatchScore, &(subst->match), sizeof(int32_t), 0, cudaMemcpyHostToDevice));
	CHECKCUDAERROR(cudaMemcpyToSymbol(_cudaMismatchScore, &(subst->mismatch), sizeof(int32_t), 0, cudaMemcpyHostToDevice));
	return;
}

