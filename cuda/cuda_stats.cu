/*  cuda_stats.cu -- CUDA GPU-accelerated statistics computation for BAM data.

    Copyright (C) 2025 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include "cuda_stats.cuh"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cooperative_groups.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <cub/cub.cuh>

namespace cg = cooperative_groups;

// CUDA kernels for statistics computation

__global__ void cuda_basic_stats_kernel(const cuda_bam_record_t *records, int num_records,
                                       cuda_basic_stats_t *stats) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    
    // Shared memory for block-level reduction
    __shared__ uint64_t shared_mapped[CUDA_STATS_BLOCK_SIZE];
    __shared__ uint64_t shared_unmapped[CUDA_STATS_BLOCK_SIZE];
    __shared__ uint64_t shared_total_length[CUDA_STATS_BLOCK_SIZE];
    __shared__ uint64_t shared_primary[CUDA_STATS_BLOCK_SIZE];
    __shared__ uint64_t shared_supplementary[CUDA_STATS_BLOCK_SIZE];
    
    int local_tid = threadIdx.x;
    
    // Initialize shared memory
    shared_mapped[local_tid] = 0;
    shared_unmapped[local_tid] = 0;
    shared_total_length[local_tid] = 0;
    shared_primary[local_tid] = 0;
    shared_supplementary[local_tid] = 0;
    
    // Process records
    for (int i = tid; i < num_records; i += stride) {
        const cuda_bam_record_t *record = &records[i];
        
        if (!(record->flag & 0x4)) {  // BAM_FUNMAP
            shared_mapped[local_tid]++;
            shared_total_length[local_tid] += record->seq_len;
            
            if (!(record->flag & 0x100)) {  // BAM_FSECONDARY
                shared_primary[local_tid]++;
            }
            if (record->flag & 0x800) {  // BAM_FSUPPLEMENTARY
                shared_supplementary[local_tid]++;
            }
        } else {
            shared_unmapped[local_tid]++;
        }
    }
    
    __syncthreads();
    
    // Block-level reduction
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (local_tid < s) {
            shared_mapped[local_tid] += shared_mapped[local_tid + s];
            shared_unmapped[local_tid] += shared_unmapped[local_tid + s];
            shared_total_length[local_tid] += shared_total_length[local_tid + s];
            shared_primary[local_tid] += shared_primary[local_tid + s];
            shared_supplementary[local_tid] += shared_supplementary[local_tid + s];
        }
        __syncthreads();
    }
    
    // Update global statistics
    if (local_tid == 0) {
        atomicAdd(&stats->mapped_reads, shared_mapped[0]);
        atomicAdd(&stats->unmapped_reads, shared_unmapped[0]);
        atomicAdd(&stats->total_length, shared_total_length[0]);
        atomicAdd(&stats->primary_mapped, shared_primary[0]);
        atomicAdd(&stats->supplementary, shared_supplementary[0]);
    }
}

__global__ void cuda_quality_stats_kernel(const cuda_bam_record_t *records, int num_records,
                                         cuda_quality_stats_t *stats) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    
    __shared__ uint64_t shared_quality_sum[CUDA_STATS_BLOCK_SIZE];
    __shared__ uint32_t shared_mapq_hist[CUDA_MAX_QUALITY_SCORE + 1];
    
    int local_tid = threadIdx.x;
    shared_quality_sum[local_tid] = 0;
    
    // Initialize MAPQ histogram in shared memory
    for (int i = local_tid; i <= CUDA_MAX_QUALITY_SCORE; i += blockDim.x) {
        shared_mapq_hist[i] = 0;
    }
    __syncthreads();
    
    // Process records
    for (int i = tid; i < num_records; i += stride) {
        const cuda_bam_record_t *record = &records[i];
        
        if (!(record->flag & 0x4)) {  // Only mapped reads
            shared_quality_sum[local_tid] += record->mapq;
            
            if (record->mapq <= CUDA_MAX_QUALITY_SCORE) {
                atomicAdd(&shared_mapq_hist[record->mapq], 1);
            }
        }
    }
    
    __syncthreads();
    
    // Reduce quality sum
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (local_tid < s) {
            shared_quality_sum[local_tid] += shared_quality_sum[local_tid + s];
        }
        __syncthreads();
    }
    
    // Update global statistics
    if (local_tid == 0) {
        atomicAdd(&stats->total_mapq, shared_quality_sum[0]);
    }
    
    // Update global MAPQ histogram
    for (int i = local_tid; i <= CUDA_MAX_QUALITY_SCORE; i += blockDim.x) {
        atomicAdd(&stats->mapq_histogram[i], shared_mapq_hist[i]);
    }
}

__global__ void cuda_gc_content_kernel(const cuda_bam_record_t *records, int num_records,
                                      cuda_gc_content_stats_t *stats) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    
    __shared__ uint32_t shared_gc_bins[CUDA_GC_BINS];
    
    // Initialize shared memory
    for (int i = threadIdx.x; i < CUDA_GC_BINS; i += blockDim.x) {
        shared_gc_bins[i] = 0;
    }
    __syncthreads();
    
    // Process records
    for (int i = tid; i < num_records; i += stride) {
        const cuda_bam_record_t *record = &records[i];
        
        if (!(record->flag & 0x4) && record->seq_len > 0) {  // Only mapped reads with sequence
            int gc_count = 0;
            
            // Count G and C bases
            for (int j = 0; j < record->seq_len; j++) {
                uint8_t base = (record->seq[j >> 1] >> ((~j & 1) << 2)) & 0xf;
                if (base == 2 || base == 4) {  // G or C in 4-bit encoding
                    gc_count++;
                }
            }
            
            // Calculate GC percentage and update histogram
            int gc_percent = (gc_count * 100) / record->seq_len;
            if (gc_percent < CUDA_GC_BINS) {
                atomicAdd(&shared_gc_bins[gc_percent], 1);
            }
        }
    }
    
    __syncthreads();
    
    // Update global GC histogram
    for (int i = threadIdx.x; i < CUDA_GC_BINS; i += blockDim.x) {
        atomicAdd(&stats->gc_histogram[i], shared_gc_bins[i]);
    }
}

__global__ void cuda_insert_size_kernel(const cuda_bam_record_t *records, int num_records,
                                       cuda_insert_size_stats_t *stats) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    
    __shared__ uint64_t shared_total_isize[CUDA_STATS_BLOCK_SIZE];
    __shared__ uint32_t shared_pair_count[CUDA_STATS_BLOCK_SIZE];
    
    int local_tid = threadIdx.x;
    shared_total_isize[local_tid] = 0;
    shared_pair_count[local_tid] = 0;
    
    // Process records
    for (int i = tid; i < num_records; i += stride) {
        const cuda_bam_record_t *record = &records[i];
        
        // Only process properly paired reads with positive insert size
        if ((record->flag & 0x2) && !(record->flag & 0x4) && record->isize > 0) {
            shared_total_isize[local_tid] += abs(record->isize);
            shared_pair_count[local_tid]++;
            
            // Update min/max atomically (simplified)
            atomicMin(&stats->min_insert_size, abs(record->isize));
            atomicMax(&stats->max_insert_size, abs(record->isize));
        }
    }
    
    __syncthreads();
    
    // Block-level reduction
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (local_tid < s) {
            shared_total_isize[local_tid] += shared_total_isize[local_tid + s];
            shared_pair_count[local_tid] += shared_pair_count[local_tid + s];
        }
        __syncthreads();
    }
    
    // Update global statistics
    if (local_tid == 0) {
        atomicAdd(&stats->total_insert_size, shared_total_isize[0]);
        atomicAdd(&stats->pair_count, shared_pair_count[0]);
    }
}

// Host functions

extern "C" int cuda_compute_basic_stats(cuda_context_t *ctx, const cuda_bam_record_t *records,
                                       int num_records, cuda_basic_stats_t *stats) {
    if (!ctx || !records || !stats || num_records <= 0) {
        return -1;
    }
    
    // Initialize stats on device
    CUDA_CHECK(cudaMemset(stats, 0, sizeof(cuda_basic_stats_t)));
    
    // Calculate optimal launch parameters
    int block_size = CUDA_STATS_BLOCK_SIZE;
    int grid_size = cuda_get_optimal_grid_size(num_records, block_size);
    
    // Launch kernel
    cudaStream_t stream = cuda_get_stream(ctx, 0);
    cuda_basic_stats_kernel<<<grid_size, block_size, 0, stream>>>(records, num_records, stats);
    
    CUDA_CHECK_KERNEL();
    CUDA_CHECK(cudaStreamSynchronize(stream));
    
    return 0;
}

extern "C" int cuda_compute_quality_stats(cuda_context_t *ctx, const cuda_bam_record_t *records,
                                         int num_records, cuda_quality_stats_t *stats) {
    if (!ctx || !records || !stats || num_records <= 0) {
        return -1;
    }
    
    // Initialize stats on device
    CUDA_CHECK(cudaMemset(stats, 0, sizeof(cuda_quality_stats_t)));
    
    // Calculate optimal launch parameters
    int block_size = CUDA_STATS_BLOCK_SIZE;
    int grid_size = cuda_get_optimal_grid_size(num_records, block_size);
    
    // Launch kernel
    cudaStream_t stream = cuda_get_stream(ctx, 1);
    cuda_quality_stats_kernel<<<grid_size, block_size, 0, stream>>>(records, num_records, stats);
    
    CUDA_CHECK_KERNEL();
    CUDA_CHECK(cudaStreamSynchronize(stream));
    
    return 0;
}

extern "C" int cuda_compute_gc_content(cuda_context_t *ctx, const cuda_bam_record_t *records,
                                      int num_records, cuda_gc_content_stats_t *stats) {
    if (!ctx || !records || !stats || num_records <= 0) {
        return -1;
    }
    
    // Initialize stats on device
    CUDA_CHECK(cudaMemset(stats, 0, sizeof(cuda_gc_content_stats_t)));
    
    // Calculate optimal launch parameters
    int block_size = CUDA_STATS_BLOCK_SIZE;
    int grid_size = cuda_get_optimal_grid_size(num_records, block_size);
    
    // Launch kernel
    cudaStream_t stream = cuda_get_stream(ctx, 2);
    cuda_gc_content_kernel<<<grid_size, block_size, 0, stream>>>(records, num_records, stats);
    
    CUDA_CHECK_KERNEL();
    CUDA_CHECK(cudaStreamSynchronize(stream));
    
    return 0;
}

extern "C" int cuda_compute_insert_size_stats(cuda_context_t *ctx, const cuda_bam_record_t *records,
                                             int num_records, cuda_insert_size_stats_t *stats) {
    if (!ctx || !records || !stats || num_records <= 0) {
        return -1;
    }
    
    // Initialize stats on device (with proper initial values for min/max)
    CUDA_CHECK(cudaMemset(stats, 0, sizeof(cuda_insert_size_stats_t)));
    
    // Set initial min/max values
    uint32_t initial_min = UINT32_MAX;
    uint32_t initial_max = 0;
    CUDA_CHECK(cudaMemcpy(&stats->min_insert_size, &initial_min, sizeof(uint32_t), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(&stats->max_insert_size, &initial_max, sizeof(uint32_t), cudaMemcpyHostToDevice));
    
    // Calculate optimal launch parameters
    int block_size = CUDA_STATS_BLOCK_SIZE;
    int grid_size = cuda_get_optimal_grid_size(num_records, block_size);
    
    // Launch kernel
    cudaStream_t stream = cuda_get_stream(ctx, 3);
    cuda_insert_size_kernel<<<grid_size, block_size, 0, stream>>>(records, num_records, stats);
    
    CUDA_CHECK_KERNEL();
    CUDA_CHECK(cudaStreamSynchronize(stream));
    
    return 0;
}

extern "C" int cuda_compute_comprehensive_stats(cuda_context_t *ctx, const cuda_bam_record_t *records,
                                               int num_records, cuda_comprehensive_stats_t *comp_stats) {
    if (!ctx || !records || !comp_stats || num_records <= 0) {
        return -1;
    }
    
    // Compute all statistics components in parallel using different streams
    int result = 0;
    
    result |= cuda_compute_basic_stats(ctx, records, num_records, &comp_stats->basic);
    result |= cuda_compute_quality_stats(ctx, records, num_records, &comp_stats->quality);
    result |= cuda_compute_gc_content(ctx, records, num_records, &comp_stats->gc_content);
    result |= cuda_compute_insert_size_stats(ctx, records, num_records, &comp_stats->insert_size);
    
    // Synchronize all streams
    cuda_synchronize_all_streams(ctx);
    
    return result;
}

extern "C" int cuda_copy_stats_to_host(const cuda_comprehensive_stats_t *device_stats,
                                      cuda_comprehensive_stats_t *host_stats) {
    if (!device_stats || !host_stats) {
        return -1;
    }
    
    CUDA_CHECK(cudaMemcpy(host_stats, device_stats, sizeof(cuda_comprehensive_stats_t), cudaMemcpyDeviceToHost));
    return 0;
}

extern "C" void cuda_print_stats_summary(const cuda_comprehensive_stats_t *stats) {
    if (!stats) return;
    
    printf("\n=== CUDA GPU Statistics Summary ===\n");
    printf("Total reads: %lu\n", stats->basic.mapped_reads + stats->basic.unmapped_reads);
    printf("Mapped reads: %lu (%.2f%%)\n", stats->basic.mapped_reads,
           (double)stats->basic.mapped_reads / (stats->basic.mapped_reads + stats->basic.unmapped_reads) * 100.0);
    printf("Unmapped reads: %lu (%.2f%%)\n", stats->basic.unmapped_reads,
           (double)stats->basic.unmapped_reads / (stats->basic.mapped_reads + stats->basic.unmapped_reads) * 100.0);
    printf("Primary mapped: %lu\n", stats->basic.primary_mapped);
    printf("Supplementary: %lu\n", stats->basic.supplementary);
    printf("Total sequence length: %lu\n", stats->basic.total_length);
    
    if (stats->basic.mapped_reads > 0) {
        printf("Average mapping quality: %.2f\n", (double)stats->quality.total_mapq / stats->basic.mapped_reads);
    }
    
    if (stats->insert_size.pair_count > 0) {
        printf("Average insert size: %.2f\n", (double)stats->insert_size.total_insert_size / stats->insert_size.pair_count);
        printf("Insert size range: %u - %u\n", stats->insert_size.min_insert_size, stats->insert_size.max_insert_size);
    }
    
    printf("===================================\n\n");
}
