/*  cuda_pileup.cu -- CUDA GPU-accelerated pileup processing implementation.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <cooperative_groups.h>
#include <cub/cub.cuh>
#include "cuda_pileup.cuh"
#include "cuda_config.h"

namespace cg = cooperative_groups;

// CUDA kernels for pileup processing

__device__ uint8_t cuda_encode_base_device(char base) {
    switch (base) {
        case 'A': case 'a': return CUDA_BASE_A;
        case 'C': case 'c': return CUDA_BASE_C;
        case 'G': case 'g': return CUDA_BASE_G;
        case 'T': case 't': return CUDA_BASE_T;
        case 'N': case 'n': return CUDA_BASE_N;
        default: return CUDA_BASE_N;
    }
}

__device__ char cuda_decode_base_device(uint8_t encoded_base) {
    switch (encoded_base) {
        case CUDA_BASE_A: return 'A';
        case CUDA_BASE_C: return 'C';
        case CUDA_BASE_G: return 'G';
        case CUDA_BASE_T: return 'T';
        case CUDA_BASE_N: return 'N';
        default: return 'N';
    }
}

// Kernel to process CIGAR operations and extract bases at specific positions
__global__ void cuda_pileup_extract_bases_kernel(
    const cuda_bam_record_t *records,
    size_t num_records,
    int32_t target_pos,
    uint32_t *base_counts,
    uint32_t *qual_sum,
    uint16_t *depth) {
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Shared memory for block-level reduction
    __shared__ uint32_t s_base_counts[CUDA_BASE_COUNT * CUDA_PILEUP_BLOCK_SIZE / 32];
    __shared__ uint32_t s_qual_sum[CUDA_PILEUP_BLOCK_SIZE / 32];
    __shared__ uint16_t s_depth[CUDA_PILEUP_BLOCK_SIZE / 32];
    
    // Initialize shared memory
    int warp_id = threadIdx.x / 32;
    int lane_id = threadIdx.x % 32;
    
    if (lane_id == 0) {
        for (int i = 0; i < CUDA_BASE_COUNT; i++) {
            s_base_counts[warp_id * CUDA_BASE_COUNT + i] = 0;
        }
        s_qual_sum[warp_id] = 0;
        s_depth[warp_id] = 0;
    }
    __syncthreads();
    
    if (tid >= num_records) return;
    
    const cuda_bam_record_t *record = &records[tid];
    
    // Check if this record overlaps the target position
    if (record->pos > target_pos || record->end_pos <= target_pos) {
        return;
    }
    
    // Parse CIGAR to find the base at target_pos
    int32_t ref_pos = record->pos;
    int32_t seq_pos = 0;
    bool found_base = false;
    uint8_t base = CUDA_BASE_N;
    uint8_t qual = 0;
    
    for (int i = 0; i < record->cigar_len && !found_base; i++) {
        uint32_t cigar_op = record->cigar_ops[i];
        uint32_t op_len = cigar_op >> 4;
        uint32_t op_type = cigar_op & 0xF;
        
        switch (op_type) {
            case 0: // M - alignment match
            case 7: // = - sequence match
            case 8: // X - sequence mismatch
                if (ref_pos <= target_pos && target_pos < ref_pos + op_len) {
                    int offset = target_pos - ref_pos;
                    if (seq_pos + offset < record->seq_len) {
                        base = record->seq[seq_pos + offset];
                        qual = record->qual[seq_pos + offset];
                        found_base = true;
                    }
                }
                ref_pos += op_len;
                seq_pos += op_len;
                break;
                
            case 1: // I - insertion
                seq_pos += op_len;
                break;
                
            case 2: // D - deletion
            case 3: // N - skipped region
                if (ref_pos <= target_pos && target_pos < ref_pos + op_len) {
                    base = CUDA_BASE_GAP;
                    qual = 0;
                    found_base = true;
                }
                ref_pos += op_len;
                break;
                
            case 4: // S - soft clipping
            case 5: // H - hard clipping
                seq_pos += op_len;
                break;
        }
    }
    
    if (found_base && base < CUDA_BASE_COUNT) {
        // Use warp-level primitives for efficient reduction
        auto warp = cg::tiled_partition<32>(cg::this_thread_block());
        
        // Count bases
        for (int i = 0; i < CUDA_BASE_COUNT; i++) {
            uint32_t count = (base == i) ? 1 : 0;
            uint32_t warp_sum = cg::reduce(warp, count, cg::plus<uint32_t>());
            if (warp.thread_rank() == 0) {
                atomicAdd(&s_base_counts[warp_id * CUDA_BASE_COUNT + i], warp_sum);
            }
        }
        
        // Sum qualities
        uint32_t warp_qual_sum = cg::reduce(warp, (uint32_t)qual, cg::plus<uint32_t>());
        if (warp.thread_rank() == 0) {
            atomicAdd(&s_qual_sum[warp_id], warp_qual_sum);
            atomicAdd(&s_depth[warp_id], 1);
        }
    }
    
    __syncthreads();
    
    // Block-level reduction to global memory
    if (threadIdx.x < CUDA_BASE_COUNT) {
        uint32_t total = 0;
        for (int i = 0; i < blockDim.x / 32; i++) {
            total += s_base_counts[i * CUDA_BASE_COUNT + threadIdx.x];
        }
        atomicAdd(&base_counts[threadIdx.x], total);
    }
    
    if (threadIdx.x == 0) {
        uint32_t total_qual = 0;
        uint16_t total_depth = 0;
        for (int i = 0; i < blockDim.x / 32; i++) {
            total_qual += s_qual_sum[i];
            total_depth += s_depth[i];
        }
        atomicAdd(qual_sum, total_qual);
        atomicAdd(depth, total_depth);
    }
}

// Kernel for processing multiple positions in parallel
__global__ void cuda_pileup_process_region_kernel(
    const cuda_bam_record_t *records,
    size_t num_records,
    int32_t start_pos,
    int32_t end_pos,
    cuda_pileup_position_t *pileup_data) {
    
    int pos_idx = blockIdx.x;
    int32_t current_pos = start_pos + pos_idx;
    
    if (current_pos >= end_pos) return;
    
    cuda_pileup_position_t *pos_data = &pileup_data[pos_idx];
    pos_data->pos = current_pos;
    pos_data->tid = records[0].tid; // Assume same chromosome
    
    // Initialize counters
    for (int i = 0; i < CUDA_BASE_COUNT; i++) {
        pos_data->base_counts[i] = 0;
    }
    pos_data->qual_sum = 0;
    pos_data->depth = 0;
    pos_data->mapping_qual_sum = 0;
    pos_data->strand_counts[0] = 0;
    pos_data->strand_counts[1] = 0;
    
    // Process all records for this position
    for (int tid = threadIdx.x; tid < num_records; tid += blockDim.x) {
        const cuda_bam_record_t *record = &records[tid];
        
        // Check if record overlaps current position
        if (record->pos > current_pos || record->end_pos <= current_pos) {
            continue;
        }
        
        // Extract base at this position (simplified CIGAR parsing)
        int32_t ref_pos = record->pos;
        int32_t seq_pos = 0;
        bool found = false;
        uint8_t base = CUDA_BASE_N;
        uint8_t qual = 0;
        
        // Simplified: assume match operations only for demo
        if (current_pos >= record->pos && current_pos < record->end_pos) {
            int offset = current_pos - record->pos;
            if (offset < record->seq_len) {
                base = record->seq[offset];
                qual = record->qual[offset];
                found = true;
            }
        }
        
        if (found && base < CUDA_BASE_COUNT) {
            atomicAdd(&pos_data->base_counts[base], 1);
            atomicAdd(&pos_data->qual_sum, qual);
            atomicAdd(&pos_data->depth, 1);
            atomicAdd(&pos_data->mapping_qual_sum, record->mapq);
            
            // Count strand
            bool reverse = (record->flag & 0x10) != 0;
            atomicAdd(&pos_data->strand_counts[reverse ? 1 : 0], 1);
        }
    }
}

// Kernel for consensus calling
__global__ void cuda_consensus_calling_kernel(
    const cuda_pileup_position_t *pileup_data,
    size_t pileup_count,
    cuda_consensus_result_t *consensus_results,
    uint16_t min_depth,
    float min_consensus_fraction) {
    
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= pileup_count) return;
    
    const cuda_pileup_position_t *pos = &pileup_data[tid];
    cuda_consensus_result_t *result = &consensus_results[tid];
    
    result->pos = pos->pos;
    result->depth = pos->depth;
    
    if (pos->depth < min_depth) {
        result->consensus_base = 'N';
        result->quality = 0;
        result->confidence = 0.0f;
        return;
    }
    
    // Find most frequent base
    uint32_t max_count = 0;
    uint8_t consensus_base_idx = CUDA_BASE_N;
    
    for (int i = 0; i < 4; i++) { // Only consider A, C, G, T
        if (pos->base_counts[i] > max_count) {
            max_count = pos->base_counts[i];
            consensus_base_idx = i;
        }
    }
    
    // Calculate consensus metrics
    float consensus_fraction = (float)max_count / pos->depth;
    
    if (consensus_fraction >= min_consensus_fraction) {
        result->consensus_base = cuda_decode_base_device(consensus_base_idx);
        result->confidence = consensus_fraction;
        
        // Simplified quality calculation
        result->quality = (uint8_t)min(60.0f, -10.0f * log10f(1.0f - consensus_fraction + 1e-6f));
    } else {
        result->consensus_base = 'N';
        result->quality = 0;
        result->confidence = 0.0f;
    }
}

// Context management functions
int cuda_pileup_init_context(cuda_pileup_context_t *ctx, 
                             size_t max_records, 
                             size_t max_positions) {
    memset(ctx, 0, sizeof(cuda_pileup_context_t));
    
    ctx->max_records = max_records;
    ctx->max_positions = max_positions;
    
    // Allocate device memory
    CUDA_CHECK(cudaMalloc(&ctx->d_records, max_records * sizeof(cuda_bam_record_t)));
    CUDA_CHECK(cudaMalloc(&ctx->d_pileup, max_positions * sizeof(cuda_pileup_position_t)));
    CUDA_CHECK(cudaMalloc(&ctx->d_consensus, max_positions * sizeof(cuda_consensus_result_t)));
    
    // Create CUDA stream
    CUDA_CHECK(cudaStreamCreate(&ctx->stream));
    
    // Set default parameters
    ctx->min_base_quality = 20;
    ctx->min_mapping_quality = 20;
    ctx->min_depth = 3;
    ctx->min_consensus_fraction = 0.75f;
    
    return 0;
}

void cuda_pileup_cleanup_context(cuda_pileup_context_t *ctx) {
    if (ctx->d_records) {
        cudaFree(ctx->d_records);
        ctx->d_records = nullptr;
    }
    
    if (ctx->d_pileup) {
        cudaFree(ctx->d_pileup);
        ctx->d_pileup = nullptr;
    }
    
    if (ctx->d_consensus) {
        cudaFree(ctx->d_consensus);
        ctx->d_consensus = nullptr;
    }
    
    if (ctx->stream) {
        cudaStreamDestroy(ctx->stream);
        ctx->stream = 0;
    }
    
    memset(ctx, 0, sizeof(cuda_pileup_context_t));
}

// Main pileup processing function
int cuda_pileup_process_region(cuda_pileup_context_t *ctx,
                              const cuda_bam_record_t *records,
                              size_t num_records,
                              int32_t start_pos,
                              int32_t end_pos,
                              cuda_pileup_position_t *output_pileup,
                              size_t *output_count) {
    
    if (num_records > ctx->max_records) {
        fprintf(stderr, "Too many records: %zu > %zu\n", num_records, ctx->max_records);
        return -1;
    }
    
    size_t num_positions = end_pos - start_pos;
    if (num_positions > ctx->max_positions) {
        fprintf(stderr, "Too many positions: %zu > %zu\n", num_positions, ctx->max_positions);
        return -1;
    }
    
    // Copy records to device
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_records, records, 
                               num_records * sizeof(cuda_bam_record_t),
                               cudaMemcpyHostToDevice, ctx->stream));
    
    // Process region
    int block_size = CUDA_PILEUP_BLOCK_SIZE;
    int grid_size = num_positions;
    
    cuda_pileup_process_region_kernel<<<grid_size, block_size, 0, ctx->stream>>>(
        ctx->d_records, num_records, start_pos, end_pos, ctx->d_pileup);
    
    CUDA_CHECK_KERNEL();
    
    // Copy results back to host
    CUDA_CHECK(cudaMemcpyAsync(output_pileup, ctx->d_pileup,
                               num_positions * sizeof(cuda_pileup_position_t),
                               cudaMemcpyDeviceToHost, ctx->stream));
    
    CUDA_CHECK(cudaStreamSynchronize(ctx->stream));
    
    *output_count = num_positions;
    return 0;
}

// Consensus generation function
int cuda_pileup_generate_consensus(cuda_pileup_context_t *ctx,
                                  const cuda_pileup_position_t *pileup_data,
                                  size_t pileup_count,
                                  cuda_consensus_result_t *consensus_results,
                                  size_t *consensus_count) {
    
    if (pileup_count > ctx->max_positions) {
        fprintf(stderr, "Too many positions for consensus: %zu > %zu\n", 
                pileup_count, ctx->max_positions);
        return -1;
    }
    
    // Copy pileup data to device
    CUDA_CHECK(cudaMemcpyAsync(ctx->d_pileup, pileup_data,
                               pileup_count * sizeof(cuda_pileup_position_t),
                               cudaMemcpyHostToDevice, ctx->stream));
    
    // Launch consensus calling kernel
    int block_size = 256;
    int grid_size = (pileup_count + block_size - 1) / block_size;
    
    cuda_consensus_calling_kernel<<<grid_size, block_size, 0, ctx->stream>>>(
        ctx->d_pileup, pileup_count, ctx->d_consensus,
        ctx->min_depth, ctx->min_consensus_fraction);
    
    CUDA_CHECK_KERNEL();
    
    // Copy results back
    CUDA_CHECK(cudaMemcpyAsync(consensus_results, ctx->d_consensus,
                               pileup_count * sizeof(cuda_consensus_result_t),
                               cudaMemcpyDeviceToHost, ctx->stream));
    
    CUDA_CHECK(cudaStreamSynchronize(ctx->stream));
    
    *consensus_count = pileup_count;
    return 0;
}

// Utility functions
uint8_t cuda_encode_base(char base) {
    return cuda_encode_base_device(base);
}

char cuda_decode_base(uint8_t encoded_base) {
    return cuda_decode_base_device(encoded_base);
}

int cuda_encode_sequence_to_gpu(const char *sequence, 
                               size_t length, 
                               uint8_t *gpu_sequence) {
    for (size_t i = 0; i < length; i++) {
        gpu_sequence[i] = cuda_encode_base(sequence[i]);
    }
    return 0;
}

int cuda_decode_sequence_from_gpu(const uint8_t *gpu_sequence, 
                                 size_t length, 
                                 char *sequence) {
    for (size_t i = 0; i < length; i++) {
        sequence[i] = cuda_decode_base(gpu_sequence[i]);
    }
    sequence[length] = '\0';
    return 0;
}

bool cuda_pileup_validate_input(const cuda_bam_record_t *records, size_t count) {
    for (size_t i = 0; i < count; i++) {
        const cuda_bam_record_t *record = &records[i];
        
        if (record->pos < 0 || record->end_pos <= record->pos) {
            fprintf(stderr, "Invalid position range for record %zu: %d-%d\n", 
                    i, record->pos, record->end_pos);
            return false;
        }
        
        if (record->seq_len == 0 || !record->seq || !record->qual) {
            fprintf(stderr, "Invalid sequence data for record %zu\n", i);
            return false;
        }
    }
    
    return true;
}

int cuda_pileup_check_memory_requirements(size_t num_records, size_t num_positions) {
    size_t total_memory = 0;
    
    total_memory += num_records * sizeof(cuda_bam_record_t);
    total_memory += num_positions * sizeof(cuda_pileup_position_t);
    total_memory += num_positions * sizeof(cuda_consensus_result_t);
    
    size_t free_mem, total_mem;
    CUDA_CHECK(cudaMemGetInfo(&free_mem, &total_mem));
    
    if (total_memory > free_mem * 0.8) { // Use only 80% of available memory
        fprintf(stderr, "Insufficient GPU memory: need %zu MB, have %zu MB available\n",
                total_memory / (1024 * 1024), free_mem / (1024 * 1024));
        return -1;
    }
    
    return 0;
}
