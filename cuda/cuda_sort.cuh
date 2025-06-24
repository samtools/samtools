/*  cuda_sort.cuh -- CUDA GPU-accelerated sorting algorithms for BAM data.

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

#ifndef CUDA_SORT_CUH
#define CUDA_SORT_CUH

#include <cuda_runtime.h>
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/execution_policy.h>
#include <stdint.h>
#include <stdbool.h>
#include "cuda_config.h"

#ifdef __cplusplus
extern "C" {
#endif

// GPU-optimized BAM record key structure for sorting
typedef struct {
    int32_t tid;        // Target ID (chromosome)
    int64_t pos;        // Position
    uint64_t hash;      // Minimizer hash for advanced sorting
    uint32_t idx;       // Original index
    bool reverse;       // Strand orientation
} cuda_bam_sort_key_t;

// Sorting mode enumeration
typedef enum {
    CUDA_SORT_COORDINATE,
    CUDA_SORT_QUERY_NAME,
    CUDA_SORT_MINIMIZER_HASH,
    CUDA_SORT_TEMPLATE_COORDINATE
} cuda_sort_mode_t;

// CUDA sorting context
typedef struct {
    cuda_bam_sort_key_t *d_keys;     // Device keys array
    uint32_t *d_indices;             // Device indices array
    size_t allocated_size;           // Currently allocated size
    cuda_sort_mode_t sort_mode;      // Current sorting mode
    cudaStream_t stream;             // CUDA stream for async operations
} cuda_sort_context_t;

// Function declarations
int cuda_sort_init_context(cuda_sort_context_t *ctx, size_t initial_size, cuda_sort_mode_t mode);
void cuda_sort_cleanup_context(cuda_sort_context_t *ctx);

// Main sorting functions
int cuda_sort_bam_records(cuda_sort_context_t *ctx, 
                         cuda_bam_sort_key_t *keys, 
                         uint32_t *indices, 
                         size_t count);

int cuda_sort_coordinate_parallel(cuda_bam_sort_key_t *keys, 
                                 uint32_t *indices, 
                                 size_t count, 
                                 cudaStream_t stream);

int cuda_sort_minimizer_hash_parallel(cuda_bam_sort_key_t *keys, 
                                     uint32_t *indices, 
                                     size_t count, 
                                     cudaStream_t stream);

// Utility functions
int cuda_generate_sort_keys(const void *bam_records, 
                           cuda_bam_sort_key_t *keys, 
                           size_t count, 
                           cuda_sort_mode_t mode);

int cuda_apply_sort_permutation(void *bam_records, 
                               const uint32_t *indices, 
                               size_t count, 
                               size_t record_size);

// Merge operations for large datasets
int cuda_merge_sorted_blocks(cuda_sort_context_t *ctx,
                           cuda_bam_sort_key_t **key_blocks,
                           uint32_t **index_blocks,
                           size_t *block_sizes,
                           int num_blocks,
                           cuda_bam_sort_key_t *output_keys,
                           uint32_t *output_indices);

// Memory management helpers
int cuda_sort_resize_buffers(cuda_sort_context_t *ctx, size_t new_size);
void cuda_sort_prefetch_data(void *data, size_t size, int device);

// Performance monitoring
typedef struct {
    float sort_time_ms;
    float memory_transfer_time_ms;
    size_t peak_memory_usage;
    int gpu_utilization_percent;
} cuda_sort_performance_t;

int cuda_sort_get_performance_stats(cuda_sort_context_t *ctx, cuda_sort_performance_t *stats);

#ifdef __cplusplus
}
#endif

#endif // CUDA_SORT_CUH
