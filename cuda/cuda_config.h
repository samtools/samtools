/*  cuda_config.h -- CUDA GPU acceleration configuration and utilities.

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

#ifndef CUDA_CONFIG_H
#define CUDA_CONFIG_H

#ifdef __cplusplus
extern "C" {
#endif

#include <cuda_runtime.h>
#include <stdint.h>
#include <stdbool.h>

// CUDA configuration constants
#define CUDA_MAX_BLOCKS 65536
#define CUDA_MAX_THREADS_PER_BLOCK 1024
#define CUDA_WARP_SIZE 32
#define CUDA_MAX_SHARED_MEMORY 49152  // 48KB shared memory per block

// Optimal block and grid sizes for different operations
#define CUDA_SORT_BLOCK_SIZE 256
#define CUDA_PILEUP_BLOCK_SIZE 512
#define CUDA_STATS_BLOCK_SIZE 256
#define CUDA_HASH_BLOCK_SIZE 256

// Memory management
#define CUDA_MEMORY_ALIGNMENT 256
#define CUDA_DEFAULT_STREAM_COUNT 4

// Error handling macros
#define CUDA_CHECK(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, \
                    cudaGetErrorString(error)); \
            exit(EXIT_FAILURE); \
        } \
    } while(0)

#define CUDA_CHECK_KERNEL() \
    do { \
        cudaError_t error = cudaGetLastError(); \
        if (error != cudaSuccess) { \
            fprintf(stderr, "CUDA kernel error at %s:%d - %s\n", __FILE__, __LINE__, \
                    cudaGetErrorString(error)); \
            exit(EXIT_FAILURE); \
        } \
        CUDA_CHECK(cudaDeviceSynchronize()); \
    } while(0)

// GPU device information structure
typedef struct {
    int device_count;
    int current_device;
    size_t global_memory;
    size_t shared_memory_per_block;
    int multiprocessor_count;
    int max_threads_per_block;
    int max_blocks_per_grid;
    int warp_size;
    int compute_capability_major;
    int compute_capability_minor;
    bool unified_memory_support;
    bool async_engine_count;
} cuda_device_info_t;

// CUDA context structure for managing resources
typedef struct {
    cuda_device_info_t device_info;
    cudaStream_t *streams;
    int stream_count;
    size_t total_memory_allocated;
    bool initialized;
} cuda_context_t;

// Function declarations
int cuda_init_context(cuda_context_t *ctx, int device_id);
void cuda_cleanup_context(cuda_context_t *ctx);
int cuda_get_optimal_block_size(int data_size, int max_block_size);
int cuda_get_optimal_grid_size(int data_size, int block_size);
void cuda_print_device_info(const cuda_device_info_t *info);
bool cuda_is_available(void);
int cuda_select_best_device(void);

// Memory management functions
void* cuda_malloc_managed(size_t size);
void cuda_free_managed(void *ptr);
void* cuda_malloc_host(size_t size);
void cuda_free_host(void *ptr);
void* cuda_malloc_device(size_t size);
void cuda_free_device(void *ptr);

// Stream management
cudaStream_t cuda_get_stream(cuda_context_t *ctx, int stream_id);
void cuda_synchronize_stream(cudaStream_t stream);
void cuda_synchronize_all_streams(cuda_context_t *ctx);

#ifdef __cplusplus
}
#endif

#endif // CUDA_CONFIG_H
