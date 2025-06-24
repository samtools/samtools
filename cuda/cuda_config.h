/*  cuda_config.h -- CUDA GPU acceleration configuration and utilities.

    Copyright (C) 2025 Akshay Dedaniya <Dedaniya08@hotmail.com>

    This module provides centralized CUDA configuration, device management,
    and memory allocation utilities for GPU-accelerated samtools operations.
    
    Key features:
    - Device detection and optimal device selection
    - Memory management with alignment optimization
    - Stream management for concurrent operations
    - Error handling macros with detailed diagnostics

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

// =============================================================================
// CUDA Hardware and Runtime Constants
// =============================================================================

// Maximum CUDA hardware limitations (safely conservative values)
#define CUDA_MAX_BLOCKS 65536
#define CUDA_MAX_THREADS_PER_BLOCK 1024
#define CUDA_WARP_SIZE 32
#define CUDA_MAX_SHARED_MEMORY 49152  // 48KB shared memory per block

// Optimal block sizes for different samtools operations (empirically tuned)
#define CUDA_SORT_BLOCK_SIZE 256      // Optimal for memory coalescing in sort
#define CUDA_PILEUP_BLOCK_SIZE 512    // Balanced computation/memory for pileup
#define CUDA_STATS_BLOCK_SIZE 256     // Good occupancy for statistics
#define CUDA_HASH_BLOCK_SIZE 256      // Efficient for hash operations

// Memory alignment and management constants
#define CUDA_MEMORY_ALIGNMENT 256     // Align to 256-byte boundaries for optimal access
#define CUDA_DEFAULT_STREAM_COUNT 4   // Number of concurrent streams for overlap

// =============================================================================
// Error Handling Macros
// =============================================================================

/**
 * Check CUDA runtime API calls and exit on error with diagnostic information.
 * Use this for CUDA runtime functions like cudaMalloc, cudaMemcpy, etc.
 */
#define CUDA_CHECK(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, \
                    cudaGetErrorString(error)); \
            exit(EXIT_FAILURE); \
        } \
    } while(0)

/**
 * Check for kernel launch errors and synchronize device.
 * Use this after kernel launches to catch both launch and execution errors.
 */
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

// =============================================================================
// CUDA Data Structures
// =============================================================================

/**
 * Comprehensive GPU device information structure.
 * Contains all relevant hardware properties for optimization decisions.
 */
typedef struct {
    int device_count;                  // Total number of CUDA devices
    int current_device;                // Currently selected device ID
    size_t global_memory;              // Total global memory in bytes
    size_t shared_memory_per_block;    // Shared memory per block in bytes
    int multiprocessor_count;          // Number of streaming multiprocessors
    int max_threads_per_block;         // Maximum threads per block
    int max_blocks_per_grid;           // Maximum blocks per grid dimension
    int warp_size;                     // Warp size (typically 32)
    int compute_capability_major;      // Compute capability major version
    int compute_capability_minor;      // Compute capability minor version
    bool unified_memory_support;       // Support for unified memory
    bool async_engine_count;           // Number of asynchronous engines
} cuda_device_info_t;

/**
 * CUDA execution context for managing GPU resources.
 * Encapsulates device state, streams, and memory tracking.
 */
typedef struct {
    cuda_device_info_t device_info;    // Device hardware information
    cudaStream_t *streams;             // Array of CUDA streams
    int stream_count;                  // Number of allocated streams
    size_t total_memory_allocated;     // Total GPU memory allocated (bytes)
    bool initialized;                  // Context initialization status
} cuda_context_t;

// =============================================================================
// CUDA Context Management Functions
// =============================================================================

/**
 * Initialize CUDA context with specified device.
 * @param ctx Pointer to context structure to initialize
 * @param device_id GPU device ID to use (-1 for auto-select)
 * @return 0 on success, negative on error
 */
int cuda_init_context(cuda_context_t *ctx, int device_id);

/**
 * Clean up and destroy CUDA context, freeing all resources.
 * @param ctx Pointer to context to cleanup
 */
void cuda_cleanup_context(cuda_context_t *ctx);

/**
 * Check if CUDA is available and functional.
 * @return true if CUDA is available, false otherwise
 */
bool cuda_is_available(void);

/**
 * Automatically select the best available GPU device.
 * @return device ID of best device, or -1 if none available
 */
int cuda_select_best_device(void);

/**
 * Print detailed device information for diagnostics.
 * @param info Pointer to device info structure
 */
void cuda_print_device_info(const cuda_device_info_t *info);

// =============================================================================
// CUDA Execution Configuration Functions
// =============================================================================

/**
 * Calculate optimal block size for given data size and constraints.
 * @param data_size Number of elements to process
 * @param max_block_size Maximum allowed block size
 * @return recommended block size
 */
int cuda_get_optimal_block_size(int data_size, int max_block_size);

/**
 * Calculate optimal grid size for given data and block size.
 * @param data_size Number of elements to process
 * @param block_size Threads per block
 * @return recommended grid size (number of blocks)
 */
int cuda_get_optimal_grid_size(int data_size, int block_size);

// =============================================================================
// CUDA Memory Management Functions
// =============================================================================

/**
 * Allocate unified memory accessible from both host and device.
 * @param size Bytes to allocate
 * @return pointer to allocated memory, or NULL on failure
 */
void* cuda_malloc_managed(size_t size);

/**
 * Free unified memory allocated with cuda_malloc_managed.
 * @param ptr Pointer to memory to free
 */
void cuda_free_managed(void *ptr);

/**
 * Allocate pinned host memory for fast transfers.
 * @param size Bytes to allocate
 * @return pointer to allocated memory, or NULL on failure
 */
void* cuda_malloc_host(size_t size);

/**
 * Free pinned host memory allocated with cuda_malloc_host.
 * @param ptr Pointer to memory to free
 */
void cuda_free_host(void *ptr);

/**
 * Allocate device memory on GPU.
 * @param size Bytes to allocate
 * @return pointer to allocated memory, or NULL on failure
 */
void* cuda_malloc_device(size_t size);

/**
 * Free device memory allocated with cuda_malloc_device.
 * @param ptr Pointer to memory to free
 */
void cuda_free_device(void *ptr);

// =============================================================================
// CUDA Stream Management Functions
// =============================================================================

/**
 * Get a specific stream from the context.
 * @param ctx Pointer to initialized context
 * @param stream_id Stream index (0 to stream_count-1)
 * @return CUDA stream handle
 */
cudaStream_t cuda_get_stream(cuda_context_t *ctx, int stream_id);

/**
 * Synchronize a specific stream (wait for completion).
 * @param stream CUDA stream to synchronize
 */
void cuda_synchronize_stream(cudaStream_t stream);

/**
 * Synchronize all streams in the context.
 * @param ctx Pointer to context containing streams
 */
void cuda_synchronize_all_streams(cuda_context_t *ctx);

#ifdef __cplusplus
}
#endif

#endif // CUDA_CONFIG_H
