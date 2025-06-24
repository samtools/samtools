/*  cuda_config.cu -- CUDA GPU acceleration configuration and utilities implementation.

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
#include "cuda_config.h"

// Global CUDA context
static cuda_context_t g_cuda_context = {0};

bool cuda_is_available(void) {
    int device_count = 0;
    cudaError_t error = cudaGetDeviceCount(&device_count);
    
    if (error != cudaSuccess) {
        // Common error cases
        if (error == cudaErrorNoDevice) {
            return false; // No CUDA-capable devices
        } else if (error == cudaErrorInsufficientDriver) {
            fprintf(stderr, "[CUDA] Driver version is insufficient for CUDA runtime version\n");
            return false;
        } else {
            fprintf(stderr, "[CUDA] Error checking device count: %s\n", cudaGetErrorString(error));
            return false;
        }
    }
    
    return (device_count > 0);
}

// New function to auto-detect and display system capabilities
int cuda_auto_detect_system(void) {
    printf("[CUDA] Auto-detecting CUDA system capabilities...\n");
    
    // Check CUDA runtime version
    int runtime_version = 0;
    cudaError_t error = cudaRuntimeGetVersion(&runtime_version);
    if (error == cudaSuccess) {
        printf("[CUDA] Runtime version: %d.%d\n", 
               runtime_version / 1000, (runtime_version % 100) / 10);
    }
    
    // Check driver version
    int driver_version = 0;
    error = cudaDriverGetVersion(&driver_version);
    if (error == cudaSuccess) {
        printf("[CUDA] Driver version: %d.%d\n", 
               driver_version / 1000, (driver_version % 100) / 10);
    }
    
    // Check device count
    int device_count = 0;
    error = cudaGetDeviceCount(&device_count);
    if (error != cudaSuccess) {
        printf("[CUDA] Error: %s\n", cudaGetErrorString(error));
        return -1;
    }
    
    if (device_count == 0) {
        printf("[CUDA] No CUDA-capable devices found\n");
        return -1;
    }
    
    printf("[CUDA] Found %d CUDA device(s):\n", device_count);
    
    // Display information for each device
    for (int i = 0; i < device_count; i++) {
        cudaDeviceProp prop;
        error = cudaGetDeviceProperties(&prop, i);
        if (error != cudaSuccess) {
            printf("[CUDA] Error getting properties for device %d: %s\n", 
                   i, cudaGetErrorString(error));
            continue;
        }
        
        printf("[CUDA] Device %d: %s\n", i, prop.name);
        printf("[CUDA]   Compute capability: %d.%d\n", prop.major, prop.minor);
        printf("[CUDA]   Global memory: %.2f GB\n", 
               (double)prop.totalGlobalMem / (1024*1024*1024));
        printf("[CUDA]   Multiprocessors: %d\n", prop.multiProcessorCount);
        printf("[CUDA]   Max threads per block: %d\n", prop.maxThreadsPerBlock);
        printf("[CUDA]   Warp size: %d\n", prop.warpSize);
        printf("[CUDA]   Memory clock rate: %.2f GHz\n", 
               prop.memoryClockRate * 1e-6);
        printf("[CUDA]   Memory bus width: %d bits\n", prop.memoryBusWidth);
        
        // Check if device supports unified memory
        if (prop.managedMemory) {
            printf("[CUDA]   Unified memory: Yes\n");
        } else {
            printf("[CUDA]   Unified memory: No\n");
        }
        
        // Check concurrent kernels
        if (prop.concurrentKernels) {
            printf("[CUDA]   Concurrent kernels: Yes\n");
        } else {
            printf("[CUDA]   Concurrent kernels: No\n");
        }
        
        printf("\n");
    }
    
    return device_count;
}

int cuda_select_best_device(void) {
    int device_count = 0;
    cudaError_t error = cudaGetDeviceCount(&device_count);
    
    if (error != cudaSuccess) {
        fprintf(stderr, "[CUDA] Error getting device count: %s\n", cudaGetErrorString(error));
        return -1;
    }
    
    if (device_count == 0) {
        fprintf(stderr, "[CUDA] No CUDA devices found\n");
        return -1;
    }
    
    int best_device = 0;
    size_t max_memory = 0;
    int max_compute_capability = 0;
    int max_multiprocessors = 0;
    
    printf("[CUDA] Selecting best device from %d available device(s):\n", device_count);
    
    for (int i = 0; i < device_count; i++) {
        cudaDeviceProp prop;
        error = cudaGetDeviceProperties(&prop, i);
        if (error != cudaSuccess) {
            fprintf(stderr, "[CUDA] Error getting properties for device %d: %s\n", 
                   i, cudaGetErrorString(error));
            continue;
        }
        
        int compute_capability = prop.major * 10 + prop.minor;
        
        printf("[CUDA] Device %d: %s (CC %d.%d, %.2f GB, %d SMs)\n", 
               i, prop.name, prop.major, prop.minor,
               (double)prop.totalGlobalMem / (1024*1024*1024),
               prop.multiProcessorCount);
        
        // Selection criteria (in order of priority):
        // 1. Compute capability (must be >= 5.0 for modern features)
        // 2. Number of multiprocessors (more parallelism)
        // 3. Global memory size (larger datasets)
        bool is_better = false;
        
        if (compute_capability < 50) {
            printf("[CUDA]   Skipping device %d (compute capability < 5.0)\n", i);
            continue;
        }
        
        if (compute_capability > max_compute_capability) {
            is_better = true;
        } else if (compute_capability == max_compute_capability) {
            if (prop.multiProcessorCount > max_multiprocessors) {
                is_better = true;
            } else if (prop.multiProcessorCount == max_multiprocessors && 
                       prop.totalGlobalMem > max_memory) {
                is_better = true;
            }
        }
        
        if (is_better) {
            best_device = i;
            max_memory = prop.totalGlobalMem;
            max_compute_capability = compute_capability;
            max_multiprocessors = prop.multiProcessorCount;
        }
    }
    
    printf("[CUDA] Selected device %d as best device\n", best_device);
    return best_device;
}

int cuda_init_context(cuda_context_t *ctx, int device_id) {
    if (ctx->initialized) {
        return 0; // Already initialized
    }
    
    memset(ctx, 0, sizeof(cuda_context_t));
    
    if (!cuda_is_available()) {
        fprintf(stderr, "CUDA is not available on this system\n");
        return -1;
    }
    
    if (device_id < 0) {
        device_id = cuda_select_best_device();
        if (device_id < 0) {
            return -1;
        }
    }
    
    CUDA_CHECK(cudaSetDevice(device_id));
    ctx->device_info.current_device = device_id;
    
    // Get device properties
    cudaDeviceProp prop;
    CUDA_CHECK(cudaGetDeviceProperties(&prop, device_id));
    
    ctx->device_info.global_memory = prop.totalGlobalMem;
    ctx->device_info.shared_memory_per_block = prop.sharedMemPerBlock;
    ctx->device_info.multiprocessor_count = prop.multiProcessorCount;
    ctx->device_info.max_threads_per_block = prop.maxThreadsPerBlock;
    ctx->device_info.max_blocks_per_grid = prop.maxGridSize[0];
    ctx->device_info.warp_size = prop.warpSize;
    ctx->device_info.compute_capability_major = prop.major;
    ctx->device_info.compute_capability_minor = prop.minor;
    ctx->device_info.unified_memory_support = prop.unifiedAddressing;
    ctx->device_info.async_engine_count = prop.asyncEngineCount;
    
    // Create CUDA streams
    ctx->stream_count = CUDA_DEFAULT_STREAM_COUNT;
    ctx->streams = (cudaStream_t*)malloc(ctx->stream_count * sizeof(cudaStream_t));
    if (!ctx->streams) {
        fprintf(stderr, "Failed to allocate memory for CUDA streams\n");
        return -1;
    }
    
    for (int i = 0; i < ctx->stream_count; i++) {
        CUDA_CHECK(cudaStreamCreate(&ctx->streams[i]));
    }
    
    ctx->initialized = true;
    printf("CUDA initialized successfully on device %d\n", device_id);
    return 0;
}

void cuda_cleanup_context(cuda_context_t *ctx) {
    if (!ctx->initialized) {
        return;
    }
    
    // Synchronize all streams before cleanup
    cuda_synchronize_all_streams(ctx);
    
    // Destroy streams
    if (ctx->streams) {
        for (int i = 0; i < ctx->stream_count; i++) {
            cudaStreamDestroy(ctx->streams[i]);
        }
        free(ctx->streams);
        ctx->streams = NULL;
    }
    
    // Reset device
    cudaDeviceReset();
    
    memset(ctx, 0, sizeof(cuda_context_t));
    printf("CUDA context cleaned up successfully\n");
}

int cuda_get_optimal_block_size(int data_size, int max_block_size) {
    int block_size = 32; // Start with minimum warp size
    
    while (block_size < max_block_size && block_size < data_size) {
        block_size *= 2;
    }
    
    return (block_size > max_block_size) ? max_block_size : block_size;
}

int cuda_get_optimal_grid_size(int data_size, int block_size) {
    int grid_size = (data_size + block_size - 1) / block_size;
    return (grid_size > CUDA_MAX_BLOCKS) ? CUDA_MAX_BLOCKS : grid_size;
}

void cuda_print_device_info(const cuda_device_info_t *info) {
    printf("CUDA Device Information:\n");
    printf("  Current Device: %d\n", info->current_device);
    printf("  Global Memory: %.2f GB\n", info->global_memory / (1024.0 * 1024.0 * 1024.0));
    printf("  Shared Memory per Block: %zu KB\n", info->shared_memory_per_block / 1024);
    printf("  Multiprocessor Count: %d\n", info->multiprocessor_count);
    printf("  Max Threads per Block: %d\n", info->max_threads_per_block);
    printf("  Max Blocks per Grid: %d\n", info->max_blocks_per_grid);
    printf("  Warp Size: %d\n", info->warp_size);
    printf("  Compute Capability: %d.%d\n", info->compute_capability_major, info->compute_capability_minor);
    printf("  Unified Memory Support: %s\n", info->unified_memory_support ? "Yes" : "No");
    printf("  Async Engine Count: %s\n", info->async_engine_count ? "Yes" : "No");
}

// Memory management functions
void* cuda_malloc_managed(size_t size) {
    void *ptr;
    CUDA_CHECK(cudaMallocManaged(&ptr, size));
    return ptr;
}

void cuda_free_managed(void *ptr) {
    if (ptr) {
        CUDA_CHECK(cudaFree(ptr));
    }
}

void* cuda_malloc_host(size_t size) {
    void *ptr;
    CUDA_CHECK(cudaMallocHost(&ptr, size));
    return ptr;
}

void cuda_free_host(void *ptr) {
    if (ptr) {
        CUDA_CHECK(cudaFreeHost(ptr));
    }
}

void* cuda_malloc_device(size_t size) {
    void *ptr;
    CUDA_CHECK(cudaMalloc(&ptr, size));
    return ptr;
}

void cuda_free_device(void *ptr) {
    if (ptr) {
        CUDA_CHECK(cudaFree(ptr));
    }
}

// Stream management
cudaStream_t cuda_get_stream(cuda_context_t *ctx, int stream_id) {
    if (!ctx->initialized || stream_id >= ctx->stream_count) {
        return 0; // Default stream
    }
    return ctx->streams[stream_id];
}

void cuda_synchronize_stream(cudaStream_t stream) {
    CUDA_CHECK(cudaStreamSynchronize(stream));
}

void cuda_synchronize_all_streams(cuda_context_t *ctx) {
    if (!ctx->initialized) {
        return;
    }
    
    for (int i = 0; i < ctx->stream_count; i++) {
        cuda_synchronize_stream(ctx->streams[i]);
    }
}

// Global context access functions
cuda_context_t* cuda_get_global_context(void) {
    return &g_cuda_context;
}

int cuda_init_global_context(int device_id) {
    return cuda_init_context(&g_cuda_context, device_id);
}

void cuda_cleanup_global_context(void) {
    cuda_cleanup_context(&g_cuda_context);
}
