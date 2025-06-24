/*  test_cuda_config.c -- Unit tests for CUDA configuration and context management.

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
#include <assert.h>
#include <string.h>

#ifdef ENABLE_CUDA
#include "../cuda/cuda_config.h"

// Test CUDA device detection
void test_cuda_device_detection() {
    printf("Testing CUDA device detection...\n");
    
    bool available = cuda_is_available();
    printf("CUDA available: %s\n", available ? "yes" : "no");
    
    if (available) {
        int device_id = cuda_select_best_device();
        printf("Best device ID: %d\n", device_id);
        assert(device_id >= 0);
    }
    
    printf("✓ CUDA device detection test passed\n");
}

// Test CUDA context initialization
void test_cuda_context_init() {
    printf("Testing CUDA context initialization...\n");
    
    if (!cuda_is_available()) {
        printf("⚠ Skipping context test - no CUDA devices available\n");
        return;
    }
    
    cuda_context_t ctx;
    int device_id = cuda_select_best_device();
    
    int result = cuda_init_context(&ctx, device_id);
    if (result == 0) {
        printf("✓ Context initialized successfully\n");
        
        // Test context properties
        assert(ctx.initialized == true);
        assert(ctx.device_info.device_count > 0);
        assert(ctx.device_info.current_device == device_id);
        assert(ctx.streams != NULL);
        assert(ctx.stream_count > 0);
        
        // Print device information
        cuda_print_device_info(&ctx.device_info);
        
        // Cleanup
        cuda_cleanup_context(&ctx);
        printf("✓ Context cleaned up successfully\n");
    } else {
        printf("✗ Context initialization failed\n");
        assert(false);
    }
    
    printf("✓ CUDA context initialization test passed\n");
}

// Test memory management functions
void test_cuda_memory_management() {
    printf("Testing CUDA memory management...\n");
    
    if (!cuda_is_available()) {
        printf("⚠ Skipping memory test - no CUDA devices available\n");
        return;
    }
    
    // Test managed memory allocation
    size_t test_size = 1024 * 1024; // 1MB
    void *managed_ptr = cuda_malloc_managed(test_size);
    
    if (managed_ptr != NULL) {
        printf("✓ Managed memory allocated: %zu bytes\n", test_size);
        
        // Test memory access (write/read)
        memset(managed_ptr, 0x42, test_size);
        unsigned char *test_ptr = (unsigned char*)managed_ptr;
        assert(test_ptr[0] == 0x42);
        assert(test_ptr[test_size-1] == 0x42);
        
        cuda_free_managed(managed_ptr);
        printf("✓ Managed memory freed\n");
    } else {
        printf("✗ Managed memory allocation failed\n");
        assert(false);
    }
    
    // Test device memory allocation
    void *device_ptr = cuda_malloc_device(test_size);
    if (device_ptr != NULL) {
        printf("✓ Device memory allocated: %zu bytes\n", test_size);
        cuda_free_device(device_ptr);
        printf("✓ Device memory freed\n");
    } else {
        printf("✗ Device memory allocation failed\n");
        assert(false);
    }
    
    // Test host memory allocation
    void *host_ptr = cuda_malloc_host(test_size);
    if (host_ptr != NULL) {
        printf("✓ Host memory allocated: %zu bytes\n", test_size);
        
        // Test memory access
        memset(host_ptr, 0x24, test_size);
        unsigned char *test_host_ptr = (unsigned char*)host_ptr;
        assert(test_host_ptr[0] == 0x24);
        
        cuda_free_host(host_ptr);
        printf("✓ Host memory freed\n");
    } else {
        printf("✗ Host memory allocation failed\n");
        assert(false);
    }
    
    printf("✓ CUDA memory management test passed\n");
}

// Test utility functions
void test_cuda_utilities() {
    printf("Testing CUDA utility functions...\n");
    
    // Test optimal block size calculation
    int block_size = cuda_get_optimal_block_size(10000, 256);
    assert(block_size > 0);
    assert(block_size <= 256);
    printf("✓ Optimal block size for 10000 elements: %d\n", block_size);
    
    // Test optimal grid size calculation
    int grid_size = cuda_get_optimal_grid_size(10000, block_size);
    assert(grid_size > 0);
    printf("✓ Optimal grid size: %d\n", grid_size);
    
    // Verify that grid_size * block_size covers all elements
    assert(grid_size * block_size >= 10000);
    
    printf("✓ CUDA utility functions test passed\n");
}

// Test stream management
void test_cuda_streams() {
    printf("Testing CUDA stream management...\n");
    
    if (!cuda_is_available()) {
        printf("⚠ Skipping stream test - no CUDA devices available\n");
        return;
    }
    
    cuda_context_t ctx;
    int device_id = cuda_select_best_device();
    
    int result = cuda_init_context(&ctx, device_id);
    if (result != 0) {
        printf("✗ Failed to initialize context for stream test\n");
        assert(false);
    }
    
    // Test stream access
    for (int i = 0; i < ctx.stream_count; i++) {
        cudaStream_t stream = cuda_get_stream(&ctx, i);
        assert(stream != NULL);
        printf("✓ Stream %d accessed successfully\n", i);
    }
    
    // Test stream synchronization
    cudaStream_t stream0 = cuda_get_stream(&ctx, 0);
    cuda_synchronize_stream(stream0);
    printf("✓ Single stream synchronized\n");
    
    cuda_synchronize_all_streams(&ctx);
    printf("✓ All streams synchronized\n");
    
    cuda_cleanup_context(&ctx);
    printf("✓ CUDA stream management test passed\n");
}

// Main test runner
int main() {
    printf("=== CUDA Configuration Unit Tests ===\n\n");
    
    // Run tests
    test_cuda_device_detection();
    printf("\n");
    
    test_cuda_context_init();
    printf("\n");
    
    test_cuda_memory_management();
    printf("\n");
    
    test_cuda_utilities();
    printf("\n");
    
    test_cuda_streams();
    printf("\n");
    
    printf("=== All CUDA Configuration Tests Passed ===\n");
    return 0;
}

#else // !ENABLE_CUDA

// Stub main function when CUDA is disabled
int main() {
    printf("CUDA support is not enabled. Skipping CUDA configuration tests.\n");
    return 0;
}

#endif // ENABLE_CUDA
