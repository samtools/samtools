/*  samtools_cuda.c -- CUDA GPU acceleration integration for samtools.

    Copyright (C) 2025 Akshay Dedaniya <Dedaniya08@hotmail.com>

    This module provides the main CUDA integration layer for samtools,
    managing global GPU context and providing high-level GPU acceleration
    functions for sorting, pileup generation, and statistics computation.
    
    Features:
    - Global CUDA context management
    - Automatic device selection and initialization  
    - Graceful fallback to CPU when GPU unavailable
    - Resource cleanup and error handling

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

#include "samtools.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef ENABLE_CUDA

// =============================================================================
// Global CUDA State Management
// =============================================================================

// Global CUDA context shared across all samtools operations
cuda_context_t *global_cuda_context = NULL;

/**
 * Initialize CUDA subsystem for samtools.
 * 
 * This function:
 * 1. Checks for CUDA availability
 * 2. Selects optimal GPU device
 * 3. Initializes global context with streams
 * 4. Provides diagnostic information
 * 
 * @return 0 on success, -1 on failure (falls back to CPU)
 */
int samtools_cuda_init(void) {
    // Prevent double initialization
    if (global_cuda_context != NULL) {
        return 0; // Already initialized successfully
    }
    
    // Check CUDA availability first
    if (!cuda_is_available()) {
        fprintf(stderr, "[CUDA] No CUDA-capable devices found. Falling back to CPU-only mode.\n");
        return -1;
    }
    
    // Allocate global context structure
    global_cuda_context = (cuda_context_t*)malloc(sizeof(cuda_context_t));
    if (!global_cuda_context) {
        fprintf(stderr, "[CUDA] Failed to allocate memory for CUDA context.\n");
        return -1;
    }
    
    // Auto-select the best available GPU device
    int device_id = cuda_select_best_device();
    if (device_id < 0) {
        fprintf(stderr, "[CUDA] Failed to select suitable CUDA device.\n");
        free(global_cuda_context);
        global_cuda_context = NULL;
        return -1;
    }
    
    // Initialize the CUDA context with selected device
    if (cuda_init_context(global_cuda_context, device_id) != 0) {
        fprintf(stderr, "[CUDA] Failed to initialize CUDA context on device %d.\n", device_id);
        free(global_cuda_context);
        global_cuda_context = NULL;
        return -1;
    }
    
    // Print device information for user feedback
    printf("[CUDA] GPU acceleration enabled.\n");
    cuda_print_device_info(&global_cuda_context->device_info);
    
    return 0;
}

/**
 * Clean up CUDA resources and reset global state.
 * 
 * This function:
 * 1. Synchronizes all pending GPU operations
 * 2. Frees CUDA context and streams
 * 3. Releases global context memory
 * 4. Resets global state for clean shutdown
 */
void samtools_cuda_cleanup(void) {
    // Only cleanup if context exists and is initialized
    if (global_cuda_context != NULL) {
        cuda_cleanup_context(global_cuda_context);
        free(global_cuda_context);
        global_cuda_context = NULL;
        printf("[CUDA] GPU acceleration cleaned up.\n");
    }
}

/**
 * Check if CUDA is currently enabled and functional.
 * 
 * @return true if CUDA context is initialized and ready for use
 */
bool samtools_cuda_enabled(void) {
    return global_cuda_context != NULL && global_cuda_context->initialized;
}

#else // !ENABLE_CUDA

// =============================================================================
// CUDA Disabled - Stub Functions
// =============================================================================

/**
 * Stub function when CUDA support is not compiled in.
 * Always returns failure to indicate CPU-only mode.
 */
int samtools_cuda_init(void) {
    return -1; // CUDA not available
}

/**
 * Stub cleanup function when CUDA is disabled.
 * No operation needed in CPU-only mode.
 */
void samtools_cuda_cleanup(void) {
    // Nothing to do - CPU-only mode
}

/**
 * Stub function to check CUDA availability when disabled.
 * Always returns false in CPU-only builds.
 */
bool samtools_cuda_enabled(void) {
    return false; // CUDA not compiled in
}

#endif // ENABLE_CUDA
