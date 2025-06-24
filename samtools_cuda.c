/*  samtools_cuda.c -- CUDA GPU acceleration integration for samtools.

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

#include "samtools.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef ENABLE_CUDA

// Global CUDA context for samtools
cuda_context_t *global_cuda_context = NULL;

int samtools_cuda_init(void) {
    if (global_cuda_context != NULL) {
        return 0; // Already initialized
    }
    
    if (!cuda_is_available()) {
        fprintf(stderr, "[CUDA] No CUDA-capable devices found. Falling back to CPU-only mode.\n");
        return -1;
    }
    
    // Allocate and initialize CUDA context
    global_cuda_context = (cuda_context_t*)malloc(sizeof(cuda_context_t));
    if (!global_cuda_context) {
        fprintf(stderr, "[CUDA] Failed to allocate memory for CUDA context.\n");
        return -1;
    }
    
    // Select the best available device
    int device_id = cuda_select_best_device();
    if (device_id < 0) {
        fprintf(stderr, "[CUDA] Failed to select CUDA device.\n");
        free(global_cuda_context);
        global_cuda_context = NULL;
        return -1;
    }
    
    // Initialize the CUDA context
    if (cuda_init_context(global_cuda_context, device_id) != 0) {
        fprintf(stderr, "[CUDA] Failed to initialize CUDA context.\n");
        free(global_cuda_context);
        global_cuda_context = NULL;
        return -1;
    }
    
    // Print device information
    printf("[CUDA] GPU acceleration enabled.\n");
    cuda_print_device_info(&global_cuda_context->device_info);
    
    return 0;
}

void samtools_cuda_cleanup(void) {
    if (global_cuda_context != NULL) {
        cuda_cleanup_context(global_cuda_context);
        free(global_cuda_context);
        global_cuda_context = NULL;
        printf("[CUDA] GPU acceleration cleaned up.\n");
    }
}

bool samtools_cuda_enabled(void) {
    return global_cuda_context != NULL && global_cuda_context->initialized;
}

#else // !ENABLE_CUDA

// Stub functions when CUDA is disabled
int samtools_cuda_init(void) {
    return -1; // CUDA not available
}

void samtools_cuda_cleanup(void) {
    // Nothing to do
}

bool samtools_cuda_enabled(void) {
    return false;
}

#endif // ENABLE_CUDA
