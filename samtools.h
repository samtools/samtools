/*  samtools.h -- utility routines and GPU acceleration interface.

    Copyright (C) 2013-2015, 2019, 2023 Genome Research Ltd.
    CUDA GPU acceleration (C) 2025 Akshay Dedaniya <Dedaniya08@hotmail.com>

    Author: Petr Danecek <pd3@sanger.ac.uk>

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

#ifndef SAMTOOLS_H
#define SAMTOOLS_H

#include "htslib/hts_defs.h"
#include "htslib/sam.h"
#include "sam_utils.h"

// =============================================================================
// Core Samtools Functions
// =============================================================================

/**
 * Get the current samtools version string.
 * @return version string (e.g., "1.0.0")
 */
const char *samtools_version(void);

// =============================================================================
// BAM Sanitizer Options and Functions
// =============================================================================

/* BAM sanitizer bit flags for various corrections */
#define FIX_POS     2      // Fix invalid positions
#define FIX_MQUAL   4      // Fix mapping quality issues  
#define FIX_UNMAP   8      // Fix unmapped read flags
#define FIX_CIGAR   16     // Fix CIGAR string issues
#define FIX_AUX     32     // Fix auxiliary field issues
#define FIX_CIGDUP  64     // Fix duplicate CIGAR operations
#define FIX_CIGARX  128    // Fix extended CIGAR operations

// Default fixes for position-sorted data (most common issues)
#define FIX_ON (FIX_MQUAL|FIX_UNMAP|FIX_CIGAR|FIX_AUX|FIX_CIGDUP)
#define FIX_ALL 127        // Apply all available fixes

/**
 * Parse sanitizer options from comma-separated string.
 * Valid keywords: "pos", "mqual", "unmap", "cigar", "cigdup", "cigarx", "aux"
 * @param str options string
 * @return bit flags for bam_sanitize()
 */
int bam_sanitize_options(const char *str);

/**
 * Sanitize a BAM record using specified fix flags.
 * @param h SAM header
 * @param b BAM record to fix
 * @param flags FIX_* bit flags
 * @return 0 on success, <0 on failure
 */
int bam_sanitize(sam_hdr_t *h, bam1_t *b, int flags);

// =============================================================================
// CUDA GPU Acceleration Support
// =============================================================================

#ifdef ENABLE_CUDA

#include "cuda/cuda_config.h"
#include "cuda/cuda_sort.cuh"
#include "cuda/cuda_pileup.cuh"
#include "cuda/cuda_stats.cuh"

// Global CUDA context shared across all samtools operations
extern cuda_context_t *global_cuda_context;

/**
 * Initialize CUDA GPU acceleration for samtools.
 * Call once at program startup before using GPU features.
 * 
 * @return 0 on success, -1 if CUDA unavailable (falls back to CPU)
 */
int samtools_cuda_init(void);

/**
 * Clean up CUDA resources and shutdown GPU acceleration.
 * Call once at program shutdown to free all GPU resources.
 */
void samtools_cuda_cleanup(void);

/**
 * Check if CUDA acceleration is currently enabled and functional.
 * 
 * @return true if GPU acceleration is available, false for CPU-only mode
 */
bool samtools_cuda_enabled(void);

#else // !ENABLE_CUDA

// Stub functions for CPU-only builds
static inline int samtools_cuda_init(void) { return -1; }
static inline void samtools_cuda_cleanup(void) { }
static inline bool samtools_cuda_enabled(void) { return false; }

#endif // ENABLE_CUDA

#endif // SAMTOOLS_H
