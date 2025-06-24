/*  cuda_pileup.cuh -- CUDA GPU-accelerated pileup processing for BAM data.

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

#ifndef CUDA_PILEUP_CUH
#define CUDA_PILEUP_CUH

#include <cuda_runtime.h>
#include <stdint.h>
#include <stdbool.h>
#include "cuda_config.h"

#ifdef __cplusplus
extern "C" {
#endif

// GPU-optimized pileup data structures
#define CUDA_MAX_PILEUP_DEPTH 10000
#define CUDA_PILEUP_BLOCK_SIZE 512
#define CUDA_BASE_COUNT 16  // A, C, G, T, N and other ambiguous bases

// Base encoding for GPU processing
typedef enum {
    CUDA_BASE_A = 0,
    CUDA_BASE_C = 1,
    CUDA_BASE_G = 2,
    CUDA_BASE_T = 3,
    CUDA_BASE_N = 4,
    CUDA_BASE_GAP = 15
} cuda_base_t;

// Compact BAM alignment representation for GPU processing
typedef struct {
    int32_t tid;           // Target sequence ID
    int32_t pos;           // 0-based position
    int32_t end_pos;       // End position
    uint16_t flag;         // BAM flags
    uint8_t mapq;          // Mapping quality
    uint16_t cigar_len;    // Number of CIGAR operations
    uint32_t *cigar_ops;   // CIGAR operations (device pointer)
    uint8_t *seq;          // Sequence data (device pointer)
    uint8_t *qual;         // Quality scores (device pointer)
    uint32_t seq_len;      // Sequence length
} cuda_bam_record_t;

// Pileup position data
typedef struct {
    int32_t tid;           // Target sequence ID
    int32_t pos;           // Position
    uint16_t depth;        // Read depth at this position
    uint32_t base_counts[CUDA_BASE_COUNT]; // Base counts A, C, G, T, N...
    uint32_t qual_sum;     // Sum of quality scores
    uint16_t mapping_qual_sum; // Sum of mapping qualities
    uint8_t strand_counts[2];  // Forward/reverse strand counts
} cuda_pileup_position_t;

// Consensus calling results
typedef struct {
    int32_t pos;           // Position
    char consensus_base;   // Called consensus base
    uint8_t quality;       // Consensus quality score
    uint16_t depth;        // Total depth
    float confidence;      // Confidence score (0.0-1.0)
} cuda_consensus_result_t;

// CUDA pileup context
typedef struct {
    cuda_bam_record_t *d_records;       // Device BAM records
    cuda_pileup_position_t *d_pileup;   // Device pileup data
    cuda_consensus_result_t *d_consensus; // Device consensus results
    
    size_t max_records;                 // Maximum number of records
    size_t max_positions;               // Maximum number of positions
    
    cudaStream_t stream;                // CUDA stream
    
    // Configuration parameters
    uint8_t min_base_quality;           // Minimum base quality
    uint8_t min_mapping_quality;        // Minimum mapping quality
    uint16_t min_depth;                 // Minimum depth for consensus
    float min_consensus_fraction;       // Minimum fraction for consensus
} cuda_pileup_context_t;

// Function declarations

// Context management
int cuda_pileup_init_context(cuda_pileup_context_t *ctx, 
                             size_t max_records, 
                             size_t max_positions);
void cuda_pileup_cleanup_context(cuda_pileup_context_t *ctx);

// Main pileup processing functions
int cuda_pileup_process_region(cuda_pileup_context_t *ctx,
                              const cuda_bam_record_t *records,
                              size_t num_records,
                              int32_t start_pos,
                              int32_t end_pos,
                              cuda_pileup_position_t *output_pileup,
                              size_t *output_count);

int cuda_pileup_generate_consensus(cuda_pileup_context_t *ctx,
                                  const cuda_pileup_position_t *pileup_data,
                                  size_t pileup_count,
                                  cuda_consensus_result_t *consensus_results,
                                  size_t *consensus_count);

// Specialized pileup operations
int cuda_pileup_count_bases_parallel(const cuda_bam_record_t *records,
                                    size_t num_records,
                                    int32_t position,
                                    uint32_t *base_counts,
                                    cudaStream_t stream);

int cuda_pileup_calculate_coverage(const cuda_bam_record_t *records,
                                  size_t num_records,
                                  int32_t start_pos,
                                  int32_t end_pos,
                                  uint16_t *coverage_array,
                                  cudaStream_t stream);

// Quality filtering and processing
int cuda_pileup_filter_by_quality(cuda_bam_record_t *records,
                                 size_t num_records,
                                 uint8_t min_base_qual,
                                 uint8_t min_map_qual,
                                 bool *filter_mask,
                                 cudaStream_t stream);

// Variant calling support
int cuda_pileup_call_variants(cuda_pileup_context_t *ctx,
                             const cuda_pileup_position_t *pileup_data,
                             size_t pileup_count,
                             cuda_consensus_result_t *variant_calls,
                             size_t *variant_count);

// Statistics computation
typedef struct {
    double mean_depth;
    double median_depth;
    uint32_t max_depth;
    uint32_t min_depth;
    uint64_t total_bases;
    uint32_t covered_positions;
    double gc_content;
} cuda_pileup_stats_t;

int cuda_pileup_compute_stats(const cuda_pileup_position_t *pileup_data,
                             size_t pileup_count,
                             cuda_pileup_stats_t *stats,
                             cudaStream_t stream);

// Memory management and optimization
int cuda_pileup_optimize_memory_layout(cuda_bam_record_t *records, size_t count);
void cuda_pileup_prefetch_data(cuda_pileup_context_t *ctx, int device);

// Utility functions
int cuda_encode_sequence_to_gpu(const char *sequence, 
                               size_t length, 
                               uint8_t *gpu_sequence);

int cuda_decode_sequence_from_gpu(const uint8_t *gpu_sequence, 
                                 size_t length, 
                                 char *sequence);

char cuda_decode_base(uint8_t encoded_base);
uint8_t cuda_encode_base(char base);

// Error handling and validation
bool cuda_pileup_validate_input(const cuda_bam_record_t *records, size_t count);
int cuda_pileup_check_memory_requirements(size_t num_records, size_t num_positions);

#ifdef __cplusplus
}
#endif

#endif // CUDA_PILEUP_CUH
