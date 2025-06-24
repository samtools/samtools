/*  cuda_stats.cuh -- CUDA GPU-accelerated statistics computation for BAM data.

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

#ifndef CUDA_STATS_CUH
#define CUDA_STATS_CUH

#include <cuda_runtime.h>
#include <stdint.h>
#include <stdbool.h>
#include "cuda_config.h"

#ifdef __cplusplus
extern "C" {
#endif

// Statistical data structures optimized for GPU processing
#define CUDA_STATS_BLOCK_SIZE 256
#define CUDA_MAX_READ_LENGTH 1000
#define CUDA_MAX_QUALITY_SCORE 100
#define CUDA_GC_BINS 101  // 0-100% GC content

// Compact BAM record for statistics processing
typedef struct {
    uint32_t flag;             // BAM flags
    int32_t tid;               // Target sequence ID
    int32_t pos;               // Position
    int32_t mpos;              // Mate position
    int32_t isize;             // Insert size
    uint8_t mapq;              // Mapping quality
    uint16_t seq_len;          // Sequence length
    uint8_t *seq;              // Sequence data (device pointer)
    uint8_t *qual;             // Quality scores (device pointer)
} cuda_bam_stats_record_t;

// Comprehensive statistics structure
typedef struct {
    // Basic counts
    uint64_t total_reads;
    uint64_t mapped_reads;
    uint64_t properly_paired_reads;
    uint64_t singleton_reads;
    uint64_t duplicate_reads;
    uint64_t secondary_reads;
    uint64_t supplementary_reads;
    uint64_t failed_qc_reads;
    
    // Quality statistics
    uint64_t total_bases;
    uint64_t q20_bases;        // Bases with quality >= 20
    uint64_t q30_bases;        // Bases with quality >= 30
    double mean_quality;
    
    // Length statistics
    uint32_t min_read_length;
    uint32_t max_read_length;
    double mean_read_length;
    uint64_t read_length_histogram[CUDA_MAX_READ_LENGTH];
    
    // Insert size statistics (for paired reads)
    int32_t min_insert_size;
    int32_t max_insert_size;
    double mean_insert_size;
    double std_insert_size;
    
    // GC content statistics
    double gc_content;
    uint32_t gc_histogram[CUDA_GC_BINS];
    
    // Coverage statistics
    double mean_coverage;
    uint32_t max_coverage;
    uint64_t covered_bases;
    
    // Mapping quality statistics
    uint32_t mapq_histogram[256];
    double mean_mapq;
    
    // Error rate estimation
    double mismatch_rate;
    double indel_rate;
} cuda_bam_statistics_t;

// Per-chromosome statistics
typedef struct {
    int32_t tid;               // Chromosome ID
    uint64_t mapped_reads;     // Reads mapped to this chromosome
    uint64_t total_bases;      // Total bases mapped
    double mean_coverage;      // Mean coverage
    double gc_content;         // GC content
} cuda_chromosome_stats_t;

// CUDA statistics processing context
typedef struct {
    cuda_bam_stats_record_t *d_records;    // Device BAM records
    cuda_bam_statistics_t *d_stats;        // Device statistics
    cuda_chromosome_stats_t *d_chr_stats;  // Device per-chromosome stats
    
    uint64_t *d_temp_storage;              // Temporary storage for reductions
    size_t temp_storage_bytes;             // Size of temporary storage
    
    size_t max_records;                    // Maximum number of records
    size_t max_chromosomes;                // Maximum number of chromosomes
    
    cudaStream_t stream;                   // CUDA stream
    
    // Processing parameters
    uint8_t min_quality_threshold;         // Minimum quality for Q20/Q30 counts
    bool include_secondary;                // Include secondary alignments
    bool include_supplementary;            // Include supplementary alignments
} cuda_stats_context_t;

// Function declarations

// Context management
int cuda_stats_init_context(cuda_stats_context_t *ctx, 
                            size_t max_records, 
                            size_t max_chromosomes);
void cuda_stats_cleanup_context(cuda_stats_context_t *ctx);

// Main statistics computation functions
int cuda_stats_compute_comprehensive(cuda_stats_context_t *ctx,
                                    const cuda_bam_stats_record_t *records,
                                    size_t num_records,
                                    cuda_bam_statistics_t *output_stats);

int cuda_stats_compute_per_chromosome(cuda_stats_context_t *ctx,
                                     const cuda_bam_stats_record_t *records,
                                     size_t num_records,
                                     cuda_chromosome_stats_t *chr_stats,
                                     size_t *num_chromosomes);

// Specialized statistics functions
int cuda_stats_compute_gc_content(const cuda_bam_stats_record_t *records,
                                 size_t num_records,
                                 double *gc_content,
                                 uint32_t *gc_histogram,
                                 cudaStream_t stream);

int cuda_stats_compute_quality_metrics(const cuda_bam_stats_record_t *records,
                                      size_t num_records,
                                      uint8_t quality_threshold,
                                      uint64_t *q_bases,
                                      double *mean_quality,
                                      cudaStream_t stream);

int cuda_stats_compute_insert_size_distribution(const cuda_bam_stats_record_t *records,
                                               size_t num_records,
                                               int32_t *min_isize,
                                               int32_t *max_isize,
                                               double *mean_isize,
                                               double *std_isize,
                                               cudaStream_t stream);

int cuda_stats_compute_coverage_stats(const cuda_bam_stats_record_t *records,
                                     size_t num_records,
                                     int32_t ref_length,
                                     double *mean_coverage,
                                     uint32_t *max_coverage,
                                     uint64_t *covered_bases,
                                     cudaStream_t stream);

// Flag-based statistics
int cuda_stats_count_flags(const cuda_bam_stats_record_t *records,
                          size_t num_records,
                          uint64_t *flag_counts,
                          cudaStream_t stream);

// Length distribution analysis
int cuda_stats_compute_length_distribution(const cuda_bam_stats_record_t *records,
                                          size_t num_records,
                                          uint32_t *min_length,
                                          uint32_t *max_length,
                                          double *mean_length,
                                          uint64_t *length_histogram,
                                          cudaStream_t stream);

// Mapping quality analysis
int cuda_stats_compute_mapq_distribution(const cuda_bam_stats_record_t *records,
                                        size_t num_records,
                                        uint32_t *mapq_histogram,
                                        double *mean_mapq,
                                        cudaStream_t stream);

// Error rate estimation
typedef struct {
    uint64_t matches;
    uint64_t mismatches;
    uint64_t insertions;
    uint64_t deletions;
    double mismatch_rate;
    double indel_rate;
} cuda_error_stats_t;

int cuda_stats_estimate_error_rates(const cuda_bam_stats_record_t *records,
                                   size_t num_records,
                                   const char *reference_sequence,
                                   size_t ref_length,
                                   cuda_error_stats_t *error_stats,
                                   cudaStream_t stream);

// Parallel processing utilities
int cuda_stats_merge_partial_results(const cuda_bam_statistics_t *partial_stats,
                                    size_t num_partials,
                                    cuda_bam_statistics_t *merged_stats,
                                    cudaStream_t stream);

// Memory optimization
int cuda_stats_optimize_memory_access(cuda_bam_stats_record_t *records, 
                                     size_t count);

// Performance monitoring
typedef struct {
    float computation_time_ms;
    float memory_transfer_time_ms;
    size_t peak_memory_usage;
    float throughput_mreads_per_sec;
} cuda_stats_performance_t;

int cuda_stats_get_performance_metrics(cuda_stats_context_t *ctx,
                                      cuda_stats_performance_t *perf);

// Utility functions for data conversion
int cuda_stats_convert_bam_record(const void *bam_record,
                                 cuda_bam_stats_record_t *cuda_record);

int cuda_stats_batch_convert_records(const void **bam_records,
                                    size_t num_records,
                                    cuda_bam_stats_record_t *cuda_records);

// Output formatting
int cuda_stats_format_output(const cuda_bam_statistics_t *stats,
                            const cuda_chromosome_stats_t *chr_stats,
                            size_t num_chromosomes,
                            char *output_buffer,
                            size_t buffer_size);

// Validation and debugging
bool cuda_stats_validate_input(const cuda_bam_stats_record_t *records, size_t count);
int cuda_stats_debug_print_summary(const cuda_bam_statistics_t *stats);

#ifdef __cplusplus
}
#endif

#endif // CUDA_STATS_CUH
