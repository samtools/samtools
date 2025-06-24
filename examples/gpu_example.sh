#!/bin/bash

# gpu_example.sh - Example usage of GPU-accelerated samtools
# 
# This script demonstrates how to use the CUDA GPU acceleration features
# in samtools for common bioinformatics workflows.

set -e

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMTOOLS="${SCRIPT_DIR}/../samtools"
INPUT_BAM="${1:-example.bam}"
OUTPUT_DIR="${2:-gpu_results}"

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

usage() {
    cat << EOF
Usage: $0 [INPUT_BAM] [OUTPUT_DIR]

Example usage of GPU-accelerated samtools operations.

ARGUMENTS:
    INPUT_BAM    Input BAM file (default: example.bam)
    OUTPUT_DIR   Output directory (default: gpu_results)

EXAMPLES:
    $0                              # Use default files
    $0 sample.bam results          # Process sample.bam, output to results/
    $0 /path/to/large.bam output   # Process large file with GPU acceleration

REQUIREMENTS:
    - NVIDIA GPU with CUDA support
    - Samtools compiled with --enable-cuda
    - Input BAM file

EOF
}

check_requirements() {
    log_info "Checking requirements..."
    
    # Check if samtools exists and has CUDA support
    if [[ ! -x "$SAMTOOLS" ]]; then
        echo "Error: Samtools not found at $SAMTOOLS"
        echo "Please build samtools with CUDA support first."
        exit 1
    fi
    
    if ! $SAMTOOLS --version | grep -q "cuda=yes"; then
        echo "Error: Samtools was not compiled with CUDA support."
        echo "Please rebuild with: ./configure --enable-cuda && make"
        exit 1
    fi
    
    # Check for NVIDIA GPU
    if ! command -v nvidia-smi &> /dev/null; then
        log_warning "nvidia-smi not found. GPU acceleration may not be available."
    else
        log_info "GPU status:"
        nvidia-smi --query-gpu=name,memory.total,memory.used,utilization.gpu --format=csv,noheader,nounits | head -1
    fi
    
    # Check input file
    if [[ ! -f "$INPUT_BAM" ]]; then
        echo "Error: Input BAM file '$INPUT_BAM' not found."
        echo "Please provide a valid BAM file or create an example file."
        exit 1
    fi
    
    log_success "Requirements check passed"
}

create_output_dir() {
    mkdir -p "$OUTPUT_DIR"
    log_info "Output directory: $OUTPUT_DIR"
}

demonstrate_gpu_sorting() {
    log_info "=== GPU-Accelerated Sorting Demo ==="
    
    local sorted_bam="$OUTPUT_DIR/sorted_gpu.bam"
    local sorted_bam_cpu="$OUTPUT_DIR/sorted_cpu.bam"
    
    # GPU sorting
    log_info "Running GPU-accelerated coordinate sorting..."
    time $SAMTOOLS sort --gpu "$INPUT_BAM" -o "$sorted_bam"
    log_success "GPU sorting completed: $sorted_bam"
    
    # CPU sorting for comparison
    log_info "Running CPU sorting for comparison..."
    time $SAMTOOLS sort "$INPUT_BAM" -o "$sorted_bam_cpu"
    log_success "CPU sorting completed: $sorted_bam_cpu"
    
    # Verify results are identical
    log_info "Verifying GPU and CPU results are identical..."
    if cmp -s <($SAMTOOLS view "$sorted_bam") <($SAMTOOLS view "$sorted_bam_cpu"); then
        log_success "GPU and CPU sorting results are identical"
        rm "$sorted_bam_cpu"  # Remove duplicate
    else
        log_warning "GPU and CPU sorting results differ - this should not happen!"
    fi
    
    # Index the sorted file
    log_info "Creating index for sorted BAM..."
    $SAMTOOLS index "$sorted_bam"
    log_success "Index created: ${sorted_bam}.bai"
}

demonstrate_gpu_statistics() {
    log_info "=== GPU-Accelerated Statistics Demo ==="
    
    local stats_gpu="$OUTPUT_DIR/stats_gpu.txt"
    local stats_cpu="$OUTPUT_DIR/stats_cpu.txt"
    
    # GPU statistics
    log_info "Computing statistics with GPU acceleration..."
    time $SAMTOOLS stats --gpu "$INPUT_BAM" > "$stats_gpu"
    log_success "GPU statistics completed: $stats_gpu"
    
    # CPU statistics for comparison
    log_info "Computing statistics with CPU for comparison..."
    time $SAMTOOLS stats "$INPUT_BAM" > "$stats_cpu"
    log_success "CPU statistics completed: $stats_cpu"
    
    # Compare key metrics
    log_info "Comparing key statistics..."
    
    local gpu_reads=$(grep "^SN" "$stats_gpu" | grep "raw total sequences:" | awk '{print $4}')
    local cpu_reads=$(grep "^SN" "$stats_cpu" | grep "raw total sequences:" | awk '{print $4}')
    
    echo "Total reads - GPU: $gpu_reads, CPU: $cpu_reads"
    
    if [[ "$gpu_reads" == "$cpu_reads" ]]; then
        log_success "GPU and CPU statistics match"
        rm "$stats_cpu"  # Remove duplicate
    else
        log_warning "GPU and CPU statistics differ - investigating differences"
        echo "Differences found in $stats_gpu vs $stats_cpu"
    fi
    
    # Show summary statistics
    log_info "Statistics summary:"
    grep "^SN" "$stats_gpu" | head -10
}

demonstrate_advanced_workflows() {
    log_info "=== Advanced GPU Workflow Demo ==="
    
    local input_file="$INPUT_BAM"
    
    # If we have a sorted BAM from previous step, use it
    if [[ -f "$OUTPUT_DIR/sorted_gpu.bam" ]]; then
        input_file="$OUTPUT_DIR/sorted_gpu.bam"
        log_info "Using previously sorted BAM: $input_file"
    fi
    
    # Combined workflow: sort + index + stats
    log_info "Running combined GPU workflow..."
    
    local workflow_output="$OUTPUT_DIR/workflow_result.bam"
    local workflow_stats="$OUTPUT_DIR/workflow_stats.txt"
    
    # Sort with GPU
    $SAMTOOLS sort --gpu "$INPUT_BAM" -o "$workflow_output"
    
    # Index
    $SAMTOOLS index "$workflow_output"
    
    # Generate comprehensive statistics with GPU
    $SAMTOOLS stats --gpu "$workflow_output" > "$workflow_stats"
    
    log_success "Combined workflow completed"
    log_info "Results:"
    echo "  - Sorted BAM: $workflow_output"
    echo "  - Index: ${workflow_output}.bai"
    echo "  - Statistics: $workflow_stats"
}

benchmark_performance() {
    log_info "=== Performance Benchmark ==="
    
    if [[ ! -f "$INPUT_BAM" ]]; then
        log_warning "Skipping benchmark - no input file"
        return
    fi
    
    local benchmark_file="$OUTPUT_DIR/benchmark_results.txt"
    
    echo "Performance Benchmark Results - $(date)" > "$benchmark_file"
    echo "========================================" >> "$benchmark_file"
    echo "Input file: $INPUT_BAM" >> "$benchmark_file"
    echo "File size: $(du -h "$INPUT_BAM" | cut -f1)" >> "$benchmark_file"
    echo "" >> "$benchmark_file"
    
    # Benchmark sorting
    echo "Sorting Performance:" >> "$benchmark_file"
    
    log_info "Benchmarking GPU sorting..."
    local gpu_sort_time=$( (time $SAMTOOLS sort --gpu "$INPUT_BAM" -o "$OUTPUT_DIR/bench_gpu.bam") 2>&1 | grep real | awk '{print $2}')
    echo "GPU Sort Time: $gpu_sort_time" >> "$benchmark_file"
    
    log_info "Benchmarking CPU sorting..."
    local cpu_sort_time=$( (time $SAMTOOLS sort "$INPUT_BAM" -o "$OUTPUT_DIR/bench_cpu.bam") 2>&1 | grep real | awk '{print $2}')
    echo "CPU Sort Time: $cpu_sort_time" >> "$benchmark_file"
    
    # Benchmark statistics
    echo "" >> "$benchmark_file"
    echo "Statistics Performance:" >> "$benchmark_file"
    
    log_info "Benchmarking GPU statistics..."
    local gpu_stats_time=$( (time $SAMTOOLS stats --gpu "$INPUT_BAM" > "$OUTPUT_DIR/bench_gpu_stats.txt") 2>&1 | grep real | awk '{print $2}')
    echo "GPU Stats Time: $gpu_stats_time" >> "$benchmark_file"
    
    log_info "Benchmarking CPU statistics..."
    local cpu_stats_time=$( (time $SAMTOOLS stats "$INPUT_BAM" > "$OUTPUT_DIR/bench_cpu_stats.txt") 2>&1 | grep real | awk '{print $2}')
    echo "CPU Stats Time: $cpu_stats_time" >> "$benchmark_file"
    
    echo "" >> "$benchmark_file"
    echo "System Information:" >> "$benchmark_file"
    if command -v nvidia-smi &> /dev/null; then
        echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader,nounits | head -1)" >> "$benchmark_file"
        echo "GPU Memory: $(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits | head -1) MB" >> "$benchmark_file"
    fi
    echo "CPU: $(nproc) cores" >> "$benchmark_file"
    echo "RAM: $(free -h | grep Mem | awk '{print $2}')" >> "$benchmark_file"
    
    log_success "Benchmark results saved to: $benchmark_file"
    
    # Cleanup benchmark files
    rm -f "$OUTPUT_DIR"/bench_*.bam "$OUTPUT_DIR"/bench_*.txt
    
    # Display results
    log_info "Performance summary:"
    cat "$benchmark_file"
}

cleanup_demo() {
    log_info "Demo completed. Results available in: $OUTPUT_DIR"
    
    # List all output files
    echo "Generated files:"
    ls -lh "$OUTPUT_DIR"
    
    if command -v nvidia-smi &> /dev/null; then
        echo ""
        log_info "Final GPU status:"
        nvidia-smi --query-gpu=name,memory.used,utilization.gpu --format=csv,noheader,nounits | head -1
    fi
}

main() {
    # Parse arguments
    if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
        usage
        exit 0
    fi
    
    log_info "GPU-Accelerated Samtools Demo"
    log_info "============================="
    
    check_requirements
    create_output_dir
    
    demonstrate_gpu_sorting
    echo ""
    
    demonstrate_gpu_statistics
    echo ""
    
    demonstrate_advanced_workflows
    echo ""
    
    benchmark_performance
    echo ""
    
    cleanup_demo
}

# Run the main function
main "$@"
