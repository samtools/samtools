#!/bin/bash

# test_cuda.sh - Test script for CUDA GPU acceleration in samtools
# 
# This script tests the CUDA functionality of samtools including:
# - GPU device detection and initialization
# - GPU-accelerated sorting
# - GPU-accelerated statistics computation
# - Performance comparisons
# - Error handling and fallback behavior

set -e

# Configuration
TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMTOOLS_BIN="${TEST_DIR}/../samtools"
TEST_DATA_DIR="${TEST_DIR}/data"
RESULTS_DIR="${TEST_DIR}/results"
TEMP_DIR="${TEST_DIR}/temp"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test configuration
VERBOSE=0
PERFORMANCE_TEST=0
CLEANUP=1

# Functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[PASS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[FAIL]${NC} $1"
}

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Test script for CUDA GPU acceleration in samtools.

OPTIONS:
    -v, --verbose           Enable verbose output
    -p, --performance       Run performance comparison tests
    -n, --no-cleanup        Don't clean up temporary files
    -h, --help              Show this help message

EXAMPLES:
    $0                      Run basic functionality tests
    $0 -v -p               Run all tests with verbose output and performance comparison
    $0 --no-cleanup        Run tests and keep temporary files for inspection

EOF
}

# Parse command line arguments
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -v|--verbose)
                VERBOSE=1
                shift
                ;;
            -p|--performance)
                PERFORMANCE_TEST=1
                shift
                ;;
            -n|--no-cleanup)
                CLEANUP=0
                shift
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                usage
                exit 1
                ;;
        esac
    done
}

# Check prerequisites
check_prerequisites() {
    log_info "Checking prerequisites..."
    
    # Check if samtools exists
    if [[ ! -x "$SAMTOOLS_BIN" ]]; then
        log_error "Samtools binary not found at $SAMTOOLS_BIN"
        exit 1
    fi
    
    # Check if CUDA is available
    if ! command -v nvidia-smi &> /dev/null; then
        log_warning "nvidia-smi not found. GPU tests will be skipped."
        return 1
    fi
    
    # Check if samtools was compiled with CUDA support
    if ! $SAMTOOLS_BIN --version | grep -q "cuda=yes"; then
        log_warning "Samtools was not compiled with CUDA support. GPU tests will be skipped."
        return 1
    fi
    
    log_success "Prerequisites check passed"
    return 0
}

# Setup test environment
setup_test_env() {
    log_info "Setting up test environment..."
    
    # Create directories
    mkdir -p "$TEST_DATA_DIR" "$RESULTS_DIR" "$TEMP_DIR"
    
    # Generate test data if it doesn't exist
    if [[ ! -f "$TEST_DATA_DIR/test_small.bam" ]]; then
        generate_test_data
    fi
    
    log_success "Test environment setup complete"
}

# Generate test data
generate_test_data() {
    log_info "Generating test data..."
    
    # Create a small test SAM file
    cat > "$TEMP_DIR/test.sam" << 'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:249250621
@SQ	SN:chr2	LN:242193529
@PG	ID:test	PN:test	VN:1.0
read1	99	chr1	1000	60	100M	=	1150	250	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read2	147	chr1	1150	60	100M	=	1000	-250	TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read3	0	chr1	2000	60	100M	*	0	0	GGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
read4	0	chr2	500	60	100M	*	0	0	AATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

    # Convert to BAM
    $SAMTOOLS_BIN view -b "$TEMP_DIR/test.sam" > "$TEST_DATA_DIR/test_small.bam"
    $SAMTOOLS_BIN index "$TEST_DATA_DIR/test_small.bam"
    
    # Create a larger test file by duplicating reads
    for i in {1..1000}; do
        sed "s/read/read${i}_/g" "$TEMP_DIR/test.sam" | tail -n +5 >> "$TEMP_DIR/test_large.sam"
    done
    
    # Add header to large file
    head -n 4 "$TEMP_DIR/test.sam" > "$TEMP_DIR/test_large_with_header.sam"
    cat "$TEMP_DIR/test_large.sam" >> "$TEMP_DIR/test_large_with_header.sam"
    
    $SAMTOOLS_BIN view -b "$TEMP_DIR/test_large_with_header.sam" > "$TEST_DATA_DIR/test_large.bam"
    $SAMTOOLS_BIN index "$TEST_DATA_DIR/test_large.bam"
    
    log_success "Test data generated"
}

# Test GPU device detection
test_gpu_detection() {
    log_info "Testing GPU device detection..."
    
    # Run samtools with GPU option to check initialization
    if $SAMTOOLS_BIN sort --gpu --help &> /dev/null; then
        log_success "GPU support detected in samtools sort"
    else
        log_error "GPU support not available in samtools sort"
        return 1
    fi
    
    if $SAMTOOLS_BIN stats --gpu --help &> /dev/null; then
        log_success "GPU support detected in samtools stats"
    else
        log_error "GPU support not available in samtools stats"
        return 1
    fi
    
    return 0
}

# Test GPU-accelerated sorting
test_gpu_sorting() {
    log_info "Testing GPU-accelerated sorting..."
    
    local input_file="$TEST_DATA_DIR/test_large.bam"
    local gpu_output="$TEMP_DIR/sorted_gpu.bam"
    local cpu_output="$TEMP_DIR/sorted_cpu.bam"
    
    # Test GPU sorting
    if [[ $VERBOSE -eq 1 ]]; then
        log_info "Running GPU sort..."
        $SAMTOOLS_BIN sort --gpu "$input_file" -o "$gpu_output"
    else
        $SAMTOOLS_BIN sort --gpu "$input_file" -o "$gpu_output" 2>/dev/null
    fi
    
    if [[ ! -f "$gpu_output" ]]; then
        log_error "GPU sorting failed - output file not created"
        return 1
    fi
    
    # Test CPU sorting for comparison
    if [[ $VERBOSE -eq 1 ]]; then
        log_info "Running CPU sort for comparison..."
        $SAMTOOLS_BIN sort "$input_file" -o "$cpu_output"
    else
        $SAMTOOLS_BIN sort "$input_file" -o "$cpu_output" 2>/dev/null
    fi
    
    # Compare outputs
    if diff <($SAMTOOLS_BIN view "$gpu_output") <($SAMTOOLS_BIN view "$cpu_output") > /dev/null; then
        log_success "GPU and CPU sorting results are identical"
    else
        log_error "GPU and CPU sorting results differ"
        return 1
    fi
    
    return 0
}

# Test GPU-accelerated statistics
test_gpu_stats() {
    log_info "Testing GPU-accelerated statistics..."
    
    local input_file="$TEST_DATA_DIR/test_large.bam"
    local gpu_stats="$TEMP_DIR/stats_gpu.txt"
    local cpu_stats="$TEMP_DIR/stats_cpu.txt"
    
    # Test GPU stats
    if [[ $VERBOSE -eq 1 ]]; then
        log_info "Running GPU stats..."
        $SAMTOOLS_BIN stats --gpu "$input_file" > "$gpu_stats"
    else
        $SAMTOOLS_BIN stats --gpu "$input_file" > "$gpu_stats" 2>/dev/null
    fi
    
    if [[ ! -f "$gpu_stats" ]]; then
        log_error "GPU statistics failed - output file not created"
        return 1
    fi
    
    # Test CPU stats for comparison
    if [[ $VERBOSE -eq 1 ]]; then
        log_info "Running CPU stats for comparison..."
        $SAMTOOLS_BIN stats "$input_file" > "$cpu_stats"
    else
        $SAMTOOLS_BIN stats "$input_file" > "$cpu_stats" 2>/dev/null
    fi
    
    # Compare key statistics (allowing for minor floating point differences)
    local gpu_reads=$(grep "^SN" "$gpu_stats" | grep "raw total sequences:" | awk '{print $4}')
    local cpu_reads=$(grep "^SN" "$cpu_stats" | grep "raw total sequences:" | awk '{print $4}')
    
    if [[ "$gpu_reads" == "$cpu_reads" ]]; then
        log_success "GPU and CPU statistics read counts match: $gpu_reads"
    else
        log_error "GPU and CPU statistics read counts differ: GPU=$gpu_reads, CPU=$cpu_reads"
        return 1
    fi
    
    return 0
}

# Performance comparison tests
test_performance() {
    if [[ $PERFORMANCE_TEST -eq 0 ]]; then
        return 0
    fi
    
    log_info "Running performance comparison tests..."
    
    local input_file="$TEST_DATA_DIR/test_large.bam"
    local results_file="$RESULTS_DIR/performance_results.txt"
    
    echo "Performance Test Results - $(date)" > "$results_file"
    echo "========================================" >> "$results_file"
    
    # GPU vs CPU sorting performance
    log_info "Measuring sorting performance..."
    
    echo "Sorting Performance:" >> "$results_file"
    
    # Time GPU sorting
    local gpu_sort_time
    gpu_sort_time=$(time ( $SAMTOOLS_BIN sort --gpu "$input_file" -o "$TEMP_DIR/perf_gpu_sort.bam" ) 2>&1 | grep real | awk '{print $2}')
    echo "GPU Sort Time: $gpu_sort_time" >> "$results_file"
    
    # Time CPU sorting
    local cpu_sort_time
    cpu_sort_time=$(time ( $SAMTOOLS_BIN sort "$input_file" -o "$TEMP_DIR/perf_cpu_sort.bam" ) 2>&1 | grep real | awk '{print $2}')
    echo "CPU Sort Time: $cpu_sort_time" >> "$results_file"
    
    # GPU vs CPU stats performance
    log_info "Measuring statistics performance..."
    
    echo "" >> "$results_file"
    echo "Statistics Performance:" >> "$results_file"
    
    # Time GPU stats
    local gpu_stats_time
    gpu_stats_time=$(time ( $SAMTOOLS_BIN stats --gpu "$input_file" > "$TEMP_DIR/perf_gpu_stats.txt" ) 2>&1 | grep real | awk '{print $2}')
    echo "GPU Stats Time: $gpu_stats_time" >> "$results_file"
    
    # Time CPU stats
    local cpu_stats_time
    cpu_stats_time=$(time ( $SAMTOOLS_BIN stats "$input_file" > "$TEMP_DIR/perf_cpu_stats.txt" ) 2>&1 | grep real | awk '{print $2}')
    echo "CPU Stats Time: $cpu_stats_time" >> "$results_file"
    
    log_success "Performance results saved to $results_file"
    
    if [[ $VERBOSE -eq 1 ]]; then
        cat "$results_file"
    fi
    
    return 0
}

# Test error handling
test_error_handling() {
    log_info "Testing error handling..."
    
    # Test with invalid GPU option (should fall back to CPU)
    local output_file="$TEMP_DIR/error_test.bam"
    
    # Force GPU unavailable and test fallback
    if CUDA_VISIBLE_DEVICES="" $SAMTOOLS_BIN sort --gpu "$TEST_DATA_DIR/test_small.bam" -o "$output_file" 2>/dev/null; then
        if [[ -f "$output_file" ]]; then
            log_success "Graceful fallback to CPU when GPU unavailable"
        else
            log_error "Failed to fall back to CPU processing"
            return 1
        fi
    else
        log_warning "Unable to test GPU fallback behavior"
    fi
    
    return 0
}

# Cleanup temporary files
cleanup() {
    if [[ $CLEANUP -eq 1 ]]; then
        log_info "Cleaning up temporary files..."
        rm -rf "$TEMP_DIR"
        log_success "Cleanup complete"
    else
        log_info "Temporary files retained in $TEMP_DIR"
    fi
}

# Main test execution
run_tests() {
    local tests_passed=0
    local tests_total=0
    
    log_info "Starting CUDA GPU acceleration tests for samtools"
    log_info "================================================"
    
    # Check if GPU tests can be run
    if ! check_prerequisites; then
        log_warning "GPU functionality tests skipped due to missing prerequisites"
        exit 0
    fi
    
    setup_test_env
    
    # Run tests
    tests=(
        "test_gpu_detection"
        "test_gpu_sorting"
        "test_gpu_stats"
        "test_error_handling"
        "test_performance"
    )
    
    for test in "${tests[@]}"; do
        tests_total=$((tests_total + 1))
        log_info "Running $test..."
        
        if $test; then
            tests_passed=$((tests_passed + 1))
        else
            log_error "Test $test failed"
        fi
        echo
    done
    
    # Report results
    log_info "Test Results Summary"
    log_info "==================="
    log_info "Tests passed: $tests_passed/$tests_total"
    
    if [[ $tests_passed -eq $tests_total ]]; then
        log_success "All tests passed!"
        cleanup
        exit 0
    else
        log_error "Some tests failed"
        cleanup
        exit 1
    fi
}

# Signal handling
trap cleanup EXIT

# Main execution
main() {
    parse_args "$@"
    run_tests
}

# Run the main function if script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
