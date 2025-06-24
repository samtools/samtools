#!/bin/bash
# =============================================================================
# CUDA GPU Detection and Capability Assessment Script for Samtools
# =============================================================================
#
# This script provides comprehensive CUDA environment detection and helps users
# determine if their system is capable of running samtools with GPU acceleration.
#
# Features:
# - Hardware detection (GPU devices, memory, compute capability)  
# - Software detection (CUDA toolkit, drivers, compilers)
# - Performance assessment and recommendations
# - Clear diagnostic output with next steps
#
# Author: Akshay Dedaniya <Dedaniya08@hotmail.com>
# Copyright (C) 2025 - MIT License

echo "==================================================================="
echo "Samtools CUDA GPU Detection & Capability Assessment"
echo "==================================================================="
echo ""

# =============================================================================
# Output Formatting and Utility Functions
# =============================================================================

# Color codes for enhanced readability
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Function to print status with consistent formatting
print_status() {
    local status=$1
    local message=$2
    case "$status" in
        "OK")
            echo -e "${GREEN}[✓]${NC} $message" ;;
        "WARN")
            echo -e "${YELLOW}[!]${NC} $message" ;;
        "ERROR")
            echo -e "${RED}[✗]${NC} $message" ;;
        "INFO")
            echo -e "${BLUE}[i]${NC} $message" ;;
        "STEP")
            echo -e "${CYAN}[→]${NC} $message" ;;
        *)
            echo -e "${NC}    $message" ;;
    esac
}

# Check 1: NVIDIA GPU Hardware Detection
echo "1. Checking for NVIDIA GPU hardware..."
echo ""

gpu_found=0

# Method 1: nvidia-smi
if command -v nvidia-smi >/dev/null 2>&1; then
    print_status "OK" "nvidia-smi found"
    
    gpu_count=$(nvidia-smi -L 2>/dev/null | wc -l)
    if [ "$gpu_count" -gt 0 ]; then
        print_status "OK" "Found $gpu_count NVIDIA GPU(s):"
        nvidia-smi -L 2>/dev/null | sed 's/^/    /'
        gpu_found=1
        
        # Show GPU details
        echo ""
        print_status "INFO" "GPU details:"
        nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv,noheader,nounits 2>/dev/null | while read line; do
            echo "    $line"
        done
    else
        print_status "ERROR" "nvidia-smi found but no GPUs detected"
    fi
else
    print_status "WARN" "nvidia-smi not found"
fi

# Method 2: Check /proc/driver/nvidia
echo ""
if [ -d "/proc/driver/nvidia" ]; then
    print_status "OK" "NVIDIA kernel driver loaded"
    if [ -r "/proc/driver/nvidia/version" ]; then
        driver_version=$(cat /proc/driver/nvidia/version 2>/dev/null | head -1)
        print_status "INFO" "Driver: $driver_version"
        gpu_found=1
    fi
else
    print_status "WARN" "NVIDIA kernel driver not found in /proc/driver/nvidia"
fi

# Method 3: lspci
echo ""
if command -v lspci >/dev/null 2>&1; then
    nvidia_devices=$(lspci 2>/dev/null | grep -i "nvidia.*\(vga\|3d\|display\)" | wc -l)
    if [ "$nvidia_devices" -gt 0 ]; then
        print_status "OK" "Found $nvidia_devices NVIDIA device(s) via lspci:"
        lspci 2>/dev/null | grep -i "nvidia.*\(vga\|3d\|display\)" | sed 's/^/    /'
        gpu_found=1
    else
        print_status "WARN" "No NVIDIA devices found via lspci"
    fi
else
    print_status "WARN" "lspci command not available"
fi

echo ""
echo "==================================================================="

# Check 2: CUDA Software Stack
echo "2. Checking CUDA software stack..."
echo ""

cuda_available=0

# Check for nvcc
if command -v nvcc >/dev/null 2>&1; then
    print_status "OK" "CUDA compiler (nvcc) found"
    nvcc_version=$(nvcc --version 2>/dev/null | grep "release" | sed 's/.*release //' | sed 's/,.*//')
    print_status "INFO" "NVCC version: $nvcc_version"
    
    # Check CUDA installation path
    cuda_path=$(dirname $(dirname $(which nvcc)))
    print_status "INFO" "CUDA installation: $cuda_path"
    
    # Check for CUDA headers
    if [ -f "$cuda_path/include/cuda_runtime.h" ]; then
        print_status "OK" "CUDA headers found"
    else
        print_status "ERROR" "CUDA headers not found"
    fi
    
    # Check for CUDA libraries
    if [ -d "$cuda_path/lib64" ]; then
        lib_path="$cuda_path/lib64"
    elif [ -d "$cuda_path/lib" ]; then
        lib_path="$cuda_path/lib"
    else
        lib_path=""
    fi
    
    if [ -n "$lib_path" ] && [ -f "$lib_path/libcudart.so" ]; then
        print_status "OK" "CUDA runtime library found"
        cuda_available=1
    else
        print_status "ERROR" "CUDA runtime library not found"
    fi
else
    print_status "ERROR" "CUDA compiler (nvcc) not found"
    print_status "INFO" "Install CUDA toolkit from: https://developer.nvidia.com/cuda-downloads"
fi

echo ""
echo "==================================================================="

# Check 3: Test CUDA functionality (if both GPU and CUDA are available)
if [ "$gpu_found" -eq 1 ] && [ "$cuda_available" -eq 1 ]; then
    echo "3. Testing CUDA functionality..."
    echo ""
    
    # Create test program
    cat > /tmp/cuda_test.cu << 'EOF'
#include <cuda_runtime.h>
#include <stdio.h>

int main() {
    int deviceCount = 0;
    cudaError_t error = cudaGetDeviceCount(&deviceCount);
    
    if (error != cudaSuccess) {
        printf("CUDA Error: %s\n", cudaGetErrorString(error));
        return 1;
    }
    
    if (deviceCount == 0) {
        printf("No CUDA devices found\n");
        return 1;
    }
    
    printf("Found %d CUDA-capable device(s):\n", deviceCount);
    
    for (int i = 0; i < deviceCount; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        
        printf("Device %d: %s\n", i, prop.name);
        printf("  Compute capability: %d.%d\n", prop.major, prop.minor);
        printf("  Global memory: %.2f GB\n", (double)prop.totalGlobalMem / (1024*1024*1024));
        printf("  Multiprocessors: %d\n", prop.multiProcessorCount);
        printf("  Max threads per block: %d\n", prop.maxThreadsPerBlock);
    }
    
    return 0;
}
EOF
    
    # Try to compile and run
    if nvcc /tmp/cuda_test.cu -o /tmp/cuda_test >/dev/null 2>&1; then
        print_status "OK" "CUDA test program compiled successfully"
        
        if /tmp/cuda_test 2>/dev/null; then
            print_status "OK" "CUDA test program executed successfully"
            echo ""
            /tmp/cuda_test 2>/dev/null | sed 's/^/    /'
        else
            print_status "ERROR" "CUDA test program failed to execute"
            echo ""
            print_status "INFO" "Possible issues:"
            echo "    - NVIDIA driver not properly loaded"
            echo "    - Driver/CUDA version mismatch"
            echo "    - Insufficient permissions"
            echo "    - Try: sudo modprobe nvidia"
        fi
    else
        print_status "ERROR" "CUDA test program failed to compile"
    fi
    
    # Cleanup
    rm -f /tmp/cuda_test.cu /tmp/cuda_test
    
else
    echo "3. Skipping CUDA functionality test"
    echo ""
    if [ "$gpu_found" -eq 0 ]; then
        print_status "WARN" "No NVIDIA GPU detected"
    fi
    if [ "$cuda_available" -eq 0 ]; then
        print_status "WARN" "CUDA software not available"
    fi
fi

echo ""
echo "==================================================================="
echo "Summary and Recommendations:"
echo "==================================================================="
echo ""

if [ "$gpu_found" -eq 1 ] && [ "$cuda_available" -eq 1 ]; then
    print_status "OK" "System is ready for CUDA acceleration!"
    echo ""
    echo "To build samtools with CUDA support:"
    echo "  ./configure --enable-cuda"
    echo "  make"
    echo ""
elif [ "$gpu_found" -eq 1 ] && [ "$cuda_available" -eq 0 ]; then
    print_status "WARN" "NVIDIA GPU found but CUDA toolkit missing"
    echo ""
    echo "To enable CUDA acceleration:"
    echo "  1. Install CUDA toolkit: https://developer.nvidia.com/cuda-downloads"
    echo "  2. Add nvcc to your PATH"
    echo "  3. Run: ./configure --enable-cuda"
    echo ""
elif [ "$gpu_found" -eq 0 ]; then
    print_status "INFO" "No NVIDIA GPU detected"
    echo ""
    echo "This system will use CPU-only acceleration."
    echo "For GPU acceleration, you need an NVIDIA GPU."
    echo ""
    echo "To build samtools (CPU-only):"
    echo "  ./configure"
    echo "  make"
    echo ""
else
    print_status "WARN" "Unexpected configuration"
    echo ""
    echo "Please check your NVIDIA GPU and CUDA installation."
    echo ""
fi

echo "For more information, see:"
echo "  - README_CUDA.md"
echo "  - CUDA_IMPLEMENTATION.md"
echo ""
