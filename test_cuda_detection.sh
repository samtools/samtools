#!/bin/bash

# test_cuda_detection.sh - Test CUDA auto-detection functionality
# This script tests the CUDA auto-detection features that would be used
# by the configure script and runtime system.

echo "=== CUDA Auto-Detection Test ==="
echo ""

# Test 1: Check for CUDA compiler
echo "1. Testing CUDA compiler detection..."
if command -v nvcc >/dev/null 2>&1; then
    echo "   ✓ nvcc found: $(which nvcc)"
    nvcc_version=$(nvcc --version 2>/dev/null | grep "release" | sed 's/.*release //' | sed 's/,.*//')
    if [ -n "$nvcc_version" ]; then
        echo "   ✓ NVCC version: $nvcc_version"
    fi
    
    # Get CUDA installation path
    nvcc_path=$(which nvcc)
    cuda_root=$(dirname $(dirname $nvcc_path))
    echo "   ✓ CUDA root path: $cuda_root"
    
    # Check for includes and libraries
    if [ -d "$cuda_root/include" ]; then
        echo "   ✓ CUDA headers found: $cuda_root/include"
    else
        echo "   ✗ CUDA headers not found in $cuda_root/include"
    fi
    
    if [ -d "$cuda_root/lib64" ]; then
        echo "   ✓ CUDA libraries found: $cuda_root/lib64"
        lib_path="$cuda_root/lib64"
    elif [ -d "$cuda_root/lib" ]; then
        echo "   ✓ CUDA libraries found: $cuda_root/lib"
        lib_path="$cuda_root/lib"
    else
        echo "   ✗ CUDA libraries not found"
        lib_path=""
    fi
else
    echo "   ✗ nvcc not found"
    echo "   → CUDA compiler not available"
fi

echo ""

# Test 2: Check for nvidia-smi
echo "2. Testing GPU detection tool..."
if command -v nvidia-smi >/dev/null 2>&1; then
    echo "   ✓ nvidia-smi found: $(which nvidia-smi)"
    
    # Get driver version
    driver_version=$(nvidia-smi --query-gpu=driver_version --format=csv,noheader,nounits 2>/dev/null | head -1)
    if [ -n "$driver_version" ]; then
        echo "   ✓ Driver version: $driver_version"
    fi
    
    # Get GPU count
    gpu_count=$(nvidia-smi --list-gpus 2>/dev/null | wc -l)
    echo "   ✓ GPU count: $gpu_count"
    
    if [ "$gpu_count" -gt 0 ]; then
        echo "   ✓ GPU devices:"
        nvidia-smi --query-gpu=name,compute_cap,memory.total --format=csv,noheader 2>/dev/null | \
        while IFS=',' read -r name compute_cap memory; do
            name=$(echo "$name" | xargs)
            compute_cap=$(echo "$compute_cap" | xargs)
            memory=$(echo "$memory" | xargs)
            echo "     - $name (CC: $compute_cap, Memory: $memory)"
        done
    fi
else
    echo "   ✗ nvidia-smi not found"
    echo "   → GPU detection tool not available"
fi

echo ""

# Test 3: Simulate CUDA runtime test
echo "3. Testing CUDA runtime availability simulation..."
cat > test_cuda_runtime.c << 'EOF'
#include <stdio.h>

// Simulate what configure.ac would test
int main() {
    printf("Simulating CUDA runtime test...\n");
    
    // This would normally be:
    // int count;
    // cudaError_t err = cudaGetDeviceCount(&count);
    // return (err == cudaSuccess && count > 0) ? 0 : 1;
    
    printf("Would check: cudaGetDeviceCount(&count)\n");
    printf("Would return: success if count > 0\n");
    
    return 0;  // Simulate success for demonstration
}
EOF

echo "   ✓ Created test program: test_cuda_runtime.c"
if gcc test_cuda_runtime.c -o test_cuda_runtime 2>/dev/null; then
    echo "   ✓ Test program compiles successfully"
    ./test_cuda_runtime
    echo "   ✓ Test program runs successfully"
    rm -f test_cuda_runtime test_cuda_runtime.c
else
    echo "   ✗ Test program compilation failed"
    rm -f test_cuda_runtime.c
fi

echo ""

# Test 4: Check configure.ac CUDA detection logic
echo "4. Testing configure.ac CUDA detection logic..."
if [ -f "configure.ac" ]; then
    echo "   ✓ configure.ac found"
    
    # Check if CUDA support is in configure.ac
    if grep -q "CUDA" configure.ac; then
        echo "   ✓ CUDA support detected in configure.ac"
        
        # Show key CUDA-related sections
        echo "   → Key CUDA detection features:"
        grep -n "AC_ARG_ENABLE.*cuda" configure.ac | head -1 | sed 's/^/     /'
        grep -n "AC_PATH_PROG.*NVCC" configure.ac | head -1 | sed 's/^/     /'
        grep -n "AC_PATH_PROG.*NVIDIA_SMI" configure.ac | head -1 | sed 's/^/     /'
        grep -n "Auto-detect CUDA" configure.ac | head -1 | sed 's/^/     /'
    else
        echo "   ✗ No CUDA support found in configure.ac"
    fi
else
    echo "   ✗ configure.ac not found"
fi

echo ""

# Test 5: Check Makefile CUDA integration
echo "5. Testing Makefile CUDA integration..."
if [ -f "Makefile" ]; then
    echo "   ✓ Makefile found"
    
    if grep -q "CUDA" Makefile; then
        echo "   ✓ CUDA integration detected in Makefile"
        
        # Show CUDA variables
        echo "   → CUDA variables in Makefile:"
        grep "^CUDA_" Makefile | sed 's/^/     /'
        grep "ENABLE_CUDA" Makefile | sed 's/^/     /'
    else
        echo "   ✗ No CUDA integration found in Makefile"
    fi
else
    echo "   ✗ Makefile not found"
fi

echo ""

# Test 6: Check if config.mk.in has CUDA support
echo "6. Testing config.mk.in CUDA template..."
if [ -f "config.mk.in" ]; then
    echo "   ✓ config.mk.in found"
    
    if grep -q "CUDA" config.mk.in; then
        echo "   ✓ CUDA template variables found"
        echo "   → CUDA template variables:"
        grep "@.*CUDA.*@\|^CUDA_\|^NVCC" config.mk.in | sed 's/^/     /'
    else
        echo "   ✗ No CUDA template variables found"
    fi
else
    echo "   ✗ config.mk.in not found"
fi

echo ""

# Test 7: Demonstrate what would happen with CUDA
echo "7. Simulating CUDA auto-setup process..."
echo "   → If CUDA were available, the auto-setup would:"
echo "     1. Detect nvcc compiler automatically"
echo "     2. Find CUDA installation path from nvcc location"
echo "     3. Verify CUDA headers in \$CUDA_ROOT/include"
echo "     4. Verify CUDA libraries in \$CUDA_ROOT/lib64 or \$CUDA_ROOT/lib"
echo "     5. Test compile and run a simple CUDA program"
echo "     6. Auto-detect GPU architecture using nvidia-smi"
echo "     7. Set appropriate compiler flags and library paths"
echo "     8. Enable CUDA compilation in Makefile"

echo ""
echo "=== Test Summary ==="
if command -v nvcc >/dev/null 2>&1 && command -v nvidia-smi >/dev/null 2>&1; then
    echo "✓ CUDA auto-detection would work - all tools available"
else
    echo "ℹ CUDA auto-detection demonstrated - tools not available on this system"
    echo "  → Install CUDA toolkit to enable actual GPU acceleration"
fi

echo ""
echo "To test with actual CUDA:"
echo "1. Install CUDA toolkit from https://developer.nvidia.com/cuda-toolkit"
echo "2. Run: ./configure --enable-cuda"
echo "3. Run: make"
echo "4. Test: ./samtools sort --gpu input.bam"
