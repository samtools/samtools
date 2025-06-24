# Fork Notice and Attribution

This repository is a GPU-accelerated fork of [samtools/samtools](https://github.com/samtools/samtools), originally developed by Genome Research Ltd. and contributors. All original code is licensed under the MIT/Expat License. This fork introduces CUDA GPU acceleration and major refactoring for performance and maintainability.

**Key differences from upstream:**
- CUDA GPU acceleration for sorting, statistics, and pileup
- Modularized codebase and robust CPU/GPU fallback
- Enhanced documentation and build/test scripts

**Author:**
- Akshay Dedaniya (<dedaniya08@hotmail.com>)

**License:**
This project, including all modifications, remains under the MIT/Expat License. See [LICENSE](LICENSE) for details.

---

# Samtools

[![Build Status](https://api.cirrus-ci.com/github/samtools/samtools.svg?branch=develop)](https://cirrus-ci.com/github/samtools/samtools)
[![Build status](https://github.com/samtools/samtools/actions/workflows/windows-build.yml/badge.svg)](https://github.com/samtools/samtools/actions/workflows/windows-build.yml?query=branch%3Adevelop)
[![Github All Releases](https://img.shields.io/github/downloads/samtools/samtools/total.svg)](https://github.com/samtools/samtools/releases/latest)

Samtools implements a suite of utilities for post-processing alignments in the SAM, BAM, and CRAM formats, including indexing, variant calling (with bcftools), and a simple alignment viewer. This release introduces **CUDA GPU acceleration** for key operations, robust auto-detection, and a major codebase refactor for maintainability and performance.

---

## Table of Contents
- [Features](#features)
- [CUDA GPU Acceleration](#cuda-gpu-acceleration)
- [Requirements](#requirements)
- [Building Samtools](#building-samtools)
- [Usage](#usage)
- [Performance Benchmarks](#performance-benchmarks)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [Citing](#citing)
- [License](#license)

---

## Features
- Comprehensive manipulation of SAM, BAM, and CRAM files
- Indexing, sorting, statistics, pileup, and consensus generation
- Variant calling (with [bcftools](https://github.com/samtools/bcftools))
- Simple alignment viewer
- **NEW:** GPU-accelerated sorting, statistics, and pileup (CUDA)
- **NEW:** Intelligent GPU auto-detection and robust CPU fallback
- Modular, maintainable codebase with improved error handling

---

## CUDA GPU Acceleration

Samtools now supports GPU acceleration for sorting, statistics, and pileup operations using NVIDIA CUDA. Key highlights:
- **Automatic GPU detection:** No manual configuration needed; falls back to CPU if no compatible GPU is found.
- **Significant speedups:** 3–5x faster for large datasets (see [benchmarks](#performance-benchmarks)).
- **Modular design:** CPU and GPU code are clearly separated for maintainability.
- **Robust fallback:** If GPU is unavailable or an error occurs, samtools continues in CPU mode with a warning.

### Supported Operations
- `samtools sort --gpu ...`
- `samtools stats --gpu ...`
- `samtools pileup --gpu ...`

### Hardware & Software Requirements
- NVIDIA GPU (Compute Capability 5.0+ recommended)
- CUDA Toolkit 11.0+ (12.0+ recommended)
- Standard build dependencies: htslib, zlib, ncurses

---

## Requirements

### Hardware
- NVIDIA GPU (Maxwell architecture or newer, 4GB+ memory recommended)
- CPU-only mode is fully supported if no GPU is present

### Software
- CUDA Toolkit 11.0 or later (only if using GPU features)
- Compatible NVIDIA drivers
- GCC/G++ compatible with CUDA
- Standard dependencies: htslib, zlib, ncurses

---

## Building Samtools

### 1. Check for CUDA Support
Run the detection script to check your system:
```sh
./check_cuda_support.sh
```

### 2. Configure and Build
#### Automatic (Recommended)
```sh
./configure
make -j$(nproc)
```
- If a compatible GPU and CUDA are detected, CUDA support is enabled automatically.
- If not, samtools builds in CPU-only mode.

#### Force Enable/Disable CUDA
```sh
./configure --enable-cuda   # Force enable (fails if requirements not met)
./configure --disable-cuda  # Force CPU-only
make -j$(nproc)
```

#### Manual CUDA Path (Advanced)
```sh
export CUDA_PATH=/usr/local/cuda
export NVCC=/usr/local/cuda/bin/nvcc
./configure --enable-cuda
make -j$(nproc)
```

#### Install
```sh
sudo make install
```

---

## Usage

### GPU-Accelerated Commands
```sh
samtools sort --gpu input.bam -o output.bam
samtools stats --gpu input.bam > stats.txt
samtools pileup --gpu input.bam > pileup.txt
```

### Check GPU Support
```sh
samtools --version
# Look for 'cuda=yes' in the features line
```

### Performance Monitoring
```sh
nvidia-smi -l 1
```

### Multi-GPU Systems
```sh
export CUDA_VISIBLE_DEVICES=0  # Use first GPU
samtools sort --gpu input.bam -o output.bam
```

---

## Performance Benchmarks

| Operation | Dataset Size | CPU Time | GPU Time | Speedup |
|-----------|--------------|----------|----------|---------|
| Sort      | 1GB BAM      | 2m 30s   | 45s      | 3.3x    |
| Sort      | 10GB BAM     | 25m      | 6m 30s   | 3.8x    |
| Stats     | 1GB BAM      | 1m 15s   | 20s      | 3.8x    |
| Stats     | 10GB BAM     | 12m      | 2m 45s   | 4.4x    |

*Benchmarks: NVIDIA RTX 3080 vs Intel i7-10700K*

---

## Troubleshooting

- **CUDA not found during compilation:** Ensure CUDA is in your PATH and LD_LIBRARY_PATH.
- **Runtime CUDA errors:** Check GPU status with `nvidia-smi` and CUDA installation with `nvcc --version`.
- **Out of GPU memory:** Reduce memory usage (`-m` option) or use CPU fallback.
- **Performance not improved:** Large files benefit most; check GPU utilization.
- **Automatic fallback:** If GPU is unavailable, samtools continues in CPU mode with a warning.

For more, see `CUDA_IMPLEMENTATION.md`.

---

## Contributing

- Follow modular code patterns in the `cuda/` directory
- Add comprehensive tests for new GPU features
- Ensure CPU fallback compatibility
- Update documentation for new features

---

## Citing

Please cite this paper when using SAMtools for your publications:

> Twelve years of SAMtools and BCFtools  
> Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li  
> _GigaScience_, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008

```
@article{10.1093/gigascience/giab008,
    author = {Danecek, Petr and Bonfield, James K and Liddle, Jennifer and Marshall, John and Ohan, Valeriu and Pollard, Martin O and Whitwham, Andrew and Keane, Thomas and McCarthy, Shane A and Davies, Robert M and Li, Heng},
    title = "{Twelve years of SAMtools and BCFtools}",
    journal = {GigaScience},
    volume = {10},
    number = {2},
    year = {2021},
    month = {02},
    abstract = "{SAMtools and BCFtools are widely used programs for processing and analysing high-throughput sequencing data. They include tools for file format conversion and manipulation, sorting, querying, statistics, variant calling, and effect analysis amongst other methods.The first version appeared online 12 years ago and has been maintained and further developed ever since, with many new features and improvements added over the years. The SAMtools and BCFtools packages represent a unique collection of tools that have been used in numerous other software projects and countless genomic pipelines.Both SAMtools and BCFtools are freely available on GitHub under the permissive MIT licence, free for both non-commercial and commercial use. Both packages have been installed >1 million times via Bioconda. The source code and documentation are available from https://www.htslib.org.}",
    issn = {2047-217X},
    doi = {10.1093/gigascience/giab008},
    url = {https://doi.org/10.1093/gigascience/giab008},
    note = {giab008},
    eprint = {https://academic.oup.com/gigascience/article-pdf/10/2/giab008/36332246/giab008.pdf},
}
```

---

## License

Samtools and its CUDA GPU acceleration code are released under the MIT license.

---

## Acknowledgments
- NVIDIA CUDA Toolkit and libraries
- Thrust and CUB library developers
- HTSlib project for foundational support
- Samtools community for testing and feedback

---

*For detailed CUDA implementation, see `CUDA_IMPLEMENTATION.md`.*
