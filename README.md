samtools
========

[![Build Status](https://api.cirrus-ci.com/github/samtools/samtools.svg?branch=develop)](https://cirrus-ci.com/github/samtools/samtools)
[![Build status](https://github.com/samtools/samtools/actions/workflows/windows-build.yml/badge.svg)](https://github.com/samtools/samtools/actions/workflows/windows-build.yml?query=branch%3Adevelop)
[![Github All Releases](https://img.shields.io/github/downloads/samtools/samtools/total.svg)](https://github.com/samtools/samtools/releases/latest)

This is the official development repository for samtools.

The original samtools package has been split into three separate
but tightly coordinated projects:
- [htslib](https://github.com/samtools/htslib): C-library for handling high-throughput sequencing data
- samtools: mpileup and other tools for handling SAM, BAM, CRAM
- [bcftools](https://github.com/samtools/bcftools): calling and other tools for handling VCF, BCF

See also http://github.com/samtools/

### Building Samtools

See [INSTALL](INSTALL) for complete details.
[Release tarballs][download] contain generated files that have not been
committed to this repository, so building the code from a Git repository
requires extra steps:

```sh
autoheader            # Build config.h.in (this may generate a warning about
                      # AC_CONFIG_SUBDIRS - please ignore it).
autoconf -Wno-syntax  # Generate the configure script
./configure           # Needed for choosing optional functionality
make
make install
```

By default, this will build against an HTSlib source tree in `../htslib`.
You can alter this to a source tree elsewhere or to a previously-installed
HTSlib by configuring with `--with-htslib=DIR`.

[download]: http://www.htslib.org/download/

### Citing

Please cite this paper when using SAMtools for your publications.

> Twelve years of SAMtools and BCFtools </br>
> Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li </br>
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
    abstract = "{SAMtools and BCFtools are widely used programs for processing and analysing high-throughput sequencing data. They include tools for file format conversion and manipulation, sorting, querying, statistics, variant calling, and effect analysis amongst other methods.The first version appeared online 12 years ago and has been maintained and further developed ever since, with many new features and improvements added over the years. The SAMtools and BCFtools packages represent a unique collection of tools that have been used in numerous other software projects and countless genomic pipelines.Both SAMtools and BCFtools are freely available on GitHub under the permissive MIT licence, free for both non-commercial and commercial use. Both packages have been installed \\&gt;1 million times via Bioconda. The source code and documentation are available from https://www.htslib.org.}",
    issn = {2047-217X},
    doi = {10.1093/gigascience/giab008},
    url = {https://doi.org/10.1093/gigascience/giab008},
    note = {giab008},
    eprint = {https://academic.oup.com/gigascience/article-pdf/10/2/giab008/36332246/giab008.pdf},
}
```
