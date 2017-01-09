samtools
========

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
requires some extra steps:

```sh
aclocal -I ../autoconf-archive/m4   # Cache configure subroutines
autoconf       # Generate the configure script, if needed
./configure    # Optional, needed for choosing optional functionality
make
make install
```

By default, this will build against an HTSlib source tree in `../htslib`.
You can alter this to a source tree elsewhere or to a previously-installed
HTSlib by configuring with `--with-htslib=DIR`.

The configure script uses several macros from the [autoconf-archive],
so you will need to point `aclocal` at a copy of the archive in order to
regenerate the script.

[autoconf-archive]: http://www.gnu.org/software/autoconf-archive/
[download]: http://www.htslib.org/download/
