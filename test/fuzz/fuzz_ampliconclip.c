/*  test/fuzz/fuzz_ampliconclip.c -- Fuzz driver for samtools ampliconclip

    Copyright (C) 2026 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

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

#include "fuzz_common.h"

int amplicon_clip_main(int argc, char *argv[]);

/*
Fuzzes samtools ampliconclip.  See fuzz_common.h for the input layout
(file byte, region bytes, extra payload) and build/run details.

The primer BED is the extra payload, so every input is self-contained.  All secondary outputs (--rejects-file, --primer-counts, -f)
go to /dev/null.

Link bam_ampliconclip.o.  Output is "-o -".  No runtime files.

Tool bytes:
 c[0] bit  0    --hard-clip (else --soft-clip)
      bit  1    --both-ends
      bit  2    --strand
      bit  3    --clipped
      bit  4    --fail
      bit  5    --no-excluded
      bit  6    --original
      bit  7    --keep-tag
 c[1] bits 0-1  --filter-len
      bits 2-3  --fail-len
      bits 4-5  --unmap-len
      bits 6-7  --tolerance
 c[2] bit  0    unused
      bit  1    --rejects-file /dev/null
      bit  2    --primer-counts /dev/null
      bit  3    -f /dev/null
      bit  4    -u
      bits 5-6  -O
      bit  7    --no-PG
*/

static void args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("ampliconclip");
    ARG2("-o", "-");
    ARG2("-b", f->extra_raw_fn);

    ARG((c[0] & 0x01) ? "--hard-clip" : "--soft-clip");
    if (c[0] & 0x02)   ARG("--both-ends");
    if (c[0] & 0x04)   ARG("--strand");
    if (c[0] & 0x08)   ARG("--clipped");
    if (c[0] & 0x10)   ARG("--fail");
    if (c[0] & 0x20)   ARG("--no-excluded");
    if (c[0] & 0x40)   ARG("--original");
    if (c[0] & 0x80)   ARG("--keep-tag");

    CHOOSE("--filter-len", c[1] & 3,        NULL, "0", "10", "1000");
    CHOOSE("--fail-len",   (c[1] >> 2) & 3, NULL, "0", "10", "1000");
    CHOOSE("--unmap-len",  (c[1] >> 4) & 3, NULL, "-1", "0", "10");
    CHOOSE("--tolerance",  (c[1] >> 6) & 3, NULL, "0", "1", "1000");

    if (c[2] & 0x02)   ARG2("--rejects-file", "/dev/null");
    if (c[2] & 0x04)   ARG2("--primer-counts", "/dev/null");
    if (c[2] & 0x08)   ARG2("-f", "/dev/null");
    if (c[2] & 0x10)   ARG("-u");
    CHOOSE("-O", (c[2] >> 5) & 3, NULL, "sam", "bam", "cram,no_ref");
    if (c[2] & 0x80)   ARG("--no-PG");

    ENDOPTS();
    ARG(f->in);
}

static const fuzz_step steps[] = {{args, amplicon_clip_main}};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 3, FUZZ_PARSE, FUZZ_EXTRA,
                               steps, 1);
}
