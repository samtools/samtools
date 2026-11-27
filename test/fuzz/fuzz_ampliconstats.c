/*  test/fuzz/fuzz_ampliconstats.c -- Fuzz driver for samtools ampliconstats

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

int main_ampliconstats(int argc, char *argv[]);

/*
Fuzzes samtools ampliconstats.  See fuzz_common.h for the input layout
(file byte, region bytes, extra payload) and build/run details.

The primer BED is the extra payload, so every input is self-contained.

Link amplicon_stats.o.  Output goes to stdout.  No runtime files.

Tool bytes:
 c[0] bits 0-1  -f
      bits 2-3  -F
      bits 4-5  -a
      bits 6-7  -l
 c[1] bits 0-1  -b
      bits 2-3  -c
      bits 4-5  -D
      bits 6-7  -d
 c[2] bits 0-1  -m
      bit  2    -s
      bit  3    -S
      bits 4-5  -t
      bit  6    unused
      bit  7    a second copy of the input
*/

static void args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("ampliconstats");

    CHOOSE("-f", c[0] & 3,        NULL, "0x1", "PAIRED,PROPER_PAIR", "0xfff");
    CHOOSE("-F", (c[0] >> 2) & 3, NULL, "0", "0xb04", "UNMAP");
    CHOOSE("-a", (c[0] >> 4) & 3, NULL, "1", "10", "1000");
    CHOOSE("-l", (c[0] >> 6) & 3, NULL, "0", "100", "10000");

    CHOOSE("-b", c[1] & 3,        NULL, "1", "10", "0");
    CHOOSE("-c", (c[1] >> 2) & 3, NULL, "0", "1", "100");
    CHOOSE("-D", (c[1] >> 4) & 3, NULL, "0", "0.01", "1");
    CHOOSE("-d", (c[1] >> 6) & 3, NULL, "1", "0,1,10", "1,10,100,1000,10000");

    CHOOSE("-m", c[2] & 3,        NULL, "0", "5", "-1");
    if (c[2] & 0x04)   ARG("-s");
    if (c[2] & 0x08)   ARG("-S");
    CHOOSE("-t", (c[2] >> 4) & 3, NULL, "0", "5", "-5");

    ENDOPTS();
    ARG(f->extra_raw_fn);
    ARG(f->in);
    if (c[2] & 0x80)   ARG(f->in);
}

static const fuzz_step steps[] = {{args, main_ampliconstats}};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 3, FUZZ_PARSE, FUZZ_EXTRA,
                               steps, 1);
}
