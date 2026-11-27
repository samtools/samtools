/*  test/fuzz/fuzz_flagstat.c -- Fuzz driver for samtools flagstat

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

int bam_flagstat(int argc, char *argv[]);

/*
Fuzzes samtools flagstat.  See fuzz_common.h for build and run details.

Link bam_stat.o.  No runtime files.

Control bytes:
 c[0] bits 0-1  -O default|json|tsv
      bits 2-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("flagstat");
    CHOOSE("-O", c[0] & 3, NULL, "default", "json", "tsv");
    ENDOPTS();
    ARG(in);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 1, FUZZ_PARSE, args, bam_flagstat);
}
