/*  test/fuzz/fuzz_split.c -- Fuzz driver for samtools split

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

int main_split(int argc, char *argv[]);

/*
Fuzzes samtools split.  See fuzz_common.h for the input layout (file byte,
region bytes, extra payload) and build/run details.

Output filenames.  split names its outputs with a format string whose
expansions include %! (the @RG ID or -d tag value) and %* (the input's
basename), both of which come from the data.  The default with -d is
"%*_%!.%.".  So by default this fuzzer uses "-f /dev/null", which has no
expansions: every output, including -u, goes to /dev/null.  This keeps the
output paths unfuzzable, at the cost of not exercising the expansion code.

Building with -DFUZZ_SPLIT_EXPANSIONS instead uses "<fuzz_prefix>_%#%%.%.",
whose only data-dependent part is %#, a decimal number (the @RG index or
tag-value count), plus %% and %.; files land under fuzz_prefix and are
removed after each run.  Only enable this if you accept digit-only,
data-dependent filenames.

The extra payload serves as the -h header file for -u output.

Link bam_split.o.  No runtime files.

Tool bytes:
 c[0] bits 0-1  -d: none, RG, BX, NM
      bit  2    -u /dev/null
      bit  3    -h <extra payload as a header file>
      bits 4-5  -M
      bit  6    -v
      bit  7    --no-PG
 c[1] bits 0-1  --output-fmt
      bits 2-3  -p (only significant with FUZZ_SPLIT_EXPANSIONS)
      bits 4-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("split");

#ifdef FUZZ_SPLIT_EXPANSIONS
    static char fmt[96];
    snprintf(fmt, sizeof(fmt), "%s_%%#%%%%.%%.", fuzz_prefix);
    ARG2("-f", fmt);
#else
    ARG2("-f", "/dev/null");
#endif

    CHOOSE("-d", c[0] & 3, NULL, "RG", "BX", "NM");
    if (c[0] & 0x04)   ARG2("-u", "/dev/null");
    if (c[0] & 0x08)   ARG2("-h", f->extra_raw_fn);
    CHOOSE("-M", (c[0] >> 4) & 3, NULL, "1", "0", "-1");
    if (c[0] & 0x40)   ARG("-v");
    if (c[0] & 0x80)   ARG("--no-PG");

    CHOOSE("--output-fmt", c[1] & 3, NULL, "sam", "bam", "cram,no_ref");
    CHOOSE("-p", (c[1] >> 2) & 3, NULL, "0", "3", "20");

    ENDOPTS();
    ARG(f->fn);
}

static const fuzz_step steps[] = {{args, main_split}};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 2, FUZZ_PARSE, FUZZ_EXTRA,
                               steps, 1);
}
