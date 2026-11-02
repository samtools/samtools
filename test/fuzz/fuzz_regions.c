/*  test/fuzz/fuzz_regions.c -- Fuzz driver for samtools view / depth / mpileup / consensus / stats region queries

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

int main_samview(int argc, char *argv[]);
int main_depth(int argc, char *argv[]);
int bam_mpileup(int argc, char *argv[]);
int main_consensus(int argc, char *argv[]);
int main_stats(int argc, char *argv[]);

/*
Fuzzes the indexed region-query paths of view, depth, mpileup, consensus
and stats: region arguments, -r, the multi-region iterator (view -M),
fetch-pairs (view -P), index-backed BED filtering, and explicit index
locations (-X and "file##idx##index").  The general options of these
commands are covered by their own fuzzers; this one concentrates on regions.

See fuzz_common.h for the input layout (file byte, region bytes) and
build/run details.  Region strings are built from the reference names in
the fuzzed header, so most of them hit real references.

Link sam_view.o, bam2depth.o, bam_plcmd.o, consensus.o, stats.o and their
dependencies.  Output goes to stdout.  No runtime files.

Inputs must be coordinate sorted, with short @SQ lines (for -a).

Tool byte:
 c[0] bits 0-2  command:
                  0 view REGION...
                  1 view -M REGION...
                  2 view -L <generated BED> (multi-region via the index)
                  3 view -P REGION... (fetch pairs; implies -M)
                  4 depth -r REGION, or -b BED
                  5 mpileup -r REGION, or -l BED
                  6 consensus -r REGION
                  7 stats FILE REGION...
      bit  3    give the index explicitly with -X (except consensus)
      bit  4    view: -c;  depth/mpileup/consensus: -a
      bit  5    view: -b;  depth/mpileup: use the BED file, not -r;
                consensus: -f pileup;  stats: -t <generated BED>
      bit  6    depth/mpileup: a second copy of the input
      bit  7    view: --unmap;  depth/mpileup/consensus: second -a (-aa)
*/

static void args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    int cmd = c[0] & 7, x = (c[0] & 0x08) && cmd != 6;
    int b4 = c[0] & 0x10, b5 = c[0] & 0x20, b6 = c[0] & 0x40, b7 = c[0] & 0x80;
    int i, positional_regions = 0;

    switch (cmd) {
    case 0: case 1: case 2: case 3:
        ARG("view");
        if (cmd == 1)  ARG("-M");
        if (cmd == 2)  ARG2("-L", f->bed_fn);
        if (cmd == 3)  ARG("-P");
        if (b4)        ARG("-c");
        if (b5)        ARG("-b");
        if (b7)        ARG("--unmap");
        positional_regions = cmd != 2;
        break;

    case 4: case 5:
        ARG(cmd == 4 ? "depth" : "mpileup");
        if (b4)        ARG("-a");
        if (b7)        ARG("-a");
        if (b5)        ARG2(cmd == 4 ? "-b" : "-l", f->bed_fn);
        else           ARG_OPTDATA("-r", f->reg[0].s);
        break;

    case 6:
        ARG("consensus");
        if (b4)        ARG("-a");
        if (b7)        ARG("-a");
        if (b5)        ARG2("-f", "pileup");
        ARG_OPTDATA("-r", f->reg[0].s);
        break;

    case 7:
        ARG("stats");
        if (b5)        ARG2("-t", f->bed_fn);
        positional_regions = 1;
        break;
    }

    if (x)  ARG("-X");
    ENDOPTS();

    // Inputs: with -X the alignment file(s), then the index file(s)
    int two = (cmd == 4 || cmd == 5) && b6;
    ARG(x ? f->fn : f->in);
    if (two)  ARG(x ? f->fn : f->in);
    if (x) {
        ARG(f->idx_path);
        if (two)  ARG(f->idx_path);
    }

    if (positional_regions)
        for (i = 0; i < f->nreg; i++)
            ARG_DATA(f->reg[i].s);
}

static const fuzz_step steps[1] = {{args, NULL}};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    // The command is chosen by the tool byte, so pick the main function
    // to match before handing over to the driver.
    static fuzz_main_fn mains[8] = {
        main_samview, main_samview, main_samview, main_samview,
        main_depth, bam_mpileup, main_consensus, main_stats
    };
    fuzz_step s = steps[0];
    s.main_fn = size > 1 ? mains[data[1] & 7] : main_samview;
    return fuzz_one_file_input(data, size, 1, FUZZ_SORTED | FUZZ_SHORT_REFS,
                               1, &s, 1);
}
