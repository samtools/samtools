/*  test/fuzz/fuzz_mpileup.c -- Fuzz driver for samtools mpileup

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

int bam_mpileup(int argc, char *argv[]);

/*
Fuzzes samtools mpileup.  See fuzz_common.h for build and run details.

Link bam_plcmd.o and its dependencies.  Output goes to stdout.  Runtime
files: ref.fa (+ .fai) for -f, regions.bed for -l, rg_exclude.txt (read
group IDs, one per line) for -G.

-r and -X need an index, and -b a file of filenames, so aren't fuzzed.
One control bit passes the same input twice, for multi-sample pileup.

Inputs must be coordinate sorted with short @SQ lines (for -a/-aa).

Control bytes:
 c[0] bit  0    -6
      bit  1    -A
      bit  2    -B
      bits 3-4  -C
      bits 5-6  -d
      bit  7    -E
 c[1] bit  0    -f ref.fa
      bit  1    -G rg_exclude.txt
      bit  2    -l regions.bed
      bits 3-4  -q
      bits 5-6  -Q
      bit  7    -R
 c[2] bits 0-1  --rf
      bits 2-3  --ff
      bit  4    -x
      bit  5    -a
      bit  6    -a (again; -aa)
      bit  7    -O
 c[3] bit  0    --output-BP-5
      bit  1    -s
      bit  2    --output-QNAME
      bits 3-4  --output-extra: none; fields and tags from c[5] and c[6];
                the same with an empty element; the same plus names that
                are neither fields nor two-character tags
      bits 5-6  --output-sep
      bit  7    --output-empty
 c[4] bit  0,1  --no-output-ins (once or twice)
      bit  2,3  --no-output-del (once or twice)
      bit  4    --no-output-ends
      bit  5    --reverse-del
      bit  6    -M
      bit  7    --no-output-ins-mods
 c[5] bit  0    a second copy of the input
      bits 1-7  --output-extra fields: QNAME FLAG RNAME POS MAPQ RNEXT PNEXT
 c[6] bits 0-7  --output-extra tags: RG(Z) NM(i) MD(Z) BX(Z) XB(B) XH(H)
                XA(A) XF(f), covering every aux type in the seed data
*/

// Builds the --output-extra list into buf from the control bytes
static void output_extra(char *buf, size_t len, const uint8_t *c) {
    static const char *fields[7] = {
        "QNAME", "FLAG", "RNAME", "POS", "MAPQ", "RNEXT", "PNEXT"
    };
    static const char *tags[8] = {
        "RG", "NM", "MD", "BX", "XB", "XH", "XA", "XF"
    };
    const char *sel[20];
    int i, n = 0, mode = (c[3] >> 3) & 3;

    for (i = 0; i < 7; i++)
        if (c[5] & (2 << i))
            sel[n++] = fields[i];
    for (i = 0; i < 8; i++)
        if (c[6] & (1 << i))
            sel[n++] = tags[i];
    if (mode == 2)
        sel[n++] = "";      // empty element, e.g. "FLAG,,NM"
    if (mode == 3) {
        sel[n++] = "XYZ";   // too long to be a tag
        sel[n++] = "Q";     // too short
    }

    // Rotate the order using the spare bits, as the order is significant
    int rot = n ? (c[6] ^ c[5]) % n : 0;
    buf[0] = '\0';
    for (i = 0; i < n; i++) {
        if (i)
            strncat(buf, ",", len - strlen(buf) - 1);
        strncat(buf, sel[(i + rot) % n], len - strlen(buf) - 1);
    }
}

static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("mpileup");

    if (c[0] & 0x01)   ARG("-6");
    if (c[0] & 0x02)   ARG("-A");
    if (c[0] & 0x04)   ARG("-B");
    CHOOSE("-C", (c[0] >> 3) & 3, NULL, "0", "50", "1000");
    CHOOSE("-d", (c[0] >> 5) & 3, NULL, "0", "1", "8000");
    if (c[0] & 0x80)   ARG("-E");

    if (c[1] & 0x01)   ARG2("-f", "ref.fa");
    if (c[1] & 0x02)   ARG2("-G", "rg_exclude.txt");
    if (c[1] & 0x04)   ARG2("-l", "regions.bed");
    CHOOSE("-q", (c[1] >> 3) & 3, NULL, "0", "20", "255");
    CHOOSE("-Q", (c[1] >> 5) & 3, NULL, "0", "13", "99");
    if (c[1] & 0x80)   ARG("-R");

    CHOOSE("--rf", c[2] & 3,        NULL, "READ2,REVERSE", "0x900", "0");
    CHOOSE("--ff", (c[2] >> 2) & 3, NULL, "0", "QCFAIL,DUP", "0xfff");
    if (c[2] & 0x10)   ARG("-x");
    if (c[2] & 0x20)   ARG("-a");
    if (c[2] & 0x40)   ARG("-a");
    if (c[2] & 0x80)   ARG("-O");

    if (c[3] & 0x01)   ARG("--output-BP-5");
    if (c[3] & 0x02)   ARG("-s");
    if (c[3] & 0x04)   ARG("--output-QNAME");
    if ((c[3] >> 3) & 3) {
        // Derived from control bytes, so marked as data; an empty list
        // (no bits set) is a deliberate edge case.
        static char extra[200];
        output_extra(extra, sizeof(extra), c);
        ARG_OPTDATA("--output-extra", extra);
    }
    CHOOSE("--output-sep",   (c[3] >> 5) & 3, NULL, ";", "\t", "|");
    if (c[3] & 0x80)   ARG2("--output-empty", "-");

    if (c[4] & 0x01)   ARG("--no-output-ins");
    if (c[4] & 0x02)   ARG("--no-output-ins");
    if (c[4] & 0x04)   ARG("--no-output-del");
    if (c[4] & 0x08)   ARG("--no-output-del");
    if (c[4] & 0x10)   ARG("--no-output-ends");
    if (c[4] & 0x20)   ARG("--reverse-del");
    if (c[4] & 0x40)   ARG("-M");
    if (c[4] & 0x80)   ARG("--no-output-ins-mods");

    ENDOPTS();
    ARG(in);
    if (c[5] & 0x01)   ARG(in);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 7, FUZZ_SORTED | FUZZ_SHORT_REFS,
                          args, bam_mpileup);
}
