/*  test/fuzz/fuzz_markdup.c -- Fuzz driver for samtools markdup

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

int bam_markdup(int argc, char *argv[]);

/*
Fuzzes samtools markdup

See fuzz_common.h for build and run details, including FUZZ_EXTRACT for
reproducing a crash outside the harness.

The leading control bytes of each input select the command line options.
The remainder is the alignment data, passed to samtools as a data: URL.

Files:
  - Output always goes to "-" (stdout), discarded by -close_fd_mask=1.
  - Stats via -f go to /dev/null, so the writing code runs but no file
    is created.
  - markdup uses temporary files internally, which can't be avoided.  We
    always pass -T fuzz_prefix (_fz.<pid>) so they land in the current directory under
    a known prefix, and remove any leftovers before and after each run
    (they survive if samtools exit()s mid-run).  Run from a scratch dir.

No other runtime files are needed.

The corpus must be coordinate sorted and carry the MC and ms tags that
fixmate -m adds; inputs without them only reach markdup's error paths.
Seeds should include Illumina-style read names (for -d optical duplicates
and --read-coords), UMIs as an 8th name field (--barcode-name), BX tags,
several read groups, existing dup flags and dt/do tags (for -c), and
supplementary/secondary reads of duplicates (for -S).

No validity of result is assumed.  We're simply checking for crashes or
illegal program behaviour.
*/

// --read-coords variants: an Illumina-style name with a lane:tile prefix
// to compare (txy), bare x:y at the start of the name (xy), and a single
// capture, which is too few.
#define RC_ILLUMINA "([^:]+:[0-9]+:[^:]+:[0-9]+:[0-9]+):([0-9]+):([0-9]+)"
#define RC_XY       "^([[:digit:]]+):([[:digit:]]+)"
#define RC_ONE      "([0-9]+)"

// --barcode-rgx: the 8th colon-separated field of the read name
#define BC_RGX "[0-9A-Za-z]+:[0-9]+:[0-9A-Za-z]+:[0-9]+:[0-9]+:[0-9]+:[0-9]+:([!-?A-~]+)"

/*
 * Control bytes:
 *
 * c[0] bit  0    -r
 *      bit  1    -S
 *      bit  2    -c
 *      bit  3    -t
 *      bit  4    --duplicate-count
 *      bit  5    -m s
 *      bit  6    --include-fails
 *      bit  7    --no-multi-dup
 * c[1] bits 0-1  -d
 *      bits 2-3  -l
 *      bit  4    -s
 *      bit  5    --json
 *      bit  6    -f /dev/null
 *      bit  7    -u
 * c[2] bits 0-1  --read-coords
 *      bits 2-4  --coords-order
 *      bit  5    --use-read-groups
 *      bit  6    --no-PG
 *      bit  7    unused
 * c[3] bits 0-1  barcode source: --barcode-tag BX, --barcode-name,
 *                or --barcode-rgx
 *      bit  2    --move-umi-to-tag
 *      bits 3-4  --umi-separator
 *      bits 5-7  unused
 */
static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("markdup");
    ARG2("-T", fuzz_prefix);

    if (c[0] & 0x01)   ARG("-r");
    if (c[0] & 0x02)   ARG("-S");
    if (c[0] & 0x04)   ARG("-c");
    if (c[0] & 0x08)   ARG("-t");
    if (c[0] & 0x10)   ARG("--duplicate-count");
    if (c[0] & 0x20) ARG2("-m", "s");
    if (c[0] & 0x40)   ARG("--include-fails");
    if (c[0] & 0x80)   ARG("--no-multi-dup");

    CHOOSE("-d", (c[1]) & 3,    NULL, "1", "100", "2500");
    CHOOSE("-l", (c[1] >> 2) & 3, NULL, "1", "10",  "100000");
    if (c[1] & 0x10)   ARG("-s");
    if (c[1] & 0x20)   ARG("--json");
    if (c[1] & 0x40) ARG2("-f", "/dev/null");
    if (c[1] & 0x80)   ARG("-u");

    CHOOSE("--read-coords", (c[2]) & 3, NULL, RC_ILLUMINA, RC_XY, RC_ONE);
    static const char *order[8] = {
        NULL, "txy", "tyx", "xyt", "yxt", "xty", "xy", "yx"
    };
    if (order[(c[2] >> 2) & 7]) ARG2("--coords-order", order[(c[2] >> 2) & 7]);
    if (c[2] & 0x20)   ARG("--use-read-groups");
    if (c[2] & 0x40)   ARG("--no-PG");

    switch (c[3] & 3) {
    case 1: ARG2("--barcode-tag", "BX"); break;
    case 2: ARG("--barcode-name"); break;
    case 3: ARG2("--barcode-rgx", BC_RGX); break;
    }
    if (c[3] & 0x04)   ARG("--move-umi-to-tag");
    CHOOSE("--umi-separator", (c[3] >> 3) & 3, NULL, ":", "_", "+");

    // Threads make runs non-deterministic.
    //ARG("-@"); ARG("2");

    ENDOPTS();
    ARG(in);
    ARG("-");   // output to stdout; never write a file

}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 4, FUZZ_SORTED, args, bam_markdup);
}
