/*  cram_size.c -- produces summary of the size of each cram data-series

    Copyright (C) 2023 Genome Research Ltd.

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

// TODO: add range query.  Eg the ability to look at size for "*" only
// (unmapped), or in a specific region such as a centromere.

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>

#include "htslib/bgzf.h"
#include "htslib/sam.h"
#include "htslib/cram.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "samtools.h"
#include "sam_opts.h"
#include "htslib/hfile.h"

/*----------------------------------------------------------------------
 * Compression method handling
 */

// A numeric version of the cram_method_details struct.
// We expand the myriad of struct field combinations into a single
// enumerated type so we can index and accumulate statistics for
// purposes of reporting.
//
// These expanded numeric values have no definition within CRAM itself
// and never occur within the file format.
enum comp_expanded {
    //----
    // Copies from htslib cram_block_method enum
    COMP_RAW    = CRAM_COMP_RAW,
    COMP_GZIP   = CRAM_COMP_GZIP,
    COMP_BZIP2  = CRAM_COMP_BZIP2,
    COMP_LZMA   = CRAM_COMP_LZMA,
    COMP_RANS8  = CRAM_COMP_RANS4x8,
    COMP_RANS16 = CRAM_COMP_RANSNx16,
    COMP_ARITH  = CRAM_COMP_ARITH,
    COMP_FQZ    = CRAM_COMP_FQZ,
    COMP_TOK3   = CRAM_COMP_TOK3,

    //----
    // Localised variants.

    // Gzip
    COMP_GZIP_1,
    COMP_GZIP_9,

    // Bzip2
    COMP_BZIP2_1,
    COMP_BZIP2_2,
    COMP_BZIP2_3,
    COMP_BZIP2_4,
    COMP_BZIP2_5,
    COMP_BZIP2_6,
    COMP_BZIP2_7,
    COMP_BZIP2_8,
    COMP_BZIP2_9,

    // rans 4x8
    COMP_RANS4x8_O0,
    COMP_RANS4x8_O1,

    // rans Nx16.  Note order here is to enable selection via bit-fields
    // bit 0: O0/O1
    // bit 1: RLE
    // bit 2: PACK
    // bit 3: 32x16
    COMP_RANS4x16_O0,
    COMP_RANS4x16_O1,
    COMP_RANS4x16_O0R,   // +RLE
    COMP_RANS4x16_O1R,
    COMP_RANS4x16_O0P,   // +PACK
    COMP_RANS4x16_O1P,
    COMP_RANS4x16_O0PR,  // +PACK+RLE
    COMP_RANS4x16_O1PR,
    COMP_RANS32x16_O0,   // SIMD variants
    COMP_RANS32x16_O1,
    COMP_RANS32x16_O0R,  // +RLE
    COMP_RANS32x16_O1R,
    COMP_RANS32x16_O0P,  // +PACK
    COMP_RANS32x16_O1P,
    COMP_RANS32x16_O0PR, // +PACK+RLE
    COMP_RANS32x16_O1PR,
    COMP_RANSNx16_STRIPE,
    COMP_RANSNx16_CAT,

    // Arith
    COMP_ARITH_O0,
    COMP_ARITH_O1,
    COMP_ARITH_O0R,   // +RLE
    COMP_ARITH_O1R,
    COMP_ARITH_O0P,   // +PACK
    COMP_ARITH_O1P,
    COMP_ARITH_O0PR,  // +PACK+RLE
    COMP_ARITH_O1PR,
    COMP_ARITH_STRIPE,
    COMP_ARITH_CAT,   // no entropy encoder
    COMP_ARITH_EXT,   // external entropy encode

    // Nake tokeniser
    COMP_TOK3_RANS,
    COMP_TOK3_ARITH,

    // To mark maximum size
    COMP_MAX,
};

static enum comp_expanded comp_method2expanded(cram_method_details *cm) {
    switch (cm->method) {
    case CRAM_COMP_GZIP:
        switch (cm->level) {
        case 1:  return COMP_GZIP_1;
        case 9:  return COMP_GZIP_9;
        default: return COMP_GZIP;
        }
        break;

    case CRAM_COMP_BZIP2:
        if (cm->level >= 1 && cm->level <= 9)
            return COMP_BZIP2_1 + cm->level-1;
        else
            return COMP_BZIP2;
        break;

    case CRAM_COMP_RANS4x8:
        return cm->order ? COMP_RANS4x8_O1 : COMP_RANS4x8_O0;

    case CRAM_COMP_RANSNx16: {
        // 8 4x16, 8 32x16 and 2 stripe/cat
        if (cm->stripe) return COMP_RANSNx16_STRIPE;
        if (cm->cat)    return COMP_RANSNx16_CAT;
        int c = COMP_RANS4x16_O0;
        c += 1*cm->order;
        c += 2*cm->rle;
        c += 4*cm->pack;
        c += 8*(cm->Nway==32);
        return c;
    }

    case CRAM_COMP_ARITH: {
        // 8 4x16, 8 32x16 and 2 stripe/cat
        if (cm->stripe) return COMP_ARITH_STRIPE;
        if (cm->cat)    return COMP_ARITH_CAT;
        if (cm->ext)    return COMP_ARITH_EXT;
        int c = COMP_ARITH_O0;
        c += 1*cm->order;
        c += 2*cm->rle;
        c += 4*cm->pack;
        return c;
    }

    case CRAM_COMP_TOK3:
        return cm->level < 10
            ? COMP_TOK3_RANS
            : COMP_TOK3_ARITH;

    default:
        // Any unspecialised method
        return (enum comp_expanded)cm->method;
    }
}

// Short form of cram_block_method_int type
static char comp_method2char[COMP_MAX] =
    ".gblr0afn"           // standard CRAM methods
    "_G"                  // gzip
    "bbbbbbbbB"           // bzip2
    "rR"                  // rans4x8
    "010101014545454582"  // ransNx16
    "aAaAaAaAaaa"         // arith
    "nN";                 // tok3

// Long form of cram_block_method_int type
static char *comp_method2str[COMP_MAX] = {
    // Standard CRAM methods
    "raw", "gzip", "bzip2", "lzma", "r4x8", "rNx16",
    "arith", "fqzcomp", "tok3",

    // custom gzip
    "gzip-min", "gzip-max",

    // custom bzip2
    "bzip2-1", "bzip2-2", "bzip2-3", "bzip2-4", "bzip2-5",
    "bzip2-6", "bzip2-7", "bzip2-8", "bzip2-9",

    // rANS 4x8
    "r4x8-o0", "r4x8-o1",

    // rANS 4x16
    "r4x16-o0",   "r4x16-o1",

    "r4x16-o0R",  "r4x16-o1R",
    "r4x16-o0P",  "r4x16-o1P",
    "r4x16-o0PR", "r4x16-o1PR",
    "r32x16-o0",  "r32x16-o1",
    "r32x16-o0R", "r32x16-o1R",
    "r32x16-o0P", "r32x16-o1P",
    "r32x16-o0PR","r32x16-o1PR",
    "rNx16-xo0",  "rNx16-cat",

    // Arith
    "arith-o0",   "arith-o1",
    "arith-o0R",  "arith-o1R",
    "arith-o0P",  "arith-o1P",
    "arith-o0PR", "arith-o1PR",
    "arith-stripe", "arith-cat", "arith-ext",

    // Name tokeniser
    "tok3-rans", "tok3-arith",
};

/*----------------------------------------------------------------------
 * Manipulation and sorting of Block Content-ID arrays and hashes
 */

typedef struct {
    int64_t csize[COMP_MAX];
    int64_t usize[COMP_MAX];
} cusize_t;

static int64_t total_csize(cusize_t *cu) {
    int i;
    int64_t tot = 0;
    for (i = 0; i < COMP_MAX; i++)
        tot += cu->csize[i];
    return tot;
}

static int64_t total_usize(cusize_t *cu) {
    int i;
    int64_t tot = 0;
    for (i = 0; i < COMP_MAX; i++)
        tot += cu->usize[i];
    return tot;
}

// cusize_t array and sorting by compressed size
static cusize_t *sort_cusize_global; // avoids a messy extra data type
static int sort_cusize_compar(const void *i1, const void *i2) {
    int64_t n = sort_cusize_global->csize[*(const int *)i2] -
                sort_cusize_global->csize[*(const int *)i1];
    return n > 0 ? 1 : (n < 0 ? -1 : *(const int *)i1 - *(const int *)i2);
}

// Sort a cusize array by size of used method.
// Returns cu->csize[comp] indices in descending size, as static mem
static int *sort_cusize(cusize_t *cu) {
    static int idx[COMP_MAX];
    int i;
    for (i = 0; i < COMP_MAX; i++)
        idx[i] = i;
    sort_cusize_global = cu;
    qsort(idx, COMP_MAX, sizeof(*idx), sort_cusize_compar);

    return idx;
}

// Hash table of cusize_t and sorting by key (content-id)
KHASH_MAP_INIT_INT(cu, cusize_t)

/* Sort by hash key. Global due to rubbish qsort API, but it's simple. */
static khash_t(cu) *global_cu_hash = NULL;
static int cu_compar(const void *i1, const void *i2) {
    return kh_key(global_cu_hash, *(const int *)i1) -
           kh_key(global_cu_hash, *(const int *)i2);
}

/*----------------------------------------------------------------------
 * Main cram_size reporting and aggregation
 */
static off_t report_size(FILE *outfp, int verbose, int ref_seq_blk,
                         khash_t(cu) *cu_size, cram_cid2ds_t *cid2ds) {
    if (!cu_size || !cid2ds)
        return -1;

    khiter_t k;
    off_t tot_size = 0;

    fprintf(outfp, "#   Content_ID  Uncomp.size    Comp.size   Ratio Method%.*s  Data_series\n", verbose ? 4 : 0, "    ");
    int *sorted_blocks = malloc(kh_end(cu_size)*sizeof(int));
    if (!sorted_blocks)
        return -1;
    int nblocks = 0;
    for (k = kh_begin(cu_size); k != kh_end(cu_size); k++) {
        if (!kh_exist(cu_size, k))
            continue;
        sorted_blocks[nblocks++] = k;
    }
    global_cu_hash = cu_size;
    qsort(sorted_blocks, nblocks, sizeof(int), cu_compar);

    int i;
    for (i = 0; i < nblocks; i++) {
        k = sorted_blocks[i];

        if (verbose) {
            // FULL output
            int *comp_idx = sort_cusize(&kh_value(cu_size, k));
            int first_line = 1, c, j;
            for (c = 0; c < COMP_MAX; c++) {
                int comp = comp_idx[c];
                if (!kh_value(cu_size, k).csize[comp] && c)
                    break;

                if (!first_line)
                    fprintf(outfp, "\n");
                first_line = 0;

                if ((int)kh_key(cu_size, k) < 0)
                    fprintf(outfp, "BLOCK %8s", "CORE");
                else
                    fprintf(outfp, "BLOCK %8d", kh_key(cu_size, k));

                fprintf(outfp, " %12"PRId64" %12"PRId64,
                        kh_value(cu_size, k).usize[comp],
                        kh_value(cu_size, k).csize[comp]);
                double f = (100.0*(kh_value(cu_size, k).csize[comp]+.0001)) /
                    (kh_value(cu_size, k).usize[comp]+.0001);
                if (f > 999)
                    fprintf(outfp, "   >999%% %-11s", comp_method2str[comp]);
                else
                    fprintf(outfp, " %6.2f%% %-11s",f, comp_method2str[comp]);

                int n, *dsa = cram_cid2ds_query(cid2ds, kh_key(cu_size, k), &n);
                for (j = 0; j < n; j++) {
                    int d = dsa[j];
                    if (d > 65535)
                        fprintf(outfp, " %c%c%c", d>>16, (d>>8)&0xff, d&0xff);
                    else
                        fprintf(outfp, " %c%c", (d>>8)&0xff, d&0xff);
                }
            }
        } else {
            // aggregate by compression type.
            int64_t csize = total_csize(&kh_value(cu_size, k));
            int64_t usize = total_usize(&kh_value(cu_size, k));
            int *comp_idx = sort_cusize(&kh_value(cu_size, k));

            char cstr[COMP_MAX+1] = {0};
            int cidx = 0, c;
            for (c = 0; c < COMP_MAX; c++) {
                if (!kh_value(cu_size, k).csize[comp_idx[c]])
                    break;
                cstr[cidx++] = comp_method2char[comp_idx[c]];
            }
            if (!*cstr) *cstr = '.';

            if ((int)kh_key(cu_size, k) < 0)
                fprintf(outfp, "BLOCK %8s", "CORE");
            else
                fprintf(outfp, "BLOCK %8d", kh_key(cu_size, k));
            fprintf(outfp, " %12"PRId64" %12"PRId64, usize, csize);
            double f = 100*(csize+.0001)/(usize+.0001);
            if (f > 999)
                fprintf(outfp, "   >999%% %-7s", cstr);
            else
                fprintf(outfp, " %6.2f%% %-7s", f, cstr);

            int n, j, *dsa = cram_cid2ds_query(cid2ds, kh_key(cu_size, k), &n);
            for (j = 0; j < n; j++) {
                int d = dsa[j];
                if (d > 65535)
                    fprintf(outfp, " %c%c%c", d>>16, (d>>8)&0xff, d&0xff);
                else
                    fprintf(outfp, " %c%c", (d>>8)&0xff, d&0xff);
            }
        }

        if ((int)kh_key(cu_size, k) >= 0 &&
            (int)kh_key(cu_size, k) == ref_seq_blk) {
            fprintf(outfp, " embedded_ref");
        }
        fprintf(outfp, "\n");

        tot_size += total_csize(&kh_value(cu_size, k));
    }

    free(sorted_blocks);

    return tot_size;
}

/* Main processing loop */
static int cram_size(hFILE *hf_in, samFile *in, sam_hdr_t *h, FILE *outfp,
                     int verbose, int encodings) {
    cram_fd *in_c;
    cram_container *c = NULL;
    cram_block *blk = NULL;
    cram_block_slice_hdr *shdr = NULL;
    khiter_t k;
    int ret;
    cram_cid2ds_t *cid2ds = NULL;
    khash_t(cu) *cu_size = kh_init(cu);
    int ref_seq_blk_used = -1;
    int64_t nseqs = 0, nbases = 0, ncont = 0, nslice = 0;

    if (!in->is_cram) {
        print_error("cram_size", "Input is not a CRAM file");
        goto err;
    }
    in_c = in->fp.cram; // low level htslib abuse?
    while ((c = cram_read_container(in_c))) {
        if (cram_container_is_empty(in_c)) {
            cram_block *blk;
            // Container compression header
            if (!(blk = cram_read_block(in_c)))
                goto err;
            cram_free_block(blk);
            cram_free_container(c);
            c = NULL; blk = NULL;
            continue;
        }

        nseqs  += cram_container_get_num_records(c);
        nbases += cram_container_get_num_bases(c);

        // Container compression header
        int32_t num_slices;
        if (!(blk = cram_read_block(in_c)))
            goto err;

        // Decode compression header...
        cram_block_compression_hdr *chdr;
        chdr = cram_decode_compression_header(in_c, blk);

        if (encodings) {
            kstring_t ks = KS_INITIALIZE;
            if (cram_describe_encodings(chdr, &ks) < 0)
                goto err;

            fprintf(outfp, "Container encodings\n%s\n", ks_str(&ks));

            ks_free(&ks);
        }

        cid2ds = cram_update_cid2ds_map(chdr, cid2ds);

        cram_free_block(blk);
        blk = NULL;

        cram_free_compression_header(chdr);

        // Container num_blocks can be invalid, due to a bug.
        // Instead we iterate in slice context instead.
        (void)cram_container_get_landmarks(c, &num_slices);
        ncont++;
        nslice += num_slices;

        int i, j;
        for (i = 0; i < num_slices; i++) {
            // Slice header
            if (!(blk = cram_read_block(in_c)))
                goto err;
            if (!(shdr = cram_decode_slice_header(in_c, blk)))
                goto err;
            cram_free_block(blk);
            blk = NULL;

            int ref_seq_blk = cram_slice_hdr_get_embed_ref_id(shdr);
            int num_blocks = cram_slice_hdr_get_num_blocks(shdr);

            // Embedded reference.  Check it's consistent (if used this is
            // an almost guaranteed certainty, so we take the easy route).
            if (ref_seq_blk >= 0) {
                if (ref_seq_blk_used == -1)
                    ref_seq_blk_used = ref_seq_blk;
                else if (ref_seq_blk_used != ref_seq_blk)
                    fprintf(stderr, "Embedded reference is not consistently using the same Content-Id.\n"
                            "Reported figures for reference will be invalid.\n");
            }

            // Slice data blocks
            for (j = 0; j < num_blocks; j++) {
                // read and discard, unless it's the ref-ID block
                if (!(blk = cram_read_block(in_c)))
                    goto err;

                int32_t csize = cram_block_get_comp_size(blk);
                int32_t usize = cram_block_get_uncomp_size(blk);
                int cid = cram_block_get_content_id(blk);
                enum cram_block_method method = cram_block_get_method(blk);

                // Expand comp to the internal sub-formats, eg
                // rANS order-0/1, PACK+RLE, etc.
                cram_method_details *cm;
                cm = cram_expand_method(cram_block_get_data(blk),
                                        cram_block_get_comp_size(blk),
                                        method);
                if (!cm)
                    goto err;
                enum comp_expanded comp
                    = comp_method2expanded(cm);
                free(cm);

                k = kh_put(cu, cu_size, cid, &ret);
                if (ret < 0)
                    goto err;
                if (ret == 0) {
                    kh_value(cu_size, k).csize[comp] += csize;
                    kh_value(cu_size, k).usize[comp] += usize;
                } else {
                    memset(&kh_value(cu_size, k), 0, sizeof(cusize_t));
                    kh_value(cu_size, k).csize[comp]  = csize;
                    kh_value(cu_size, k).usize[comp]  = usize;
                }

                cram_free_block(blk);
                blk = NULL;
            }
            cram_free_slice_header(shdr);
            shdr = NULL;
        }

        cram_free_container(c);
        c = NULL;
    }

    off_t tot_size = report_size(outfp, verbose, ref_seq_blk_used,
                                 cu_size, cid2ds);
    if (tot_size < 0)
        goto err;

    kh_destroy(cu, cu_size);
    cram_cid2ds_free(cid2ds);

    off_t end = htell(hf_in);

    fprintf(outfp, "\n");
    fprintf(outfp, "Number of containers  %18"PRId64"\n", ncont);
    fprintf(outfp, "Number of slices      %18"PRId64"\n", nslice);
    fprintf(outfp, "Number of sequences   %18"PRId64"\n", nseqs);
    fprintf(outfp, "Number of bases       %18"PRId64"\n", nbases);
    fprintf(outfp, "Total file size       %18"PRId64"\n", end);
    fprintf(outfp, "Format overhead size  %18"PRId64"\n", end - tot_size);

    return 0;

 err:
    // Report anyway so we can get stats on partial files, but be
    // sure to error too.
    report_size(outfp, verbose, ref_seq_blk_used, cu_size, cid2ds);

    print_error("cram_size", "Failed in decoding CRAM file");
    if (blk)
        cram_free_block(blk);
    if (shdr)
        cram_free_slice_header(shdr);
    if (c)
        cram_free_container(c);
    if (cid2ds)
        cram_cid2ds_free(cid2ds);

    return -1;
}

/* main() for cram_size */
int main_cram_size(int argc, char *argv[]) {
    int c, usage = 0, verbose = 0, encodings = 0;
    sam_hdr_t *h = 0;
    hFILE *hf_in = NULL;
    samFile *in = NULL;
    sam_global_args ga;
    FILE *outfp = stdout;

    static const struct option lopts[] = {
        {"output", required_argument, NULL, 'o'},
        {"verbose",  no_argument, NULL, 'v'},
        {"encodings", no_argument, NULL, 'e'},
        SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', '-', '-'),
        { NULL, 0, NULL, 0 }
    };

    sam_global_args_init(&ga);

    while ((c = getopt_long(argc, argv, "vo:e", lopts, NULL)) >= 0) {
        switch (c) {
        case 'o':
            if (!(outfp = fopen(optarg, "w"))) {
                perror(optarg);
                goto err;
            }
            break;

        case 'v':
            verbose++;
            break;

        case 'e':
            encodings++;
            break;

        default:
            if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
        case '?': usage=1; break;
        }
    }

    if ((optind == argc && isatty(0)) || usage) {
        printf("Usage: samtools cram_size [-ve] [-o out.size] [in.cram]\n");
        return 0;
    }

    char *fn = optind < argc ? argv[optind] : "-";

    // We want access to in->fp.cram->fp, but this is an opaque struct so we
    // can't get that.  However we opened with hopen and then reopen as
    // CRAM with hts_hopen, which will swallow the initial hFILE and take
    // owenership of it.  Hence we now know in->fp.cram->fp.
    if (!(hf_in = hopen(fn, "r"))) {
        print_error_errno("cram_size", "failed to open file '%s'", fn);
        return 1;
    }
    if (!(in = hts_hopen(hf_in, fn, "r"))) {
        print_error_errno("cram_size", "failed to open file '%s'", fn);
        goto err;
    }

    if (!(h = sam_hdr_read(in)))
        goto err;

    int ret = cram_size(hf_in, in, h, outfp, verbose, encodings);
    sam_hdr_destroy(h);
    sam_close(in);
    if (outfp != stdout)
        fclose(outfp);

    return ret ? 1 : 0;

 err:
    if (in)
        sam_close(in);
    if (h)
        sam_hdr_destroy(h);

    return 1;
}
