/*  bedcov.c -- bedcov subcommand.

    Copyright (C) 2012 Broad Institute.
    Copyright (C) 2013-2014, 2018-2022 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

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

#include <config.h>

#include <zlib.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"
#include "samtools.h"
#include "sam_opts.h"

#include "htslib/kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

#define DEFAULT_DEPTH 64000

typedef struct {
    htsFile *fp;
    sam_hdr_t *header;
    hts_itr_t *iter;
    int min_mapQ;
    uint32_t flags;  // read filtering flags
    int64_t rcnt;
} aux_t;

static int read_bam(void *data, bam1_t *b)
{
    aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int ret;
    while (1)
    {
        ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->header, b);
        if ( ret<0 ) break;
        if ( b->core.flag & aux->flags ) continue;
        if ( (int)b->core.qual < aux->min_mapQ ) continue;
        break;
    }
    return ret;
}

static int incr_rcnt(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    aux_t *aux = (aux_t *)data;
    aux->rcnt++;
    return 0;
}

int main_bedcov(int argc, char *argv[])
{
    gzFile fp;
    kstring_t str;
    kstream_t *ks;
    hts_idx_t **idx;
    aux_t **aux;
    int *n_plp, dret, i, j, m, n, c, ret, status = 0, min_mapQ = 0;
    int skip_DN = 0, do_rcount = 0;
    int64_t *cnt, *pcov = NULL;
    const bam_pileup1_t **plp;
    int usage = 0, has_index_file = 0;
    uint32_t flags = (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP);
    int tflags = 0, min_depth = -1;

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        {"min-MQ", required_argument, NULL, 'Q'},
        {"min-mq", required_argument, NULL, 'Q'},
        SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', 0, '-'),
        { NULL, 0, NULL, 0 }
    };

    while ((c = getopt_long(argc, argv, "Q:Xg:G:jd:c", lopts, NULL)) >= 0) {
        switch (c) {
        case 'Q': min_mapQ = atoi(optarg); break;
        case 'X': has_index_file = 1; break;
        case 'c': do_rcount = 1; break;
        case 'g':
            tflags = bam_str2flag(optarg);
            if (tflags < 0 || tflags > ((BAM_FSUPPLEMENTARY << 1) - 1)) {
                print_error("bedcov", "Flag value \"%s\" is not supported", optarg);
                return 1;
            }
            flags &= ~tflags;
            break;
        case 'G':
            tflags = bam_str2flag(optarg);
            if (tflags < 0 || tflags > ((BAM_FSUPPLEMENTARY << 1) - 1)) {
                print_error("bedcov", "Flag value \"%s\" is not supported", optarg);
                return 1;
            }
            flags |= tflags;
            break;
        case 'j': skip_DN = 1; break;
        case 'd': min_depth = atoi(optarg); break;
        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
                  /* else fall-through */
        case '?': usage = 1; break;
        }
        if (usage) break;
    }
    if (usage || optind + 2 > argc) {
        fprintf(stderr, "Usage: samtools bedcov [options] <in.bed> <in1.bam> [...]\n\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "      -Q, --min-MQ <int>  mapping quality threshold [0]\n");
        fprintf(stderr, "      -X                  use customized index files\n");
        fprintf(stderr, "      -g <flags>          remove the specified flags from the set used to filter out reads\n");
        fprintf(stderr, "      -G <flags>          add the specified flags to the set used to filter out reads\n"
                        "                          The default set is UNMAP,SECONDARY,QCFAIL,DUP or 0x704\n");
        fprintf(stderr, "      -j                  do not include deletions (D) and ref skips (N) in bedcov computation\n");
        fprintf(stderr, "      -d <int>            depth threshold. Number of reference bases with coverage above and"
                        "                          including this value will be displayed in a separate column\n");
        fprintf(stderr, "      -c                  add an additional column showing read count\n");
        sam_global_opt_help(stderr, "-.--.--.");
        return 1;
    }
    if (has_index_file) {
        if ((argc - optind - 1) % 2 != 0) { // Calculate # of input BAM files
            fprintf(stderr, "ERROR: odd number of filenames detected! Each BAM file should have an index file\n");
            return 1;
        }
        n = (argc - optind - 1) / 2;
    } else {
        n = argc - optind - 1;
    }

    memset(&str, 0, sizeof(kstring_t));
    aux = calloc(n, sizeof(aux_t*));
    idx = calloc(n, sizeof(hts_idx_t*));
    for (i = 0; i < n; ++i) {
        aux[i] = calloc(1, sizeof(aux_t));
        aux[i]->min_mapQ = min_mapQ;
        aux[i]->fp = sam_open_format(argv[i+optind+1], "r", &ga.in);
        if (aux[i]->fp) {
            // If index filename has not been specfied, look in BAM folder
            if (has_index_file) {
                idx[i] = sam_index_load2(aux[i]->fp, argv[i+optind+1], argv[i+optind+n+1]);
            } else {
                idx[i] = sam_index_load(aux[i]->fp, argv[i+optind+1]);
            }
        }
        if (aux[i]->fp == 0 || idx[i] == 0) {
            fprintf(stderr, "ERROR: fail to open index BAM file '%s'\n", argv[i+optind+1]);
            return 2;
        }
        // TODO bgzf_set_cache_size(aux[i]->fp, 20);
        aux[i]->header = sam_hdr_read(aux[i]->fp);
        if (aux[i]->header == NULL) {
            fprintf(stderr, "ERROR: failed to read header for '%s'\n",
                    argv[i+optind+1]);
            return 2;
        }
        aux[i]->flags = flags;
    }
    cnt = calloc(n, sizeof(*cnt));

    if (min_depth >= 0) pcov = calloc(n, sizeof(*pcov));
    if (!cnt || (min_depth >= 0 && !pcov)) {
        print_error_errno("bedcov", "failed to allocate memory");
        return 2;
    }

    fp = gzopen(argv[optind], "rb");
    if (fp == NULL) {
        print_error_errno("bedcov", "can't open BED file '%s'", argv[optind]);
        return 2;
    }
    ks = ks_init(fp);
    n_plp = calloc(n, sizeof(int));
    plp = calloc(n, sizeof(bam_pileup1_t*));
    while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
        char *p, *q;
        int tid, pos, num = 0;
        int64_t beg = 0, end = 0;
        bam_mplp_t mplp;

        if (str.l == 0 || *str.s == '#') continue; /* empty or comment line */
        /* Track and browser lines.  Also look for a trailing *space* in
           case someone has badly-chosen a chromosome name (it would
           be followed by a tab in that case). */
        if (strncmp(str.s, "track ", 6) == 0) continue;
        if (strncmp(str.s, "browser ", 8) == 0) continue;
        for (p = q = str.s; *p && !isspace(*p); ++p);
        if (*p == 0) goto bed_error;
        char c = *p;
        *p = 0; tid = bam_name2id(aux[0]->header, q); *p = c;
        if (tid < 0) goto bed_error;
        num = sscanf(p + 1, "%"SCNd64" %"SCNd64, &beg, &end);
        if (num < 2 || end < beg) goto bed_error;

        for (i = 0; i < n; ++i) {
            if (aux[i]->iter) hts_itr_destroy(aux[i]->iter);
            aux[i]->iter = sam_itr_queryi(idx[i], tid, beg, end);
            aux[i]->rcnt = 0;
        }

        mplp = bam_mplp_init(n, read_bam, (void**)aux);
        if (min_depth > DEFAULT_DEPTH)
            bam_mplp_set_maxcnt(mplp, min_depth);
        else
            bam_mplp_set_maxcnt(mplp, DEFAULT_DEPTH);

        memset(cnt, 0, sizeof(*cnt) * n);
        if (min_depth >= 0) memset(pcov, 0, sizeof(*pcov) * n);

        if (do_rcount)
            bam_mplp_constructor(mplp, incr_rcnt);

        while ((ret = bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)) > 0)
            if (pos >= beg && pos < end) {
                for (i = 0; i < n; ++i) {
                    m = 0;
                    if (skip_DN || min_depth >= 0) {
                        for (j = 0; j < n_plp[i]; ++j) {
                            const bam_pileup1_t *pi = plp[i] + j;
                            if (pi->is_del || pi->is_refskip) ++m;
                        }
                    }
                    int pd = n_plp[i] - m;
                    cnt[i] += pd;
                    if (min_depth >= 0 && pd >= min_depth) pcov[i]++;
                }
            }

        if (ret < 0) {
            print_error("bedcov", "error reading from input file");
            status = 2;
            bam_mplp_destroy(mplp);
            break;
        }

        for (i = 0; i < n; ++i) {
            kputc('\t', &str);
            kputl(cnt[i], &str);
        }
        if (min_depth >= 0) {
            for (i = 0; i < n; ++i) {
                kputc('\t', &str);
                kputl(pcov[i], &str);
            }
        }
        if (do_rcount) {
            for (i = 0; i < n; ++i) {
                kputc('\t', &str);
                kputl(aux[i]->rcnt, &str);
            }
        }
        puts(str.s);
        bam_mplp_destroy(mplp);
        continue;

bed_error:
        fprintf(stderr, "Errors in BED line '%s'\n", str.s);
        status = 2;
    }
    free(n_plp); free(plp);
    ks_destroy(ks);
    gzclose(fp);

    free(cnt);
    free(pcov);
    for (i = 0; i < n; ++i) {
        if (aux[i]->iter) hts_itr_destroy(aux[i]->iter);
        hts_idx_destroy(idx[i]);
        sam_hdr_destroy(aux[i]->header);
        sam_close(aux[i]->fp);
        free(aux[i]);
    }
    free(aux); free(idx);
    free(str.s);
    sam_global_args_free(&ga);
    return status;
}
