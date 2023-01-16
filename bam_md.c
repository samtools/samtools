/*  bam_md.c -- calmd subcommand.

    Copyright (C) 2009-2011, 2014-2015, 2019-2020, 2022 Genome Research Ltd.
    Portions copyright (C) 2009-2011 Broad Institute.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <errno.h>
#include <assert.h>
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/thread_pool.h"
#include "sam_opts.h"
#include "samtools.h"

#define USE_EQUAL 1
#define DROP_TAG  2
#define BIN_QUAL  4
#define UPDATE_NM 8
#define UPDATE_MD 16
#define HASH_QNM  32

typedef struct cached_ref_entry {
    char *ref;
    hts_pos_t len;
} cached_ref_entry;

typedef struct ref_cache {
    cached_ref_entry *refs;
    char *last_ref;
    hts_pos_t last_len;
    int nref;
    int last_tid;
} ref_cache;

int bam_aux_drop_other(bam1_t *b, uint8_t *s);

static int bam_fillmd1_core(const char *ref_name, bam1_t *b, char *ref,
                            hts_pos_t ref_len, int flag, int max_nm,
                            int quiet_mode, uint32_t *skipped)
{
    uint8_t *seq = bam_get_seq(b);
    uint32_t *cigar = bam_get_cigar(b);
    bam1_core_t *c = &b->core;
    int i, qpos, matched = 0;
    hts_pos_t rpos;
    kstring_t str = KS_INITIALIZE;
    int32_t old_nm_i = -1, nm = 0;
    uint32_t err = 0;

    if (c->l_qseq == 0) {
        if (!quiet_mode) {
            if (ref_name) {
                fprintf(stderr, "[bam_fillmd1] no sequence in alignment "
                        "record for '%s' at %s:%"PRIhts_pos", skipped\n",
                        bam_get_qname(b), ref_name, c->pos + 1);
            } else {
                fprintf(stderr, "[bam_fillmd1] no sequence in alignment "
                        "record for '%s', skipped", bam_get_qname(b));
            }
        }
        if (skipped) (*skipped)++;
        return 0;
    }

    for (i = qpos = 0, rpos = c->pos; i < c->n_cigar; ++i) {
        int j, oplen = cigar[i]>>4, op = cigar[i]&0xf;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (j = 0; j < oplen; ++j) {
                int c1, c2, z = qpos + j;
                if (rpos+j >= ref_len || z >= c->l_qseq || ref[rpos+j] == '\0')
                    break; // out of bounds
                c1 = bam_seqi(seq, z);
                c2 = seq_nt16_table[(uint8_t)ref[rpos+j]];
                if ((c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0) { // a match
                    if (flag&USE_EQUAL) seq[z/2] &= (z&1)? 0xf0 : 0x0f;
                    ++matched;
                } else {
                    err |= kputw(matched, &str) < 0;
                    err |= kputc(toupper(ref[rpos+j]), &str) < 0;
                    matched = 0; ++nm;
                }
            }
            if (j < oplen) break;
            rpos += oplen; qpos += oplen;
        } else if (op == BAM_CDEL) {
            err |= kputw(matched, &str) < 0;
            err |= kputc('^', &str) < 0;
            for (j = 0; j < oplen; ++j) {
                if (rpos+j >= ref_len || ref[rpos+j] == '\0') break;
                err |= kputc(toupper(ref[rpos+j]), &str) < 0;
            }
            matched = 0;
            rpos += j; nm += j;
            if (j < oplen) break;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            qpos += oplen;
            if (op == BAM_CINS) nm += oplen;
        } else if (op == BAM_CREF_SKIP) {
            rpos += oplen;
        }
    }
    err |= kputw(matched, &str) < 0;
    if (err) {
        print_error_errno("calmd", "Couldn't build new MD string");
        goto fail;
    }
    // apply max_nm
    if (max_nm > 0 && nm >= max_nm) {
        for (i = qpos = 0, rpos = c->pos; i < c->n_cigar; ++i) {
            int j, oplen = cigar[i]>>4, op = cigar[i]&0xf;
            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                for (j = 0; j < oplen; ++j) {
                    int c1, c2, z = qpos + j;
                    if (rpos+j >= ref_len || z >= c->l_qseq || ref[rpos+j] == '\0')
                        break; // out of bounds
                    c1 = bam_seqi(seq, z);
                    c2 = seq_nt16_table[(uint8_t)ref[rpos+j]];
                    if ((c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0) { // a match
                        seq[z/2] |= (z&1)? 0x0f : 0xf0;
                        bam_get_qual(b)[z] = 0;
                    }
                }
                if (j < oplen) break;
                rpos += oplen; qpos += oplen;
            } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) rpos += oplen;
            else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) qpos += oplen;
        }
    }
    // update NM
    if ((flag & UPDATE_NM) && !(c->flag & BAM_FUNMAP)) {
        uint8_t *old_nm = bam_aux_get(b, "NM");
        if (old_nm) old_nm_i = bam_aux2i(old_nm);
        if (!old_nm) {
            if (bam_aux_append(b, "NM", 'i', 4, (uint8_t*)&nm) < 0)
                goto aux_fail;
        }
        else if (nm != old_nm_i) {
            if (!quiet_mode) {
                fprintf(stderr, "[bam_fillmd1] different NM for read '%s': %d -> %d\n", bam_get_qname(b), old_nm_i, nm);
            }
            if (bam_aux_del(b, old_nm) < 0) goto aux_fail;
            if (bam_aux_append(b, "NM", 'i', 4, (uint8_t*)&nm) < 0)
                goto aux_fail;
        }
    }
    // update MD
    if ((flag & UPDATE_MD) && !(c->flag & BAM_FUNMAP)) {
        uint8_t *old_md = bam_aux_get(b, "MD");
        if (!old_md) {
            if (bam_aux_append(b, "MD", 'Z', str.l + 1, (uint8_t*)str.s) < 0)
                goto aux_fail;
        } else {
            int is_diff = 0;
            if (strlen((char*)old_md+1) == str.l) {
                for (i = 0; i < str.l; ++i)
                    if (toupper(old_md[i+1]) != toupper(str.s[i]))
                        break;
                if (i < str.l) is_diff = 1;
            } else is_diff = 1;
            if (is_diff) {
                if (!quiet_mode) {
                    fprintf(stderr, "[bam_fillmd1] different MD for read '%s': '%s' -> '%s'\n", bam_get_qname(b), old_md+1, str.s);
                }
                if (bam_aux_del(b, old_md) < 0) goto aux_fail;
                if (bam_aux_append(b, "MD", 'Z', str.l + 1, (uint8_t*)str.s) < 0)
                    goto aux_fail;
            }
        }
    }

    // drop all tags but RG
    if (flag&DROP_TAG) {
        uint8_t *q = bam_aux_get(b, "RG");
        bam_aux_drop_other(b, q);
    }
    // reduce the resolution of base quality
    if (flag&BIN_QUAL) {
        uint8_t *qual = bam_get_qual(b);
        for (i = 0; i < b->core.l_qseq; ++i)
            if (qual[i] >= 3) qual[i] = qual[i]/10*10 + 7;
    }

    free(str.s);
    return 0;

 aux_fail:
    if (errno == ENOMEM) {
        print_error("calmd", "Couldn't add aux tag (too long)");
    } else if (errno == EINVAL) {
        print_error("calmd", "Corrupt aux data");
    } else {
        print_error_errno("calmd", "Couldn't add aux tag");
    }
 fail:
    free(str.s);
    return -1;
}

int bam_fillmd1(bam1_t *b, char *ref, int flag, int quiet_mode)
{
    return bam_fillmd1_core(NULL, b, ref, INT_MAX, flag, 0, quiet_mode, NULL);
}

// Get a new reference sequence.
// For position-sorted inputs, the previous reference should never be
// needed again and can be discarded to save memory.  For other orderings,
// references are stored in a cache in case they're required in the future.
// The caching mode is turned on if the requested  tid is less than the last
// one used, indicating the file ordering doesn't match the sequence dictionary.
static int get_ref(faidx_t *fai, sam_hdr_t *header, ref_cache *cache,
                   int tid, char **ref_out, const char **ref_name_out,
                   hts_pos_t *len_out)
{
    char *ref = NULL;
    const char *ref_name;
    hts_pos_t len = 0;

    // This should only be called when tid changes
    assert(tid != cache->last_tid);

    // Array lookup, should be fast
    ref_name = sam_hdr_tid2name(header, tid);
    *ref_name_out = ref_name;

    // Return a cached entry, if available
    if (cache->refs && tid >= 0 && tid < cache->nref
        && cache->refs[tid].ref) {
        assert(cache->last_ref == NULL);
        *ref_out = cache->refs[tid].ref;
        *len_out = cache->refs[tid].len;
        cache->last_tid = tid;
        return 0;
    }

    // Try to get the reference
    if (ref_name)
        ref = fai_fetch64(fai, ref_name, &len);

    if (!ref) {
        // Historically, calmd doesn't worry too much about missing refs
        *ref_out = NULL;
        *len_out = 0;
        return 0;
    }

    if (!cache->refs && cache->last_tid > tid) {
        // Going backwards throught the list of tids implies
        // a non-position-ordered file, so turn on caching mode
        cache->nref = sam_hdr_nref(header);
        if (cache->nref < 0) {
            print_error("calmd", "couldn't get number of refs from header");
            return -1;
        }
        if (cache->nref > 0) {
            cache->refs = calloc(cache->nref, sizeof(cache->refs[0]));
            if (!cache->refs) {
                print_error_errno("calmd",
                                  "couldn't allocate reference cache");
                return -1;
            }
            // Add the reference we already have as the first entry
            if (cache->last_tid >= 0 && cache->last_tid < cache->nref) {
                cache->refs[cache->last_tid].ref = cache->last_ref;
                cache->refs[cache->last_tid].len = cache->last_len;
            } else {
                free(cache->last_ref);
            }
            cache->last_ref = NULL;
        }
    }

    if (cache->refs) {
        assert(cache->last_ref == NULL);  // Shouldn't be set when caching
        // Add the new reference to the cache
        if (tid >= 0 && tid < cache->nref) {
            cache->refs[tid].ref = ref;
            cache->refs[tid].len = len;
        }
    } else {
        // Streaming mode - free the last ref and replace it with this one
        free(cache->last_ref);
        cache->last_ref = ref;
        cache->last_len = len;
    }

    *ref_out = ref;
    *len_out = len;
    cache->last_tid = tid;
    return 0;
}

static void refs_destroy(ref_cache *cache) {
    if (cache->refs) {
        int i;
        assert(cache->last_ref == NULL);
        for (i = 0; i < cache->nref; i++)
            free(cache->refs[i].ref);
        free(cache->refs);
    } else {
        free(cache->last_ref);
    }
}

int calmd_usage() {
    fprintf(stderr,
"Usage: samtools calmd [-eubrAESQ] <aln.bam> <ref.fasta>\n"
"Options:\n"
"  -e       change identical bases to '='\n"
"  -u       uncompressed BAM output (for piping)\n"
"  -b       compressed BAM output\n"
"  -S       ignored (input format is auto-detected)\n"
"  -A       modify the quality string\n"
"  -Q       use quiet mode to output less debug info to stdout\n"
"  -r       compute the BQ tag (without -A) or cap baseQ by BAQ (with -A)\n"
"  -E       extended BAQ for better sensitivity but lower specificity\n"
"  --no-PG  do not add a PG line\n");

    sam_global_opt_help(stderr, "-....@-.");
    return 1;
}

int bam_fillmd(int argc, char *argv[])
{
    int c, flt_flag, ret, is_bam_out, is_uncompressed, max_nm, is_realn, capQ, baq_flag, quiet_mode, no_pg = 0;
    hts_pos_t len = 0;
    htsThreadPool p = {NULL, 0};
    samFile *fp = NULL, *fpout = NULL;
    sam_hdr_t *header = NULL;
    faidx_t *fai = NULL;
    char *ref = NULL, mode_w[8], *ref_file, *arg_list = NULL;
    ref_cache refs = { NULL, NULL, 0, 0, -2 };
    const char *ref_name = NULL;
    bam1_t *b = NULL;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    uint32_t skipped = 0;

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 0, 0, 0,'@'),
        {"no-PG", no_argument, NULL, 1},
        { NULL, 0, NULL, 0 }
    };

    flt_flag = UPDATE_NM | UPDATE_MD;
    is_bam_out = is_uncompressed = is_realn = max_nm = capQ = baq_flag = quiet_mode = 0;
    strcpy(mode_w, "w");
    while ((c = getopt_long(argc, argv, "EqQreuNhbSC:n:Ad@:", lopts, NULL)) >= 0) {
        switch (c) {
        case 'r': is_realn = 1; break;
        case 'e': flt_flag |= USE_EQUAL; break;
        case 'd': flt_flag |= DROP_TAG; break;
        case 'q': flt_flag |= BIN_QUAL; break;
        case 'h': flt_flag |= HASH_QNM; break;
        case 'N': flt_flag &= ~(UPDATE_MD|UPDATE_NM); break;
        case 'b': is_bam_out = 1; break;
        case 'u': is_uncompressed = is_bam_out = 1; break;
        case 'S': break;
        case 'n': max_nm = atoi(optarg); break;
        case 'C': capQ = atoi(optarg); break;
        case 'A': baq_flag |= 1; break;
        case 'E': baq_flag |= 2; break;
        case 'Q': quiet_mode = 1; break;
        case 1: no_pg = 1; break;
        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            fprintf(stderr, "[bam_fillmd] unrecognized option '-%c'\n\n", c);
            /* else fall-through */
        case '?': return calmd_usage();
        }
    }
    if (is_bam_out) strcat(mode_w, "b");
    else strcat(mode_w, "h");
    if (is_uncompressed) strcat(mode_w, "0");
    if (optind + (ga.reference == NULL) >= argc)
        return calmd_usage();
    fp = sam_open_format(argv[optind], "r", &ga.in);
    if (fp == NULL) {
        print_error_errno("calmd", "Failed to open input file '%s'", argv[optind]);
        return 1;
    }

    if (!no_pg && !(arg_list = stringify_argv(argc+1, argv-1))) {
        print_error("calmd", "failed to create arg_list");
        return 1;
    }

    header = sam_hdr_read(fp);
    if (header == NULL || sam_hdr_nref(header) == 0) {
        fprintf(stderr, "[bam_fillmd] input SAM does not have header. Abort!\n");
        goto fail;
    }

    fpout = sam_open_format("-", mode_w, &ga.out);
    if (fpout == NULL) {
        print_error_errno("calmd", "Failed to open output");
        goto fail;
    }
    if (!no_pg && sam_hdr_add_pg(header, "samtools",
                                 "VN", samtools_version(),
                                 arg_list ? "CL": NULL,
                                 arg_list ? arg_list : NULL,
                                 NULL)) {
        print_error("calmd", "failed to add PG line to header");
        goto fail;
    }
    if (sam_hdr_write(fpout, header) < 0) {
        print_error_errno("calmd", "Failed to write sam header");
        goto fail;
    }

    if (ga.nthreads > 0) {
        if (!(p.pool = hts_tpool_init(ga.nthreads))) {
            fprintf(stderr, "Error creating thread pool\n");
            goto fail;
        }
        hts_set_opt(fp,    HTS_OPT_THREAD_POOL, &p);
        hts_set_opt(fpout, HTS_OPT_THREAD_POOL, &p);
    }

    ref_file = argc > optind + 1 ? argv[optind+1] : ga.reference;
    fai = fai_load(ref_file);

    if (!fai) {
        print_error_errno("calmd", "Failed to open reference file '%s'", ref_file);
        goto fail;
    }

    b = bam_init1();
    if (!b) {
        fprintf(stderr, "[bam_fillmd] Failed to allocate bam struct\n");
        goto fail;
    }
    while ((ret = sam_read1(fp, header, b)) >= 0) {
        if (b->core.tid >= 0) {
            if (refs.last_tid != b->core.tid) {
                if (get_ref(fai, header, &refs, b->core.tid,
                            &ref, &ref_name, &len) < 0) {
                    goto fail;
                }
                if (ref == 0) { // FIXME: Should this always be fatal?
                    fprintf(stderr, "[bam_fillmd] fail to find sequence '%s' in the reference.\n",
                            ref_name ? ref_name : "(unknown)");
                    if (is_realn || capQ > 10) goto fail; // Would otherwise crash
                }
            }
            if (is_realn) {
                if (sam_prob_realn(b, ref, len, baq_flag) < -3) {
                    print_error_errno("calmd", "BAQ alignment failed");
                    goto fail;
                }
            }
            if (capQ > 10) {
                int q = sam_cap_mapq(b, ref, len, capQ);
                if (b->core.qual > q) b->core.qual = q;
            }
            if (ref) {
                if (bam_fillmd1_core(ref_name, b, ref, len, flt_flag, max_nm,
                                     quiet_mode, &skipped) < 0)
                    goto fail;
            }
        }
        if (sam_write1(fpout, header, b) < 0) {
            print_error_errno("calmd", "failed to write to output file");
            goto fail;
        }
    }
    if (ret < -1) {
        fprintf(stderr, "[bam_fillmd] Error reading input.\n");
        goto fail;
    }

    if (skipped) {
        fprintf(stderr, "[calmd] Warning: %"PRIu32" records skipped due "
                "to no query sequence\n",
                skipped);
    }

    bam_destroy1(b);
    sam_hdr_destroy(header);

    free(arg_list);
    refs_destroy(&refs);
    fai_destroy(fai);
    sam_close(fp);
    if (sam_close(fpout) < 0) {
        fprintf(stderr, "[bam_fillmd] error when closing output file\n");
        return 1;
    }
    if (p.pool) hts_tpool_destroy(p.pool);

    return 0;

 fail:
    free(arg_list);
    refs_destroy(&refs);
    if (b) bam_destroy1(b);
    if (header) sam_hdr_destroy(header);
    if (fai) fai_destroy(fai);
    if (fp) sam_close(fp);
    if (fpout) sam_close(fpout);
    if (p.pool) hts_tpool_destroy(p.pool);

    return 1;
}
