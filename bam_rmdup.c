/*  bam_rmdup.c -- duplicate read detection.

    Copyright (C) 2009, 2015 Genome Research Ltd.
    Portions copyright (C) 2009 Broad Institute.

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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <unistd.h>
#include "htslib/sam.h"
#include "sam_opts.h"
#include "bam.h" // for bam_get_library

typedef bam1_t *bam1_p;

#include "htslib/khash.h"
KHASH_SET_INIT_STR(name)
KHASH_MAP_INIT_INT64(pos, bam1_p)

#define BUFFER_SIZE 0x40000

typedef struct {
    uint64_t n_checked, n_removed;
    khash_t(pos) *best_hash;
} lib_aux_t;
KHASH_MAP_INIT_STR(lib, lib_aux_t)

typedef struct {
    int n, max;
    bam1_t **a;
} tmp_stack_t;

static inline void stack_insert(tmp_stack_t *stack, bam1_t *b)
{
    if (stack->n == stack->max) {
        stack->max = stack->max? stack->max<<1 : 0x10000;
        stack->a = (bam1_t**)realloc(stack->a, sizeof(bam1_t*) * stack->max);
    }
    stack->a[stack->n++] = b;
}

static inline void dump_best(tmp_stack_t *stack, samFile *out, bam_hdr_t *hdr)
{
    int i;
    for (i = 0; i != stack->n; ++i) {
        sam_write1(out, hdr, stack->a[i]);
        bam_destroy1(stack->a[i]);
    }
    stack->n = 0;
}

static void clear_del_set(khash_t(name) *del_set)
{
    khint_t k;
    for (k = kh_begin(del_set); k < kh_end(del_set); ++k)
        if (kh_exist(del_set, k))
            free((char*)kh_key(del_set, k));
    kh_clear(name, del_set);
}

static lib_aux_t *get_aux(khash_t(lib) *aux, const char *lib)
{
    khint_t k = kh_get(lib, aux, lib);
    if (k == kh_end(aux)) {
        int ret;
        char *p = strdup(lib);
        lib_aux_t *q;
        k = kh_put(lib, aux, p, &ret);
        q = &kh_val(aux, k);
        q->n_checked = q->n_removed = 0;
        q->best_hash = kh_init(pos);
        return q;
    } else return &kh_val(aux, k);
}

static void clear_best(khash_t(lib) *aux, int max)
{
    khint_t k;
    for (k = kh_begin(aux); k != kh_end(aux); ++k) {
        if (kh_exist(aux, k)) {
            lib_aux_t *q = &kh_val(aux, k);
            if (kh_size(q->best_hash) >= max)
                kh_clear(pos, q->best_hash);
        }
    }
}

static inline int sum_qual(const bam1_t *b)
{
    int i, q;
    uint8_t *qual = bam_get_qual(b);
    for (i = q = 0; i < b->core.l_qseq; ++i) q += qual[i];
    return q;
}

void bam_rmdup_core(samFile *in, bam_hdr_t *hdr, samFile *out)
{
    bam1_t *b;
    int last_tid = -1, last_pos = -1;
    tmp_stack_t stack;
    khint_t k;
    khash_t(lib) *aux;
    khash_t(name) *del_set;

    aux = kh_init(lib);
    del_set = kh_init(name);
    b = bam_init1();
    memset(&stack, 0, sizeof(tmp_stack_t));

    kh_resize(name, del_set, 4 * BUFFER_SIZE);
    while (sam_read1(in, hdr, b) >= 0) {
        bam1_core_t *c = &b->core;
        if (c->tid != last_tid || last_pos != c->pos) {
            dump_best(&stack, out, hdr); // write the result
            clear_best(aux, BUFFER_SIZE);
            if (c->tid != last_tid) {
                clear_best(aux, 0);
                if (kh_size(del_set)) { // check
                    fprintf(stderr, "[bam_rmdup_core] %llu unmatched pairs\n", (long long)kh_size(del_set));
                    clear_del_set(del_set);
                }
                if ((int)c->tid == -1) { // append unmapped reads
                    sam_write1(out, hdr, b);
                    while (sam_read1(in, hdr, b) >= 0) sam_write1(out, hdr, b);
                    break;
                }
                last_tid = c->tid;
                fprintf(stderr, "[bam_rmdup_core] processing reference %s...\n", hdr->target_name[c->tid]);
            }
        }
        if (!(c->flag&BAM_FPAIRED) || (c->flag&(BAM_FUNMAP|BAM_FMUNMAP)) || (c->mtid >= 0 && c->tid != c->mtid)) {
            sam_write1(out, hdr, b);
        } else if (c->isize > 0) { // paired, head
            uint64_t key = (uint64_t)c->pos<<32 | c->isize;
            const char *lib;
            lib_aux_t *q;
            int ret;
            lib = bam_get_library(hdr, b);
            q = lib? get_aux(aux, lib) : get_aux(aux, "\t");
            ++q->n_checked;
            k = kh_put(pos, q->best_hash, key, &ret);
            if (ret == 0) { // found in best_hash
                bam1_t *p = kh_val(q->best_hash, k);
                ++q->n_removed;
                if (sum_qual(p) < sum_qual(b)) { // the current alignment is better; this can be accelerated in principle
                    kh_put(name, del_set, strdup(bam_get_qname(p)), &ret); // p will be removed
                    bam_copy1(p, b); // replaced as b
                } else kh_put(name, del_set, strdup(bam_get_qname(b)), &ret); // b will be removed
                if (ret == 0)
                    fprintf(stderr, "[bam_rmdup_core] inconsistent BAM file for pair '%s'. Continue anyway.\n", bam_get_qname(b));
            } else { // not found in best_hash
                kh_val(q->best_hash, k) = bam_dup1(b);
                stack_insert(&stack, kh_val(q->best_hash, k));
            }
        } else { // paired, tail
            k = kh_get(name, del_set, bam_get_qname(b));
            if (k != kh_end(del_set)) {
                free((char*)kh_key(del_set, k));
                kh_del(name, del_set, k);
            } else sam_write1(out, hdr, b);
        }
        last_pos = c->pos;
    }

    for (k = kh_begin(aux); k != kh_end(aux); ++k) {
        if (kh_exist(aux, k)) {
            lib_aux_t *q = &kh_val(aux, k);
            dump_best(&stack, out, hdr);
            fprintf(stderr, "[bam_rmdup_core] %lld / %lld = %.4lf in library '%s'\n", (long long)q->n_removed,
                    (long long)q->n_checked, (double)q->n_removed/q->n_checked, kh_key(aux, k));
            kh_destroy(pos, q->best_hash);
            free((char*)kh_key(aux, k));
        }
    }
    kh_destroy(lib, aux);

    clear_del_set(del_set);
    kh_destroy(name, del_set);
    free(stack.a);
    bam_destroy1(b);
}

void bam_rmdupse_core(samFile *in, bam_hdr_t *hdr, samFile *out, int force_se);

static int rmdup_usage(void) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:  samtools rmdup [-sS] <input.srt.bam> <output.bam>\n\n");
    fprintf(stderr, "Option: -s    rmdup for SE reads\n");
    fprintf(stderr, "        -S    treat PE reads as SE in rmdup (force -s)\n");

    sam_global_opt_help(stderr, "-....");
    return 1;
}

int bam_rmdup(int argc, char *argv[])
{
    int c, is_se = 0, force_se = 0;
    samFile *in, *out;
    bam_hdr_t *header;
    char wmode[3] = {'w', 'b', 0};
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 0, 0, 0),
        { NULL, 0, NULL, 0 }
    };

    while ((c = getopt_long(argc, argv, "sS", lopts, NULL)) >= 0) {
        switch (c) {
        case 's': is_se = 1; break;
        case 'S': force_se = is_se = 1; break;
        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
        case '?': return rmdup_usage();
        }
    }
    if (optind + 2 > argc)
        return rmdup_usage();

    in = sam_open_format(argv[optind], "r", &ga.in);
    header = sam_hdr_read(in);
    if (header == NULL || header->n_targets == 0) {
        fprintf(stderr, "[bam_rmdup] input SAM does not have header. Abort!\n");
        return 1;
    }

    sam_open_mode(wmode+1, argv[optind+1], NULL);
    out = sam_open_format(argv[optind+1], wmode, &ga.out);
    if (in == 0 || out == 0) {
        fprintf(stderr, "[bam_rmdup] fail to read/write input files\n");
        return 1;
    }
    sam_hdr_write(out, header);

    if (is_se) bam_rmdupse_core(in, header, out, force_se);
    else bam_rmdup_core(in, header, out);
    bam_hdr_destroy(header);
    sam_close(in); sam_close(out);
    return 0;
}
