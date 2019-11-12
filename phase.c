/*  phase.c -- phase subcommand.

    Copyright (C) 2011 Broad Institute.
    Copyright (C) 2013-2016, 2019 Genome Research Ltd.

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
#include <unistd.h>
#include <stdint.h>
#include <math.h>
#include <zlib.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "sam_opts.h"
#include "samtools.h"
#include "htslib/hts_os.h"

#include "htslib/kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

#define MAX_VARS 256
#define FLIP_PENALTY 2
#define FLIP_THRES 4
#define MASK_THRES 3

#define FLAG_FIX_CHIMERA 0x1
#define FLAG_LIST_EXCL   0x4
#define FLAG_DROP_AMBI   0x8

typedef struct {
    // configurations, initialized in the main function
    int flag, k, min_baseQ, min_varLOD, max_depth, no_pg;
    // other global variables
    int vpos_shift;
    samFile* fp;
    sam_hdr_t* fp_hdr;
    char *pre, *arg_list;
    char *out_name[3];
    samFile* out[3];
    sam_hdr_t* out_hdr[3];
    // alignment queue
    int n, m;
    bam1_t **b;
} phaseg_t;

typedef struct {
    int8_t seq[MAX_VARS]; // TODO: change to dynamic memory allocation!
    int vpos, beg, end;
    uint32_t vlen:16, single:1, flip:1, phase:1, phased:1, ambig:1;
    uint32_t in:16, out:16; // in-phase and out-phase
} frag_t, *frag_p;

#define rseq_lt(a,b) ((a)->vpos < (b)->vpos)

#include "htslib/khash.h"
KHASH_SET_INIT_INT64(set64)
KHASH_MAP_INIT_INT64(64, frag_t)

typedef khash_t(64) nseq_t;

#include "htslib/ksort.h"
KSORT_INIT(rseq, frag_p, rseq_lt)

static inline uint64_t X31_hash_string(const char *s)
{
    uint64_t h = *s;
    if (h) for (++s ; *s; ++s) h = (h << 5) - h + *s;
    return h;
}

static void count1(int l, const uint8_t *seq, int *cnt)
{
    int i, j, n_ambi;
    uint32_t z, x;
    if (seq[l-1] == 0) return; // do nothing is the last base is ambiguous
    for (i = n_ambi = 0; i < l; ++i) // collect ambiguous bases
        if (seq[i] == 0) ++n_ambi;
    if (l - n_ambi <= 1) return; // only one SNP
    for (x = 0; x < 1u<<n_ambi; ++x) { // count
        for (i = j = 0, z = 0; i < l; ++i) {
            int c;
            if (seq[i]) c = seq[i] - 1;
            else {
                c = x>>j&1;
                ++j;
            }
            z = z<<1 | c;
        }
        ++cnt[z];
    }
}

static int **count_all(int l, int vpos, nseq_t *hash)
{
    khint_t k;
    int i, j, **cnt;
    uint8_t *seq;
    seq = calloc(l, 1);
    cnt = calloc(vpos, sizeof(int*));
    for (i = 0; i < vpos; ++i) cnt[i] = calloc(1<<l, sizeof(int));
    for (k = 0; k < kh_end(hash); ++k) {
        if (kh_exist(hash, k)) {
            frag_t *f = &kh_val(hash, k);
            if (f->vpos >= vpos || f->single) continue; // out of region; or singleton
            if (f->vlen == 1) { // such reads should be flagged as deleted previously if everything is right
                f->single = 1;
                continue;
            }
            for (j = 1; j < f->vlen; ++j) {
                for (i = 0; i < l; ++i)
                    seq[i] = j < l - 1 - i? 0 : f->seq[j - (l - 1 - i)];
                count1(l, seq, cnt[f->vpos + j]);
            }
        }
    }
    free(seq);
    return cnt;
}

// phasing
static int8_t *dynaprog(int l, int vpos, int **w)
{
    int *f[2], *curr, *prev, max, i;
    int8_t **b, *h = 0;
    uint32_t x, z = 1u<<(l-1), mask = (1u<<l) - 1;
    f[0] = calloc(z, sizeof(int));
    f[1] = calloc(z, sizeof(int));
    b = calloc(vpos, sizeof(int8_t*));
    prev = f[0]; curr = f[1];
    // fill the backtrack matrix
    for (i = 0; i < vpos; ++i) {
        int *wi = w[i], *tmp;
        int8_t *bi;
        bi = b[i] = calloc(z, 1);
        /* In the following, x is the current state, which is the
         * lexicographically smaller local haplotype. xc is the complement of
         * x, or the larger local haplotype; y0 and y1 are the two predecessors
         * of x. */
        for (x = 0; x < z; ++x) { // x0 is the smaller
            uint32_t y0, y1, xc;
            int c0, c1;
            xc = ~x&mask; y0 = x>>1; y1 = xc>>1;
            c0 = prev[y0] + wi[x] + wi[xc];
            c1 = prev[y1] + wi[x] + wi[xc];
            if (c0 > c1) bi[x] = 0, curr[x] = c0;
            else bi[x] = 1, curr[x] = c1;
        }
        tmp = prev; prev = curr; curr = tmp; // swap
    }
    { // backtrack
        uint32_t max_x = 0;
        int which = 0;
        h = calloc(vpos, 1);
        for (x = 0, max = 0, max_x = 0; x < z; ++x)
            if (prev[x] > max) max = prev[x], max_x = x;
        for (i = vpos - 1, x = max_x; i >= 0; --i) {
            h[i] = which? (~x&1) : (x&1);
            which = b[i][x]? !which : which;
            x = b[i][x]? (~x&mask)>>1 : x>>1;
        }
    }
    // free
    for (i = 0; i < vpos; ++i) free(b[i]);
    free(f[0]); free(f[1]); free(b);
    return h;
}

// phase each fragment
static uint64_t *fragphase(int vpos, const int8_t *path, nseq_t *hash, int flip)
{
    khint_t k;
    uint64_t *pcnt;
    uint32_t *left, *rght, max;
    left = rght = 0; max = 0;
    pcnt = calloc(vpos, 8);
    for (k = 0; k < kh_end(hash); ++k) {
        if (kh_exist(hash, k)) {
            int i, c[2];
            frag_t *f = &kh_val(hash, k);
            if (f->vpos >= vpos) continue;
            // get the phase
            c[0] = c[1] = 0;
            for (i = 0; i < f->vlen; ++i) {
                if (f->seq[i] == 0) continue;
                ++c[f->seq[i] == path[f->vpos + i] + 1? 0 : 1];
            }
            f->phase = c[0] > c[1]? 0 : 1;
            f->in = c[f->phase]; f->out = c[1 - f->phase];
            f->phased = f->in == f->out? 0 : 1;
            f->ambig = (f->in && f->out && f->out < 3 && f->in <= f->out + 1)? 1 : 0;
            // fix chimera
            f->flip = 0;
            if (flip && c[0] >= 3 && c[1] >= 3) {
                int sum[2], m, mi, md;
                if (f->vlen > max) { // enlarge the array
                    max = f->vlen;
                    kroundup32(max);
                    left = realloc(left, max * 4);
                    rght = realloc(rght, max * 4);
                }
                for (i = 0, sum[0] = sum[1] = 0; i < f->vlen; ++i) { // get left counts
                    if (f->seq[i]) {
                        int c = f->phase? 2 - f->seq[i] : f->seq[i] - 1;
                        ++sum[c == path[f->vpos + i]? 0 : 1];
                    }
                    left[i] = sum[1]<<16 | sum[0];
                }
                for (i = f->vlen - 1, sum[0] = sum[1] = 0; i >= 0; --i) { // get right counts
                    if (f->seq[i]) {
                        int c = f->phase? 2 - f->seq[i] : f->seq[i] - 1;
                        ++sum[c == path[f->vpos + i]? 0 : 1];
                    }
                    rght[i] = sum[1]<<16 | sum[0];
                }
                // find the best flip point
                for (i = m = 0, mi = -1, md = -1; i < f->vlen - 1; ++i) {
                    int a[2];
                    a[0] = (left[i]&0xffff) + (rght[i+1]>>16&0xffff) - (rght[i+1]&0xffff) * FLIP_PENALTY;
                    a[1] = (left[i]>>16&0xffff) + (rght[i+1]&0xffff) - (rght[i+1]>>16&0xffff) * FLIP_PENALTY;
                    if (a[0] > a[1]) {
                        if (a[0] > m) m = a[0], md = 0, mi = i;
                    } else {
                        if (a[1] > m) m = a[1], md = 1, mi = i;
                    }
                }
                if (m - c[0] >= FLIP_THRES && m - c[1] >= FLIP_THRES) { // then flip
                    f->flip = 1;
                    if (md == 0) { // flip the tail
                        for (i = mi + 1; i < f->vlen; ++i)
                            if (f->seq[i] == 1) f->seq[i] = 2;
                            else if (f->seq[i] == 2) f->seq[i] = 1;
                    } else { // flip the head
                        for (i = 0; i <= mi; ++i)
                            if (f->seq[i] == 1) f->seq[i] = 2;
                            else if (f->seq[i] == 2) f->seq[i] = 1;
                    }
                }
            }
            // update pcnt[]
            if (!f->single) {
                for (i = 0; i < f->vlen; ++i) {
                    int c;
                    if (f->seq[i] == 0) continue;
                    c = f->phase? 2 - f->seq[i] : f->seq[i] - 1;
                    if (c == path[f->vpos + i]) {
                        if (f->phase == 0) ++pcnt[f->vpos + i];
                        else pcnt[f->vpos + i] += 1ull<<32;
                    } else {
                        if (f->phase == 0) pcnt[f->vpos + i] += 1<<16;
                        else pcnt[f->vpos + i] += 1ull<<48;
                    }
                }
            }
        }
    }
    free(left); free(rght);
    return pcnt;
}

static uint64_t *genmask(int vpos, const uint64_t *pcnt, int *_n)
{
    int i, max = 0, max_i = -1, m = 0, n = 0, beg = 0, score = 0;
    uint64_t *list = 0;
    for (i = 0; i < vpos; ++i) {
        uint64_t x = pcnt[i];
        int c[4], pre = score, s;
        c[0] = x&0xffff; c[1] = x>>16&0xffff; c[2] = x>>32&0xffff; c[3] = x>>48&0xffff;
        s = (c[1] + c[3] == 0)? -(c[0] + c[2]) : (c[1] + c[3] - 1);
        if (c[3] > c[2]) s += c[3] - c[2];
        if (c[1] > c[0]) s += c[1] - c[0];
        score += s;
        if (score < 0) score = 0;
        if (pre == 0 && score > 0) beg = i; // change from zero to non-zero
        if ((i == vpos - 1 || score == 0) && max >= MASK_THRES) {
            if (n == m) {
                m = m? m<<1 : 4;
                list = realloc(list, m * 8);
            }
            list[n++] = (uint64_t)beg<<32 | max_i;
            i = max_i; // reset i to max_i
            score = 0;
        } else if (score > max) max = score, max_i = i;
        if (score == 0) max = 0;
    }
    *_n = n;
    return list;
}

// trim heading and tailing ambiguous bases; mark deleted and remove sequence
static int clean_seqs(int vpos, nseq_t *hash)
{
    khint_t k;
    int ret = 0;
    for (k = 0; k < kh_end(hash); ++k) {
        if (kh_exist(hash, k)) {
            frag_t *f = &kh_val(hash, k);
            int beg, end, i;
            if (f->vpos >= vpos) {
                ret = 1;
                continue;
            }
            for (i = 0; i < f->vlen; ++i)
                if (f->seq[i] != 0) break;
            beg = i;
            for (i = f->vlen - 1; i >= 0; --i)
                if (f->seq[i] != 0) break;
            end = i + 1;
            if (end - beg <= 0) kh_del(64, hash, k);
            else {
                if (beg != 0) memmove(f->seq, f->seq + beg, end - beg);
                f->vpos += beg; f->vlen = end - beg;
                f->single = f->vlen == 1? 1 : 0;
            }
        }
    }
    return ret;
}

static int dump_aln(phaseg_t *g, int min_pos, const nseq_t *hash)
{
    int i, is_flip, drop_ambi;
    drop_ambi = g->flag & FLAG_DROP_AMBI;
    is_flip = (drand48() < 0.5);
    for (i = 0; i < g->n; ++i) {
        int end, which;
        uint64_t key;
        khint_t k;
        bam1_t *b = g->b[i];
        key = X31_hash_string(bam_get_qname(b));
        end = bam_endpos(b);
        if (end > min_pos) break;
        k = kh_get(64, hash, key);
        if (k == kh_end(hash)) which = 3;
        else {
            frag_t *f = &kh_val(hash, k);
            if (f->ambig) which = drop_ambi? 2 : 3;
            else if (f->phased && f->flip) which = 2;
            else if (f->phased == 0) which = 3;
            else { // phased and not flipped
                char c = 'Y';
                which = f->phase;
                bam_aux_append(b, "ZP", 'A', 1, (uint8_t*)&c);
            }
            if (which < 2 && is_flip) which = 1 - which; // increase the randomness
        }
        if (which == 3) which = (drand48() < 0.5);
        if (sam_write1(g->out[which], g->out_hdr[which], b) < 0) {
            print_error_errno("phase", "error writing to '%s'", g->out_name[which]);
            return -1;
        }
        bam_destroy1(b);
        g->b[i] = 0;
    }
    memmove(g->b, g->b + i, (g->n - i) * sizeof(void*));
    g->n -= i;
    return 0;
}

static int phase(phaseg_t *g, const char *chr, int vpos, uint64_t *cns, nseq_t *hash)
{
    int i, j, n_seqs = kh_size(hash), n_masked = 0, min_pos;
    khint_t k;
    frag_t **seqs;
    int8_t *path, *sitemask;
    uint64_t *pcnt, *regmask;

    if (vpos == 0) return 0;
    i = clean_seqs(vpos, hash); // i is true if hash has an element with its vpos >= vpos
    min_pos = i? cns[vpos]>>32 : 0x7fffffff;
    if (vpos == 1) {
        printf("PS\t%s\t%d\t%d\n", chr, (int)(cns[0]>>32) + 1, (int)(cns[0]>>32) + 1);
        printf("M0\t%s\t%d\t%d\t%c\t%c\t%d\t0\t0\t0\t0\n//\n", chr, (int)(cns[0]>>32) + 1, (int)(cns[0]>>32) + 1,
            "ACGTX"[cns[0]&3], "ACGTX"[cns[0]>>16&3], g->vpos_shift + 1);
        for (k = 0; k < kh_end(hash); ++k) {
            if (kh_exist(hash, k)) {
                frag_t *f = &kh_val(hash, k);
                if (f->vpos) continue;
                f->flip = 0;
                if (f->seq[0] == 0) f->phased = 0;
                else f->phased = 1, f->phase = f->seq[0] - 1;
            }
        }
        if (dump_aln(g, min_pos, hash) < 0) return -1;
        ++g->vpos_shift;
        return 1;
    }
    { // phase
        int **cnt;
        uint64_t *mask;
        printf("PS\t%s\t%d\t%d\n", chr, (int)(cns[0]>>32) + 1, (int)(cns[vpos-1]>>32) + 1);
        sitemask = calloc(vpos, 1);
        cnt = count_all(g->k, vpos, hash);
        path = dynaprog(g->k, vpos, cnt);
        for (i = 0; i < vpos; ++i) free(cnt[i]);
        free(cnt);
        pcnt = fragphase(vpos, path, hash, 0); // do not fix chimeras when masking
        mask = genmask(vpos, pcnt, &n_masked);
        regmask = calloc(n_masked, 8);
        for (i = 0; i < n_masked; ++i) {
            regmask[i] = cns[mask[i]>>32]>>32<<32 | cns[(uint32_t)mask[i]]>>32;
            for (j = mask[i]>>32; j <= (int32_t)mask[i]; ++j)
                sitemask[j] = 1;
        }
        free(mask);
        if (g->flag & FLAG_FIX_CHIMERA) {
            free(pcnt);
            pcnt = fragphase(vpos, path, hash, 1);
        }
    }
    for (i = 0; i < n_masked; ++i)
        printf("FL\t%s\t%d\t%d\n", chr, (int)(regmask[i]>>32) + 1, (int)regmask[i] + 1);
    for (i = 0; i < vpos; ++i) {
        uint64_t x = pcnt[i];
        int8_t c[2];
        c[0] = (cns[i]&0xffff)>>2 == 0? 4 : (cns[i]&3);
        c[1] = (cns[i]>>16&0xffff)>>2 == 0? 4 : (cns[i]>>16&3);
        printf("M%d\t%s\t%d\t%d\t%c\t%c\t%d\t%d\t%d\t%d\t%d\n", sitemask[i]+1, chr, (int)(cns[0]>>32) + 1, (int)(cns[i]>>32) + 1, "ACGTX"[c[path[i]]], "ACGTX"[c[1-path[i]]],
            i + g->vpos_shift + 1, (int)(x&0xffff), (int)(x>>16&0xffff), (int)(x>>32&0xffff), (int)(x>>48&0xffff));
    }
    free(path); free(pcnt); free(regmask); free(sitemask);
    seqs = calloc(n_seqs, sizeof(frag_t*));
    for (k = 0, i = 0; k < kh_end(hash); ++k)
        if (kh_exist(hash, k) && kh_val(hash, k).vpos < vpos && !kh_val(hash, k).single)
            seqs[i++] = &kh_val(hash, k);
    n_seqs = i;
    ks_introsort_rseq(n_seqs, seqs);
    for (i = 0; i < n_seqs; ++i) {
        frag_t *f = seqs[i];
        printf("EV\t0\t%s\t%d\t40\t%dM\t*\t0\t0\t", chr, f->vpos + 1 + g->vpos_shift, f->vlen);
        for (j = 0; j < f->vlen; ++j) {
            uint32_t c = cns[f->vpos + j];
            if (f->seq[j] == 0) putchar('N');
            else putchar("ACGT"[f->seq[j] == 1? (c&3) : (c>>16&3)]);
        }
        printf("\t*\tYP:i:%d\tYF:i:%d\tYI:i:%d\tYO:i:%d\tYS:i:%d\n", f->phase, f->flip, f->in, f->out, f->beg+1);
    }
    free(seqs);
    printf("//\n");
    fflush(stdout);
    g->vpos_shift += vpos;
    if (dump_aln(g, min_pos, hash) < 0) return -1;
    return vpos;
}

static void update_vpos(int vpos, nseq_t *hash)
{
    khint_t k;
    for (k = 0; k < kh_end(hash); ++k) {
        if (kh_exist(hash, k)) {
            frag_t *f = &kh_val(hash, k);
            if (f->vpos < vpos) kh_del(64, hash, k); // TODO: if frag_t::seq is allocated dynamically, free it
            else f->vpos -= vpos;
        }
    }
}

static nseq_t *shrink_hash(nseq_t *hash) // TODO: to implement
{
    return hash;
}

static int readaln(void *data, bam1_t *b)
{
    phaseg_t *g = (phaseg_t*)data;
    int ret;
    while (1)
    {
        ret = sam_read1(g->fp, g->fp_hdr, b);
        if (ret < 0) break;
        if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
        if ( g->pre ) {
            if (g->n == g->m) {
                g->m = g->m? g->m<<1 : 16;
                g->b = realloc(g->b, g->m * sizeof(bam1_t*));
            }
            g->b[g->n++] = bam_dup1(b);
        }
        break;
    }
    return ret;
}

static khash_t(set64) *loadpos(const char *fn, sam_hdr_t *h)
{
    gzFile fp;
    kstream_t *ks;
    int ret, dret;
    kstring_t *str;
    khash_t(set64) *hash;

    fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
    if (fp == NULL) {
        print_error_errno("phase", "Couldn't open site file '%s'", fn);
        return NULL;
    }

    hash = kh_init(set64);
    str = calloc(1, sizeof(kstring_t));

    ks = ks_init(fp);
    while (ks_getuntil(ks, 0, str, &dret) >= 0) {
        int tid = bam_name2id(h, str->s);
        if (tid >= 0 && dret != '\n') {
            if (ks_getuntil(ks, 0, str, &dret) >= 0) {
                uint64_t x = (uint64_t)tid<<32 | (atoi(str->s) - 1);
                kh_put(set64, hash, x, &ret);
            } else break;
        }
        if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
        if (dret < 0) break;
    }
    ks_destroy(ks);
    gzclose(fp);
    free(str->s); free(str);
    return hash;
}

static int gl2cns(float q[16])
{
    int i, j, min_ij;
    float min, min2;
    min = min2 = 1e30; min_ij = -1;
    for (i = 0; i < 4; ++i) {
        for (j = i; j < 4; ++j) {
            if (q[i<<2|j] < min) min_ij = i<<2|j, min2 = min, min = q[i<<2|j];
            else if (q[i<<2|j] < min2) min2 = q[i<<2|j];
        }
    }
    return (min_ij>>2&3) == (min_ij&3)? 0 : 1<<18 | (min_ij>>2&3)<<16 | (min_ij&3) | (int)(min2 - min + .499) << 2;
}

static int start_output(phaseg_t *g, int c, const char *middle, const htsFormat *fmt)
{
    kstring_t s = { 0, 0, NULL };
    ksprintf(&s, "%s.%s.%s", g->pre, middle, hts_format_file_extension(fmt));
    g->out_name[c] = ks_release(&s);
    g->out[c] = sam_open_format(g->out_name[c], "wb", fmt);
    if (! g->out[c]) {
        print_error_errno("phase", "Failed to open output file '%s'", g->out_name[c]);
        return -1;
    }

    g->out_hdr[c] = sam_hdr_dup(g->fp_hdr);
    if (!g->no_pg && sam_hdr_add_pg(g->out_hdr[c], "samtools",
                                    "VN", samtools_version(),
                                    g->arg_list ? "CL": NULL,
                                    g->arg_list ? g->arg_list : NULL,
                                    NULL)) {
        print_error("phase", "failed to add PG line to header");
        return -1;
    }
    if (sam_hdr_write(g->out[c], g->out_hdr[c]) < 0) {
        print_error_errno("phase", "Failed to write header for '%s'", g->out_name[c]);
        return -1;
    }

    return 0;
}

int main_phase(int argc, char *argv[])
{
    int c, tid, pos, vpos = 0, n, lasttid = -1, max_vpos = 0, usage = 0;
    const bam_pileup1_t *plp;
    bam_plp_t iter;
    nseq_t *seqs;
    uint64_t *cns = 0;
    phaseg_t g;
    char *fn_list = 0;
    khash_t(set64) *set = 0;
    errmod_t *em;
    uint16_t *bases;

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 0, 0, 0, '-'),
        {"no-PG", no_argument, NULL, 1},
        { NULL, 0, NULL, 0 }
    };

    // FIXME Leaks galore in the case of error returns

    memset(&g, 0, sizeof(phaseg_t));
    g.flag = FLAG_FIX_CHIMERA;
    g.min_varLOD = 37; g.k = 13; g.min_baseQ = 13; g.max_depth = 256;
    while ((c = getopt_long(argc, argv, "Q:eFq:k:b:l:D:A", lopts, NULL)) >= 0) {
        switch (c) {
            case 'D': g.max_depth = atoi(optarg); break;
            case 'q': g.min_varLOD = atoi(optarg); break;
            case 'Q': g.min_baseQ = atoi(optarg); break;
            case 'k': g.k = atoi(optarg); break;
            case 'F': g.flag &= ~FLAG_FIX_CHIMERA; break;
            case 'e': g.flag |= FLAG_LIST_EXCL; break;
            case 'A': g.flag |= FLAG_DROP_AMBI; break;
            case 'b': g.pre = strdup(optarg); break;
            case 'l': fn_list = strdup(optarg); break;
            case 1: g.no_pg = 1; break;
            default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
                      /* else fall-through */
            case '?': usage=1; break;
        }
        if (usage) break;
    }
    if (usage || argc == optind) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   samtools phase [options] <in.bam>\n\n");
        fprintf(stderr, "Options: -k INT    block length [%d]\n", g.k);
        fprintf(stderr, "         -b STR    prefix of BAMs to output [null]\n");
        fprintf(stderr, "         -q INT    min het phred-LOD [%d]\n", g.min_varLOD);
        fprintf(stderr, "         -Q INT    min base quality in het calling [%d]\n", g.min_baseQ);
        fprintf(stderr, "         -D INT    max read depth [%d]\n", g.max_depth);
//      fprintf(stderr, "         -l FILE   list of sites to phase [null]\n");
        fprintf(stderr, "         -F        do not attempt to fix chimeras\n");
        fprintf(stderr, "         -A        drop reads with ambiguous phase\n");
        fprintf(stderr, "         --no-PG   do not add a PG line\n");
//      fprintf(stderr, "         -e        do not discover SNPs (effective with -l)\n");
        fprintf(stderr, "\n");

        sam_global_opt_help(stderr, "-....--.");

        return 1;
    }
    g.fp = sam_open_format(argv[optind], "r", &ga.in);
    if (!g.fp) {
        print_error_errno("phase", "Couldn't open '%s'", argv[optind]);
        return 1;
    }
    g.fp_hdr = sam_hdr_read(g.fp);
    if (g.fp_hdr == NULL) {
        fprintf(stderr, "[%s] Failed to read header for '%s'\n",
                __func__, argv[optind]);
        return 1;
    }
    if (!g.no_pg && !(g.arg_list = stringify_argv(argc+1, argv-1))) {
        print_error("phase", "failed to create arg_list");
        return 1;
    }
    if (fn_list) { // read the list of sites to phase
        set = loadpos(fn_list, g.fp_hdr);
        if (set == NULL) return 1;
        free(fn_list);
    } else g.flag &= ~FLAG_LIST_EXCL;
    if (g.pre) { // open BAMs to write
        if (ga.out.format == unknown_format)
            ga.out.format = bam; // default via "wb".

        // Open each output file g.out[0..2], dupping and writing the header
        if (start_output(&g, 0, "0", &ga.out) < 0 ||
            start_output(&g, 1, "1", &ga.out) < 0 ||
            start_output(&g, 2, "chimera", &ga.out) < 0) return 1;
    }

    iter = bam_plp_init(readaln, &g);
    g.vpos_shift = 0;
    seqs = kh_init(64);
    em = errmod_init(1. - 0.83);
    bases = calloc(g.max_depth, 2);
    printf("CC\n");
    printf("CC\tDescriptions:\nCC\n");
    printf("CC\t  CC      comments\n");
    printf("CC\t  PS      start of a phase set\n");
    printf("CC\t  FL      filtered region\n");
    printf("CC\t  M[012]  markers; 0 for singletons, 1 for phased and 2 for filtered\n");
    printf("CC\t  EV      supporting reads; SAM format\n");
    printf("CC\t  //      end of a phase set\nCC\n");
    printf("CC\tFormats of PS, FL and M[012] lines (1-based coordinates):\nCC\n");
    printf("CC\t  PS  chr  phaseSetStart  phaseSetEnd\n");
    printf("CC\t  FL  chr  filterStart    filterEnd\n");
    printf("CC\t  M?  chr  PS  pos  allele0  allele1  hetIndex  #supports0  #errors0  #supp1  #err1\n");
    printf("CC\nCC\n");
    fflush(stdout);
    while ((plp = bam_plp_auto(iter, &tid, &pos, &n)) != 0) {
        int i, k, c, tmp, dophase = 1, in_set = 0;
        float q[16];
        if (tid < 0) break;
        if (tid != lasttid) { // change of chromosome
            g.vpos_shift = 0;
            if (lasttid >= 0) {
                seqs = shrink_hash(seqs);
                if (phase(&g, sam_hdr_tid2name(g.fp_hdr, lasttid),
                          vpos, cns, seqs) < 0) {
                    return 1;
                }
                update_vpos(0x7fffffff, seqs);
            }
            lasttid = tid;
            vpos = 0;
        }
        if (set && kh_get(set64, set, (uint64_t)tid<<32 | pos) != kh_end(set)) in_set = 1;
        if (n > g.max_depth) continue; // do not proceed if the depth is too high
        // fill the bases array and check if there is a variant
        for (i = k = 0; i < n; ++i) {
            const bam_pileup1_t *p = plp + i;
            uint8_t *seq;
            int q, baseQ, b;
            if (p->is_del || p->is_refskip) continue;
            baseQ = bam_get_qual(p->b)[p->qpos];
            if (baseQ < g.min_baseQ) continue;
            seq = bam_get_seq(p->b);
            b = seq_nt16_int[bam_seqi(seq, p->qpos)];
            if (b > 3) continue;
            q = baseQ < p->b->core.qual? baseQ : p->b->core.qual;
            if (q < 4) q = 4;
            if (q > 63) q = 63;
            bases[k++] = q<<5 | (int)bam_is_rev(p->b)<<4 | b;
        }
        if (k == 0) continue;
        errmod_cal(em, k, 4, bases, q); // compute genotype likelihood
        c = gl2cns(q); // get the consensus
        // tell if to proceed
        if (set && (g.flag&FLAG_LIST_EXCL) && !in_set) continue; // not in the list
        if (!in_set && (c&0xffff)>>2 < g.min_varLOD) continue; // not a variant
        // add the variant
        if (vpos == max_vpos) {
            max_vpos = max_vpos? max_vpos<<1 : 128;
            cns = realloc(cns, max_vpos * 8);
        }
        cns[vpos] = (uint64_t)pos<<32 | c;
        for (i = 0; i < n; ++i) {
            const bam_pileup1_t *p = plp + i;
            uint64_t key;
            khint_t k;
            uint8_t *seq = bam_get_seq(p->b);
            frag_t *f;
            if (p->is_del || p->is_refskip) continue;
            if (p->b->core.qual == 0) continue;
            // get the base code
            c = seq_nt16_int[bam_seqi(seq, p->qpos)];
            if (c == (cns[vpos]&3)) c = 1;
            else if (c == (cns[vpos]>>16&3)) c = 2;
            else c = 0;
            // write to seqs
            key = X31_hash_string(bam_get_qname(p->b));
            k = kh_put(64, seqs, key, &tmp);
            f = &kh_val(seqs, k);
            if (tmp == 0) { // present in the hash table
                if (vpos - f->vpos + 1 < MAX_VARS) {
                    f->vlen = vpos - f->vpos + 1;
                    f->seq[f->vlen-1] = c;
                    f->end = bam_endpos(p->b);
                }
                dophase = 0;
            } else { // absent
                memset(f->seq, 0, MAX_VARS);
                f->beg = p->b->core.pos;
                f->end = bam_endpos(p->b);
                f->vpos = vpos, f->vlen = 1, f->seq[0] = c, f->single = f->phased = f->flip = f->ambig = 0;
            }
        }
        if (dophase) {
            seqs = shrink_hash(seqs);
            if (phase(&g, sam_hdr_tid2name(g.fp_hdr, tid), vpos, cns, seqs) < 0) {
                return 1;
            }
            update_vpos(vpos, seqs);
            cns[0] = cns[vpos];
            vpos = 0;
        }
        ++vpos;
    }
    if (tid >= 0) {
        if (phase(&g, sam_hdr_tid2name(g.fp_hdr, tid), vpos, cns, seqs) < 0) {
            return 1;
        }
    }
    sam_hdr_destroy(g.fp_hdr);
    bam_plp_destroy(iter);
    sam_close(g.fp);
    kh_destroy(64, seqs);
    kh_destroy(set64, set);
    free(cns);
    errmod_destroy(em);
    free(bases);
    if (g.pre) {
        int res = 0;
        for (c = 0; c <= 2; ++c) {
            if (sam_close(g.out[c]) < 0) {
                fprintf(stderr, "[%s] error on closing '%s'\n",
                        __func__, g.out_name[c]);
                res = 1;
            }
            sam_hdr_destroy(g.out_hdr[c]);
            free(g.out_name[c]);
        }
        free(g.pre); free(g.b);
        if (res) return 1;
    }
    free(g.arg_list);
    sam_global_args_free(&ga);
    return 0;
}
