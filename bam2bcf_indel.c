/*  bam2bcf_indel.c -- indel caller.

    Copyright (C) 2010, 2011 Broad Institute.
    Copyright (C) 2012-2014 Genome Research Ltd.

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

#include <assert.h>
#include <ctype.h>
#include <string.h>
#include "htslib/sam.h"
#include "bam2bcf.h"
#include "kprobaln.h"
#include "htslib/khash.h"
KHASH_SET_INIT_STR(rg)

#include "htslib/ksort.h"
KSORT_INIT_GENERIC(uint32_t)

#define MINUS_CONST 0x10000000
#define INDEL_WINDOW_SIZE 50

void *bcf_call_add_rg(void *_hash, const char *hdtext, const char *list)
{
    const char *s, *p, *q, *r, *t;
    khash_t(rg) *hash;
    if (list == 0 || hdtext == 0) return _hash;
    if (_hash == 0) _hash = kh_init(rg);
    hash = (khash_t(rg)*)_hash;
    if ((s = strstr(hdtext, "@RG\t")) == 0) return hash;
    do {
        t = strstr(s + 4, "@RG\t"); // the next @RG
        if ((p = strstr(s, "\tID:")) != 0) p += 4;
        if ((q = strstr(s, "\tPL:")) != 0) q += 4;
        if (p && q && (t == 0 || (p < t && q < t))) { // ID and PL are both present
            int lp, lq;
            char *x;
            for (r = p; *r && *r != '\t' && *r != '\n'; ++r) { }
            lp = r - p;
            for (r = q; *r && *r != '\t' && *r != '\n'; ++r) { }
            lq = r - q;
            x = calloc((lp > lq? lp : lq) + 1, 1);
            for (r = q; *r && *r != '\t' && *r != '\n'; ++r) x[r-q] = *r;
            if (strstr(list, x)) { // insert ID to the hash table
                khint_t k;
                int ret;
                for (r = p; *r && *r != '\t' && *r != '\n'; ++r) x[r-p] = *r;
                x[r-p] = 0;
                k = kh_get(rg, hash, x);
                if (k == kh_end(hash)) k = kh_put(rg, hash, x, &ret);
                else free(x);
            } else free(x);
        }
        s = t;
    } while (s);
    return hash;
}

void bcf_call_del_rghash(void *_hash)
{
    khint_t k;
    khash_t(rg) *hash = (khash_t(rg)*)_hash;
    if (hash == 0) return;
    for (k = kh_begin(hash); k < kh_end(hash); ++k)
        if (kh_exist(hash, k))
            free((char*)kh_key(hash, k));
    kh_destroy(rg, hash);
}

static int tpos2qpos(const bam1_core_t *c, const uint32_t *cigar, int32_t tpos, int is_left, int32_t *_tpos)
{
    int k, x = c->pos, y = 0, last_y = 0;
    *_tpos = c->pos;
    for (k = 0; k < c->n_cigar; ++k) {
        int op = cigar[k] & BAM_CIGAR_MASK;
        int l = cigar[k] >> BAM_CIGAR_SHIFT;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (c->pos > tpos) return y;
            if (x + l > tpos) {
                *_tpos = tpos;
                return y + (tpos - x);
            }
            x += l; y += l;
            last_y = y;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) y += l;
        else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            if (x + l > tpos) {
                *_tpos = is_left? x : x + l;
                return y;
            }
            x += l;
        }
    }
    *_tpos = x;
    return last_y;
}
// FIXME: check if the inserted sequence is consistent with the homopolymer run
// l is the relative gap length and l_run is the length of the homopolymer on the reference
static inline int est_seqQ(const bcf_callaux_t *bca, int l, int l_run)
{
    int q, qh;
    q = bca->openQ + bca->extQ * (abs(l) - 1);
    qh = l_run >= 3? (int)(bca->tandemQ * (double)abs(l) / l_run + .499) : 1000;
    return q < qh? q : qh;
}

static inline int est_indelreg(int pos, const char *ref, int l, char *ins4)
{
    int i, j, max = 0, max_i = pos, score = 0;
    l = abs(l);
    for (i = pos + 1, j = 0; ref[i]; ++i, ++j) {
        if (ins4) score += (toupper(ref[i]) != "ACGTN"[(int)ins4[j%l]])? -10 : 1;
        else score += (toupper(ref[i]) != toupper(ref[pos+1+j%l]))? -10 : 1;
        if (score < 0) break;
        if (max < score) max = score, max_i = i;
    }
    return max_i - pos;
}

/*
    notes:
        - n .. number of samples
        - the routine sets bam_pileup1_t.aux of each read as follows:
            - 6: unused
            - 6: the call; index to bcf_callaux_t.indel_types   .. (aux>>16)&0x3f
            - 8: estimated sequence quality                     .. (aux>>8)&0xff
            - 8: indel quality                                  .. aux&0xff
 */
int bcf_call_gap_prep(int n, int *n_plp, bam_pileup1_t **plp, int pos, bcf_callaux_t *bca, const char *ref,
                      const void *rghash)
{
    int i, s, j, k, t, n_types, *types, max_rd_len, left, right, max_ins, *score1, *score2, max_ref2;
    int N, K, l_run, ref_type, n_alt;
    char *inscns = 0, *ref2, *query, **ref_sample;
    khash_t(rg) *hash = (khash_t(rg)*)rghash;
    if (ref == 0 || bca == 0) return -1;
    // mark filtered reads
    if (rghash) {
        N = 0;
        for (s = N = 0; s < n; ++s) {
            for (i = 0; i < n_plp[s]; ++i) {
                bam_pileup1_t *p = plp[s] + i;
                const uint8_t *rg = bam_aux_get(p->b, "RG");
                p->aux = 1; // filtered by default
                if (rg) {
                    khint_t k = kh_get(rg, hash, (const char*)(rg + 1));
                    if (k != kh_end(hash)) p->aux = 0, ++N; // not filtered
                }
            }
        }
        if (N == 0) return -1; // no reads left
    }
    // determine if there is a gap
    for (s = N = 0; s < n; ++s) {
        for (i = 0; i < n_plp[s]; ++i)
            if (plp[s][i].indel != 0) break;
        if (i < n_plp[s]) break;
    }
    if (s == n) return -1; // there is no indel at this position.
    for (s = N = 0; s < n; ++s) N += n_plp[s]; // N is the total number of reads
    { // find out how many types of indels are present
        bca->max_support = bca->max_frac = 0;
        int m, n_alt = 0, n_tot = 0, indel_support_ok = 0;
        uint32_t *aux;
        aux = calloc(N + 1, 4);
        m = max_rd_len = 0;
        aux[m++] = MINUS_CONST; // zero indel is always a type
        for (s = 0; s < n; ++s) {
            int na = 0, nt = 0;
            for (i = 0; i < n_plp[s]; ++i) {
                const bam_pileup1_t *p = plp[s] + i;
                if (rghash == 0 || p->aux == 0) {
                    ++nt;
                    if (p->indel != 0) {
                        ++na;
                        aux[m++] = MINUS_CONST + p->indel;
                    }
                }
                j = bam_cigar2qlen(p->b->core.n_cigar, bam_get_cigar(p->b));
                if (j > max_rd_len) max_rd_len = j;
            }
            double frac = (double)na/nt;
            if ( !indel_support_ok && na >= bca->min_support && frac >= bca->min_frac )
                indel_support_ok = 1;
            if ( na > bca->max_support && frac > 0 ) bca->max_support = na, bca->max_frac = frac;
            n_alt += na;
            n_tot += nt;
        }
        // To prevent long stretches of N's to be mistaken for indels (sometimes thousands of bases),
        //  check the number of N's in the sequence and skip places where half or more reference bases are Ns.
        int nN=0; for (i=pos; i-pos<max_rd_len && ref[i]; i++) if ( ref[i]=='N' ) nN++;
        if ( nN*2>(i-pos) ) { free(aux); return -1; }

        ks_introsort(uint32_t, m, aux);
        // squeeze out identical types
        for (i = 1, n_types = 1; i < m; ++i)
            if (aux[i] != aux[i-1]) ++n_types;
        // Taking totals makes it hard to call rare indels
        if ( !bca->per_sample_flt )
            indel_support_ok = ( (double)n_alt / n_tot < bca->min_frac || n_alt < bca->min_support ) ? 0 : 1;
        if ( n_types == 1 || !indel_support_ok ) { // then skip
            free(aux); return -1;
        }
        if (n_types >= 64) {
            free(aux);
            // TODO revisit how/whether to control printing this warning
            if (hts_verbose >= 2)
                fprintf(stderr, "[%s] excessive INDEL alleles at position %d. Skip the position.\n", __func__, pos + 1);
            return -1;
        }
        types = (int*)calloc(n_types, sizeof(int));
        t = 0;
        types[t++] = aux[0] - MINUS_CONST;
        for (i = 1; i < m; ++i)
            if (aux[i] != aux[i-1])
                types[t++] = aux[i] - MINUS_CONST;
        free(aux);
        for (t = 0; t < n_types; ++t)
            if (types[t] == 0) break;
        ref_type = t; // the index of the reference type (0)
    }
    { // calculate left and right boundary
        left = pos > INDEL_WINDOW_SIZE? pos - INDEL_WINDOW_SIZE : 0;
        right = pos + INDEL_WINDOW_SIZE;
        if (types[0] < 0) right -= types[0];
        // in case the alignments stand out the reference
        for (i = pos; i < right; ++i)
            if (ref[i] == 0) break;
        right = i;
    }
    /* The following block fixes a long-existing flaw in the INDEL
     * calling model: the interference of nearby SNPs. However, it also
     * reduces the power because sometimes, substitutions caused by
     * indels are not distinguishable from true mutations. Multiple
     * sequence realignment helps to increase the power.
     *
     * Masks mismatches present in at least 70% of the reads with 'N'.
     */
    { // construct per-sample consensus
        int L = right - left + 1, max_i, max2_i;
        uint32_t *cns, max, max2;
        char *ref0, *r;
        ref_sample = calloc(n, sizeof(char*));
        cns = calloc(L, 4);
        ref0 = calloc(L, 1);
        for (i = 0; i < right - left; ++i)
            ref0[i] = seq_nt16_table[(int)ref[i+left]];
        for (s = 0; s < n; ++s) {
            r = ref_sample[s] = calloc(L, 1);
            memset(cns, 0, sizeof(int) * L);
            // collect ref and non-ref counts
            for (i = 0; i < n_plp[s]; ++i) {
                bam_pileup1_t *p = plp[s] + i;
                bam1_t *b = p->b;
                uint32_t *cigar = bam_get_cigar(b);
                uint8_t *seq = bam_get_seq(b);
                int x = b->core.pos, y = 0;
                for (k = 0; k < b->core.n_cigar; ++k) {
                    int op = cigar[k]&0xf;
                    int j, l = cigar[k]>>4;
                    if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                        for (j = 0; j < l; ++j)
                            if (x + j >= left && x + j < right)
                                cns[x+j-left] += (bam_seqi(seq, y+j) == ref0[x+j-left])? 1 : 0x10000;
                        x += l; y += l;
                    } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) x += l;
                    else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) y += l;
                }
            }
            // determine the consensus
            for (i = 0; i < right - left; ++i) r[i] = ref0[i];
            max = max2 = 0; max_i = max2_i = -1;
            for (i = 0; i < right - left; ++i) {
                if (cns[i]>>16 >= max>>16) max2 = max, max2_i = max_i, max = cns[i], max_i = i;
                else if (cns[i]>>16 >= max2>>16) max2 = cns[i], max2_i = i;
            }
            if ((double)(max&0xffff) / ((max&0xffff) + (max>>16)) >= 0.7) max_i = -1;
            if ((double)(max2&0xffff) / ((max2&0xffff) + (max2>>16)) >= 0.7) max2_i = -1;
            if (max_i >= 0) r[max_i] = 15;
            if (max2_i >= 0) r[max2_i] = 15;
            //for (i = 0; i < right - left; ++i) fputc("=ACMGRSVTWYHKDBN"[(int)r[i]], stderr); fputc('\n', stderr);
        }
        free(ref0); free(cns);
    }
    { // the length of the homopolymer run around the current position
        int c = seq_nt16_table[(int)ref[pos + 1]];
        if (c == 15) l_run = 1;
        else {
            for (i = pos + 2; ref[i]; ++i)
                if (seq_nt16_table[(int)ref[i]] != c) break;
            l_run = i;
            for (i = pos; i >= 0; --i)
                if (seq_nt16_table[(int)ref[i]] != c) break;
            l_run -= i + 1;
        }
    }
    // construct the consensus sequence
    max_ins = types[n_types - 1];   // max_ins is at least 0
    if (max_ins > 0) {
        int *inscns_aux = calloc(5 * n_types * max_ins, sizeof(int));
        // count the number of occurrences of each base at each position for each type of insertion
        for (t = 0; t < n_types; ++t) {
            if (types[t] > 0) {
                for (s = 0; s < n; ++s) {
                    for (i = 0; i < n_plp[s]; ++i) {
                        bam_pileup1_t *p = plp[s] + i;
                        if (p->indel == types[t]) {
                            uint8_t *seq = bam_get_seq(p->b);
                            for (k = 1; k <= p->indel; ++k) {
                                int c = seq_nt16_int[bam_seqi(seq, p->qpos + k)];
                                assert(c<5);
                                ++inscns_aux[(t*max_ins+(k-1))*5 + c];
                            }
                        }
                    }
                }
            }
        }
        // use the majority rule to construct the consensus
        inscns = calloc(n_types * max_ins, 1);
        for (t = 0; t < n_types; ++t) {
            for (j = 0; j < types[t]; ++j) {
                int max = 0, max_k = -1, *ia = &inscns_aux[(t*max_ins+j)*5];
                for (k = 0; k < 5; ++k)
                    if (ia[k] > max)
                        max = ia[k], max_k = k;
                inscns[t*max_ins + j] = max? max_k : 4;
                if ( max_k==4 ) { types[t] = 0; break; } // discard insertions which contain N's
            }
        }
        free(inscns_aux);
    }
    // compute the likelihood given each type of indel for each read
    max_ref2 = right - left + 2 + 2 * (max_ins > -types[0]? max_ins : -types[0]);
    ref2  = calloc(max_ref2, 1);
    query = calloc(right - left + max_rd_len + max_ins + 2, 1);
    score1 = calloc(N * n_types, sizeof(int));
    score2 = calloc(N * n_types, sizeof(int));
    bca->indelreg = 0;
    for (t = 0; t < n_types; ++t) {
        int l, ir;
        kpa_par_t apf1 = { 1e-4, 1e-2, 10 }, apf2 = { 1e-6, 1e-3, 10 };
        apf1.bw = apf2.bw = abs(types[t]) + 3;
        // compute indelreg
        if (types[t] == 0) ir = 0;
        else if (types[t] > 0) ir = est_indelreg(pos, ref, types[t], &inscns[t*max_ins]);
        else ir = est_indelreg(pos, ref, -types[t], 0);
        if (ir > bca->indelreg) bca->indelreg = ir;
//      fprintf(stderr, "%d, %d, %d\n", pos, types[t], ir);
        // realignment
        for (s = K = 0; s < n; ++s) {
            // write ref2
            for (k = 0, j = left; j <= pos; ++j)
                ref2[k++] = seq_nt16_int[(int)ref_sample[s][j-left]];
            if (types[t] <= 0) j += -types[t];
            else for (l = 0; l < types[t]; ++l)
                     ref2[k++] = inscns[t*max_ins + l];
            for (; j < right && ref[j]; ++j)
                ref2[k++] = seq_nt16_int[(int)ref_sample[s][j-left]];
            for (; k < max_ref2; ++k) ref2[k] = 4;
            if (j < right) right = j;
            // align each read to ref2
            for (i = 0; i < n_plp[s]; ++i, ++K) {
                bam_pileup1_t *p = plp[s] + i;
                int qbeg, qend, tbeg, tend, sc, kk;
                uint8_t *seq = bam_get_seq(p->b);
                uint32_t *cigar = bam_get_cigar(p->b);
                if (p->b->core.flag&4) continue; // unmapped reads
                // FIXME: the following loop should be better moved outside; nonetheless, realignment should be much slower anyway.
                for (kk = 0; kk < p->b->core.n_cigar; ++kk)
                    if ((cigar[kk]&BAM_CIGAR_MASK) == BAM_CREF_SKIP) break;
                if (kk < p->b->core.n_cigar) continue;
                // FIXME: the following skips soft clips, but using them may be more sensitive.
                // determine the start and end of sequences for alignment
                qbeg = tpos2qpos(&p->b->core, bam_get_cigar(p->b), left,  0, &tbeg);
                qend = tpos2qpos(&p->b->core, bam_get_cigar(p->b), right, 1, &tend);
                if (types[t] < 0) {
                    int l = -types[t];
                    tbeg = tbeg - l > left?  tbeg - l : left;
                }
                // write the query sequence
                for (l = qbeg; l < qend; ++l)
                    query[l - qbeg] = seq_nt16_int[bam_seqi(seq, l)];
                { // do realignment; this is the bottleneck
                    const uint8_t *qual = bam_get_qual(p->b), *bq;
                    uint8_t *qq;
                    qq = calloc(qend - qbeg, 1);
                    bq = (uint8_t*)bam_aux_get(p->b, "ZQ");
                    if (bq) ++bq; // skip type
                    for (l = qbeg; l < qend; ++l) {
                        qq[l - qbeg] = bq? qual[l] + (bq[l] - 64) : qual[l];
                        if (qq[l - qbeg] > 30) qq[l - qbeg] = 30;
                        if (qq[l - qbeg] < 7) qq[l - qbeg] = 7;
                    }
                    sc = kpa_glocal((uint8_t*)ref2 + tbeg - left, tend - tbeg + abs(types[t]),
                                    (uint8_t*)query, qend - qbeg, qq, &apf1, 0, 0);
                    l = (int)(100. * sc / (qend - qbeg) + .499); // used for adjusting indelQ below
                    if (l > 255) l = 255;
                    score1[K*n_types + t] = score2[K*n_types + t] = sc<<8 | l;
                    if (sc > 5) {
                        sc = kpa_glocal((uint8_t*)ref2 + tbeg - left, tend - tbeg + abs(types[t]),
                                        (uint8_t*)query, qend - qbeg, qq, &apf2, 0, 0);
                        l = (int)(100. * sc / (qend - qbeg) + .499);
                        if (l > 255) l = 255;
                        score2[K*n_types + t] = sc<<8 | l;
                    }
                    free(qq);
                }
/*
                for (l = 0; l < tend - tbeg + abs(types[t]); ++l)
                    fputc("ACGTN"[(int)ref2[tbeg-left+l]], stderr);
                fputc('\n', stderr);
                for (l = 0; l < qend - qbeg; ++l) fputc("ACGTN"[(int)query[l]], stderr);
                fputc('\n', stderr);
                fprintf(stderr, "pos=%d type=%d read=%d:%d name=%s qbeg=%d tbeg=%d score=%d\n", pos, types[t], s, i, bam1_qname(p->b), qbeg, tbeg, sc);
*/
            }
        }
    }
    free(ref2); free(query);
    { // compute indelQ
        int *sc, tmp, *sumq;
        sc   = alloca(n_types * sizeof(int));
        sumq = alloca(n_types * sizeof(int));
        memset(sumq, 0, sizeof(int) * n_types);
        for (s = K = 0; s < n; ++s) {
            for (i = 0; i < n_plp[s]; ++i, ++K) {
                bam_pileup1_t *p = plp[s] + i;
                int *sct = &score1[K*n_types], indelQ1, indelQ2, seqQ, indelQ;
                for (t = 0; t < n_types; ++t) sc[t] = sct[t]<<6 | t;
                for (t = 1; t < n_types; ++t) // insertion sort
                    for (j = t; j > 0 && sc[j] < sc[j-1]; --j)
                        tmp = sc[j], sc[j] = sc[j-1], sc[j-1] = tmp;
                /* errmod_cal() assumes that if the call is wrong, the
                 * likelihoods of other events are equal. This is about
                 * right for substitutions, but is not desired for
                 * indels. To reuse errmod_cal(), I have to make
                 * compromise for multi-allelic indels.
                 */
                if ((sc[0]&0x3f) == ref_type) {
                    indelQ1 = (sc[1]>>14) - (sc[0]>>14);
                    seqQ = est_seqQ(bca, types[sc[1]&0x3f], l_run);
                } else {
                    for (t = 0; t < n_types; ++t) // look for the reference type
                        if ((sc[t]&0x3f) == ref_type) break;
                    indelQ1 = (sc[t]>>14) - (sc[0]>>14);
                    seqQ = est_seqQ(bca, types[sc[0]&0x3f], l_run);
                }
                tmp = sc[0]>>6 & 0xff;
                indelQ1 = tmp > 111? 0 : (int)((1. - tmp/111.) * indelQ1 + .499); // reduce indelQ
                sct = &score2[K*n_types];
                for (t = 0; t < n_types; ++t) sc[t] = sct[t]<<6 | t;
                for (t = 1; t < n_types; ++t) // insertion sort
                    for (j = t; j > 0 && sc[j] < sc[j-1]; --j)
                        tmp = sc[j], sc[j] = sc[j-1], sc[j-1] = tmp;
                if ((sc[0]&0x3f) == ref_type) {
                    indelQ2 = (sc[1]>>14) - (sc[0]>>14);
                } else {
                    for (t = 0; t < n_types; ++t) // look for the reference type
                        if ((sc[t]&0x3f) == ref_type) break;
                    indelQ2 = (sc[t]>>14) - (sc[0]>>14);
                }
                tmp = sc[0]>>6 & 0xff;
                indelQ2 = tmp > 111? 0 : (int)((1. - tmp/111.) * indelQ2 + .499);
                // pick the smaller between indelQ1 and indelQ2
                indelQ = indelQ1 < indelQ2? indelQ1 : indelQ2;
                if (indelQ > 255) indelQ = 255;
                if (seqQ > 255) seqQ = 255;
                p->aux = (sc[0]&0x3f)<<16 | seqQ<<8 | indelQ; // use 22 bits in total
                sumq[sc[0]&0x3f] += indelQ < seqQ? indelQ : seqQ;
//              fprintf(stderr, "pos=%d read=%d:%d name=%s call=%d indelQ=%d seqQ=%d\n", pos, s, i, bam1_qname(p->b), types[sc[0]&0x3f], indelQ, seqQ);
            }
        }
        // determine bca->indel_types[] and bca->inscns
        bca->maxins = max_ins;
        bca->inscns = realloc(bca->inscns, bca->maxins * 4);
        for (t = 0; t < n_types; ++t)
            sumq[t] = sumq[t]<<6 | t;
        for (t = 1; t < n_types; ++t) // insertion sort
            for (j = t; j > 0 && sumq[j] > sumq[j-1]; --j)
                tmp = sumq[j], sumq[j] = sumq[j-1], sumq[j-1] = tmp;
        for (t = 0; t < n_types; ++t) // look for the reference type
            if ((sumq[t]&0x3f) == ref_type) break;
        if (t) { // then move the reference type to the first
            tmp = sumq[t];
            for (; t > 0; --t) sumq[t] = sumq[t-1];
            sumq[0] = tmp;
        }
        for (t = 0; t < 4; ++t) bca->indel_types[t] = B2B_INDEL_NULL;
        for (t = 0; t < 4 && t < n_types; ++t) {
            bca->indel_types[t] = types[sumq[t]&0x3f];
            memcpy(&bca->inscns[t * bca->maxins], &inscns[(sumq[t]&0x3f) * max_ins], bca->maxins);
        }
        // update p->aux
        for (s = n_alt = 0; s < n; ++s) {
            for (i = 0; i < n_plp[s]; ++i) {
                bam_pileup1_t *p = plp[s] + i;
                int x = types[p->aux>>16&0x3f];
                for (j = 0; j < 4; ++j)
                    if (x == bca->indel_types[j]) break;
                p->aux = j<<16 | (j == 4? 0 : (p->aux&0xffff));
                if ((p->aux>>16&0x3f) > 0) ++n_alt;
                //fprintf(stderr, "X pos=%d read=%d:%d name=%s call=%d type=%d seqQ=%d indelQ=%d\n", pos, s, i, bam1_qname(p->b), (p->aux>>16)&0x3f, bca->indel_types[(p->aux>>16)&0x3f], (p->aux>>8)&0xff, p->aux&0xff);
            }
        }
    }
    free(score1); free(score2);
    // free
    for (i = 0; i < n; ++i) free(ref_sample[i]);
    free(ref_sample);
    free(types); free(inscns);
    return n_alt > 0? 0 : -1;
}
