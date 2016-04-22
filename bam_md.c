/*  bam_md.c -- calmd subcommand.

    Copyright (C) 2009-2011, 2014-2015 Genome Research Ltd.
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

#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "kprobaln.h"
#include "sam_opts.h"
#include "samtools.h"

#define USE_EQUAL 1
#define DROP_TAG  2
#define BIN_QUAL  4
#define UPDATE_NM 8
#define UPDATE_MD 16
#define HASH_QNM  32

int bam_aux_drop_other(bam1_t *b, uint8_t *s);

void bam_fillmd1_core(bam1_t *b, char *ref, int ref_len, int flag, int max_nm)
{
    uint8_t *seq = bam_get_seq(b);
    uint32_t *cigar = bam_get_cigar(b);
    bam1_core_t *c = &b->core;
    int i, x, y, u = 0;
    kstring_t *str;
    int32_t old_nm_i = -1, nm = 0;

    str = (kstring_t*)calloc(1, sizeof(kstring_t));
    for (i = y = 0, x = c->pos; i < c->n_cigar; ++i) {
        int j, l = cigar[i]>>4, op = cigar[i]&0xf;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (j = 0; j < l; ++j) {
                int c1, c2, z = y + j;
                if (x+j >= ref_len || ref[x+j] == '\0') break; // out of bounds
                c1 = bam_seqi(seq, z), c2 = seq_nt16_table[(int)ref[x+j]];
                if ((c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0) { // a match
                    if (flag&USE_EQUAL) seq[z/2] &= (z&1)? 0xf0 : 0x0f;
                    ++u;
                } else {
                    kputw(u, str); kputc(ref[x+j], str);
                    u = 0; ++nm;
                }
            }
            if (j < l) break;
            x += l; y += l;
        } else if (op == BAM_CDEL) {
            kputw(u, str); kputc('^', str);
            for (j = 0; j < l; ++j) {
                if (x+j >= ref_len || ref[x+j] == '\0') break;
                kputc(ref[x+j], str);
            }
            u = 0;
            x += j; nm += j;
            if (j < l) break;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            y += l;
            if (op == BAM_CINS) nm += l;
        } else if (op == BAM_CREF_SKIP) {
            x += l;
        }
    }
    kputw(u, str);
    // apply max_nm
    if (max_nm > 0 && nm >= max_nm) {
        for (i = y = 0, x = c->pos; i < c->n_cigar; ++i) {
            int j, l = cigar[i]>>4, op = cigar[i]&0xf;
            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                for (j = 0; j < l; ++j) {
                    int c1, c2, z = y + j;
                    if (x+j >= ref_len || ref[x+j] == '\0') break; // out of bounds
                    c1 = bam_seqi(seq, z), c2 = seq_nt16_table[(int)ref[x+j]];
                    if ((c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0) { // a match
                        seq[z/2] |= (z&1)? 0x0f : 0xf0;
                        bam_get_qual(b)[z] = 0;
                    }
                }
                if (j < l) break;
                x += l; y += l;
            } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) x += l;
            else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) y += l;
        }
    }
    // update NM
    if ((flag & UPDATE_NM) && !(c->flag & BAM_FUNMAP)) {
        uint8_t *old_nm = bam_aux_get(b, "NM");
        if (old_nm) old_nm_i = bam_aux2i(old_nm);
        if (!old_nm) bam_aux_append(b, "NM", 'i', 4, (uint8_t*)&nm);
        else if (nm != old_nm_i) {
            fprintf(stderr, "[bam_fillmd1] different NM for read '%s': %d -> %d\n", bam_get_qname(b), old_nm_i, nm);
            bam_aux_del(b, old_nm);
            bam_aux_append(b, "NM", 'i', 4, (uint8_t*)&nm);
        }
    }
    // update MD
    if ((flag & UPDATE_MD) && !(c->flag & BAM_FUNMAP)) {
        uint8_t *old_md = bam_aux_get(b, "MD");
        if (!old_md) bam_aux_append(b, "MD", 'Z', str->l + 1, (uint8_t*)str->s);
        else {
            int is_diff = 0;
            if (strlen((char*)old_md+1) == str->l) {
                for (i = 0; i < str->l; ++i)
                    if (toupper(old_md[i+1]) != toupper(str->s[i]))
                        break;
                if (i < str->l) is_diff = 1;
            } else is_diff = 1;
            if (is_diff) {
                fprintf(stderr, "[bam_fillmd1] different MD for read '%s': '%s' -> '%s'\n", bam_get_qname(b), old_md+1, str->s);
                bam_aux_del(b, old_md);
                bam_aux_append(b, "MD", 'Z', str->l + 1, (uint8_t*)str->s);
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

    free(str->s); free(str);
}

void bam_fillmd1(bam1_t *b, char *ref, int flag)
{
    bam_fillmd1_core(b, ref, INT_MAX, flag, 0);
}

int bam_cap_mapQ(bam1_t *b, char *ref, int ref_len, int thres)
{
    uint8_t *seq = bam_get_seq(b), *qual = bam_get_qual(b);
    uint32_t *cigar = bam_get_cigar(b);
    bam1_core_t *c = &b->core;
    int i, x, y, mm, q, len, clip_l, clip_q;
    double t;
    if (thres < 0) thres = 40; // set the default
    mm = q = len = clip_l = clip_q = 0;
    for (i = y = 0, x = c->pos; i < c->n_cigar; ++i) {
        int j, l = cigar[i]>>4, op = cigar[i]&0xf;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (j = 0; j < l; ++j) {
                int c1, c2, z = y + j;
                if (x+j >= ref_len || ref[x+j] == '\0') break; // out of bounds
                c1 = bam_seqi(seq, z), c2 = seq_nt16_table[(int)ref[x+j]];
                if (c2 != 15 && c1 != 15 && qual[z] >= 13) { // not ambiguous
                    ++len;
                    if (c1 && c1 != c2 && qual[z] >= 13) { // mismatch
                        ++mm;
                        q += qual[z] > 33? 33 : qual[z];
                    }
                }
            }
            if (j < l) break;
            x += l; y += l; len += l;
        } else if (op == BAM_CDEL) {
            for (j = 0; j < l; ++j)
                if (x+j >= ref_len || ref[x+j] == '\0') break;
            if (j < l) break;
            x += l;
        } else if (op == BAM_CSOFT_CLIP) {
            for (j = 0; j < l; ++j) clip_q += qual[y+j];
            clip_l += l;
            y += l;
        } else if (op == BAM_CHARD_CLIP) {
            clip_q += 13 * l;
            clip_l += l;
        } else if (op == BAM_CINS) y += l;
        else if (op == BAM_CREF_SKIP) x += l;
    }
    for (i = 0, t = 1; i < mm; ++i)
        t *= (double)len / (i+1);
    t = q - 4.343 * log(t) + clip_q / 5.;
    if (t > thres) return -1;
    if (t < 0) t = 0;
    t = sqrt((thres - t) / thres) * thres;
//  fprintf(stderr, "%s %lf %d\n", bam_get_qname(b), t, q);
    return (int)(t + .499);
}

int bam_prob_realn_core(bam1_t *b, const char *ref, int ref_len, int flag)
{
    int k, i, bw, x, y, yb, ye, xb, xe, apply_baq = flag&1, extend_baq = flag>>1&1, redo_baq = flag&4;
    uint32_t *cigar = bam_get_cigar(b);
    bam1_core_t *c = &b->core;
    kpa_par_t conf = kpa_par_def;
    uint8_t *bq = 0, *zq = 0, *qual = bam_get_qual(b);
    if ((c->flag & BAM_FUNMAP) || b->core.l_qseq == 0 || qual[0] == (uint8_t)-1)
        return -1; // do nothing

    // test if BQ or ZQ is present
    if ((bq = bam_aux_get(b, "BQ")) != 0) ++bq;
    if ((zq = bam_aux_get(b, "ZQ")) != 0 && *zq == 'Z') ++zq;
    if (bq && redo_baq)
    {
        bam_aux_del(b, bq-1);
        bq = 0;
    }
    if (bq && zq) { // remove the ZQ tag
        bam_aux_del(b, zq-1);
        zq = 0;
    }
    if (bq || zq) {
        if ((apply_baq && zq) || (!apply_baq && bq)) return -3; // in both cases, do nothing
        if (bq && apply_baq) { // then convert BQ to ZQ
            for (i = 0; i < c->l_qseq; ++i)
                qual[i] = qual[i] + 64 < bq[i]? 0 : qual[i] - ((int)bq[i] - 64);
            *(bq - 3) = 'Z';
        } else if (zq && !apply_baq) { // then convert ZQ to BQ
            for (i = 0; i < c->l_qseq; ++i)
                qual[i] += (int)zq[i] - 64;
            *(zq - 3) = 'B';
        }
        return 0;
    }
    // find the start and end of the alignment
    x = c->pos, y = 0, yb = ye = xb = xe = -1;
    for (k = 0; k < c->n_cigar; ++k) {
        int op, l;
        op = cigar[k]&0xf; l = cigar[k]>>4;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (yb < 0) yb = y;
            if (xb < 0) xb = x;
            ye = y + l; xe = x + l;
            x += l; y += l;
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) y += l;
        else if (op == BAM_CDEL) x += l;
        else if (op == BAM_CREF_SKIP) return -1; // do nothing if there is a reference skip
    }
    // set bandwidth and the start and the end
    bw = 7;
    if (abs((xe - xb) - (ye - yb)) > bw)
        bw = abs((xe - xb) - (ye - yb)) + 3;
    conf.bw = bw;
    xb -= yb + bw/2; if (xb < 0) xb = 0;
    xe += c->l_qseq - ye + bw/2;
    if (xe - xb - c->l_qseq > bw)
        xb += (xe - xb - c->l_qseq - bw) / 2, xe -= (xe - xb - c->l_qseq - bw) / 2;
    { // glocal
        uint8_t *s, *r, *q, *seq = bam_get_seq(b), *bq;
        int *state;
        bq = calloc(c->l_qseq + 1, 1);
        memcpy(bq, qual, c->l_qseq);
        s = calloc(c->l_qseq, 1);
        for (i = 0; i < c->l_qseq; ++i) s[i] = seq_nt16_int[bam_seqi(seq, i)];
        r = calloc(xe - xb, 1);
        for (i = xb; i < xe; ++i) {
            if (i >= ref_len || ref[i] == '\0') { xe = i; break; }
            r[i-xb] = seq_nt16_int[seq_nt16_table[(int)ref[i]]];
        }
        state = calloc(c->l_qseq, sizeof(int));
        q = calloc(c->l_qseq, 1);
        kpa_glocal(r, xe-xb, s, c->l_qseq, qual, &conf, state, q);
        if (!extend_baq) { // in this block, bq[] is capped by base quality qual[]
            for (k = 0, x = c->pos, y = 0; k < c->n_cigar; ++k) {
                int op = cigar[k]&0xf, l = cigar[k]>>4;
                if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                    for (i = y; i < y + l; ++i) {
                        if ((state[i]&3) != 0 || state[i]>>2 != x - xb + (i - y)) bq[i] = 0;
                        else bq[i] = bq[i] < q[i]? bq[i] : q[i];
                    }
                    x += l; y += l;
                } else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) y += l;
                else if (op == BAM_CDEL) x += l;
            }
            for (i = 0; i < c->l_qseq; ++i) bq[i] = qual[i] - bq[i] + 64; // finalize BQ
        } else { // in this block, bq[] is BAQ that can be larger than qual[] (different from the above!)
            uint8_t *left, *rght;
            left = calloc(c->l_qseq, 1); rght = calloc(c->l_qseq, 1);
            for (k = 0, x = c->pos, y = 0; k < c->n_cigar; ++k) {
                int op = cigar[k]&0xf, l = cigar[k]>>4;
                if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                    for (i = y; i < y + l; ++i)
                        bq[i] = ((state[i]&3) != 0 || state[i]>>2 != x - xb + (i - y))? 0 : q[i];
                    for (left[y] = bq[y], i = y + 1; i < y + l; ++i)
                        left[i] = bq[i] > left[i-1]? bq[i] : left[i-1];
                    for (rght[y+l-1] = bq[y+l-1], i = y + l - 2; i >= y; --i)
                        rght[i] = bq[i] > rght[i+1]? bq[i] : rght[i+1];
                    for (i = y; i < y + l; ++i)
                        bq[i] = left[i] < rght[i]? left[i] : rght[i];
                    x += l; y += l;
                } else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) y += l;
                else if (op == BAM_CDEL) x += l;
            }
            for (i = 0; i < c->l_qseq; ++i) bq[i] = 64 + (qual[i] <= bq[i]? 0 : qual[i] - bq[i]); // finalize BQ
            free(left); free(rght);
        }
        if (apply_baq) {
            for (i = 0; i < c->l_qseq; ++i) qual[i] -= bq[i] - 64; // modify qual
            bam_aux_append(b, "ZQ", 'Z', c->l_qseq + 1, bq);
        } else bam_aux_append(b, "BQ", 'Z', c->l_qseq + 1, bq);
        free(bq); free(s); free(r); free(q); free(state);
    }
    return 0;
}

int bam_prob_realn(bam1_t *b, const char *ref)
{
    return bam_prob_realn_core(b, ref, INT_MAX, 1);
}

int calmd_usage() {
    fprintf(stderr,
"Usage: samtools calmd [-eubrAES] <aln.bam> <ref.fasta>\n"
"Options:\n"
"  -e       change identical bases to '='\n"
"  -u       uncompressed BAM output (for piping)\n"
"  -b       compressed BAM output\n"
"  -S       ignored (input format is auto-detected)\n"
"  -A       modify the quality string\n"
"  -r       compute the BQ tag (without -A) or cap baseQ by BAQ (with -A)\n"
"  -E       extended BAQ for better sensitivity but lower specificity\n");

    sam_global_opt_help(stderr, "-....");
    return 1;
}

int bam_fillmd(int argc, char *argv[])
{
    int c, flt_flag, tid = -2, ret, len, is_bam_out, is_uncompressed, max_nm, is_realn, capQ, baq_flag;
    samFile *fp = NULL, *fpout = NULL;
    bam_hdr_t *header = NULL;
    faidx_t *fai = NULL;
    char *ref = NULL, mode_w[8], *ref_file;
    bam1_t *b = NULL;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 0, 0, 0),
        { NULL, 0, NULL, 0 }
    };

    flt_flag = UPDATE_NM | UPDATE_MD;
    is_bam_out = is_uncompressed = is_realn = max_nm = capQ = baq_flag = 0;
    strcpy(mode_w, "w");
    while ((c = getopt_long(argc, argv, "EqreuNhbSC:n:Ad", lopts, NULL)) >= 0) {
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

    header = sam_hdr_read(fp);
    if (header == NULL || header->n_targets == 0) {
        fprintf(stderr, "[bam_fillmd] input SAM does not have header. Abort!\n");
        goto fail;
    }

    fpout = sam_open_format("-", mode_w, &ga.out);
    if (fpout == NULL) {
        print_error_errno("calmd", "Failed to open output");
        goto fail;
    }
    if (sam_hdr_write(fpout, header) < 0) {
        print_error_errno("calmd", "Failed to write sam header");
        goto fail;
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
            if (tid != b->core.tid) {
                free(ref);
                ref = fai_fetch(fai, header->target_name[b->core.tid], &len);
                tid = b->core.tid;
                if (ref == 0) { // FIXME: Should this always be fatal?
                    fprintf(stderr, "[bam_fillmd] fail to find sequence '%s' in the reference.\n",
                            header->target_name[tid]);
                    if (is_realn || capQ > 10) goto fail; // Would otherwise crash
                }
            }
            if (is_realn) bam_prob_realn_core(b, ref, len, baq_flag);
            if (capQ > 10) {
                int q = bam_cap_mapQ(b, ref, len, capQ);
                if (b->core.qual > q) b->core.qual = q;
            }
            if (ref) bam_fillmd1_core(b, ref, len, flt_flag, max_nm);
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
    bam_destroy1(b);
    bam_hdr_destroy(header);

    free(ref);
    fai_destroy(fai);
    sam_close(fp);
    if (sam_close(fpout) < 0) {
        fprintf(stderr, "[bam_fillmd] error when closing output file\n");
        return 1;
    }
    return 0;

 fail:
    free(ref);
    if (b) bam_destroy1(b);
    if (header) bam_hdr_destroy(header);
    if (fai) fai_destroy(fai);
    if (fp) sam_close(fp);
    if (fpout) sam_close(fpout);
    return 1;
}
