/*  padding.c -- depad subcommand.

    Copyright (C) 2011, 2012 Broad Institute.
    Copyright (C) 2014 Genome Research Ltd.
    Portions copyright (C) 2012, 2013 Peter Cock, The James Hutton Institute.

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

#include <string.h>
#include <assert.h>
#include <unistd.h>
#include "htslib/kstring.h"
#include "sam_header.h"
#include "sam.h"
#include "bam.h"
#include "htslib/faidx.h"

bam_header_t *bam_header_dup(const bam_header_t *h0); /*in sam.c*/

static void replace_cigar(bam1_t *b, int n, uint32_t *cigar)
{
    if (n != b->core.n_cigar) {
        int o = b->core.l_qname + b->core.n_cigar * 4;
        if (b->data_len + (n - b->core.n_cigar) * 4 > b->m_data) {
            b->m_data = b->data_len + (n - b->core.n_cigar) * 4;
            kroundup32(b->m_data);
            b->data = (uint8_t*)realloc(b->data, b->m_data);
        }
        memmove(b->data + b->core.l_qname + n * 4, b->data + o, b->data_len - o);
        memcpy(b->data + b->core.l_qname, cigar, n * 4);
        b->data_len += (n - b->core.n_cigar) * 4;
        b->core.n_cigar = n;
    } else memcpy(b->data + b->core.l_qname, cigar, n * 4);
}

#define write_cigar(_c, _n, _m, _v) do { \
        if (_n == _m) { \
            _m = _m? _m<<1 : 4; \
            _c = (uint32_t*)realloc(_c, _m * 4); \
        } \
        _c[_n++] = (_v); \
    } while (0)

static void unpad_seq(bam1_t *b, kstring_t *s)
{
    int k, j, i;
    int length;
    uint32_t *cigar = bam1_cigar(b);
    uint8_t *seq = bam1_seq(b);
    // b->core.l_qseq gives length of the SEQ entry (including soft clips, S)
    // We need the padded length after alignment from the CIGAR (excluding
    // soft clips S, but including pads from CIGAR D operations)
    length = 0;
    for (k = 0; k < b->core.n_cigar; ++k) {
        int op, ol;
        op= bam_cigar_op(cigar[k]);
        ol = bam_cigar_oplen(cigar[k]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL)
            length += ol;
    }
    ks_resize(s, length);
    for (k = 0, s->l = 0, j = 0; k < b->core.n_cigar; ++k) {
        int op, ol;
        op = bam_cigar_op(cigar[k]);
        ol = bam_cigar_oplen(cigar[k]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (i = 0; i < ol; ++i, ++j) s->s[s->l++] = bam1_seqi(seq, j);
        } else if (op == BAM_CSOFT_CLIP) {
            j += ol;
        } else if (op == BAM_CHARD_CLIP) {
            /* do nothing */
        } else if (op == BAM_CDEL) {
            for (i = 0; i < ol; ++i) s->s[s->l++] = 0;
        } else {
            fprintf(stderr, "[depad] ERROR: Didn't expect CIGAR op %c in read %s\n", BAM_CIGAR_STR[op], bam1_qname(b));
            exit(1);
        }
    }
    assert(length == s->l);
}

int load_unpadded_ref(faidx_t *fai, char *ref_name, int ref_len, kstring_t *seq)
{
    char base;
    char *fai_ref = 0;
    int fai_ref_len = 0, k;

    fai_ref = fai_fetch(fai, ref_name, &fai_ref_len);
    if (fai_ref_len != ref_len) {
        fprintf(stderr, "[depad] ERROR: FASTA sequence %s length %i, expected %i\n", ref_name, fai_ref_len, ref_len);
        free(fai_ref);
        return -1;
    }
    ks_resize(seq, ref_len);
    seq->l = 0;
    for (k = 0; k < ref_len; ++k) {
        base = fai_ref[k];
        if (base == '-' || base == '*') {
            // Map gaps to null to match unpad_seq function
            seq->s[seq->l++] = 0;
        } else {
            int i = bam_nt16_table[(int)base];
            if (i == 0 || i==16) { // Equals maps to 0, anything unexpected to 16
                fprintf(stderr, "[depad] ERROR: Invalid character %c (ASCII %i) in FASTA sequence %s\n", base, (int)base, ref_name);
                free(fai_ref);
                return -1;
            }
            seq->s[seq->l++] = i;
        }
    }
    assert(ref_len == seq->l);
    free(fai_ref);
    return 0;
}

int get_unpadded_len(faidx_t *fai, char *ref_name, int padded_len)
{
    char base;
    char *fai_ref = 0;
    int fai_ref_len = 0, k;
    int bases=0, gaps=0;

    fai_ref = fai_fetch(fai, ref_name, &fai_ref_len);
    if (fai_ref_len != padded_len) {
        fprintf(stderr, "[depad] ERROR: FASTA sequence '%s' length %i, expected %i\n", ref_name, fai_ref_len, padded_len);
        free(fai_ref);
        return -1;
    }
    for (k = 0; k < padded_len; ++k) {
        //fprintf(stderr, "[depad] checking base %i of %i or %i\n", k+1, ref_len, strlen(fai_ref));
        base = fai_ref[k];
        if (base == '-' || base == '*') {
            gaps += 1;
        } else {
            int i = bam_nt16_table[(int)base];
            if (i == 0 || i==16) { // Equals maps to 0, anything unexpected to 16
                fprintf(stderr, "[depad] ERROR: Invalid character %c (ASCII %i) in FASTA sequence '%s'\n", base, (int)base, ref_name);
                free(fai_ref);
                return -1;
            }
            bases += 1;
        }
    }
    free(fai_ref);
    assert (padded_len == bases + gaps);
    return bases;
}

static inline int * update_posmap(int *posmap, kstring_t ref)
{
    int i, k;
    posmap = realloc(posmap, ref.m * sizeof(int));
    for (i = k = 0; i < ref.l; ++i) {
        posmap[i] = k;
        if (ref.s[i]) ++k;
    }
    return posmap;
}

int bam_pad2unpad(samfile_t *in, samfile_t *out, faidx_t *fai)
{
    bam_header_t *h = 0;
    bam1_t *b = 0;
    kstring_t r, q;
    int r_tid = -1;
    uint32_t *cigar2 = 0;
    int ret = 0, n2 = 0, m2 = 0, *posmap = 0;

    b = bam_init1();
    r.l = r.m = q.l = q.m = 0; r.s = q.s = 0;
    int read_ret;
    h = in->header;
    while ((read_ret = samread(in, b)) >= 0) { // read one alignment from `in'
        uint32_t *cigar = bam1_cigar(b);
        n2 = 0;
        if (b->core.pos == 0 && b->core.tid >= 0 && strcmp(bam1_qname(b), h->target_name[b->core.tid]) == 0) {
            // fprintf(stderr, "[depad] Found embedded reference '%s'\n", bam1_qname(b));
            r_tid = b->core.tid;
            unpad_seq(b, &r);
            if (h->target_len[r_tid] != r.l) {
                fprintf(stderr, "[depad] ERROR: (Padded) length of '%s' is %u in BAM header, but %llu in embedded reference\n", bam1_qname(b), h->target_len[r_tid], (unsigned long long)(r.l));
                return -1;
            }
            if (fai) {
                // Check the embedded reference matches the FASTA file
                if (load_unpadded_ref(fai, h->target_name[b->core.tid], h->target_len[b->core.tid], &q)) {
                    fprintf(stderr, "[depad] ERROR: Failed to load embedded reference '%s' from FASTA\n", h->target_name[b->core.tid]);
                    return -1;
                }
                assert(r.l == q.l);
                int i;
                for (i = 0; i < r.l; ++i) {
                    if (r.s[i] != q.s[i]) {
                        // Show gaps as ASCII 45
                        fprintf(stderr, "[depad] ERROR: Embedded sequence and reference FASTA don't match for %s base %i, '%c' vs '%c'\n",
                            h->target_name[b->core.tid], i+1,
                            r.s[i] ? bam_nt16_rev_table[(int)r.s[i]] : 45,
                            q.s[i] ? bam_nt16_rev_table[(int)q.s[i]] : 45);
                        return -1;
                    }
                }
            }
            write_cigar(cigar2, n2, m2, bam_cigar_gen(b->core.l_qseq, BAM_CMATCH));
            replace_cigar(b, n2, cigar2);
            posmap = update_posmap(posmap, r);
        } else if (b->core.n_cigar > 0) {
            int i, k, op;
            if (b->core.tid < 0) {
                fprintf(stderr, "[depad] ERROR: Read '%s' has CIGAR but no RNAME\n", bam1_qname(b));
                return -1;
            } else if (b->core.tid == r_tid) {
                ; // good case, reference available
                //fprintf(stderr, "[depad] Have ref '%s' for read '%s'\n", h->target_name[b->core.tid], bam1_qname(b));
            } else if (fai) {
                if (load_unpadded_ref(fai, h->target_name[b->core.tid], h->target_len[b->core.tid], &r)) {
                    fprintf(stderr, "[depad] ERROR: Failed to load '%s' from reference FASTA\n", h->target_name[b->core.tid]);
                    return -1;
                }
                posmap = update_posmap(posmap, r);
                r_tid = b->core.tid;
                // fprintf(stderr, "[depad] Loaded %s from FASTA file\n", h->target_name[b->core.tid]);
            } else {
                fprintf(stderr, "[depad] ERROR: Missing %s embedded reference sequence (and no FASTA file)\n", h->target_name[b->core.tid]);
                return -1;
            }
            unpad_seq(b, &q);
            if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) {
                write_cigar(cigar2, n2, m2, cigar[0]);
            } else if (bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP) {
                write_cigar(cigar2, n2, m2, cigar[0]);
                if (b->core.n_cigar > 2 && bam_cigar_op(cigar[1]) == BAM_CSOFT_CLIP) {
                    write_cigar(cigar2, n2, m2, cigar[1]);
                }
            }
            /* Determine CIGAR operator for each base in the aligned read */
            for (i = 0, k = b->core.pos; i < q.l; ++i, ++k)
                q.s[i] = q.s[i]? (r.s[k]? BAM_CMATCH : BAM_CINS) : (r.s[k]? BAM_CDEL : BAM_CPAD);
            /* Include any pads if starts with an insert */
            if (q.s[0] == BAM_CINS) {
                for (k = 0; k+1 < b->core.pos && !r.s[b->core.pos - k - 1]; ++k);
                if (k) write_cigar(cigar2, n2, m2, bam_cigar_gen(k, BAM_CPAD));
            }
            /* Count consecutive CIGAR operators to turn into a CIGAR string */
            for (i = k = 1, op = q.s[0]; i < q.l; ++i) {
                if (op != q.s[i]) {
                    write_cigar(cigar2, n2, m2, bam_cigar_gen(k, op));
                    op = q.s[i]; k = 1;
                } else ++k;
            }
            write_cigar(cigar2, n2, m2, bam_cigar_gen(k, op));
            if (bam_cigar_op(cigar[b->core.n_cigar-1]) == BAM_CSOFT_CLIP) {
                write_cigar(cigar2, n2, m2, cigar[b->core.n_cigar-1]);
            } else if (bam_cigar_op(cigar[b->core.n_cigar-1]) == BAM_CHARD_CLIP) {
                if (b->core.n_cigar > 2 && bam_cigar_op(cigar[b->core.n_cigar-2]) == BAM_CSOFT_CLIP) {
                    write_cigar(cigar2, n2, m2, cigar[b->core.n_cigar-2]);
                }
                write_cigar(cigar2, n2, m2, cigar[b->core.n_cigar-1]);
            }
            /* Remove redundant P operators between M/X/=/D operators, e.g. 5M2P10M -> 15M */
            int pre_op, post_op;
            for (i = 2; i < n2; ++i)
                if (bam_cigar_op(cigar2[i-1]) == BAM_CPAD) {
                    pre_op = bam_cigar_op(cigar2[i-2]);
                    post_op = bam_cigar_op(cigar2[i]);
                    /* Note don't need to check for X/= as code above will use M only */
                    if ((pre_op == BAM_CMATCH || pre_op == BAM_CDEL) && (post_op == BAM_CMATCH || post_op == BAM_CDEL)) {
                        /* This is a redundant P operator */
                        cigar2[i-1] = 0; // i.e. 0M
                        /* If had same operator either side, combine them in post_op */
                        if (pre_op == post_op) {
                            /* If CIGAR M, could treat as simple integers since BAM_CMATCH is zero*/
                            cigar2[i] = bam_cigar_gen(bam_cigar_oplen(cigar2[i-2]) + bam_cigar_oplen(cigar2[i]), post_op);
                            cigar2[i-2] = 0; // i.e. 0M
                        }
                    }
                }
            /* Remove the zero'd operators (0M) */
            for (i = k = 0; i < n2; ++i)
                if (cigar2[i]) cigar2[k++] = cigar2[i];
            n2 = k;
            replace_cigar(b, n2, cigar2);
        }
        /* Even unmapped reads can have a POS value, e.g. if their mate was mapped */
        if (b->core.pos != -1) b->core.pos = posmap[b->core.pos];
        if (b->core.mtid < 0 || b->core.mpos < 0) {
            /* Nice case, no mate to worry about*/
            // fprintf(stderr, "[depad] Read '%s' mate not mapped\n", bam1_qname(b));
            /* TODO - Warning if FLAG says mate should be mapped? */
            /* Clean up funny input where mate position is given but mate reference is missing: */
            b->core.mtid = -1;
            b->core.mpos = -1;
        } else if (b->core.mtid == b->core.tid) {
            /* Nice case, same reference */
            // fprintf(stderr, "[depad] Read '%s' mate mapped to same ref\n", bam1_qname(b));
            b->core.mpos = posmap[b->core.mpos];
        } else {
            /* Nasty case, Must load alternative posmap */
            // fprintf(stderr, "[depad] Loading reference '%s' temporarily\n", h->target_name[b->core.mtid]);
            if (!fai) {
                fprintf(stderr, "[depad] ERROR: Needed reference %s sequence for mate (and no FASTA file)\n", h->target_name[b->core.mtid]);
                return -1;
            }
            /* Temporarily load the other reference sequence */
            if (load_unpadded_ref(fai, h->target_name[b->core.mtid], h->target_len[b->core.mtid], &r)) {
                fprintf(stderr, "[depad] ERROR: Failed to load '%s' from reference FASTA\n", h->target_name[b->core.mtid]);
                return -1;
            }
            posmap = update_posmap(posmap, r);
            b->core.mpos = posmap[b->core.mpos];
            /* Restore the reference and posmap*/
            if (load_unpadded_ref(fai, h->target_name[b->core.tid], h->target_len[b->core.tid], &r)) {
                fprintf(stderr, "[depad] ERROR: Failed to load '%s' from reference FASTA\n", h->target_name[b->core.tid]);
                return -1;
            }
            posmap = update_posmap(posmap, r);
        }
        /* Most reads will have been moved so safest to always recalculate the BIN value */
        b->core.bin = bam_reg2bin(b->core.pos, bam_calend(&b->core, bam1_cigar(b)));
        samwrite(out, b);
    }
    if (read_ret < -1) {
        fprintf(stderr, "[depad] truncated file.\n");
        ret = 1;
    }
    free(r.s); free(q.s); free(posmap);
    bam_destroy1(b);
    return ret;
}

bam_header_t * fix_header(bam_header_t *old, faidx_t *fai)
{
#if 0
    int i = 0, unpadded_len = 0;
    bam_header_t *header = 0 ;

    header = bam_header_dup(old);
    for (i = 0; i < old->n_targets; ++i) {
        unpadded_len = get_unpadded_len(fai, old->target_name[i], old->target_len[i]);
        if (unpadded_len < 0) {
            fprintf(stderr, "[depad] ERROR getting unpadded length of '%s', padded length %i\n", old->target_name[i], old->target_len[i]);
        } else {
            header->target_len[i] = unpadded_len;
            //fprintf(stderr, "[depad] Recalculating '%s' length %i -> %i\n", old->target_name[i], old->target_len[i], header->target_len[i]);
        }
    }
    /* Duplicating the header allocated new buffer for header string */
    /* After modifying the @SQ lines it will only get smaller, since */
    /* the LN entries will be the same or shorter, and we'll remove */
    /* any MD entries (MD5 checksums). */
    assert(strlen(old->text) == strlen(header->text));
    assert (0==strcmp(old->text, header->text));
    const char *text;
    text = old->text;
    header->text[0] = '\0'; /* Resuse the allocated buffer */
    char * newtext = header->text;
    char * end=NULL;
    while (text[0]=='@') {
        end = strchr(text, '\n');
        assert(end != 0);
        if (text[1]=='S' && text[2]=='Q' && text[3]=='\t') {
            /* TODO - edit the @SQ line here to remove MD and fix LN. */
            /* For now just remove the @SQ line, and samtools will */
            /* automatically generate a minimal replacement with LN. */
            /* However, that discards any other tags like AS, SP, UR. */
            //fprintf(stderr, "[depad] Removing @SQ line\n");
        } else {
            /* Copy this line to the new header */
            strncat(newtext, text, end - text + 1);
        }
        text = end + 1;
    }
    assert (text[0]=='\0');
    /* Check we didn't overflow the buffer */
    assert (strlen(header->text) <= strlen(old->text));
    if (strlen(header->text) < header->l_text) {
        //fprintf(stderr, "[depad] Reallocating header buffer\n");
        assert (newtext == header->text);
        newtext = malloc(strlen(header->text) + 1);
        strcpy(newtext, header->text);
        free(header->text);
        header->text = newtext;
        header->l_text = strlen(newtext);
    }
    //fprintf(stderr, "[depad] Here is the new header (pending @SQ lines),\n\n%s\n(end)\n", header->text);
    return header;
#else
    fprintf(stderr, "Samtools-htslib: fix_header() header parsing not yet implemented\n");
    abort();
#endif
}

static int usage(int is_long_help);

int main_pad2unpad(int argc, char *argv[])
{
    samfile_t *in = 0, *out = 0;
    bam_header_t *h = 0;
    faidx_t *fai = 0;
    int c, is_bamin = 1, compress_level = -1, is_bamout = 1, is_long_help = 0;
    char in_mode[5], out_mode[5], *fn_out = 0, *fn_list = 0, *fn_ref = 0;
    int ret=0;

    /* parse command-line options */
    strcpy(in_mode, "r"); strcpy(out_mode, "w");
    while ((c = getopt(argc, argv, "Sso:u1T:?")) >= 0) {
        switch (c) {
        case 'S': is_bamin = 0; break;
        case 's': assert(compress_level == -1); is_bamout = 0; break;
        case 'o': fn_out = strdup(optarg); break;
        case 'u': assert(is_bamout == 1); compress_level = 0; break;
        case '1': assert(is_bamout == 1); compress_level = 1; break;
        case 'T': fn_ref = strdup(optarg); break;
        case '?': is_long_help = 1; break;
        default: return usage(is_long_help);
        }
    }
    if (argc == optind) return usage(is_long_help);

    if (is_bamin) strcat(in_mode, "b");
    if (is_bamout) strcat(out_mode, "b");
    strcat(out_mode, "h");
    if (compress_level >= 0) {
        char tmp[2];
        tmp[0] = compress_level + '0'; tmp[1] = '\0';
        strcat(out_mode, tmp);
    }

    // Load FASTA reference (also needed for SAM -> BAM if missing header)
    if (fn_ref) {
        fn_list = samfaipath(fn_ref);
        fai = fai_load(fn_ref);
    }
    // open file handlers
    if ((in = samopen(argv[optind], in_mode, fn_list)) == 0) {
        fprintf(stderr, "[depad] failed to open \"%s\" for reading.\n", argv[optind]);
        ret = 1;
        goto depad_end;
    }
    if (in->header == 0) {
        fprintf(stderr, "[depad] failed to read the header from \"%s\".\n", argv[optind]);
        ret = 1;
        goto depad_end;
    }
    if (in->header->text == 0 || in->header->l_text == 0) {
        fprintf(stderr, "[depad] Warning - failed to read any header text from \"%s\".\n", argv[optind]);
        assert (0 == in->header->l_text);
        assert (0 == in->header->text);
    }
    if (fn_ref) {
        h = fix_header(in->header, fai);
    } else {
        fprintf(stderr, "[depad] Warning - reference lengths will not be corrected without FASTA reference\n");
        h = in->header;
    }
    if ((out = samopen(fn_out? fn_out : "-", out_mode, h)) == 0) {
        fprintf(stderr, "[depad] failed to open \"%s\" for writing.\n", fn_out? fn_out : "standard output");
        ret = 1;
        goto depad_end;
    }

    // Do the depad
    ret = bam_pad2unpad(in, out, fai);

depad_end:
    // close files, free and return
    if (fai) fai_destroy(fai);
    if (! (in && h == in->header)) bam_header_destroy(h);
    samclose(in);
    samclose(out);
    free(fn_list); free(fn_out);
    return ret;
}

static int usage(int is_long_help)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   samtools depad <in.bam>\n\n");
    fprintf(stderr, "Options: -s       output is SAM (default is BAM)\n");
    fprintf(stderr, "         -S       input is SAM (default is BAM)\n");
    fprintf(stderr, "         -u       uncompressed BAM output (can't use with -s)\n");
    fprintf(stderr, "         -1       fast compression BAM output (can't use with -s)\n");
    fprintf(stderr, "         -T FILE  padded reference sequence file [null]\n");
    fprintf(stderr, "         -o FILE  output file name [stdout]\n");
    fprintf(stderr, "         -?       longer help\n");
    fprintf(stderr, "\n");
    if (is_long_help)
        fprintf(stderr, "Notes:\n\
\n\
  1. Requires embedded reference sequences (before the reads for that reference),\n\
     or ideally a FASTA file of the padded reference sequences (via the -T argument).\n\
\n\
  2. The input padded alignment read's CIGAR strings must not use P or I operators.\n\
\n");
    return 1;
}
