/*  bam_checksum.c -- produces checksums on SAM/BAM/CRAM/FASTA/FASTQ data

    Copyright (C) 2024 Genome Research Ltd.

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

/*
 * This is inspired by Biobambam's bamseqchksum tool written by
 * David Jackson and amended by German Tischler.
 *
 * It computes order agnostic checksums for a variety of SAM fields, allowing
 * validation that all the data is still present at different stages of an
 * analysis pipeline.  This may be useful to detect sequences which have been
 * lost by an aligner, memory corruptions flipping individual sequence bases,
 * or file format decoding errors.
 *
 * We start with something basic such as a FASTQ file, and name, seq and qual
 * checksums should still all match after aligning and sorting.
 */

/*
TODO
- Mechanisms for merging checksums together.  Eg merging two read-groups into
  a new file. This can be done by calling update_hash on each CRC in the rows
  (except "combined") to get a new "all" category which would match a merged
  file checksum.

- Query regions.  When we get differences, this can help bisect the data.
  (If unmapped it's hard, but see "samtools cat -r" for CRAM)

- Remove dead HASH_ADD code.  I don't think it's ever going to be used.
 */

#include <config.h>
#include <stdio.h>
#include <unistd.h>
#include <zlib.h>

#include <htslib/sam.h>
#include <htslib/khash.h>

#include "sam_opts.h"
#include "sam_utils.h"
#include "samtools.h"

// TODO: consider adding hts_crc32 to htslib, similar to our hts_md5
// This exposes libdeflate from htslib instead of punting the configure /
// building requirements downstream.
#ifdef HAVE_LIBDEFLATE
#  include <libdeflate.h>
#  define crc32 libdeflate_crc32
#endif

typedef struct {
    int req_flags, excl_flags;  // BAM flags filtering
    int flag_mask, rev_comp, in_order, sanitize;
    int check_pos, check_cigar, check_mate;
    char *tag_str; // X,Y,Z or "*,X,Y,Z" for negation
    char **tags;   // parsed and split tag_str
    int ntags;
    int64_t nrec;
    int verbose;   // whether to show zero count lines
    int show_qc;   // show pass & fail stats
    FILE *fp;
} opts;

// FIXME: qual+33 is a pain, but only for the benefit of compatability with
// biobambam's bamseqchksum.  It's also wrong for QUAL "*" as it triggers a
// wraparound and turning from BAM's 0xff-run to ASCII makes no sense in a
// checksum.

/* ----------------------------------------------------------------------
 * Utility functions.  Possible candidates for moving to htslib?
 */

#if 1
// Nibble at a time.  This could be sped up further.  Eg see htslib's simd.c.
// That code ought to be expanded upon and exposed from htslib.
//
// However this is still 2.4x quicker than the naive implementation below
// It's now around 8% of CPU for a NovaSeq BAM, so some optimisation is
// possible but we're at deminishing returns.
void fill_seq_qual(opts *o, bam1_t *b, uint8_t *restrict seq_buf,
                   uint8_t *restrict qual_buf) {
    // Tables mapping a pair of nibbles to a pair of ASCII bytes
    static const char code2fwdbase[512] =
        "===A=C=M=G=R=S=V=T=W=Y=H=K=D=B=N"
        "A=AAACAMAGARASAVATAWAYAHAKADABAN"
        "C=CACCCMCGCRCSCVCTCWCYCHCKCDCBCN"
        "M=MAMCMMMGMRMSMVMTMWMYMHMKMDMBMN"
        "G=GAGCGMGGGRGSGVGTGWGYGHGKGDGBGN"
        "R=RARCRMRGRRRSRVRTRWRYRHRKRDRBRN"
        "S=SASCSMSGSRSSSVSTSWSYSHSKSDSBSN"
        "V=VAVCVMVGVRVSVVVTVWVYVHVKVDVBVN"
        "T=TATCTMTGTRTSTVTTTWTYTHTKTDTBTN"
        "W=WAWCWMWGWRWSWVWTWWWYWHWKWDWBWN"
        "Y=YAYCYMYGYRYSYVYTYWYYYHYKYDYBYN"
        "H=HAHCHMHGHRHSHVHTHWHYHHHKHDHBHN"
        "K=KAKCKMKGKRKSKVKTKWKYKHKKKDKBKN"
        "D=DADCDMDGDRDSDVDTDWDYDHDKDDDBDN"
        "B=BABCBMBGBRBSBVBTBWBYBHBKBDBBBN"
        "N=NANCNMNGNRNSNVNTNWNYNHNKNDNBNN";

    static const char code2revbase[512] =
        "==T=G=K=C=Y=S=B=A=W=R=D=M=H=V=N="
        "=TTTGTKTCTYTSTBTATWTRTDTMTHTVTNT"
        "=GTGGGKGCGYGSGBGAGWGRGDGMGHGVGNG"
        "=KTKGKKKCKYKSKBKAKWKRKDKMKHKVKNK"
        "=CTCGCKCCCYCSCBCACWCRCDCMCHCVCNC"
        "=YTYGYKYCYYYSYBYAYWYRYDYMYHYVYNY"
        "=STSGSKSCSYSSSBSASWSRSDSMSHSVSNS"
        "=BTBGBKBCBYBSBBBABWBRBDBMBHBVBNB"
        "=ATAGAKACAYASABAAAWARADAMAHAVANA"
        "=WTWGWKWCWYWSWBWAWWWRWDWMWHWVWNW"
        "=RTRGRKRCRYRSRBRARWRRRDRMRHRVRNR"
        "=DTDGDKDCDYDSDBDADWDRDDDMDHDVDND"
        "=MTMGMKMCMYMSMBMAMWMRMDMMMHMVMNM"
        "=HTHGHKHCHYHSHBHAHWHRHDHMHHHVHNH"
        "=VTVGVKVCVYVSVBVAVWVRVDVMVHVVVNV"
        "=NTNGNKNCNYNSNBNANWNRNDNMNHNVNNN";

    uint8_t *seq = bam_get_seq(b);
    uint8_t *qual = bam_get_qual(b);

    if ((b->core.flag & BAM_FREVERSE) && o->rev_comp) {
        int i, j, len2 = b->core.l_qseq & ~1;
        for (i=0, j=b->core.l_qseq-1; i < len2; i+=2, j-=2) {
            memcpy(&seq_buf[j-1], &code2revbase[(size_t)seq[i>>1]*2], 2);
            qual_buf[j-0] = qual[i+0]+33;
            qual_buf[j-1] = qual[i+1]+33;
        }
        if (i < b->core.l_qseq) {
            seq_buf[j] = "=TGKCYSBAWRDMHVN"[bam_seqi(seq, i)];
            qual_buf[j] = qual[i]+33;
        }
    } else {
        int i, j, len2 = b->core.l_qseq & ~1;
        for (i = j = 0; i < len2; i+=2, j++) {
            // Note size_t cast helps gcc optimiser.
            memcpy(&seq_buf[i], &code2fwdbase[(size_t)seq[j]*2], 2);
            // Simple, but a union approach is a little faster with clang.
            qual_buf[i+0] = qual[i+0]+33;
            qual_buf[i+1] = qual[i+1]+33;
        }
        if (i < b->core.l_qseq) {
            seq_buf[i] = seq_nt16_str[bam_seqi(seq, i)];
            qual_buf[i] = qual[i]+33;
        }
    }
}

#else
// Simple version
void fill_seq_qual(opts *o, bam1_t *b, uint8_t *restrict seq_buf,
                   uint8_t *restrict qual_buf) {
    uint8_t *seq = bam_get_seq(b);
    uint8_t *qual = bam_get_qual(b);

    if ((b->core.flag & BAM_FREVERSE) && o->rev_comp) {
        for (int i=0, j=b->core.l_qseq-1; i < b->core.l_qseq; i++,j--) {
            seq_buf[j] = "=TGKCYSBAWRDMHVN"[bam_seqi(seq, i)];
            qual_buf[j] = qual[i]+33;
        }
    } else {
        for (int i = 0; i < b->core.l_qseq; i++) {
            seq_buf[i] = seq_nt16_str[bam_seqi(seq, i)];
            qual_buf[i] = qual[i]+33;
        }
    }
}
#endif


/* ----------------------------------------------------------------------
 * Checksum aggregation
 */

/*
 * The hash is multiplicative within a finite field, modulo PRIME.
 * We need to avoid zeros, and the data type has to be large enough to ensure
 * no wraparound happens (other than the intended modulo).
 */
#define PRIME ((1u<<31)-1)
#ifdef HASH_ADD
uint64_t update_hash(uint64_t hash, uint32_t crc) {
    return (hash + crc) % PRIME;
}
#else
uint64_t update_hash(uint64_t hash, uint32_t crc) {
    crc &= PRIME;
    if (crc == 0 || crc == PRIME)
        crc = 1;

    return (hash * crc) % PRIME;
}
#endif

typedef struct {
    uint64_t seq[3];   // flag + seq
    uint64_t name[3];  // name + flag + seq
    uint64_t qual[3];  // flag + seq + qual
    uint64_t aux[3];   // flag + seq + aux
    uint64_t pos[3];   // flag + seq + chr/pos
    uint64_t cigar[3]; // flag + seq + cigar
    uint64_t mate[3];  // flag + seq + rnext/pnext/tlen
    uint64_t count[3];
} sums_t;

typedef struct {
    uint32_t seq;
    uint32_t name;
    uint32_t qual;
    uint32_t aux;
    uint32_t pos;
    uint32_t cigar;
    uint32_t mate;
} crcs_t;

KHASH_MAP_INIT_STR(chk, sums_t)

void
sums_update(int qcfail, sums_t *h32, const crcs_t *c, opts *o) {
    uint32_t count_crc = 0;
    if (o->in_order) {
        uint8_t c[8];
        u64_to_le(h32->count[0], c);
        count_crc = crc32(0L, c, 8);
    }

    for (int i = 0, j = 0; i < 1 + o->show_qc; i++, j++) {
        h32->seq[j]  = update_hash(h32->seq[j],  count_crc ^ c->seq);
        h32->name[j] = update_hash(h32->name[j], count_crc ^ c->name);
        h32->qual[j] = update_hash(h32->qual[j], count_crc ^ c->qual);
        h32->aux[j]  = update_hash(h32->aux[j],  count_crc ^ c->aux);
        h32->pos[j]  = update_hash(h32->pos[j],  count_crc ^ c->pos);
        h32->cigar[j]= update_hash(h32->cigar[j],count_crc ^ c->cigar);
        h32->mate[j] = update_hash(h32->mate[j], count_crc ^ c->mate);
        h32->count[j]++;
        if (qcfail)
            j++;  // 0/1 for pass, 0/2 for fail
    }
}

void sums_report(opts *o, sums_t *h32, const char *set) {
    for (int i = 0; i < 1 + 2*o->show_qc; i++) {
        uint64_t hc = 1;
        char *pass[] = {"all", "pass", "fail"};
        if (!o->verbose && !h32->count[i])
            continue;

        fprintf(o->fp, "%-10s %-4s %10"PRIu64"  %08"PRIx64"  %08"PRIx64
                "  %08"PRIx64"  %08"PRIx64, set, pass[i], h32->count[i],
                h32->seq[i], h32->name[i], h32->qual[i], h32->aux[i]);
        if (o->check_pos)
            fprintf(o->fp, "  %08"PRIx64, h32->pos[i]);
        if (o->check_cigar)
            fprintf(o->fp, "  %08"PRIx64, h32->cigar[i]);
        if (o->check_mate)
            fprintf(o->fp, "  %08"PRIx64, h32->mate[i]);

        // Merge all
        hc = update_hash(hc, h32->count[i]>>32);
        hc = update_hash(hc, h32->count[i] & 0xffffffff);
        hc = update_hash(hc, h32->seq[i]);
        hc = update_hash(hc, h32->name[i]);
        hc = update_hash(hc, h32->seq[i]);
        hc = update_hash(hc, h32->aux[i]);
        if (o->check_pos)
            hc = update_hash(hc, h32->pos[i]);
        if (o->check_cigar)
            hc = update_hash(hc, h32->cigar[i]);
        if (o->check_mate)
            hc = update_hash(hc, h32->mate[i]);

        fprintf(o->fp, "  %08"PRIx64"\n", hc);
    }
}

void sums_init(sums_t *h32) {
    for (int i = 0; i < 3; i++) {
        h32->seq[i]   = 1;
        h32->name[i]  = 1;
        h32->qual[i]  = 1;
        h32->aux[i]   = 1;
        h32->pos[i]   = 1;
        h32->cigar[i] = 1;
        h32->mate[i]  = 1;
        h32->count[i] = 0;
    }
}

/* ----------------------------------------------------------------------
 * Main checksumming algorithm
 */

/*
 * Canonicalised integer tags.
 * We can store CcSsIi for unsigned and signed char, short and integer.
 * (This can also happen for B arrays, but we don't yet canonicalise these.)
 *
 * Unfortunately some BAMs have degenerate encs, eg XAs\000\001 for XA:s:1.
 * Also CRAM's computed NM can change, so NM:i:0 could be NMc0 or NMC0.
 *
 * Rules: unsigned if >= 0
 *        smallest encoding necessary
 *
 * Returns a tag pointer (possibly local static, or original ptr),
 *         plus rewrites *tag_len if needed.
 */
uint8_t *canonical_tag(uint8_t *tag, size_t *tag_len) {
    switch (tag[2]) {
        static uint8_t ct[7], code;
        int64_t val;

    case 'C': case 'c':
    case 'S': case 's':
    case 'I': case 'i':
        val = bam_aux2i(tag+2);
        if (val >= 0) {
            if (val <= 255) code = 'C';
            else if (val <= 65535) code = 'S';
            else code = 'i';
        } else {
            if (val >= -128 && val <= 127) code = 'c';
            else if (val >= -32768 && val <= 32767) code = 's';
            else code = 'i';
        }
        if (code == tag[2])
            return tag;

        // Otherwise rewrite it;
        ct[0] = tag[0];
        ct[1] = tag[1];
        ct[2] = code;
        switch (code) {
        case 'C': case 'c':
            ct[3] = val;
            *tag_len = 4;
            break;

        case 'S': case 's':
            // Don't care about sign as it's defined anyway
            u16_to_le(val, tag+3);
            *tag_len = 5;
            break;

        case 'I': case 'i':
            // Don't care about sign as it's defined anyway
            u32_to_le(val, tag+3);
            *tag_len = 7;
            break;
        }
        return ct;

    default:
        return tag;
    }
}

// Qsort callback, by integer
static int tag_qsort(const void *t1, const void *t2) {
    return *(const int *)t1 - *(const int *)t2;
}

/*
 * Produces a concatenated string of aux tags in <ID><TYPE><VAL> binary
 * representation,  with the tag names and orders defined in tag_ids[],
 * checksums it, and combines it with the flag-seq CRC.
 * If *tag_str is "*" then we negate tag_ids and encode everything but those.
 * This is a bit trickier as we can no longer use the order specified and
 * instead encode in ASCII sorted order instead.
 *
 * If the read-group is found in the RG:Z: aux, this is returned in
 * the *RGZ ptr (which points to the <VAL> field.
 *
 * Returns 0 on success, updating *crc_aux,
 *        -1 on error
 */
int hash_aux(bam1_t *b, kstring_t *ks, int ntags,
             char **tag_ids,
             uint8_t **tag_ptr, size_t *tag_len,
             const char *tag_str, short (*tag_keep)[125],
             uint32_t crc_seq, uint32_t *crc_aux,
             uint8_t **RGZ) {
    size_t aux_len = bam_get_l_aux(b);
    // 1 byte minimum forces a non-NULL pointer so CRC works
    if (ks_resize(ks, aux_len+1) < 0)
        return -1;
    uint8_t *aux_ptr = (uint8_t *)ks->s;

    // Pass 1: find all tags to copy and their lengths
    uint8_t *aux = bam_aux_first(b), *aux_next;
    memset(tag_len, 0, ntags * sizeof(*tag_len));

    if (*tag_str == '*') {
        // All tags bar specific ones, in alphanumeric order.
        // Select the tags by name on pass 1, then sort by name to get
        // a canonical order, and finally concatenate tags in order.
        int *tag_id = malloc(aux_len/4 * sizeof(*tag_id));
        ntags = 0;
        while (aux) {
            if (aux[-2] == 'R' && aux[-1] == 'G' && aux[0] == 'Z' && RGZ)
                *RGZ = aux+1;
            aux_next = bam_aux_next(b, aux);
            if (tag_keep[aux[-2]-'0'][aux[-1]-'0'] == 0) {
                size_t tag_sz = aux_next
                    ? aux_next - aux
                    : b->data + b->l_data - aux + 2;
                tag_id[ntags] = (aux[-2]<<24) | (aux[-1]<<16) | ntags;
                tag_ptr[ntags] = aux-2;
                tag_len[ntags] = tag_sz;
                ntags++;
            }

            aux = aux_next;
        }

        // Sort
        qsort(tag_id, ntags, sizeof(*tag_id), tag_qsort);

        // Now we have tag_ptr2 in order of occurrence and tag_id in
        // lexicalgraphical order.  Stitch together
        for (int i = 0; i < ntags; i++) {
            int orig_pos = tag_id[i]&0xffff;
            size_t len = tag_len[orig_pos];
            uint8_t *tag = canonical_tag(tag_ptr[orig_pos], &len);
            memcpy(aux_ptr, tag, len);
            aux_ptr += len;
        }

        free(tag_id);
    } else {
        // Selected tags only, in the order requested
        while (aux) {
            if (aux[-2] == 'R' && aux[-1] == 'G' && aux[0] == 'Z' && RGZ)
                *RGZ = aux+1;
            aux_next = bam_aux_next(b, aux);
            if (!(aux[-2] >= '0' && aux[-2] <= 'z' &&
                  aux[-1] >= '0' && aux[-2] <= 'z'))
                continue; // skip illegal tag names
            int i = tag_keep[aux[-2]-'0'][aux[-1]-'0']-1;
            if (i>=0) {
                // found one
                size_t tag_sz = aux_next
                    ? aux_next - aux
                    : b->data + b->l_data - aux + 2;

                tag_ptr[i] = aux-2;
                tag_len[i] = tag_sz;
            }

            aux = aux_next;
        }

        // Pass 2: copy tags in the order we requested
        for (int i = 0; i < ntags; i++) {
            if (tag_len[i]) {
                size_t len = tag_len[i];
                uint8_t *tag = canonical_tag(tag_ptr[i], &len);
                memcpy(aux_ptr, tag, len);
                aux_ptr += len;
            }
        }
    }

    *crc_aux = crc32(crc_seq, (uint8_t *)ks->s, aux_ptr - (uint8_t *)ks->s);

    return 0;
}

int checksum(sam_global_args *ga, opts *o, char *fn) {
    samFile *fp = NULL;
    sam_hdr_t *hdr = NULL;
    bam1_t *b = bam_init1();
    char **tags = o->tags;
    int ntags = o->ntags;
    uint8_t **tag_ptr = calloc(65536, sizeof(*tag_ptr));
    size_t   *tag_len = calloc(65536, sizeof(*tag_len));
    kstring_t aux_ks  = KS_INITIALIZE;
    kstring_t seq_ks  = KS_INITIALIZE;
    kstring_t qual_ks = KS_INITIALIZE;
    khash_t(chk) *h = kh_init(chk);
    int ret = -1;
    int64_t nrec = o->nrec;

    if (!b || !tag_ptr || !tag_len || !h)
        goto err;

#ifndef HTS_LITTLE_ENDIAN
    kstring_t cigar_ks = KS_INITIALIZE;
#endif

    // A precomputed lookup table to speed up selection of tags
    short tag_keep[125][125] = {0};
    for (int i = 0; i < ntags; i++) {
        char *t = tags[i];
        if (t[0] != '*' &&
            !(t[0] >= '0' && t[0] <= 'z' &&
              t[1] >= '0' && t[1] <= 'z')) {
            fprintf(stderr, "[checksum] Illegal tag ID '%.2s'\n", t);
            goto err;
        }
        tag_keep[t[0]-'0'][t[1]-'0'] = i+1;
    }

    sums_t h32, noRG;
    sums_init(&h32);
    sums_init(&noRG);
    uint32_t crc32_start = crc32(0L, Z_NULL, 0);

    fp = sam_open_format(fn, "r", &ga->in);
    if (!fp) {
        print_error_errno("checksum", "Cannot open input file \"%s\"", fn);
        goto err;
    }

    if (ga->nthreads > 0)
        hts_set_threads(fp, ga->nthreads);

    if (!(hdr = sam_hdr_read(fp)))
        goto err;

    int r;
    while ((r = sam_read1(fp, hdr, b)) >= 0) {
        crcs_t c;

        if (b->core.flag & o->excl_flags)
            continue;

        if ((b->core.flag & o->req_flags) != o->req_flags)
            continue;

        if (o->sanitize)
            bam_sanitize(hdr, b, o->sanitize);

        // 8 bits of flag corresponding to original instrument data
        uint8_t flags = b->core.flag & o->flag_mask;

        // Copy sequence out from nibble to base, and reverse complement
        // seq / qual if required.  Qual is +33 (ASCII format) only for
        // compatibility with biobambam's bamseqchksum tool.
        // The +1 here and elsewhere is to force zero byte allocations to
        // always return a pointer rather than NULL.  This in turn prevents
        // crc32() from considering it as a reinitialisation.
        if (ks_resize(&seq_ks, b->core.l_qseq+1) < 0 ||
            ks_resize(&qual_ks, b->core.l_qseq+1) < 0)
            goto err;

        fill_seq_qual(o, b, (uint8_t *)seq_ks.s, (uint8_t *)qual_ks.s);

        // flag + seq
        uint32_t crc = crc32(crc32_start, &flags, 1);
        c.seq = crc32(crc, (uint8_t *)seq_ks.s, b->core.l_qseq);

        // name + flag + seq.
        // flag + seq + name would be faster, but bamseqchksum does this.
        // Also include single nul for compatibility too.
        crc = crc32(crc32_start, (uint8_t *)bam_get_qname(b),
                    b->core.l_qname - b->core.l_extranul);
        crc = crc32(crc, &flags, 1);
        c.name = crc32(crc, (uint8_t *)seq_ks.s, b->core.l_qseq);

        // flag + seq + qual
        c.qual = crc32(c.seq, (uint8_t *)qual_ks.s, b->core.l_qseq);

        // flag + seq + aux tags
        uint8_t *RGZ = NULL;
        if (hash_aux(b, &aux_ks, ntags, tags, tag_ptr, tag_len,
                     o->tag_str, tag_keep, c.seq, &c.aux, &RGZ) < 0)
            goto err;

        // + chr + pos
        if (o->check_pos) {
            uint8_t chr_pos[12];
            u32_to_le(b->core.tid, chr_pos);
            u64_to_le(b->core.pos, chr_pos+4);
            c.pos = crc32(c.seq, chr_pos, 12);
        }

        // + mate rnext/pnext/tlen
        if (o->check_mate) {
            uint8_t mate[4+8+8];
            u32_to_le(b->core.mtid, mate);
            u64_to_le(b->core.mpos, mate+4);
            u64_to_le(b->core.isize, mate+12);
            c.mate = crc32(c.seq, mate, 12);
        }

        // + cigar + mapq
        if (o->check_cigar) {
            uint8_t *cigar = (uint8_t *)bam_get_cigar(b);
#ifndef HTS_LITTLE_ENDIAN
            if (ks_resize(&cigar_ks, 4 * b->core.n_cigar+1) < 0)
                goto err;
            uint32_t *cig32 = bam_get_cigar(b);
            cigar = (uint8_t *)cigar_ks.s;

            for (int i = 0; i < b->core.n_cigar; i++)
                u32_to_le(cig32[i], cigar + 4*i);
#endif
            uint8_t mapq[4];
            u32_to_le(b->core.qual, mapq);
            c.cigar = crc32(c.seq, mapq, 4);
            c.cigar = crc32(c.cigar, cigar, 4 * b->core.n_cigar);
        }

        // Aggregate checksum hashes
        sums_update(b->core.flag & BAM_FQCFAIL, &h32, &c, o);

        if (RGZ) {
            sums_t *h32p;

            // create func
            int ret;
            khiter_t k = kh_get(chk, h, (char *)RGZ);
            if (k == kh_end(h)) {
                char *rgz_;
                k = kh_put(chk, h, rgz_ = strdup((char *)RGZ), &ret);
                if (ret < 0)
                    goto err;
                sums_init(&kh_value(h, k));
            }
            h32p = &kh_value(h, k);

            sums_update(b->core.flag & BAM_FQCFAIL, h32p, &c, o);
        } else {
            sums_update(b->core.flag & BAM_FQCFAIL, &noRG, &c, o);
        }

        if (nrec && --nrec == 0)
            break;
    }

    // Report hashes
    fprintf(o->fp, "# Checksum for file: %s\n", fn);
    fprintf(o->fp, "# Aux tags:          %s\n", o->tag_str);
    char *s=bam_flag2str(o->flag_mask);
    fprintf(o->fp, "# BAM flags:         %s\n", s);
    free(s);
    fprintf(o->fp, "\n# Group    QC        count  flag+seq  +name     +qual     +aux    ");
    if (o->check_pos)
        fprintf(o->fp, "  +chr/pos");
    if (o->check_cigar)
        fprintf(o->fp, "  +cigar  ");
    if (o->check_mate)
        fprintf(o->fp, "  +mate   ");
    fprintf(o->fp, "  combined\n");
    sums_report(o, &h32,  "all");
    if (o->verbose || (noRG.count[0]+noRG.count[1]))
        sums_report(o, &noRG,  "-");
    for (khiter_t k = kh_begin(h); k != kh_end(h); k++) {
        if (!kh_exist(h, k))
            continue;

        sums_report(o, &kh_value(h, k), kh_key(h, k));

    }

    if (r < -1)
        goto err;

    if (sam_close(fp) < 0) {
        fp = NULL;
        print_error_errno("checksum", "Closing input file \"%s\"", fn);
        goto err;
    }
    fp = NULL;
    ret = 0;

 err:
    if (b)   bam_destroy1(b);
    if (hdr) sam_hdr_destroy(hdr);
    if (fp)  sam_close(fp);

    free(tag_ptr);
    free(tag_len);
    ks_free(&aux_ks);
    ks_free(&seq_ks);
    ks_free(&qual_ks);
#ifndef HTS_LITTLE_ENDIAN
    ks_free(&cigar_ks);
#endif

    if (h) {
        for (khiter_t k = kh_begin(h); k != kh_end(h); k++) {
            if (!kh_exist(h, k))
                continue;

            free((char *)kh_key(h, k));
        }
        kh_destroy(chk, h);
    }

    return ret;
}

void usage_exit(FILE *fp, int ret) {
    fprintf(stderr, "Usage: samtools checksum [options] [file ...]\n");
    fprintf(stderr, "Options:\n\
  -F, --exclude-flags FLAG    Filter if any FLAGs are present [0x900]\n\
  -f, --require-flags FLAG    Filter unless all FLAGs are present [0]\n\
  -b, --flag-mask FLAG        BAM FLAGs to use in checksums [0x0c1]\n\
  -c, --no-rev-comp           Do not reverse-complement sequences [off]\n\
  -t, --tags STR[,STR]        Select tags to checksum [BC,FI,QT,RT,TC]\n\
  -O, --in-order              Use order-specific checksumming [off]\n\
  -P, --check-pos             Also checksum CHR / POS [off]\n\
  -C, --check-cigar           Also checksum MAPQ / CIGAR [off]\n\
  -M, --check_mate            Also checksum PNEXT / RNEXT / TLEN [off]\n\
  -z, --sanitize FLAGS        Perform sanity checks and fix records [off]\n\
  -N, --count INT             Stop after INT number of records [0]\n\
  -o, --output FILE           Write report to FILE [stdout]\n\
  -q, --show-qc               Also show QC pass/fail lines\n\
  -v, --verbose               Increase verbosity: show lines with 0 counts\n\
  -a, --all                   Check all: -PCMOc -b 0xfff -f0 -F0 -z all,cigarx\n");
    fprintf(fp, "\nGlobal options:\n");
    sam_global_opt_help(fp, "-.---@--");
    exit(ret);
}

int parse_tags(opts *o) {
    // Count
    int nt = 0;
    for (char *t = o->tag_str; *t; t++) {
        nt++;
        char *l = t;
        while (*t && *t != ',')
            t++;
        if (t-l != 2 && !(t-l == 1 && *l == '*')) {
            fprintf(stderr, "Bad tag string.  Should be XX,YY,ZZ syntax\n");
            return 1;
        }
        if (!*t)
            break;
    }

    // Split by tag
    o->ntags = nt;
    o->tags = calloc(nt, sizeof(*o->tags));
    if (!o->tags)
        return 1;

    nt = 0;
    for (char *t = o->tag_str; *t; t++, nt++) {
        o->tags[nt] =  t;
        while (*t && *t != ',')
            t++;
        if (!*t)
            break;
    }

    return 0;
}

/* ----------------------------------------------------------------------
 */
int main_checksum(int argc, char **argv) {
    opts opts = {
        .req_flags    = 0,
        .excl_flags   = BAM_FSECONDARY | BAM_FSUPPLEMENTARY,
        .flag_mask    = BAM_FPAIRED | BAM_FREAD1 | BAM_FREAD2,
        .rev_comp     = 1,
        .tag_str      = "BC,FI,QT,RT,TC",
        .check_pos    = 0,
        .check_cigar  = 0,
        .check_mate   = 0,
        .in_order     = 0,
        .sanitize     = 0,
        .nrec         = 0,
        .verbose      = 0,
        .show_qc      = 0,
        .fp           = stdout,
    };

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 'I', '-', '-', '.', '@'),
        {"--exclude-flags", required_argument, NULL, 'F'},
        {"--require-flags", required_argument, NULL, 'f'},
        {"--flag-mask",     required_argument, NULL, 'b'},
        {"--tags",          required_argument, NULL, 't'},
        {"--no-rev-comp",   no_argument,       NULL, 'c'},
        {"--in-order",      no_argument,       NULL, 'O'},
        {"--check-pos",     no_argument,       NULL, 'P'},
        {"--check-cigar",   no_argument,       NULL, 'C'},
        {"--check-mate",    no_argument,       NULL, 'M'},
        {"--count",         required_argument, NULL, 'N'},
        {"--sanitize",      required_argument, NULL, 'z'},
        {"--output",        required_argument, NULL, 'o'},
        {"--show-qc",       no_argument,       NULL, 'q'},
        {"--verbose",       no_argument,       NULL, 'v'},
        {"--all",           no_argument,       NULL, 'a'},
        {NULL, 0, NULL, 0}
    };

    if (argc == 1 && isatty(STDIN_FILENO))
        usage_exit(stdout, EXIT_SUCCESS);

    int c;
    while ((c = getopt_long(argc, argv, "@:f:F:t:cPCMOb:z:aN:vqo:",
                            lopts, NULL)) >= 0) {
        switch (c) {
        case 'O':
            opts.in_order = 1;
            break;
        case 'F':
            if ((opts.excl_flags = bam_str2flag(optarg)) < 0) {
                print_error("checksum", "could not parse flag %s", optarg);
                return 1;
            }
            break;
        case 'f':
            if ((opts.req_flags  = bam_str2flag(optarg)) < 0) {
                print_error("checksum", "could not parse flag %s", optarg);
                return 1;
            }
            break;
        case 'b':
            if ((opts.flag_mask  = bam_str2flag(optarg)) < 0) {
                print_error("checksum", "could not parse flag %s", optarg);
                return 1;
            }
            break;
        case 'P':
            opts.check_pos = 1;
            break;
        case 'C':
            opts.check_cigar = 1;
            break;
        case 'M':
            opts.check_mate = 1;
            break;
        case 't':
            opts.tag_str = optarg;
            break;
        case 'c':
            opts.rev_comp = 0;
            break;
        case 'N':
            opts.nrec = strtoll(optarg, NULL, 0);
            break;

        case 'v':
            opts.verbose++;
            break;
        case 'q':
            opts.show_qc=1;
            break;

        case 'z':
            if ((opts.sanitize = bam_sanitize_options(optarg)) < 0)
                return 1;
            break;

        case 'a':
            // ALL: a shorthand for a bunch of options to checksum the entire
            // file contents.  TODO: we still need tag wildcards.
            opts.req_flags = 0;
            opts.excl_flags = 0;
            opts.flag_mask = -1;
            opts.rev_comp = 0;
            opts.in_order = 1;
            opts.check_pos = 1;
            opts.check_cigar = 1;
            opts.check_mate = 1;
            opts.sanitize = FIX_ALL | FIX_CIGARX;
            opts.tag_str = "*,cF";
            break;

        case 'o':
            opts.fp = fopen(optarg, "w");
            if (!opts.fp) {
                perror(optarg);
                return 1;
            }
            break;

        default:
            if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0)
                break;
            /* else fall-through */
        case '?':
            usage_exit(stderr, EXIT_FAILURE);
        }
    }

    if (!opts.tags) {
        if (parse_tags(&opts) < 0)
            return 1;
    }

    int ret = 0;
    if (argc-optind) {
        while (optind < argc)
            ret |= checksum(&ga, &opts, argv[optind++]) < 0;
    } else {
        ret = checksum(&ga, &opts, "-") < 0;
    }

    if (opts.fp != stdout)
        ret |= fclose(opts.fp) < 0;

    free(opts.tags);

    return ret;
}
