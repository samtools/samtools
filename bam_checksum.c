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
- Support multiple read-groups, which aids spliting pooled samples and
  tracking that data isn't lost.
- Separate "all" from "pass" only (dropping QC fail)
- Mechanisms for merging checksums together.  Eg merging two read-groups into
  a new file.
- Make tags configurable
- More components so we can checksum also any combo.  Eg CIGAR, MAPQ,
  RNEXT/PNEXT/TLEN, etc.  This provides a route to detecting more types of
  data loss. (Also see bam_mate.c:bam_sanitize() function)
- Query regions.  When we get differences, this can help bisect the data.
  (If unmapped it's hard, but see "samtools cat -r" for CRAM)
 */

#include <config.h>
#include <stdio.h>
#include <unistd.h>
#include <zlib.h>

#include <htslib/sam.h>
#include "sam_opts.h"
#include "sam_utils.h"

// TODO: consider adding hts_crc32 to htslib, similar to our hts_md5
// This exposes libdeflate from htslib instead of punting the configure /
// building requirements downstream.
#ifdef HAVE_LIBDEFLATE
#  include <libdeflate.h>
#  define crc32 libdeflate_crc32
#endif

typedef struct {
    int incl_flags, req_flags, excl_flags;  // BAM flags filtering
    int rev_comp;
} opts;

// FIXME: qual+33 is a pain, but only for the benefit of compatability with
// biobambam's bamseqchksum.  It's also wrong for QUAL "*" as it triggers a
// wraparound and turning from BAM's 0xff-run to ASCII makes no sense in a
// checksum.

#if 1
// Nibble at a time.  This could be sped up further.  Eg see htslib's simd.c.
// That code ought to be expanded upon and exposed from htslib.
//
// However this is still 2.4x quicker than the naive implementation below
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
    uint64_t seq;   // flag + seq
    uint64_t name;  // name + flag + seq
    uint64_t qual;  // flag + seq + qual
    uint64_t aux;   // flag + seq + aux
} hashes;

void
update_hashes(hashes *h32, uint32_t s, uint32_t n, uint32_t q, uint32_t a) {
    h32->seq  = update_hash(h32->seq,  s);
    h32->name = update_hash(h32->name, n);
    h32->qual = update_hash(h32->qual, q);
    h32->aux  = update_hash(h32->aux,  a);
}

#ifdef HASH_ADD
#  define H32_INIT {0,0,0,0}
#else
#  define H32_INIT {1,1,1,1}
#endif

/*
 * Produces a concatenated string of aux tags in <ID><TYPE><VAL> binary
 * representation,  with the tag names and orders defined in tag_ids[],
 * checksums it, and combines it with the flag-seq CRC.
 *
 * Returns 0 on success, updating *crc_aux,
 *        -1 on error
 */
int hash_aux(bam1_t *b, kstring_t *ks, int ntags,
             const char *tag_ids[],
             uint8_t **tag_ptr, size_t *tag_len,
             short (*tag_keep)[125],
             uint32_t crc_seq, uint32_t *crc_aux) {
    size_t aux_len = bam_get_l_aux(b);
    if (ks_resize(ks, aux_len) < 0)
        return -1;
    uint8_t *aux_ptr = (uint8_t *)ks->s;

    // Pass 1: find all tags to copy and their lengths
    uint8_t *aux = bam_aux_first(b), *aux_next;
    memset(tag_len, 0, ntags * sizeof(*tag_len));
    while (aux) {
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
            memcpy(aux_ptr, tag_ptr[i], tag_len[i]);
            aux_ptr += tag_len[i];
        }
    }

    *crc_aux = aux_ptr > (uint8_t *)ks->s
        ? crc32(crc_seq, (uint8_t *)ks->s, aux_ptr - (uint8_t *)ks->s)
        : crc_seq;

    return 0;
}

int checksum(sam_global_args *ga, opts *o, char *fn) {
    samFile *fp = NULL;
    sam_hdr_t *hdr = NULL;
    bam1_t *b = bam_init1();
    static const char *tags[] = {"BC","FI","QT","RT","TC"};
    const int ntags = sizeof(tags) / sizeof(*tags);
    uint8_t **tag_ptr = calloc(ntags, sizeof(*tag_ptr));
    size_t   *tag_len = calloc(ntags, sizeof(*tag_len));
    kstring_t aux_ks  = KS_INITIALIZE;
    kstring_t seq_ks  = KS_INITIALIZE;
    kstring_t qual_ks = KS_INITIALIZE;

    if (!b || !tag_ptr || !tag_len)
        goto err;

    // A precomputed lookup table to speed up selection of tags
    short tag_keep[125][125] = {0};
    for (int i = 0; i < ntags; i++) {
        if (!(tags[i][0] >= '0' && tags[i][0] <= 'z' &&
              tags[i][1] >= '0' && tags[i][1] <= 'z')) {
            fprintf(stderr, "[checksum] Illegal tag ID '%.2s'\n", tags[i]);
            goto err;
        }
        tag_keep[tags[i][0]-'0'][tags[i][1]-'0'] = i+1;
    }

    hashes h32 = H32_INIT;

    const uint32_t crc32_start = crc32(0L, Z_NULL, 0);

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
    uint64_t count = 0;

    while ((r = sam_read1(fp, hdr, b)) >= 0) {
        // TODO: configurable filter
        if (b->core.flag & o->excl_flags)
            continue;

        // 8 bits of flag corresponding to original instrument data
        uint8_t flags = b->core.flag & (BAM_FPAIRED | BAM_FREAD1 | BAM_FREAD2);

        // Copy sequence out from nibble to base, and reverse complement
        // seq / qual if required.  Qual is +33 (ASCII format) only for
        // compatibility with biobambam's bamseqchksum tool.
        if (ks_resize(&seq_ks, b->core.l_qseq) < 0 ||
            ks_resize(&qual_ks, b->core.l_qseq) < 0)
            goto err;

        fill_seq_qual(o, b, (uint8_t *)seq_ks.s, (uint8_t *)qual_ks.s);

        // flag + seq
        uint32_t crc = crc32(crc32_start, &flags, 1);
        uint32_t crc_seq = crc32(crc, (uint8_t *)seq_ks.s, b->core.l_qseq);

        // name + flag + seq.
        // flag + seq + name would be faster, but bamseqchksum does this.
        // Also include single nul for compatibility too.
        crc = crc32(crc32_start, (uint8_t *)bam_get_qname(b),
                    b->core.l_qname - b->core.l_extranul);
        crc = crc32(crc, &flags, 1);
        uint32_t crc_name = crc32(crc, (uint8_t *)seq_ks.s, b->core.l_qseq);

        // flag + seq + qual
        uint32_t crc_qual = crc32(crc_seq, (uint8_t *)qual_ks.s,
                                  b->core.l_qseq);

        // flag + seq + aux tags
        uint32_t crc_aux;
        if (hash_aux(b, &aux_ks, ntags, tags, tag_ptr, tag_len,
                     tag_keep, crc_seq, &crc_aux) < 0)
            goto err;

        // Aggregate hashes
        update_hashes(&h32, crc_seq, crc_name, crc_qual, crc_aux);

        count++;
    }

    printf("Count          %"PRIu64"\n", count);
    printf("Flag+Seq       %08"PRIx64"\n", h32.seq);
    printf("Name+Flag+Seq  %08"PRIx64"\n", h32.name);
    printf("Flag+Seq+Qual  %08"PRIx64"\n", h32.qual);
    printf("Flag+Seq+Aux   %08"PRIx64"\n", h32.aux);
    puts("");

    if (r <= -1)
        goto err;
    if (hdr)
        sam_hdr_destroy(hdr);

    if (sam_close(fp) < 0) {
        print_error_errno("checksum", "Closing input file \"%s\"", fn);
        goto err;
    }

    free(tag_ptr);
    free(tag_len);
    ks_free(&aux_ks);
    ks_free(&seq_ks);
    ks_free(&qual_ks);

    bam_destroy1(b);
    return 0;

 err:
    if (b)   bam_destroy1(b);
    if (hdr) sam_hdr_destroy(hdr);
    if (fp)  sam_close(fp);

    free(tag_ptr);
    free(tag_len);
    ks_free(&aux_ks);
    ks_free(&seq_ks);
    ks_free(&qual_ks);

    return -1;
}

void usage_exit(FILE *fp, int ret) {
    fprintf(stderr, "Usage: samtools checksum [options] [file]\n");
    exit(ret);
}

int main_checksum(int argc, char **argv) {
    opts opts = {
        .incl_flags   = 0xffff,
        .req_flags    = 0,
        .excl_flags   = BAM_FSECONDARY | BAM_FSUPPLEMENTARY,
        .rev_comp     = 1,
    };

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 'I', '-', '-', '.', '@'),
        {"--excl-flags",    required_argument, NULL, 'F'},
        {"--exclude-flags", required_argument, NULL, 'F'},
        {"--require-flags", required_argument, NULL, 'f'},
        {"--incl-flags",    required_argument, NULL, 1},
        {"--include-flags", required_argument, NULL, 1},
        {NULL, 0, NULL, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "@:f:F:", lopts, NULL)) >= 0) {
        switch (c) {
        case 'F': opts.excl_flags = atoi(optarg); break;
        case 'f': opts.req_flags  = atoi(optarg); break;
        case  1 : opts.incl_flags = atoi(optarg); break;
        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
        case '?':
            usage_exit(stderr, EXIT_FAILURE);
        }
    }

    int ret = 0;
    if (argc-optind) {
        while (optind < argc)
            ret |= checksum(&ga, &opts, argv[optind++]) < 0;
    } else {
        ret = checksum(&ga, &opts, "-") < 0;
    }

    return ret;
}
