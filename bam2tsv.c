/*  bam2tsv.c -- tsv subcommand.

    Copyright (C) 2021 Pierre Lindenbaum
    Institut du Thorax. u1087 Nantes. France.
    @yokofakun

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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <unistd.h>
#include <inttypes.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "samtools.h"
#include "bedidx.h"
#include "sam_opts.h"
#include "htslib/thread_pool.h"
#include "htslib/khash.h"

#define WHERE do {fprintf(stderr, "[%s:%d]", __FILE__, __LINE__);} while (0)
#define DEBUG(...) do {WHERE;fprintf(stderr, __VA_ARGS__);fputc('\n', stderr);} while (0)

/** all the information for a given base */
typedef struct base_info {
    /* current record */
    bam1_t* b;
    /* 1-based reference position . Do not use hts_pos_t because can be negative due to clipping */
    int64_t ref1;
    /* reference base */
    char ref_base;
    /* 0-based position in the read */
    hts_pos_t read0;
    /* 0-based position in the unclipped read */
    int64_t unclipped_read0;
    /* cigar operator as char */
    char op_chr;
    /* cigar operator as int */
    char op_int;
    /* read base */
    char base;
    /* read qual */
    char qual;
} BaseInfo, *BaseInfoPtr;

/** program parameters */
typedef struct tsv_context {
    /* output FILE */
    FILE* out;
    /* input bam */
    htsFile* in ;
    /* SAM header */
    sam_hdr_t* header;
    /* fasta reference */
    faidx_t* fai;
    /* skip 'N' cigar operator flag */
    int skip_N;
    /* user query as a string containing some opcodes defined in 'bam2tsv.h' */
    char* query;
    /** we keep a subs-string of REFERENCE, used in function get_reference_base_at */
    /* current reference chromosome */
    int32_t ref_tid;
    /* current reference start 1-based */
    hts_pos_t ref_start;
    /* current reference end 1-based */
    hts_pos_t ref_end;
    /* current reference end sequence */
    char* ref_seq;
} TsvContext, *TsvContextPtr;

/** code copied from markdup.c but use a signed position */
static int64_t unclipped_start(bam1_t* b) {
    uint32_t* cigar = bam_get_cigar(b);
    int64_t clipped = 0;
    uint32_t i;

    for (i = 0; i < b->core.n_cigar; i++) {
        char c = bam_cigar_opchr(cigar[i]);

        if (c == 'S' || c == 'H') { // clips
            clipped += (int64_t) bam_cigar_oplen(cigar[i]);
        } else {
            break;
        }
    }
    return (int64_t) b->core.pos - clipped + 1;
}



/* get the REF base a tid:ref1 */
static int get_reference_base_at(TsvContextPtr ctx, int32_t tid, int64_t ref1, char* base) {
    /* ref1 negative, for example for unclipped part of read */
    if (ref1 < 1 || ctx->fai == NULL) return 'N';
    /* the position is not in the current interval */
    if (!(tid == ctx->ref_tid && ref1 >= ctx->ref_start && ref1 <= ctx->ref_end)) {
        int len;
        /* get len of chromosome tid */
        hts_pos_t ref_len = sam_hdr_tid2len(ctx->header, tid);
        if (ref1 > ref_len) {
            *base = 'N';
            return 0;
        }
        /* release previous memory */
        free(ctx->ref_seq);
        /* we keep a buffer extended by REF_BUFFER_SIZE */
        #define REF_BUFFER_SIZE 1000000
        ctx->ref_start = (ref1 <= REF_BUFFER_SIZE ? 1 : ref1 - REF_BUFFER_SIZE);
        ctx->ref_end = (ref1 + REF_BUFFER_SIZE >= ref_len ? ref_len : ref1 + REF_BUFFER_SIZE);
        /* fetch interval */
        ctx->ref_seq = faidx_fetch_seq(ctx->fai, sam_hdr_tid2name(ctx->header, tid), ctx->ref_start - 1, ctx->ref_end, & len);
        if (ctx->ref_seq == NULL) {
            print_error("tsv","warning. Cannot fetch reference.");
            *base = 'N';
            return -1;
        }
    }
    *base = ctx->ref_seq[ref1 - ctx->ref_start];
    return 0;
}

/** loop over ctx->query and print each column for the current aln */
static int sam2tsv_base(TsvContextPtr ctx, BaseInfoPtr aln) {
    char* p = ctx->query;
    while (*p != 0) {
        if (p != ctx->query) fputc('\t', ctx->out);
        switch (*p) {
            #define CASE_OPCODE(OPCODE, LABEL, DESC, FUN) case OPCODE: FUN;break
            #include "bam2tsv.h"
            #undef CASE_OPCODE
        }
        p++;
    }
    return fputc('\n', ctx->out) == EOF ? -1 : 0;
}

/** set aln read base at 'i' */
#define READ_BASE_AT(i) (read_bases == NULL ? 'N' : seq_nt16_str[bam_seqi(read_bases, i)])
/** set aln read qual at 'i' */
#define READ_QUAL_AT(i) (read_quals == NULL ? '*' : (read_quals[0] == 0xff ? 'B' : read_quals[i] + 33))
/** fetch aln ref base at 'i' */
#define REF_BASE_AT(i) if (get_reference_base_at(ctx, b->core.tid, i, &aln.ref_base) != 0) {\
    print_error("tsv", "Cannot fetch base.");\
    ret = -1;\
    break;\
    }
/** call sam2tsv_base */
#define INVOKE_SAM2TSV_BASE if (sam2tsv_base(ctx, &aln) != 0) {\
        ret = -1;\
        break;\
    }

static int sam2tsv_aln(TsvContextPtr ctx, bam1_t* b) {
    int ret = 0, i, j;
    // read sequence
    uint8_t* read_bases = NULL;
    // read qualities
    uint8_t* read_quals = NULL;
    // current information under each base
    BaseInfo aln;

    // skip unmapped records
    if ((b->core.flag & BAM_FMUNMAP)) return 0;
    // seq and qual
    if (b->core.l_qseq) {
        read_bases = bam_get_seq(b);
        read_quals = bam_get_qual(b);
    }
    uint32_t* cig = bam_get_cigar(b);
    // one based reference positon
    int64_t ref1 = unclipped_start(b);
    int n_cigar = b->core.n_cigar;
    // 0-based position of the read
    int64_t read0 = 0;
    int64_t unclipped_read0 = 0;
    aln.b = b;

    //loop over each cigar element
    for (i = 0; i < n_cigar; i++) {
        aln.op_chr = bam_cigar_opchr(cig[i]);
        aln.op_int = bam_cigar_op(cig[i]);
        int oplen = bam_cigar_oplen(cig[i]);
        switch (aln.op_int) {
        case BAM_CPAD:
            break;
        case BAM_CHARD_CLIP:
            for (j = 0; j < oplen; j++) {
                aln.base = 'N';
                aln.qual = '.';
                REF_BASE_AT(ref1)
                aln.ref1 = ref1;
                aln.read0 = -1;
                aln.unclipped_read0 = unclipped_read0;
                INVOKE_SAM2TSV_BASE
                ref1++;
                unclipped_read0++;
            }
            break;
        case BAM_CINS:
            for (j = 0; j < oplen; j++) {
                aln.base = READ_BASE_AT(read0);
                aln.qual = READ_QUAL_AT(read0);
                aln.ref_base = '-';
                aln.ref1 = -1;
                aln.read0 = read0;
                aln.unclipped_read0 = unclipped_read0;
                INVOKE_SAM2TSV_BASE
                read0++;
                unclipped_read0++;
            }
            break;
        case BAM_CREF_SKIP:
            if (ctx->skip_N) {
                ref1 += oplen;
                break;
            }
            // NO break here, continue
            case BAM_CDEL:
                for (j = 0; j < oplen; j++) {
                    aln.base = '-';
                    aln.qual = '-';
                    REF_BASE_AT(ref1)
                    aln.ref1 = ref1;
                    aln.read0 = -1;
                    aln.unclipped_read0 = -1;
                    INVOKE_SAM2TSV_BASE
                    ref1++;
                }
                break;
            case BAM_CSOFT_CLIP:
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                for (j = 0; j < oplen; j++) {
                    aln.base = READ_BASE_AT(read0);
                    aln.qual = READ_QUAL_AT(read0);
                    REF_BASE_AT(ref1)
                    aln.ref1 = ref1;
                    aln.read0 = read0;
                    aln.unclipped_read0 = unclipped_read0;
                    INVOKE_SAM2TSV_BASE
                    ref1++;
                    read0++;
                    unclipped_read0++;
                }
                break;

            default: {
                print_error("tsv", "Unsupported cigar operator '%c'.", aln.op_chr);
                ret = -1;
                break;
            }
        }

    }
    return ret;
}

/* loop over each sam record */
static int sam2tsv_core(TsvContextPtr param) {
    int ret = 0, r;
    bam1_t* b = bam_init1();
    while ((r = sam_read1(param->in , param->header, b)) >= 0) {
        if (sam2tsv_aln(param, b) != 0) {
            ret = -1;
            break;
        }
    }
    bam_destroy1(b);
    return ret;
}

#define DEFAULT_QUERY "NQFTOBbqRpu"

static void sam2tsv_usage(FILE* fp) {
    fprintf(fp, "Usage: samtools tsv [options] (in.bam|stdin)\n");
    fprintf(fp, "\nOptions:\n");
    fprintf(fp, "  -o|--output FILE      Write output to FILE [stdout]\n");
    fprintf(fp, "  -q|--query query      A string of column operators (see option -l) [" DEFAULT_QUERY "]\n");
    fprintf(fp, "  -N|--skip-N           Skip 'N' cigar operator (reference skip)\n");
    fprintf(fp, "  -l|--list             print available columns operators on stdout and exit.\n");
    sam_global_opt_help(fp, "-.--T@-.");
}

int main_bam2tsv(int argc, char* argv[]) {
    int c, ret = EXIT_SUCCESS;
    char* p = NULL;
    TsvContext ctx;
    char* out_fname = NULL;
    htsThreadPool pool = {
        NULL,
        0
    };
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static
    const struct option lopts[] = {
        {"list", no_argument, NULL, 'l'},
        {"query", required_argument, NULL, 'q'},
        {"skip-N", no_argument, NULL, 'N'},
        {"output", required_argument, NULL, 'o'},
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 0, 0, 'T', '@'),
        {NULL,0,NULL,0}
    };

    memset((void*)&ctx, 0, sizeof(TsvContext));
    ctx.out = stdout;

    while ((c = getopt_long(argc, argv, "Nlo:q:T:@:",
            lopts, NULL)) >= 0) {
        switch (c) {
        case 'o':
            if (strcmp(optarg,"-")!=0) out_fname = optarg;
            break;
        case 'N':
            ctx.skip_N = 1;
            break;
        case 'q':
            free(ctx.query);
            ctx.query = strdup(optarg);
            if (ctx.query==NULL) {
                print_error("tsv", "out of memory.");
                ret = EXIT_FAILURE;
                goto cleanup;
            }
            break;
        case 'l': {
            printf("#operator\tlabel\tdescription\n");
            #define CASE_OPCODE(OPCODE, LABEL, DESC, FUN) printf("%c\t%-15s\t%s\n", OPCODE, LABEL, DESC)
            #include "bam2tsv.h"
            #undef CASE_OPCODE
            return EXIT_SUCCESS;
        }
        default:
            if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
        case '?':
            sam2tsv_usage(stderr);
            exit(EXIT_FAILURE);
        }
    }

    if (optind == argc && isatty(STDIN_FILENO)) {
       sam2tsv_usage(stderr);
       exit(EXIT_FAILURE);
    }

    if (ctx.query == NULL) {
        ctx.query = strdup(DEFAULT_QUERY);
        if (ctx.query==NULL) {
            print_error("tsv", "out of memory.");
            ret = EXIT_FAILURE;
            goto cleanup;
        }
    }

    // check query
    p = ctx.query;
    while ( *p != 0) {
        switch (*p) {
            #define CASE_OPCODE(OPCODE, LABEL, DESC, FUN) case OPCODE: break
            #include "bam2tsv.h"
            #undef CASE_OPCODE
        default:
            print_error("tsv", "In query \"%s\" unknown opcode \"%c\".", ctx.query, * p);
            ret = EXIT_FAILURE;
            goto cleanup;
        }
        p++;
    }

    /* load reference index */
    if (ga.reference == NULL) {
        print_error("tsv", "undefined reference.\n");
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    char* fn_fai = fai_path(ga.reference);
    if (fn_fai == NULL) {
        print_error("tsv", "Cannot get fasta index.");
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    ctx.fai = fai_load3(ga.reference, fn_fai, NULL, FAI_CREATE);
    if (ctx.fai == NULL) {
        print_error("tsv", "Cannot load fasta index.");
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    /* open the SAM */
    ctx.in = sam_open_format(optind == argc ? "-" : argv[optind], "r", &ga.in);
    if (!ctx.in) {
        print_error_errno("tsv", "failed to open \"%s\" for input.", argv[optind]);
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    if (fn_fai && hts_set_fai_filename(ctx.in, fn_fai) != 0) {
        print_error("tsv","Failed to load reference file \"%s\".\n", fn_fai);
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    ctx.header = sam_hdr_read(ctx.in);
    if (ctx.header == NULL) {
        print_error("tsv", "cannot read SAM header.\n");
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    if (out_fname != NULL) {
        ctx.out = fopen(out_fname, "w");
        if (!ctx.out) {
            print_error_errno("tsv", "failed to open \"%s\" for writing", out_fname);
            ret = EXIT_FAILURE;
            goto cleanup;
        }
    }

    //print header
    p = ctx.query;
    while (*p != 0) {
        fputc(p == ctx.query ? '#' : '\t', ctx.out);
        switch (*p) {
            #define CASE_OPCODE(OPCODE, LABEL, DESC, FUN) case OPCODE: fputs(LABEL, ctx.out);break
            #include "bam2tsv.h"
            #undef CASE_OPCODE
        }
        p++;
    }
    fputc('\n', ctx.out);

    if (ga.nthreads > 0) {
        if (!(pool.pool = hts_tpool_init(ga.nthreads))) {
            print_error("tsv","error creating thread pool\n");
            return 1;
        }

        hts_set_opt(ctx.in, HTS_OPT_THREAD_POOL, &pool);
    }

    ret = sam2tsv_core( &ctx);

    cleanup:
        if (ctx.in != NULL) {
            sam_close(ctx.in);
        }
    if (out_fname != NULL) {
        fflush(ctx.out);
        fclose(ctx.out);
        free(out_fname);
    }
    if (ctx.header != NULL) {
        sam_hdr_destroy(ctx.header);
    }
    if (ctx.fai) fai_destroy(ctx.fai);
    if (pool.pool) hts_tpool_destroy(pool.pool);
    sam_global_args_free( &ga);
    free(ctx.query);
    return ret;
}
