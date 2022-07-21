/*  bam_fastq.c -- FASTA and FASTQ file generation

    Copyright (C) 2009-2017, 2019-2020 Genome Research Ltd.
    Portions copyright (C) 2009, 2011, 2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notices and this permission notice shall be included in
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
#include <strings.h>
#include <stdbool.h>
#include <ctype.h>
#include <assert.h>
#include <inttypes.h>
#include <unistd.h>

#include "htslib/sam.h"
#include "htslib/klist.h"
#include "htslib/kstring.h"
#include "htslib/bgzf.h"
#include "htslib/thread_pool.h"
#include "samtools.h"
#include "sam_opts.h"

#define DEFAULT_BARCODE_TAG "BC"
#define DEFAULT_QUALITY_TAG "QT"
#define INDEX_SEPARATOR "+"

int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };
static void bam2fq_usage(FILE *to, const char *command)
{
    int fq = strcasecmp("fastq", command) == 0 || strcasecmp("bam2fq", command) == 0;
    fprintf(to,
"Usage: samtools %s [options...] <in.bam>\n", command);
    fprintf(to,
"\n"
"Description:\n"
"Converts a SAM, BAM or CRAM to %s format.\n"
"\n"
"Options:\n"
"  -0 FILE      write reads designated READ_OTHER to FILE\n"
"  -1 FILE      write reads designated READ1 to FILE\n"
"  -2 FILE      write reads designated READ2 to FILE\n"
"  -o FILE      write reads designated READ1 or READ2 to FILE\n"
"               note: if a singleton file is specified with -s, only\n"
"               paired reads will be written to the -1 and -2 files.\n"
"  -f INT       only include reads with all  of the FLAGs in INT present [0]\n"       //   F&x == x
"  -F INT       only include reads with none of the FLAGS in INT present [0x900]\n"       //   F&x == 0
"  -G INT       only EXCLUDE reads with all  of the FLAGs in INT present [0]\n"       // !(F&x == x)
"  -n           don't append /1 and /2 to the read name\n"
"  -N           always append /1 and /2 to the read name\n",
    fq ? "FASTQ" : "FASTA");
    if (fq) fprintf(to,
"  -O           output quality in the OQ tag if present\n");
    fprintf(to,
"  -s FILE      write singleton reads designated READ1 or READ2 to FILE\n"
"  -t           copy RG, BC and QT tags to the %s header line\n",
    fq ? "FASTQ" : "FASTA");
    fprintf(to,
"  -T TAGLIST   copy arbitrary tags to the %s header line, '*' for all\n",
    fq ? "FASTQ" : "FASTA");
    if (fq) fprintf(to,
"  -v INT       default quality score if not given in file [1]\n"
"  -i           add Illumina Casava 1.8 format entry to header (eg 1:N:0:ATCACG)\n"
"  -c INT       compression level [0..9] to use when writing bgzf files [1]\n"
"  --i1 FILE    write first index reads to FILE\n"
"  --i2 FILE    write second index reads to FILE\n"
"  --barcode-tag TAG\n"
"               Barcode tag [" DEFAULT_BARCODE_TAG "]\n"
"  --quality-tag TAG\n"
"               Quality tag [" DEFAULT_QUALITY_TAG "]\n"
"  --index-format STR\n"
"               How to parse barcode and quality tags\n\n");
    sam_global_opt_help(to, "-.--.@-.");
    fprintf(to,
"\n"
"The files will be automatically compressed if the file names have a .gz\n"
"or .bgzf extension.  The input to this program must be collated by name.\n"
"Run 'samtools collate' or 'samtools sort -n' to achieve this.\n"
"\n"
"Reads are designated READ1 if FLAG READ1 is set and READ2 is not set.\n"
"Reads are designated READ2 if FLAG READ1 is not set and READ2 is set.\n"
"Otherwise reads are designated READ_OTHER (both flags set or both flags unset).\n"
"Run 'samtools flags' for more information on flag codes and meanings.\n");
    fprintf(to,
"\n"
"The index-format string describes how to parse the barcode and quality tags.\n"
"It is made up of 'i' or 'n' followed by a length or '*'.  For example:\n"
"   i14i8       The first 14 characters are index 1, the next 8 are index 2\n"
"   n8i14       Ignore the first 8 characters, and use the next 14 for index 1\n\n"
"If the tag contains a separator, then the numeric part can be replaced with\n"
"'*' to mean 'read until the separator or end of tag', for example:\n"
"   i*i*        Break the tag at the separator into index 1 and index 2\n"
"   n*i*        Ignore the left part of the tag until the separator,\n"
"               then use the second part of the tag as index 1\n");
    fprintf(to,
"\n"
"Examples:\n"
"To get just the paired reads in separate files, use:\n"
"   samtools %s -1 pair1.%s -2 pair2.%s -0 /dev/null -s /dev/null -n in.bam\n"
"\nTo get all non-supplementary/secondary reads in a single file, redirect\n"
"the output:\n"
"   samtools %s in.bam > all_reads.%s\n",
            command, fq ? "fq" : "fa", fq ? "fq" : "fa",
            command, fq ? "fq" : "fa");
}

typedef enum { READ_UNKNOWN = 0, READ_1 = 1, READ_2 = 2 } readpart;
typedef enum { FASTA, FASTQ } fastfile;
typedef struct bam2fq_opts {
    char *fnse;
    char *fnr[3];
    char *fn_input; // pointer to input filename in argv do not free
    bool has12, has12always, use_oq, copy_tags, illumina_tag;
    int flag_on, flag_off, flag_alloff;
    sam_global_args ga;
    fastfile filetype;
    int def_qual;
    char *barcode_tag;
    char *quality_tag;
    char *index_file[2];
    char *index_format;
    char *extra_tags;
    char compression_level;
} bam2fq_opts_t;

typedef struct bam2fq_state {
    samFile *fp;
    samFile *fpse;
    samFile *fpr[3];
    samFile *fpi[3];
    samFile *hstdout;
    sam_hdr_t *h;
    bool has12, use_oq, copy_tags, illumina_tag;
    int flag_on, flag_off, flag_alloff;
    fastfile filetype;
    int def_qual;
    char *index_sequence;
    char compression_level;
    htsThreadPool p;
} bam2fq_state_t;

static readpart which_readpart(const bam1_t *b)
{
    if ((b->core.flag & BAM_FREAD1) && !(b->core.flag & BAM_FREAD2)) {
        return READ_1;
    } else if ((b->core.flag & BAM_FREAD2) && !(b->core.flag & BAM_FREAD1)) {
        return READ_2;
    } else {
        return READ_UNKNOWN;
    }
}

static void free_opts(bam2fq_opts_t *opts)
{
    free(opts);
}

// return true if valid
static bool parse_opts(int argc, char *argv[], bam2fq_opts_t** opts_out)
{
    // Parse args
    bam2fq_opts_t* opts = calloc(1, sizeof(bam2fq_opts_t));
    opts->has12 = true;
    opts->has12always = false;
    opts->filetype = FASTQ;
    opts->def_qual = 1;
    opts->barcode_tag = NULL;
    opts->quality_tag = NULL;
    opts->index_format = NULL;
    opts->index_file[0] = NULL;
    opts->index_file[1] = NULL;
    opts->extra_tags = NULL;
    opts->compression_level = 1;
    opts->flag_off = BAM_FSECONDARY|BAM_FSUPPLEMENTARY;
    int flag_off_set = 0;

    int c;
    sam_global_args_init(&opts->ga);
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', 0, '@'),
        {"i1", required_argument, NULL, 1},
        {"I1", required_argument, NULL, 1},
        {"i2", required_argument, NULL, 2},
        {"I2", required_argument, NULL, 2},
        {"if", required_argument, NULL, 3},
        {"IF", required_argument, NULL, 3},
        {"index-format", required_argument, NULL, 3},
        {"barcode-tag", required_argument, NULL, 'b'},
        {"quality-tag", required_argument, NULL, 'q'},
        { NULL, 0, NULL, 0 }
    };
    while ((c = getopt_long(argc, argv, "0:1:2:o:f:F:G:niNOs:c:tT:v:@:",
                            lopts, NULL)) > 0) {
        switch (c) {
            case 'b': opts->barcode_tag = optarg; break;
            case 'q': opts->quality_tag = optarg; break;
            case  1 : opts->index_file[0] = optarg; break;
            case  2 : opts->index_file[1] = optarg; break;
            case  3 : opts->index_format = optarg; break;
            case '0': opts->fnr[0] = optarg; break;
            case '1': opts->fnr[1] = optarg; break;
            case '2': opts->fnr[2] = optarg; break;
            case 'o': opts->fnr[1] = optarg; opts->fnr[2] = optarg; break;
            case 'f': opts->flag_on |= strtol(optarg, 0, 0); break;
            case 'F':
                if (!flag_off_set) {
                    flag_off_set = 1;
                    opts->flag_off = 0;
                }
                opts->flag_off |= strtol(optarg, 0, 0);
                break;
            case 'G': opts->flag_alloff |= strtol(optarg, 0, 0); break;
            case 'n': opts->has12 = false; break;
            case 'N': opts->has12always = true; break;
            case 'O': opts->use_oq = true; break;
            case 's': opts->fnse = optarg; break;
            case 't': opts->copy_tags = true; break;
            case 'i': opts->illumina_tag = true; break;
            case 'c':
                opts->compression_level = atoi(optarg);
                if (opts->compression_level < 0)
                    opts->compression_level = 0;
                if (opts->compression_level > 9)
                    opts->compression_level = 9;
                break;
            case 'T': opts->extra_tags = optarg; break;
            case 'v': opts->def_qual = atoi(optarg); break;

            case '?':
                bam2fq_usage(stderr, argv[0]);
                free_opts(opts);
                return false;
            default:
                if (parse_sam_global_opt(c, optarg, lopts, &opts->ga) != 0) {
                    bam2fq_usage(stderr, argv[0]);
                    free_opts(opts);
                    return false;
                }
                break;
        }
    }

    if (opts->fnr[1] || opts->fnr[2]) opts->has12 = false;
    if (opts->has12always) opts->has12 = true;

    if (!opts->barcode_tag) opts->barcode_tag = DEFAULT_BARCODE_TAG;
    if (!opts->quality_tag) opts->quality_tag = DEFAULT_QUALITY_TAG;

    int nIndex = 0;
    if (opts->index_format) {
        char *s;
        for (s = opts->index_format; *s; s++) {
            if (*s == 'i') nIndex++;
        }
    }
    if (nIndex>2) {
        fprintf(stderr,"Invalid index format: more than 2 indexes\n");
        bam2fq_usage(stderr, argv[0]);
        free_opts(opts);
        return false;
    }

    if (opts->index_file[1] && !opts->index_file[0]) {
        fprintf(stderr, "Index one specified, but index two not given\n");
        bam2fq_usage(stderr, argv[0]);
        free_opts(opts);
        return false;
    }

    if (opts->illumina_tag && !nIndex) {
        fprintf(stderr, "You must specify an index format (--index-format) with the Illumina Casava (-i) option\n");
        bam2fq_usage(stderr, argv[0]);
        free_opts(opts);
        return false;
    }

    if (nIndex==0 && opts->index_file[0]) {
        fprintf(stderr, "index_format not specified, but index file given\n");
        bam2fq_usage(stderr, argv[0]);
        free_opts(opts);
        return false;
    }

    if (opts->def_qual < 0 || 93 < opts->def_qual) {
        fprintf(stderr, "Invalid -v default quality %i, allowed range 0 to 93\n", opts->def_qual);
        bam2fq_usage(stderr, argv[0]);
        free_opts(opts);
        return false;
    }

    const char* type_str = argv[0];
    if (strcasecmp("fastq", type_str) == 0 ||
        strcasecmp("bam2fq", type_str) == 0) {
        opts->filetype = FASTQ;
    } else if (strcasecmp("fasta", type_str) == 0) {
        opts->filetype = FASTA;
    } else {
        print_error("bam2fq", "Unrecognised type call \"%s\", this should be impossible... but you managed it!", type_str);
        bam2fq_usage(stderr, argv[0]);
        free_opts(opts);
        return false;
    }

    if (argc == optind && isatty(STDIN_FILENO)) {
        bam2fq_usage(stdout, argv[0]);
        free_opts(opts);
        return true;
    }

    if (argc - optind > 1) {
        fprintf(stderr, "Too many arguments.\n");
        bam2fq_usage(stderr, argv[0]);
        free_opts(opts);
        return false;
    }
    opts->fn_input = argc > optind ? argv[optind] : "-";
    *opts_out = opts;
    return true;
}

void set_sam_opts(samFile *fp, bam2fq_state_t *state,
                  const bam2fq_opts_t *opts) {
    if (state->has12)
        hts_set_opt(fp, FASTQ_OPT_RNUM, 1);

    if (state->illumina_tag)
        hts_set_opt(fp, FASTQ_OPT_CASAVA, 1);

    hts_set_opt(fp, FASTQ_OPT_BARCODE, opts->barcode_tag);

    if (opts->extra_tags && (*opts->extra_tags == '*' || *opts->extra_tags == '\0'))
        hts_set_opt(fp, FASTQ_OPT_AUX, NULL);
    else {
        kstring_t tag_list = {0,0};
        if (state->copy_tags)
            kputs("RG,BC,QT", &tag_list);
        if (opts->extra_tags) {
            if (tag_list.l)
                kputc(',', &tag_list);
            kputs(opts->extra_tags, &tag_list);
        }
        if (tag_list.l)
            hts_set_opt(fp, FASTQ_OPT_AUX, tag_list.s);
        ks_free(&tag_list);
    }
}

// Open a file as normal or gzipped based on filename.
// Note we always use bgzf and don't bother to attempt non-blocked
// gzip streams.  This is a departure from the old fastq code.
static samFile *sam_open_z(char *fn, char *mode, bam2fq_state_t *state) {
    char modez[6];
    strcpy(modez, mode);

    size_t l = strlen(fn);
    if ((l > 3 && strcmp(fn+l-3, ".gz") == 0) ||
        (l > 4 && strcmp(fn+l-4, ".bgz") == 0) ||
        (l > 5 && strcmp(fn+l-5, ".bgzf") == 0)) {
        char m[3] = {'z', state->compression_level+'0', '\0'};
        strcat(modez, m);
    }

    samFile *fp = sam_open(fn, modez);
    if (!fp)
        return NULL;

    if (state->p.pool)
        hts_set_thread_pool(fp, &state->p);

    return fp;
}

static bool init_state(const bam2fq_opts_t* opts, bam2fq_state_t** state_out)
{
    char *mode = opts->filetype == FASTA ? "wF" : "wf";

    bam2fq_state_t* state = calloc(1, sizeof(bam2fq_state_t));
    if (!state)
        return false;
    state->flag_on = opts->flag_on;
    state->flag_off = opts->flag_off;
    state->flag_alloff = opts->flag_alloff;
    state->has12 = opts->has12;
    state->use_oq = opts->use_oq;
    state->illumina_tag = opts->illumina_tag;
    state->copy_tags = opts->copy_tags;
    state->filetype = opts->filetype;
    state->def_qual = opts->def_qual;
    state->index_sequence = NULL;
    state->hstdout = NULL;
    state->compression_level = opts->compression_level;

    state->fp = sam_open(opts->fn_input, "r");
    if (state->fp == NULL) {
        print_error_errno("bam2fq","Cannot read file \"%s\"", opts->fn_input);
        free(state);
        return false;
    }

    state->p.pool = NULL;
    if (opts->ga.nthreads > 0) {
        if (!(state->p.pool = hts_tpool_init(opts->ga.nthreads))) {
            fprintf(stderr, "Failed to create thread pool\n");
            free(state);
            return false;
        }
        state->p.qsize = opts->ga.nthreads*2;
        hts_set_thread_pool(state->fp, &state->p);
    }

    uint32_t rf = SAM_QNAME | SAM_FLAG | SAM_SEQ | SAM_QUAL;
    if (opts->use_oq || opts->extra_tags || opts->index_file[0]) rf |= SAM_AUX;
    if (hts_set_opt(state->fp, CRAM_OPT_REQUIRED_FIELDS, rf)) {
        fprintf(stderr, "Failed to set CRAM_OPT_REQUIRED_FIELDS value\n");
        free(state);
        return false;
    }
    if (hts_set_opt(state->fp, CRAM_OPT_DECODE_MD, 0)) {
        fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
        free(state);
        return false;
    }
    if (opts->fnse) {
        if (!(state->fpse = sam_open_z(opts->fnse, mode, state))) {
            print_error_errno("bam2fq", "Cannot open singleton file \"%s\"", opts->fnse);
            free(state);
            return false;
        }
        set_sam_opts(state->fpse, state, opts);
    }

    if (opts->ga.reference) {
        if (hts_set_fai_filename(state->fp, opts->ga.reference) != 0) {
            print_error_errno("bam2fq", "cannot load reference \"%s\"", opts->ga.reference);
            free(state);
            return false;
        }
    }

    // single, read1, read2
    int i, j;
    for (i = 0; i < 3; ++i) {
        if (opts->fnr[i]) {
            for (j = 0; j < i; j++)
                if (opts->fnr[j] && strcmp(opts->fnr[j], opts->fnr[i]) == 0)
                    break;
            if (j == i) {
                if (!(state->fpr[i] = sam_open_z(opts->fnr[i], mode, state))) {
                    print_error_errno("bam2fq", "Cannot open r%d file \"%s\"",
                                      i, opts->fnr[i]);
                    free(state);
                    return false;
                }
                set_sam_opts(state->fpr[i], state, opts);
            } else {
                state->fpr[i] = state->fpr[j];
            }
        } else {
            if (!state->hstdout) {
                if (!(state->hstdout = sam_open_z("-", mode, state))) {
                    print_error_errno("bam2fq", "Cannot open STDOUT");
                    free(state);
                    return false;
                }
                set_sam_opts(state->hstdout, state, opts);
                autoflush_if_stdout(state->hstdout, "-");
            }
            state->fpr[i] = state->hstdout;
        }
    }

    // index 1, index 2
    for (i = 0; i < 2; i++) {
        state->fpi[i] = NULL;
        if (opts->index_file[i]) {
            for (j = 0; j < 3; j++)
                if (opts->fnr[j] && strcmp(opts->fnr[j], opts->index_file[i]) == 0)
                    break;
            for (j -= 3; j >= 0 && j < i; j++)
                if (opts->index_file[j] && strcmp(opts->index_file[j], opts->index_file[i]) == 0)
                    break;
            if (i == j) {
                if (!(state->fpi[i] = sam_open_z(opts->index_file[i], mode,
                                                 state))) {
                    print_error_errno("bam2fq", "Cannot open i%d file \"%s\"",
                                      i+1, opts->index_file[i]);
                    free(state);
                    return false;
                }
                set_sam_opts(state->fpi[i], state, opts);
            } else if (j < 0) {
                state->fpi[i] = state->fpr[j+3];
            } else {
                state->fpi[i] = state->fpi[j];
            }
        }
    }

    state->h = sam_hdr_read(state->fp);
    if (state->h == NULL) {
        fprintf(stderr, "Failed to read header for \"%s\"\n", opts->fn_input);
        free(state);
        return false;
    }

    *state_out = state;
    return true;
}

static bool destroy_state(const bam2fq_opts_t *opts, bam2fq_state_t *state, int* status)
{
    bool valid = true;
    sam_hdr_destroy(state->h);
    check_sam_close("bam2fq", state->fp, opts->fn_input, "file", status);
    if (state->fpse && sam_close(state->fpse) < 0) {
        print_error_errno("bam2fq", "Error closing singleton file \"%s\"", opts->fnse);
        valid = false;
    }

    int i, j;
    for (i = 0; i < 3; ++i) {
        if (state->fpr[i] != state->hstdout) {
            for (j = 0; j < i; j++)
                if (state->fpr[i] == state->fpr[j])
                    break;
            if (j == i && sam_close(state->fpr[i])) {
                print_error_errno("bam2fq", "Error closing r%d file \"%s\"", i, opts->fnr[i]);
                valid = false;
            }
        }
    }
    if (state->hstdout) {
        release_autoflush(state->hstdout);
        if (sam_close(state->hstdout) < 0) {
            print_error_errno("bam2fq", "Error closing STDOUT");
            valid = false;
        }
    }
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 3; j++)
            if (state->fpi[i] == state->fpr[j])
                break;
        for (j -= 3; j >= 0 && j < i; j++)
            if (state->fpi[i] == state->fpi[j])
                break;
        if (j == i && state->fpi[i] && sam_close(state->fpi[i]) < 0) {
            print_error_errno("bam2fq", "Error closing i%d file \"%s\"", i+1, opts->index_file[i]);
            valid = false;
        }
    }
    free(state->index_sequence);
    if (state->p.pool)
        hts_tpool_destroy(state->p.pool);
    free(state);
    return valid;
}

static inline bool filter_it_out(const bam1_t *b, const bam2fq_state_t *state)
{
    return ((b->core.flag&(state->flag_on)) != state->flag_on // or reads indicated by filter flags
        ||  (b->core.flag&(state->flag_off)) != 0
        ||  (b->core.flag&(state->flag_alloff) && (b->core.flag&(state->flag_alloff)) == state->flag_alloff));

}

int write_index_rec(samFile *fp, bam1_t *b, bam2fq_state_t *state,
                    bam2fq_opts_t* opts, char *seq, int seq_len,
                    char *qual, int qual_len) {
    if (!fp || !b || !seq_len)
        return 0;

    int ret = -1;
    bam1_t *b2 = bam_init1(); // FIXME: reuse
    if (!b2)
        return -1;

    size_t aux_len = b->data + b->l_data - bam_get_aux(b);
    if (bam_set1(b2, b->core.l_qname, bam_get_qname(b),
                 (b->core.flag | BAM_FUNMAP) & ~BAM_FREVERSE,
                 -1, -1, 0,    // refid, pos, mapq
                 0, NULL,      // cigar
                 -1, -1, 0,    // rnext, pnext, tlen
                 seq_len, seq, qual,
                 aux_len) < 0)
        goto err;

    uint8_t *q = bam_get_qual(b2);
    if (qual) {
        int i;
        for (i = 0; i < seq_len; i++)
            q[i] -= '!';
    } else {
        memset(q, opts->def_qual, seq_len);
    }

    memcpy(bam_get_aux(b2), bam_get_aux(b), aux_len);
    b2->l_data += aux_len;
    if (sam_write1(fp, state->h, b2) < 0)
        goto err;

    ret = 0;
 err:
    if (b2)
        bam_destroy1(b2);
    return ret;
}

int output_index(bam1_t *b1, bam1_t *b2, bam2fq_state_t *state,
                 bam2fq_opts_t* opts) {
    bam1_t *b = b1 ? b1 : b2;

    char *ifmt = opts->index_format;
    if (!ifmt)
        ifmt = "i*i*";

    // Get seq / qual elements
    char *bc = NULL, *qt = NULL;
    if (b1)
        bc = (char *)bam_aux_get(b1, opts->barcode_tag);
    if (b2 && !bc)
        bc = (char *)bam_aux_get(b2, opts->barcode_tag);
    if (!bc)
        return 0;
    else
        bc++; // skip Z

    if (b1)
        qt = (char *)bam_aux_get(b1, opts->quality_tag);
    if (b2 && !qt)
        qt = (char *)bam_aux_get(b2, opts->quality_tag);
    if (qt && strlen(bc) != strlen(qt)-1)
        qt = NULL;
    else if (qt)
        qt++;

    int inum = 0;
    while (inum < 2) {
        char fc = *ifmt++;
        if (!fc)
            break; // ran out of index-format

        long len, rem = 0;
        if (isdigit(*ifmt)) {
            rem = len = strtol(ifmt, &ifmt, 10);
        } else {
            ifmt++;
            len = 0;
        }

        char *bc_end = bc, *qt_end = qt;
        while (len ? *bc_end && rem-- : isalpha(*bc_end))
            bc_end++, qt_end += qt != NULL;

        switch (fc) {
        case 'n':
            // skip
            bc = bc_end + (len==0);
            if (qt)
                qt = qt_end + (len==0);
            break;

        case 'i':
            if (write_index_rec(state->fpi[inum], b, state, opts,
                                bc, bc_end-bc, qt, qt_end-qt) < 0)
                return -1;
            bc = bc_end + (len==0);
            if (qt)
                qt = qt_end + (len==0);
            inum++;
            break;

        default:
            fprintf(stderr, "Unknown index-format code\n");
            return -1;
        }
    }

    return 0;
}

static int flush_rec(bam2fq_state_t *state, bam2fq_opts_t* opts,
                     bam1_t *b[4], int score[3], int best[3],
                     int64_t *n_singletons) {
    // Paired data, with 1 or 2 ends present.
    if (score[1] > 0 && score[2] > 0) {
        // If CASAVA tag is required and barcode is only on R1,
        // copy it to R2
        if (state->illumina_tag) {
            char *tag;
            if ((tag = (char *)bam_aux_get(b[best[1]],
                                           opts->barcode_tag)))
                if (bam_aux_update_str(b[best[2]],
                                       opts->barcode_tag,
                                       strlen(tag), tag+1) < 0)
                    goto err;
            if ((tag = (char *)bam_aux_get(b[best[1]],
                                           opts->quality_tag)))
                if (bam_aux_update_str(b[best[2]],
                                       opts->quality_tag,
                                       strlen(tag), tag+1) < 0)
                    goto err;

        }
        if (sam_write1(state->fpr[1], state->h, b[best[1]]) < 0)
            goto err;
        if (sam_write1(state->fpr[2], state->h, b[best[2]]) < 0)
            goto err;

        if (output_index(b[best[1]], b[best[2]], state, opts) < 0)
            goto err;
    } else if (score[1] > 0 || score[2] > 0) {
        if (state->fpse) {
            // print whichever one exists to fpse
            if (score[1] > 0) {
                if (sam_write1(state->fpse, state->h, b[best[1]]) < 0)
                    goto err;
            } else {
                if (sam_write1(state->fpse, state->h, b[best[2]]) < 0)
                    goto err;
            }
            ++(*n_singletons);
        } else {
            if (score[1] > 0) {
                if (sam_write1(state->fpr[1], state->h, b[best[1]]) < 0)
                    goto err;
            } else {
                if (sam_write1(state->fpr[2], state->h, b[best[2]]) < 0)
                    goto err;
            }
        }

        if (output_index(score[1] > 0 ? b[best[1]] : NULL,
                         score[2] > 0 ? b[best[2]] : NULL,
                         state, opts) < 0)
            goto err;
    }

    if (score[0]) { // single ended data (neither READ1 nor READ2)
        if (sam_write1(state->fpr[0], state->h, b[best[0]]) < 0)
            goto err;

        if (output_index(b[best[0]], NULL, state, opts) < 0)
            goto err;
    }

    return 0;

 err:
    return -1;
}

static bool bam2fq_mainloop(bam2fq_state_t *state, bam2fq_opts_t* opts)
{
    int n;
    char *current_qname = NULL;
    int64_t n_reads = 0, n_singletons = 0; // Statistics
    int score[3];
    int at_eof;
    bool valid = false;
    int best[3] = {-1, -1, -1}; // map R0, R1, single to b[] indices;
                                // indexed by [readpart]
    bam1_t *b[4];               // 3 readparts, plus current record

    for (n = 0; n < 4; n++) {
        if (!(b[n] = bam_init1())) {
            perror("[bam2fq_mainloop] Malloc error for bam record buffer.");
            return false;
        }
    }

    n = 0;
    while (true) {
        int res = sam_read1(state->fp, state->h, b[n]);
        if (res < -1) {
            print_error("bam2fq", "Failed to read bam record");
            goto err;
        }
        at_eof = res < 0;

        if (!at_eof && filter_it_out(b[n], state))
            continue;
        if (!at_eof) {
            ++n_reads;

            // Handle -O option: use OQ for qual
            uint8_t *oq;
            if (state->use_oq && (oq = bam_aux_get(b[n],"OQ")) && *oq == 'Z') {
                int i, l = strlen((char *)++oq);
                uint8_t *qual = bam_get_qual(b[n]);
                for (i = 0; i < l && i < b[n]->core.l_qseq; i++)
                    qual[i] = oq[i] - '!';
            }
        }

        if (at_eof
            || !current_qname
            || (strcmp(current_qname, bam_get_qname(b[n])) != 0)) {
            // New name, so flush best examples of previous name.
            if (current_qname)
                if (flush_rec(state, opts, b, score, best, &n_singletons) < 0)
                    goto err;

            current_qname = bam_get_qname(b[n]);
            score[0] = score[1] = score[2] = 0;

            if (at_eof) { break; }
        }

        // Prefer a copy of the read that has base qualities
        int b_score = bam_get_qual(b[n])[0] != 0xff? 2 : 1;
        readpart rp = which_readpart(b[n]);
        if (score[rp] < b_score) {
            score[rp] = b_score;
            // Record b[n] slot for best copy of readpair and find a new
            // slot for next bam read
            best[rp] = n;
            int used_slot[4] = {0}, i;
            for (i = 0; i < 3; i++)
                if (best[i] >= 0)
                    used_slot[best[i]] = 1;
            for (i = 0; i < 4 && used_slot[i]; i++)
                ;
            n = i;
        }
    }

    valid = true;
 err:
    if (!valid)
        print_error_errno("bam2fq", "Error writing to FASTx files.");

    for (n = 0; n < 4; n++)
        bam_destroy1(b[n]);

    fprintf(stderr, "[M::%s] discarded %" PRId64 " singletons\n",
            __func__, n_singletons);
    fprintf(stderr, "[M::%s] processed %" PRId64 " reads\n",
            __func__, n_reads);

    return valid;
}

int main_bam2fq(int argc, char *argv[])
{
    int status = EXIT_FAILURE;
    bam2fq_opts_t* opts = NULL;
    bam2fq_state_t* state = NULL;

    bool valid = parse_opts(argc, argv, &opts);
    if (!valid || opts == NULL) return valid ? EXIT_SUCCESS : EXIT_FAILURE;

    if (!init_state(opts, &state)) goto err;

    if (!bam2fq_mainloop(state,opts)) goto err;

    if (!destroy_state(opts, state, &status)) goto err;

    status = EXIT_SUCCESS;
 err:
    sam_global_args_free(&opts->ga);
    free_opts(opts);

    return status;
}
