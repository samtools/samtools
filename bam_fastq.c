/*  bam_fastq.c -- FASTA and FASTQ file generation

    Copyright (C) 2009-2017, 2019 Genome Research Ltd.
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

#define taglist_free(p)
KLIST_INIT(ktaglist, char*, taglist_free)

#define DEFAULT_BARCODE_TAG "BC"
#define DEFAULT_QUALITY_TAG "QT"
#define INDEX_SEPARATOR "+"

int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };
static const char *copied_tags[] = { "RG", "BC", "QT", NULL };

static void bam2fq_usage(FILE *to, const char *command)
{
    int fq = strcasecmp("fastq", command) == 0 || strcasecmp("bam2fq", command) == 0;
    fprintf(to,
"Usage: samtools %s [options...] <in.bam>\n", command);
    fprintf(to,
"\n"
"Description:\n"
"Converts a SAM, BAM or CRAM into either FASTQ or FASTA format depending on the command invoked.\n"
"\n"
"Options:\n"
"  -0 FILE              write reads designated READ_OTHER to FILE\n"
"  -1 FILE              write reads designated READ1 to FILE\n"
"  -2 FILE              write reads designated READ2 to FILE\n"
"  -o FILE              write reads designated READ1 or READ2 to FILE\n"
"                       note: if a singleton file is specified with -s, only\n"
"                       paired reads will be written to the -1 and -2 files.\n"
"  -f INT               only include reads with all  of the FLAGs in INT present [0]\n"       //   F&x == x
"  -F INT               only include reads with none of the FLAGS in INT present [0x900]\n"       //   F&x == 0
"  -G INT               only EXCLUDE reads with all  of the FLAGs in INT present [0]\n"       // !(F&x == x)
"  -n                   don't append /1 and /2 to the read name\n"
"  -N                   always append /1 and /2 to the read name\n");
    if (fq) fprintf(to,
"  -O                   output quality in the OQ tag if present\n");
    fprintf(to,
"  -s FILE              write singleton reads designated READ1 or READ2 to FILE\n"
"  -t                   copy RG, BC and QT tags to the %s header line\n",
    fq ? "FASTQ" : "FASTA");
    fprintf(to,
"  -T TAGLIST           copy arbitrary tags to the %s header line\n",
    fq ? "FASTQ" : "FASTA");
    if (fq) fprintf(to,
"  -v INT               default quality score if not given in file [1]\n"
"  -i                   add Illumina Casava 1.8 format entry to header (eg 1:N:0:ATCACG)\n"
"  -c                   compression level [0..9] to use when creating gz or bgzf fastq files [1]\n"
"  --i1 FILE            write first index reads to FILE\n"
"  --i2 FILE            write second index reads to FILE\n"
"  --barcode-tag TAG    Barcode tag [default: " DEFAULT_BARCODE_TAG "]\n"
"  --quality-tag TAG    Quality tag [default: " DEFAULT_QUALITY_TAG "]\n"
"  --index-format STR   How to parse barcode and quality tags\n\n");
    sam_global_opt_help(to, "-.--.@-.");
    fprintf(to,
"\n"
"The files will be automatically compressed if the file names have a .gz or .bgzf extension.\n"
"The input to this program must be collated by name. Run 'samtools collate' or 'samtools sort -n'.\n"
"\n"
"Reads are designated READ1 if FLAG READ1 is set and READ2 is not set.\n"
"Reads are designated READ2 if FLAG READ1 is not set and READ2 is set.\n"
"Reads are designated READ_OTHER if FLAGs READ1 and READ2 are either both set\n"
"or both unset.\n"
"Run 'samtools flags' for more information on flag codes and meanings.\n");
    fprintf(to,
"\n"
"The index-format string describes how to parse the barcode and quality tags, for example:\n"
"   i14i8       the first 14 characters are index 1, the next 8 characters are index 2\n"
"   n8i14       ignore the first 8 characters, and use the next 14 characters for index 1\n"
"If the tag contains a separator, then the numeric part can be replaced with '*' to mean\n"
"'read until the separator or end of tag', for example:\n"
"   n*i*        ignore the left part of the tag until the separator, then use the second part\n"
"               of the tag as index 1\n");
    fprintf(to,
"\n"
"Examples:\n"
" To get just the paired reads in separate files, use:\n"
"   samtools %s -1 paired1.%s -2 paired2.%s -0 /dev/null -s /dev/null -n in.bam\n"
"\n To get all non-supplementary/secondary reads in a single file, redirect the output:\n"
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
    BGZF *fpse;
    BGZF *fpr[3];
    BGZF *fpi[2];
    BGZF *hstdout;
    sam_hdr_t *h;
    bool has12, use_oq, copy_tags, illumina_tag;
    int flag_on, flag_off, flag_alloff;
    fastfile filetype;
    int def_qual;
    klist_t(ktaglist) *taglist;
    char *index_sequence;
    char compression_level;
    htsThreadPool p;
} bam2fq_state_t;

/*
 * Get and decode the read from a BAM record.
 *
 * TODO: htslib really needs an interface for this.  Consider this or perhaps
 * bam_get_seq_str (current vs original orientation) and bam_get_qual_str
 * functions as string formatted equivalents to bam_get_{seq,qual}?
 */

/*
 * Reverse a string in place.
 * From http://stackoverflow.com/questions/8534274/is-the-strrev-function-not-available-in-linux.
 * Author Sumit-naik: http://stackoverflow.com/users/4590926/sumit-naik
 */
static char *reverse(char *str)
{
    int i = strlen(str)-1,j=0;
    char ch;
    while (i>j) {
        ch = str[i];
        str[i]= str[j];
        str[j] = ch;
        i--;
        j++;
    }
    return str;
}

/* return the read, reverse complemented if necessary */
static char *get_read(const bam1_t *rec)
{
    int len = rec->core.l_qseq + 1;
    char *read = calloc(1, len);
    char *seq = (char *)bam_get_seq(rec);
    int n;

    if (!read) return NULL;

    for (n=0; n < rec->core.l_qseq; n++) {
        if (rec->core.flag & BAM_FREVERSE) read[n] = seq_nt16_str[seq_comp_table[bam_seqi(seq,n)]];
        else                               read[n] = seq_nt16_str[bam_seqi(seq,n)];
    }
    if (rec->core.flag & BAM_FREVERSE) reverse(read);
    return read;
}

/*
 * get and decode the quality from a BAM record
 */
static int get_quality(const bam1_t *rec, char **qual_out)
{
    char *quality = calloc(1, rec->core.l_qseq + 1);
    char *q = (char *)bam_get_qual(rec);
    int n;

    if (!quality) return -1;

    if (*q == '\xff') {
        free(quality);
        *qual_out = NULL;
        return 0;
    }

    for (n=0; n < rec->core.l_qseq; n++) {
        quality[n] = q[n]+33;
    }
    if (rec->core.flag & BAM_FREVERSE) reverse(quality);
    *qual_out = quality;
    return 0;
}

//
// End of htslib complaints
//


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

/*
 * parse the length part from the index-format string
 */
static int getLength(char **s)
{
    int n = 0;
    while (**s) {
        if (**s == '*') { n=-1; (*s)++; break; }
        if ( !isdigit(**s)) break;
        n = n*10 + ((**s)-'0');
        (*s)++;
    }
    return n;
}

static bool copy_tag(const char *tag, const bam1_t *rec, kstring_t *linebuf)
{
    uint8_t *s = bam_aux_get(rec, tag);
    if (s) {
        char aux_type = *s;
        switch (aux_type) {
            case 'C':
            case 'S': aux_type = 'I'; break;
            case 'c':
            case 's': aux_type = 'i'; break;
            case 'd': aux_type = 'f'; break;
        }

        // Ensure space.  Need 6 chars + length of tag.  Max length of
        // i is 16, A is 21, B currently 26, Z is unknown, so
        // have to check that one later.
        if (ks_resize(linebuf, ks_len(linebuf) + 64) < 0) return false;

        kputc('\t', linebuf);
        kputsn(tag, 2, linebuf);
        kputc(':', linebuf);
        kputc(aux_type=='I'? 'i': aux_type, linebuf);
        kputc(':', linebuf);
        switch (aux_type) {
            case 'H':
            case 'Z':
                if (kputs(bam_aux2Z(s), linebuf) < 0) return false;
                break;
            case 'i': kputw(bam_aux2i(s), linebuf); break;
            case 'I': kputuw(bam_aux2i(s), linebuf); break;
            case 'A': kputc(bam_aux2A(s), linebuf); break;
            case 'f': kputd(bam_aux2f(s), linebuf); break;
            case 'B': kputs("*** Unhandled aux type ***", linebuf); return false;
            default:  kputs("*** Unknown aux type ***", linebuf); return false;
       }
    }
    return true;
}

static int insert_index_sequence_into_linebuf(char *index_sequence, kstring_t *linebuf, bam1_t *rec)
{
    if (!index_sequence) return 0;

    kstring_t new = {0,0,NULL};
    if (linebuf->s) {
        char *s = strchr(linebuf->s, '\n');
        if (s) {
            if (ks_resize(&new, linebuf->l + strlen(index_sequence) + 16) < 0)
                return -1;
            *s = 0;
            kputs(linebuf->s, &new);
            kputc(' ', &new);
            readpart readpart = which_readpart(rec);
            if (readpart == READ_1) kputc('1', &new);
            else if (readpart == READ_2) kputc('2', &new);
            else kputc('0', &new);

            kputc(':', &new);
            if (rec->core.flag & BAM_FQCFAIL) kputc('Y', &new);
            else                              kputc('N', &new);

            kputs(":0:", &new);
            kputs(index_sequence, &new);
            kputc('\n', &new);
            kputs(s+1, &new);
            free(ks_release(linebuf));
            linebuf->s = new.s; linebuf->l = new.l; linebuf->m = new.m;
        }
    }
    return 0;
}

static bool make_fq_line(const bam1_t *rec, char *seq, char *qual, kstring_t *linebuf, const bam2fq_state_t *state)
{
    int i;

    linebuf->l = 0;
    // Write read name
    if (kputc(state->filetype == FASTA? '>' : '@', linebuf) < 0) return false;
    if (kputs(bam_get_qname(rec), linebuf) < 0) return false;
    // Add the /1 /2 if requested
    if (state->has12) {
        readpart readpart = which_readpart(rec);
        if (readpart == READ_1) {
            if (kputs("/1", linebuf) < 0) return false;
        } else if (readpart == READ_2) {
            if (kputs("/2", linebuf) < 0) return false;
        }
    }
    if (state->copy_tags) {
        for (i = 0; copied_tags[i]; ++i) {
            if (!copy_tag(copied_tags[i], rec, linebuf)) {
                fprintf(stderr, "Problem copying aux tags: [%s]\n", linebuf->s);
                return false;
            }
        }
    }

    if (state->taglist->size) {
        kliter_t(ktaglist) *p;
        for (p = kl_begin(state->taglist); p != kl_end(state->taglist); p = kl_next(p)) {
            if (!copy_tag(kl_val(p), rec, linebuf)) {
                fprintf(stderr, "Problem copying aux tags: [%s]\n", linebuf->s);
                return false;
            }
        }
    }

    if (kputc('\n', linebuf) < 0) return false;
    if (kputs(seq, linebuf) < 0) return false;
    if (kputc('\n', linebuf) < 0) return false;

    if (state->filetype == FASTQ) {
        // Write quality
        if (kputs("+\n", linebuf) < 0) return false;
        if (qual && *qual) {
            if (kputs(qual, linebuf) < 0) return false;
        } else {
            int len = strlen(seq);
            if (ks_resize(linebuf, ks_len(linebuf) + len + 1) < 0) return false;
            for (i = 0; i < len; ++i) {
                kputc(33 + state->def_qual, linebuf);
            }
        }
        if (kputc('\n', linebuf) < 0) return false;
    }
    return true;
}

/*
 * Create FASTQ lines from the barcode tag using the index-format
 */
static bool tags2fq(bam1_t *rec, bam2fq_state_t *state, const bam2fq_opts_t* opts)
{
    uint8_t *p;
    char *ifmt = opts->index_format;
    char *tag = NULL;
    char *qual = NULL;
    char *sub_tag = NULL;
    char *sub_qual = NULL;
    size_t tag_len;
    int file_number = 0;
    kstring_t linebuf = { 0, 0, NULL }; // Buffer

    if (!ifmt) return true;

    // read barcode tag
    p = bam_aux_get(rec,opts->barcode_tag);
    if (p) tag = bam_aux2Z(p);

    if (!tag) return true; // there is no tag

    tag_len = strlen(tag);
    sub_tag = calloc(1, tag_len + 1);
    if (!sub_tag) goto fail;
    sub_qual = calloc(1, tag_len + 1);
    if (!sub_qual) goto fail;

    // read quality tag
    p = bam_aux_get(rec, opts->quality_tag);
    if (p) qual = bam_aux2Z(p);

    // Parse the index-format string
    while (*ifmt) {
        if (file_number > 1) break;     // shouldn't happen if we've validated paramaters correctly
        char action = *ifmt;        // should be 'i' or 'n'
        ifmt++; // skip over action
        int index_len = getLength(&ifmt);
        int n = 0;

        if (index_len < 0) {
            // read until separator
            while (isalpha(*tag)) {
                sub_tag[n] = *tag++;
                if (qual) sub_qual[n] = *qual++;
                n++;
            }
            if (*tag) { // skip separator
                tag++;
                if (qual) qual++;
            }
        } else {
            // read index_len characters
            while (index_len-- && *tag) {
                sub_tag[n] = *tag++;
                if (qual) sub_qual[n] = *qual++;
                n++;
            }
        }
        sub_tag[n] = '\0';
        sub_qual[n] = '\0';

        if (action=='i' && *sub_tag) {
            if (state->index_sequence) {
                char *new_index_sequence = realloc(state->index_sequence, strlen(state->index_sequence) + strlen(sub_tag) + 2);
                if (!new_index_sequence) goto fail;
                state->index_sequence = new_index_sequence;
                strcat(state->index_sequence, INDEX_SEPARATOR);
                strcat(state->index_sequence, sub_tag);
            } else {
                state->index_sequence = strdup(sub_tag);    // we're going to need this later...
            }
            if (!state->index_sequence) goto fail;
            if (!make_fq_line(rec, sub_tag, sub_qual, &linebuf, state)) goto fail;
            if (state->illumina_tag) {
                if (insert_index_sequence_into_linebuf(sub_tag, &linebuf, rec) < 0) {
                    goto fail;
                }
            }
            if (state->fpi[file_number]) {
                if (bgzf_write(state->fpi[file_number++], linebuf.s, linebuf.l) < 0)
                    goto fail;
            }
        }

    }

    free(sub_qual); free(sub_tag);
    free(linebuf.s);
    return true;

 fail:
    perror(__func__);
    free(sub_qual); free(sub_tag);
    free(linebuf.s);
    return false;
}

// Transform a bam1_t record into a string with the FASTQ representation of it
// @returns false for error, true for success
static bool bam1_to_fq(const bam1_t *b, kstring_t *linebuf, const bam2fq_state_t *state)
{
    int32_t qlen = b->core.l_qseq;
    assert(qlen >= 0);
    const uint8_t *oq = NULL;
    char *qual = NULL;

    char *seq = get_read(b);
    if (!seq) return false;

    if (state->use_oq) oq = bam_aux_get(b, "OQ");
    if (oq && *oq=='Z') {
        qual = strdup(bam_aux2Z(oq));
        if (!qual) goto fail;
        if (b->core.flag & BAM_FREVERSE) { // read is reverse complemented
            reverse(qual);
        }
    } else {
        if (get_quality(b, &qual) < 0) goto fail;
    }

    if (!make_fq_line(b, seq, qual, linebuf, state)) goto fail;

    free(qual);
    free(seq);
    return true;

 fail:
    free(seq);
    free(qual);
    return false;
}

static void free_opts(bam2fq_opts_t *opts)
{
    free(opts->barcode_tag);
    free(opts->quality_tag);
    free(opts->index_format);
    free(opts->extra_tags);
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
    while ((c = getopt_long(argc, argv, "0:1:2:o:f:F:G:niNOs:c:tT:v:@:", lopts, NULL)) > 0) {
        switch (c) {
            case 'b': opts->barcode_tag = strdup(optarg); break;
            case 'q': opts->quality_tag = strdup(optarg); break;
            case  1 : opts->index_file[0] = optarg; break;
            case  2 : opts->index_file[1] = optarg; break;
            case  3 : opts->index_format = strdup(optarg); break;
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
                opts->flag_off |= strtol(optarg, 0, 0); break;
            case 'G': opts->flag_alloff |= strtol(optarg, 0, 0); break;
            case 'n': opts->has12 = false; break;
            case 'N': opts->has12always = true; break;
            case 'O': opts->use_oq = true; break;
            case 's': opts->fnse = optarg; break;
            case 't': opts->copy_tags = true; break;
            case 'i': opts->illumina_tag = true; break;
            case 'c': opts->compression_level = atoi(optarg); break;
            case 'T': opts->extra_tags = strdup(optarg); break;
            case 'v': opts->def_qual = atoi(optarg); break;
            case '?': bam2fq_usage(stderr, argv[0]); free_opts(opts); return false;
            default:
                if (parse_sam_global_opt(c, optarg, lopts, &opts->ga) != 0) {
                    bam2fq_usage(stderr, argv[0]); free_opts(opts); return false;
                }
                break;
        }
    }

    if (opts->fnr[1] || opts->fnr[2]) opts->has12 = false;
    if (opts->has12always) opts->has12 = true;

    if (!opts->barcode_tag) opts->barcode_tag = strdup(DEFAULT_BARCODE_TAG);
    if (!opts->quality_tag) opts->quality_tag = strdup(DEFAULT_QUALITY_TAG);

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
    if (strcasecmp("fastq", type_str) == 0 || strcasecmp("bam2fq", type_str) == 0) {
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

static BGZF *open_fqfile(char *filename, int c, htsThreadPool *tp)
{
    char mode[4] = "w";
    size_t len = strlen(filename);

    mode[2] = 0; mode[3] = 0;
    if (len > 3 && strstr(filename + (len - 3),".gz")) {
        mode[1] = 'g'; mode[2] = c+'0';
    } else if ((len > 4 && strstr(filename + (len - 4),".bgz"))
               || (len > 5 && strstr(filename + (len - 5),".bgzf"))) {
        mode[1] = c+'0';
    } else {
        mode[1] = 'u';
    }

    BGZF *fp = bgzf_open(filename,mode);
    if (!fp)
        return fp;
    if (tp->pool && bgzf_thread_pool(fp, tp->pool, tp->qsize) < 0) {
        bgzf_close(fp);
        return NULL;
    }
    return fp;
}

static bool init_state(const bam2fq_opts_t* opts, bam2fq_state_t** state_out)
{
    bam2fq_state_t* state = calloc(1, sizeof(bam2fq_state_t));
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

    state->taglist = kl_init(ktaglist);
    if (opts->extra_tags) {
        char *save_p;
        char *s = strtok_r(opts->extra_tags, ",", &save_p);
        while (s) {
            if (strlen(s) != 2) {
                fprintf(stderr, "Parsing extra tags - '%s' is not two characters\n", s);
                free(state);
                return false;
            }
            char **et = kl_pushp(ktaglist, state->taglist);
            *et = s;
            s = strtok_r(NULL, ",", &save_p);
        }
    }

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
        state->fpse = open_fqfile(opts->fnse, state->compression_level, &state->p);
        if (state->fpse == NULL) {
            print_error_errno("bam2fq", "Cannot write to singleton file \"%s\"", opts->fnse);
            free(state);
            return false;
        }
    }

    if (opts->ga.reference) {
        if (hts_set_fai_filename(state->fp, opts->ga.reference) != 0) {
            print_error_errno("bam2fq", "cannot load reference \"%s\"", opts->ga.reference);
            free(state);
            return false;
        }
    }

    int i, j;
    for (i = 0; i < 3; ++i) {
        if (opts->fnr[i]) {
            for (j = 0; j < i; j++)
                if (opts->fnr[j] && strcmp(opts->fnr[j], opts->fnr[i]) == 0)
                    break;
            if (j == i) {
                state->fpr[i] = open_fqfile(opts->fnr[i], state->compression_level, &state->p);
                if (state->fpr[i] == NULL) {
                    print_error_errno("bam2fq", "Cannot write to r%d file \"%s\"",
                                      i, opts->fnr[i]);
                    free(state);
                    return false;
                }
            } else {
                state->fpr[i] = state->fpr[j];
            }
        } else {
            if (!state->hstdout) {
                state->hstdout = bgzf_dopen(fileno(stdout), "wu");
                if (!state->hstdout) {
                    print_error_errno("bam2fq", "Cannot open STDOUT");
                    free(state);
                    return false;
                }
            }
            state->fpr[i] = state->hstdout;
        }
    }
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
                state->fpi[i] = open_fqfile(opts->index_file[i], state->compression_level, &state->p);
                if (state->fpi[i] == NULL) {
                    print_error_errno("bam2fq", "Cannot write to i%d file \"%s\"",
                                      i+1, opts->index_file[i]);
                    free(state);
                    return false;
                }
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
    if (state->fpse && bgzf_close(state->fpse)) { print_error_errno("bam2fq", "Error closing singleton file \"%s\"", opts->fnse); valid = false; }
    int i, j;
    for (i = 0; i < 3; ++i) {
        if (state->fpr[i] != state->hstdout) {
            for (j = 0; j < i; j++)
                if (state->fpr[i] == state->fpr[j])
                    break;
            if (j == i && bgzf_close(state->fpr[i])) {
                print_error_errno("bam2fq", "Error closing r%d file \"%s\"", i, opts->fnr[i]);
                valid = false;
            }
        }
    }
    if (state->hstdout) {
        if (bgzf_close(state->hstdout)) {
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
        if (j == i && state->fpi[i] && bgzf_close(state->fpi[i])) {
            print_error_errno("bam2fq", "Error closing i%d file \"%s\"", i+1, opts->index_file[i]);
            valid = false;
        }
    }
    kl_destroy(ktaglist,state->taglist);
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

static bool bam2fq_mainloop(bam2fq_state_t *state, bam2fq_opts_t* opts)
{
    int n;
    bam1_t *records[3] = {NULL, NULL, NULL};
    char *current_qname = NULL;
    int64_t n_reads = 0, n_singletons = 0; // Statistics
    kstring_t linebuf[3] = {{0,0,NULL},{0,0,NULL},{0,0,NULL}};
    int score[3];
    int at_eof;
    bool valid = true;
    bam1_t* b = NULL;

    while (true) {
        if (!b)
            b = bam_init1();
        if (b == NULL) {
            perror("[bam2fq_mainloop] Malloc error for bam record buffer.");
            valid = false;
            break;
        }
        int res = sam_read1(state->fp, state->h, b);
        if (res < -1) {
            fprintf(stderr, "[bam2fq_mainloop] Failed to read bam record.\n");
            valid = false;
            break;
        }
        at_eof = res < 0;

        if (!at_eof && filter_it_out(b, state))
            continue;
        if (!at_eof) ++n_reads;

        if (at_eof || !current_qname || (strcmp(current_qname, bam_get_qname(b)) != 0)) {
            if (current_qname) {
                if (state->illumina_tag) {
                    for (n=0; valid && n<3; n++) {
                        if (!records[n]) continue;
                        if (insert_index_sequence_into_linebuf(state->index_sequence, &linebuf[n], records[n]) < 0) valid = false;
                    }
                    if (!valid) break;
                }
                free(state->index_sequence); state->index_sequence = NULL;
                if (score[1] > 0 && score[2] > 0) {
                    // print linebuf[1] to fpr[1], linebuf[2] to fpr[2]
                    if (bgzf_write(state->fpr[1], linebuf[1].s, linebuf[1].l) < 0) { valid = false; break; }
                    if (bgzf_write(state->fpr[2], linebuf[2].s, linebuf[2].l) < 0) { valid = false; break; }
                } else if (score[1] > 0 || score[2] > 0) {
                    if (state->fpse) {
                        // print whichever one exists to fpse
                        if (score[1] > 0) {
                            if (bgzf_write(state->fpse, linebuf[1].s, linebuf[1].l) < 0) { valid = false; break; }
                        } else {
                            if (bgzf_write(state->fpse, linebuf[2].s, linebuf[2].l) < 0) { valid = false; break; }
                        }
                        ++n_singletons;
                    } else {
                        if (score[1] > 0) {
                            if (bgzf_write(state->fpr[1], linebuf[1].s, linebuf[1].l) < 0) { valid = false; break; }
                        } else {
                            if (bgzf_write(state->fpr[2], linebuf[2].s, linebuf[2].l) < 0) { valid = false; break; }
                        }
                    }
                }
                if (score[0]) { // TODO: check this
                    // print linebuf[0] to fpr[0]
                    if (bgzf_write(state->fpr[0], linebuf[0].s, linebuf[0].l) < 0) { valid = false; break; }
                }
            }


            free(current_qname); current_qname = NULL;
            score[0] = score[1] = score[2] = 0;
            for (n=0; n < 3; n++) {
                bam_destroy1(records[n]); records[n]=NULL;
            }

            if (at_eof) { break; }

            current_qname = strdup(bam_get_qname(b));
            if (!current_qname) { valid = false; break; }
        }

        // Prefer a copy of the read that has base qualities
        int b_score = bam_get_qual(b)[0] != 0xff? 2 : 1;
        readpart rp = which_readpart(b);
        if (b_score > score[rp]) {
            if (!tags2fq(b, state, opts)) { valid = false; break; }
            if (records[rp]) bam_destroy1(records[rp]);
            records[rp] = b;
            score[rp] = b_score;
            b = NULL;
            if(!bam1_to_fq(records[rp], &linebuf[rp], state)) {
                fprintf(stderr, "[%s] Error converting read to FASTA/Q\n", __func__);
                valid = false; break;
            }
        }
    }
    if (!valid)
    {
        perror("[bam2fq_mainloop] Error writing to FASTx files.");
    }
    bam_destroy1(b);
    for (n=0; n < 3; n++) {
        bam_destroy1(records[n]);
    }
    free(current_qname);
    free(linebuf[0].s);
    free(linebuf[1].s);
    free(linebuf[2].s);
    fprintf(stderr, "[M::%s] discarded %" PRId64 " singletons\n", __func__, n_singletons);
    fprintf(stderr, "[M::%s] processed %" PRId64 " reads\n", __func__, n_reads);

    return valid;
}

int main_bam2fq(int argc, char *argv[])
{
    int status = EXIT_SUCCESS;
    bam2fq_opts_t* opts = NULL;
    bam2fq_state_t* state = NULL;

    bool valid = parse_opts(argc, argv, &opts);
    if (!valid || opts == NULL) return valid ? EXIT_SUCCESS : EXIT_FAILURE;

    if (!init_state(opts, &state)) return EXIT_FAILURE;

    if (!bam2fq_mainloop(state,opts)) status = EXIT_FAILURE;

    if (!destroy_state(opts, state, &status)) return EXIT_FAILURE;
    sam_global_args_free(&opts->ga);
    free_opts(opts);

    return status;
}
