/*  faidx.c -- faidx subcommand.

    Copyright (C) 2008, 2009, 2013, 2016, 2018-2020, 2022, 2024 Genome Research Ltd.
    Portions copyright (C) 2011 Broad Institute.

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
DEALINGS IN THE SOFTWARE.

History:

  * 2016-01-12: Pierre Lindenbaum @yokofakun : added options -o -n

*/

#include <config.h>

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdarg.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/hfile.h>
#include <htslib/kstring.h>
#include <htslib/bgzf.h>
#include <htslib/thread_pool.h>
#include "sam_opts.h"
#include "samtools.h"

// Negative indicates the same as input data
#define DEFAULT_FASTA_LINE_LEN -60

#ifndef ABS
#   define ABS(x) ((x)>=0?(x):-(x))
#endif

//new params required for output creation
typedef struct output {
    int isbgzip;                //is bgzip or uncompressed file
    FILE *fp;                   //uncompressed file pointer
    BGZF *bgzf_fp;              //bgzf file pointer
    sam_global_args *gopt;      //options
    kstring_t buffer;
} output;

static unsigned char comp_base[256] = {
  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
 16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
 32, '!', '"', '#', '$', '%', '&', '\'','(', ')', '*', '+', ',', '-', '.', '/',
'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?',
'@', 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z', '[', '\\',']', '^', '_',
'`', 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', '{', '|', '}', '~', 127,
128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255,
};

static void reverse_complement(char *str, const hts_pos_t len) {
    char c;
    hts_pos_t i = 0, j = len - 1;

    while (i <= j) {
        c = str[i];
        str[i] = comp_base[(unsigned char)str[j]];
        str[j] = comp_base[(unsigned char)c];
        i++;
        j--;
    }
}

static void reverse(char *str, const hts_pos_t len) {
    char c;
    hts_pos_t i = 0, j = len - 1;

    while (i < j) {
        c = str[i];
        str[i] = str[j];
        str[j] = c;
        i++;
        j--;
    }
}

/// wrappedwrite  - wraps the fwrite and bgzf_write
/** @param out    - pointer to data required to write output
*   @param buffer - data to write
*   @param length - data length
* returns error or length written on success
*/
static inline size_t wrappedwrite(output *out, const char *buffer, size_t length)
{
    return out->isbgzip ? bgzf_write(out->bgzf_fp, buffer, length) :
        fwrite(buffer, 1, length, out->fp);
}

static int write_line(faidx_t *faid, output *out, const char *line, const char *name,
                      const int ignore, const hts_pos_t length, const hts_pos_t seq_len) {
    int id;
    hts_pos_t beg, end;

    if (seq_len < 0) {
        fprintf(stderr, "[faidx] Failed to fetch sequence in %s\n", name);

        if (ignore && seq_len == -2) {
            return EXIT_SUCCESS;
        } else {
            return EXIT_FAILURE;
        }
    } else if (seq_len == 0) {
        fprintf(stderr, "[faidx] Zero length sequence: %s\n", name);
    } else if (fai_parse_region(faid, name, &id, &beg, &end, 0)
               && (end < HTS_POS_MAX) && (seq_len != end - beg)) {
        fprintf(stderr, "[faidx] Truncated sequence: %s\n", name);
    }

    hts_pos_t i, seq_sz = seq_len;

    for (i = 0; i < seq_sz; i += length)
    {
        hts_pos_t len = i + length < seq_sz ? length : seq_sz - i;
        if (wrappedwrite(out, line + i, len) < len ||
              wrappedwrite(out, "\n", 1) < 1) {
            print_error_errno("faidx", "failed to write output");
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}


static int write_output(faidx_t *faid, output *out, const char *name, const int ignore,
                        const hts_pos_t length, const int rev,
                        const char *pos_strand_name, const char *neg_strand_name,
                        enum fai_format_options format) {
    hts_pos_t seq_len, wrap_len = length, len = 0;
    char *seq =  NULL, *qual = NULL;
    int ret = EXIT_FAILURE;

    if (wrap_len < 0)
        wrap_len = fai_line_length(faid, name);
    if (wrap_len <= 0)
        wrap_len = HTS_POS_MAX;

    seq = fai_fetch64(faid, name, &seq_len);
    if (rev && seq_len > 0) {
        reverse_complement(seq, seq_len);
    }
    //write the name
    len = ksprintf(&out->buffer, "%c%s%s\n", format == FAI_FASTA ? '>' : '@', name, rev ? neg_strand_name : pos_strand_name);
    if (wrappedwrite(out, out->buffer.s, out->buffer.l) < len) {
        fprintf(stderr,"[faidx] Failed to write buffer\n");
        goto exit;
    }
    ks_clear(&out->buffer);
    //write bases
    if ((ret = write_line(faid, out, seq, name, ignore, wrap_len, seq_len) == EXIT_FAILURE)) {
        goto exit;
    }

    if (format == FAI_FASTQ) {
        //write quality
        qual = fai_fetchqual64(faid, name, &seq_len);
        if (rev && seq_len > 0) {
            reverse(qual, seq_len);
        }

        len = ksprintf(&out->buffer, "+\n");
        if (wrappedwrite(out, out->buffer.s, out->buffer.l) < len) {
            fprintf(stderr,"[faidx] Failed to write buffer\n");
            goto exit;
        }
        ks_clear(&out->buffer);
        if ((ret = write_line(faid, out, qual, name, ignore, wrap_len, seq_len) == EXIT_FAILURE)) {
            goto exit;
        }
    }
    ret = EXIT_SUCCESS;

exit:
    if (seq) {
        free(seq);
    }
    if (qual) {
        free(qual);
    }
    return ret;
}


static int read_regions_from_file(faidx_t *faid, hFILE *in_file, output *out, const int ignore,
                                  const hts_pos_t length, const int rev,
                                  const char *pos_strand_name,
                                  const char *neg_strand_name,
                                  enum fai_format_options format) {
    kstring_t line = {0, 0, NULL};
    int ret = EXIT_FAILURE;

    while (line.l = 0, kgetline(&line, (kgets_func *)hgets, in_file) >= 0) {
        if ((ret = write_output(faid, out, line.s, ignore, length, rev, pos_strand_name, neg_strand_name, format)) == EXIT_FAILURE) {
            break;
        }
    }

    free(line.s);

    return ret;
}

static int usage(FILE *fp, enum fai_format_options format, int exit_status)
{
    char *tool, *file_type, *index_name;

    if (format == FAI_FASTA) {
        tool = "faidx <file.fa|file.fa.gz>";
        file_type = "FASTA";
        index_name = "file.fa";
    } else {
        tool = "fqidx <file.fq|file.fq.gz>";
        file_type = "FASTQ";
        index_name = "file.fq";
    }

    fprintf(fp, "Usage: samtools %s [<reg> [...]]\n", tool);
    fprintf(fp, "Option: \n"
                "  -o, --output FILE        Write %s to file.\n"
                "  -n, --length INT         Length of %s sequence line. [60]\n"
                "  -c, --continue           Continue after trying to retrieve missing region.\n"
                "  -r, --region-file FILE   File of regions.  Format is chr:from-to. One per line.\n"
                "  -i, --reverse-complement Reverse complement sequences.\n"
                "      --mark-strand TYPE   Add strand indicator to sequence name\n"
                "                           TYPE = rc   for /rc on negative strand (default)\n"
                "                                  no   for no strand indicator\n"
                "                                  sign for (+) / (-)\n"
                "                                  custom,<pos>,<neg> for custom indicator\n"
                "      --fai-idx      FILE  name of the index file (default %s.fai).\n"
                "      --gzi-idx      FILE  name of compressed file index (default %s.gz.gzi).\n",
                file_type, file_type, index_name, index_name);


    if (format == FAI_FASTA) {
       fprintf(fp, "  -f, --fastq              File and index in FASTQ format.\n");
    }

    fprintf(fp, "  -h, --help               This message.\n");
    sam_global_opt_help(fp, "---.-@--");

    return exit_status;
}

int faidx_core(int argc, char *argv[], enum fai_format_options format)
{
    int c, ignore_error = 0, rev = 0;
    hts_pos_t line_len = DEFAULT_FASTA_LINE_LEN ;/* fasta line len */
    char* output_file = NULL; /* output file (default is stdout ) */
    char *region_file = NULL; // list of regions from file, one per line
    char *pos_strand_name = ""; // Extension to add to name for +ve strand
    char *neg_strand_name = "/rc"; // Extension to add to name for -ve strand
    char *strand_names = NULL; // Used for custom strand annotation
    char *fai_name = NULL; // specified index name
    char *gzi_name = NULL; // specified compressed index name
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    int exit_status = EXIT_FAILURE, flushed = 0;
    struct output out = { 0, stdout, NULL, &ga, KS_INITIALIZE}; //data required for output writing
    faidx_t *fai = NULL;
    hts_tpool *pool = NULL;

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', '-', '-', 0, '-', '@'),     //output format opt and thread count - long options
        { "output", required_argument,       NULL, 'o' },
        { "help",   no_argument,             NULL, 'h' },
        { "length", required_argument,       NULL, 'n' },
        { "continue", no_argument,           NULL, 'c' },
        { "region-file", required_argument,  NULL, 'r' },
        { "fastq", no_argument,              NULL, 'f' },
        { "reverse-complement", no_argument, NULL, 'i' },
        { "mark-strand", required_argument, NULL, 1000 },
        { "fai-idx", required_argument,     NULL, 1001 },
        { "gzi-idx", required_argument,     NULL, 1002 },
        { NULL, 0, NULL, 0 }
    };

    while ((c = getopt_long(argc, argv, "ho:n:cr:fi@:", lopts, NULL)) >= 0) {
        switch (c) {
            case 'o':
                {
                    output_file = optarg;
                    char *ext = strrchr(output_file, '.');
                    if (!ext) {
                        break;
                    }
                    if (!strcmp(ext, ".gz") || !strcmp(ext, ".bgz") || !strcmp(ext, ".bgzf")) {
                        out.isbgzip = 1;        //bgzip output
                    }
                    break;
                }
            case 'n': line_len = strtol(optarg, NULL, 10);
                      if (line_len < 0) {
                        fprintf(stderr,"[faidx] bad line length '%s', using default:%d\n",optarg,ABS(DEFAULT_FASTA_LINE_LEN));
                        line_len= ABS(DEFAULT_FASTA_LINE_LEN);
                      }
                      break;
            case 'c': ignore_error = 1; break;
            case 'r': region_file = optarg; break;
            case 'f': format = FAI_FASTQ; break;
            case 'i': rev = 1; break;
            case '?':
                exit_status = usage(stderr, format, EXIT_FAILURE);
                goto exit2;
                break;
            case 'h':
                exit_status = usage(stdout, format, EXIT_SUCCESS);
                goto exit2;
                break;
            case 1000:
                if (strcmp(optarg, "no") == 0) {
                    pos_strand_name = neg_strand_name = "";
                } else if (strcmp(optarg, "sign") == 0) {
                    pos_strand_name = "(+)";
                    neg_strand_name = "(-)";
                } else if (strcmp(optarg, "rc") == 0) {
                    pos_strand_name = "";
                    neg_strand_name = "/rc";
                } else if (strncmp(optarg, "custom,", 7) == 0) {
                    size_t len = strlen(optarg + 7);
                    size_t comma = strcspn(optarg + 7, ",");
                    free(strand_names);
                    strand_names = pos_strand_name = malloc(len + 2);
                    if (!strand_names) {
                        fprintf(stderr, "[faidx] Out of memory\n");
                        goto exit2;
                    }
                    neg_strand_name = pos_strand_name + comma + 1;
                    memcpy(pos_strand_name, optarg + 7, comma);
                    pos_strand_name[comma] = '\0';
                    if (comma < len)
                        memcpy(neg_strand_name, optarg + 7 + comma + 1,
                               len - comma);
                    neg_strand_name[len - comma] = '\0';
                } else {
                    fprintf(stderr, "[faidx] Unknown --mark-strand option \"%s\"\n", optarg);
                    exit_status = usage(stderr, format, EXIT_FAILURE);
                    goto exit2;
                }
                break;
            case 1001: fai_name = optarg; break;
            case 1002: gzi_name = optarg; break;
            // handle standard samtools options like thread count, compression level...
            default:
                if (parse_sam_global_opt(c, optarg, lopts, &ga)) {
                    fprintf(stderr, "[faidx] Invalid option \"%s\"\n", optarg);
                    goto exit2;
                }
                break;
        }
    }

    if ( argc==optind ) {
        exit_status = usage(stdout, format, EXIT_SUCCESS);
        goto exit2;
    }

    if (optind+1 == argc && !region_file) {
        if (output_file && !fai_name)
            fai_name = output_file;

        if (fai_build3(argv[optind], fai_name, gzi_name) != 0) {
            if (fai_name)
                fprintf(stderr, "[faidx] Could not build fai index %s", fai_name);
            else
                fprintf(stderr, "[faidx] Could not build fai index %s.fai", argv[optind]);

            if (gzi_name)
                fprintf(stderr, " or compressed index %s\n", gzi_name);
            else
                fprintf(stderr, "\n");

            goto exit2;
        }
        exit_status = EXIT_SUCCESS;
        goto exit2;
    }

    if (out.gopt->nthreads > 0) {       //setup thread pool
        if (!(pool = hts_tpool_init(out.gopt->nthreads))) {
            fprintf(stderr, "Failed to setup thread pool\n");
        }
    }

    fai = fai_load3_format(argv[optind], fai_name, gzi_name, FAI_CREATE, format);

    if (!fai) {
        if (fai_name)
            fprintf(stderr, "[faidx] Could not load fai index %s", fai_name);
        else
            fprintf(stderr, "[faidx] Could not load fai index %s.fai", argv[optind]);

        if (gzi_name)
            fprintf(stderr, " or compressed index %s\n", gzi_name);
        else
            fprintf(stderr, "\n");

        goto exit2;
    }

    if (pool) {                         //use thread pool if set
        if (fai_thread_pool(fai, pool, 0)) {
            fprintf(stderr, "Failed to set thread pool for reading\n");
        }
    }

    /** output file provided by user */
    if( output_file != NULL ) {
        if( strcmp( output_file, argv[optind] ) == 0 ) {
            fprintf(stderr,"[faidx] Same input/output : %s\n", output_file);
            goto exit2;
        }
        if (!out.isbgzip) {
            out.fp = fopen( output_file, "w" );
        } else {
            hts_opt *opts = (hts_opt *)(out.gopt->out.specific);
            char mode[13] = "w";
            int level = 4;                                      //default compression level
            while (opts) {
                if (opts->opt == HTS_OPT_COMPRESSION_LEVEL) {   //compression level
                    level = opts->val.i;
                    break;
                }
                opts = opts->next;
            }
            if (level >= 0) {
                snprintf(mode, sizeof(mode), "w%d", level);     //pass compression with mode
            }
            out.bgzf_fp = bgzf_open(output_file, mode);
        }

        if( (!out.isbgzip && out.fp == NULL) || (out.isbgzip && out.bgzf_fp == NULL)) {
            fprintf(stderr,"[faidx] Cannot open \"%s\" for writing :%s.\n", output_file, strerror(errno) );
            goto exit2;
        }
        if (out.isbgzip && pool) {                              //use thread pool if set
            if (bgzf_thread_pool(out.bgzf_fp, pool, 0)) {
                fprintf(stderr, "Failed to set thread pool for writing\n");
            }
        }
    }

    if (region_file) {
        hFILE *rf;

        if ((rf = hopen(region_file, "r"))) {
            exit_status = read_regions_from_file(fai, rf, &out, ignore_error, line_len, rev, pos_strand_name, neg_strand_name, format);

            if (hclose(rf) != 0) {
                fprintf(stderr, "[faidx] Warning: failed to close %s", region_file);
            }
            if (exit_status == EXIT_FAILURE) {
                goto exit1;
            }
        } else {
            fprintf(stderr, "[faidx] Failed to open \"%s\" for reading.\n", region_file);
            goto exit1;
        }
    }

    exit_status = EXIT_SUCCESS;
    while ( ++optind<argc && exit_status == EXIT_SUCCESS) {
        exit_status = write_output(fai, &out, argv[optind], ignore_error, line_len, rev, pos_strand_name, neg_strand_name, format);
    }

    flushed = out.isbgzip ? bgzf_flush(out.bgzf_fp) : fflush(out.fp);
    if (flushed == EOF) {
        print_error_errno("faidx", "Failed to flush output\n");
        exit_status = EXIT_FAILURE;
    }

exit1:
    if( output_file != NULL && !out.isbgzip) {
        fclose(out.fp);     //no need to check result as already flushed
    } else if( output_file != NULL && out.isbgzip) {
        if (bgzf_close(out.bgzf_fp) < 0) {
            print_error_errno("faidx", "Failed to close output\n");
            exit_status = EXIT_FAILURE;
        }
    }

exit2:
    if (strand_names) {
        free(strand_names);
    }
    if (fai) {
        fai_destroy(fai);
    }
    if (pool) {
        hts_tpool_destroy(pool);
    }
    sam_global_args_free(&ga);
    ks_free(&out.buffer);

    return exit_status;
}


int faidx_main(int argc, char *argv[]) {
    return faidx_core(argc, argv, FAI_FASTA);
}


int fqidx_main(int argc, char *argv[]) {
    return faidx_core(argc, argv, FAI_FASTQ);
}

