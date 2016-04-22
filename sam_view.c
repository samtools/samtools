/*  sam_view.c -- SAM<->BAM<->CRAM conversion.

    Copyright (C) 2009-2015 Genome Research Ltd.
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
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <getopt.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "samtools.h"
#include "sam_opts.h"
KHASH_SET_INIT_STR(rg)

typedef khash_t(rg) *rghash_t;

// This structure contains the settings for a samview run
typedef struct samview_settings {
    rghash_t rghash;
    int min_mapQ;
    int flag_on;
    int flag_off;
    int min_qlen;
    int remove_B;
    uint32_t subsam_seed;
    double subsam_frac;
    char* library;
    void* bed;
    size_t remove_aux_len;
    char** remove_aux;
} samview_settings_t;


// TODO Add declarations of these to a viable htslib or samtools header
extern const char *bam_get_library(bam_hdr_t *header, const bam1_t *b);
extern int bam_remove_B(bam1_t *b);
extern char *samfaipath(const char *fn_ref);
void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);

// Returns 0 to indicate read should be output 1 otherwise
static int process_aln(const bam_hdr_t *h, bam1_t *b, samview_settings_t* settings)
{
    if (settings->remove_B) bam_remove_B(b);
    if (settings->min_qlen > 0) {
        int k, qlen = 0;
        uint32_t *cigar = bam_get_cigar(b);
        for (k = 0; k < b->core.n_cigar; ++k)
            if ((bam_cigar_type(bam_cigar_op(cigar[k]))&1) || bam_cigar_op(cigar[k]) == BAM_CHARD_CLIP)
                qlen += bam_cigar_oplen(cigar[k]);
        if (qlen < settings->min_qlen) return 1;
    }
    if (b->core.qual < settings->min_mapQ || ((b->core.flag & settings->flag_on) != settings->flag_on) || (b->core.flag & settings->flag_off))
        return 1;
    if (settings->bed && (b->core.tid < 0 || !bed_overlap(settings->bed, h->target_name[b->core.tid], b->core.pos, bam_endpos(b))))
        return 1;
    if (settings->subsam_frac > 0.) {
        uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(bam_get_qname(b)) ^ settings->subsam_seed);
        if ((double)(k&0xffffff) / 0x1000000 >= settings->subsam_frac) return 1;
    }
    if (settings->rghash) {
        uint8_t *s = bam_aux_get(b, "RG");
        if (s) {
            khint_t k = kh_get(rg, settings->rghash, (char*)(s + 1));
            if (k == kh_end(settings->rghash)) return 1;
        }
    }
    if (settings->library) {
        const char *p = bam_get_library((bam_hdr_t*)h, b);
        if (!p || strcmp(p, settings->library) != 0) return 1;
    }
    if (settings->remove_aux_len) {
        size_t i;
        for (i = 0; i < settings->remove_aux_len; ++i) {
            uint8_t *s = bam_aux_get(b, settings->remove_aux[i]);
            if (s) {
                bam_aux_del(b, s);
            }
        }
    }
    return 0;
}

static char *drop_rg(char *hdtxt, rghash_t h, int *len)
{
    char *p = hdtxt, *q, *r, *s;
    kstring_t str;
    memset(&str, 0, sizeof(kstring_t));
    while (1) {
        int toprint = 0;
        q = strchr(p, '\n');
        if (q == 0) q = p + strlen(p);
        if (q - p < 3) break; // the line is too short; then stop
        if (strncmp(p, "@RG\t", 4) == 0) {
            int c;
            khint_t k;
            if ((r = strstr(p, "\tID:")) != 0) {
                r += 4;
                for (s = r; *s != '\0' && *s != '\n' && *s != '\t'; ++s);
                c = *s; *s = '\0';
                k = kh_get(rg, h, r);
                *s = c;
                if (k != kh_end(h)) toprint = 1;
            }
        } else toprint = 1;
        if (toprint) {
            kputsn(p, q - p, &str); kputc('\n', &str);
        }
        p = q + 1;
    }
    *len = str.l;
    return str.s;
}

static int usage(FILE *fp, int exit_status, int is_long_help);

static int add_read_group_single(const char *subcmd, samview_settings_t *settings, char *name)
{
    char *d = strdup(name);
    int ret = 0;

    if (d == NULL) goto err;

    if (settings->rghash == NULL) {
        settings->rghash = kh_init(rg);
        if (settings->rghash == NULL) goto err;
    }

    kh_put(rg, settings->rghash, d, &ret);
    if (ret == -1) goto err;
    if (ret ==  0) free(d); /* Duplicate */
    return 0;

 err:
    print_error(subcmd, "Couldn't add \"%s\" to read group list: memory exhausted?", name);
    free(d);
    return -1;
}

static int add_read_groups_file(const char *subcmd, samview_settings_t *settings, char *fn)
{
    FILE *fp;
    char buf[1024];
    int ret = 0;
    if (settings->rghash == NULL) {
        settings->rghash = kh_init(rg);
        if (settings->rghash == NULL) {
            perror(NULL);
            return -1;
        }
    }

    fp = fopen(fn, "r");
    if (fp == NULL) {
        print_error_errno(subcmd, "failed to open \"%s\" for reading", fn);
        return -1;
    }

    while (ret != -1 && !feof(fp) && fscanf(fp, "%1023s", buf) > 0) {
        char *d = strdup(buf);
        if (d != NULL) {
            kh_put(rg, settings->rghash, d, &ret);
            if (ret == 0) free(d); /* Duplicate */
        } else {
            ret = -1;
        }
    }
    if (ferror(fp)) ret = -1;
    if (ret == -1) {
        print_error_errno(subcmd, "failed to read \"%s\"", fn);
    }
    fclose(fp);
    return (ret != -1) ? 0 : -1;
}

static inline int check_sam_write1(samFile *fp, const bam_hdr_t *h, const bam1_t *b, const char *fname, int *retp)
{
    int r = sam_write1(fp, h, b);
    if (r >= 0) return r;

    if (fname) print_error_errno("view", "writing to \"%s\" failed", fname);
    else print_error_errno("view", "writing to standard output failed");

    *retp = EXIT_FAILURE;
    return r;
}

static void check_sam_close(const char *subcmd, samFile *fp, const char *fname, const char *null_fname, int *retp)
{
    int r = sam_close(fp);
    if (r >= 0) return;

    // TODO Need error infrastructure so we can print a message instead of r
    if (fname) print_error(subcmd, "error closing \"%s\": %d", fname, r);
    else print_error(subcmd, "error closing %s: %d", null_fname, r);

    *retp = EXIT_FAILURE;
}

int main_samview(int argc, char *argv[])
{
    int c, is_header = 0, is_header_only = 0, ret = 0, compress_level = -1, is_count = 0;
    int is_long_help = 0, n_threads = 0;
    int64_t count = 0;
    samFile *in = 0, *out = 0, *un_out=0;
    bam_hdr_t *header = NULL;
    char out_mode[5], out_un_mode[5], *out_format = "";
    char *fn_in = 0, *fn_out = 0, *fn_list = 0, *q, *fn_un_out = 0;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;

    samview_settings_t settings = {
        .rghash = NULL,
        .min_mapQ = 0,
        .flag_on = 0,
        .flag_off = 0,
        .min_qlen = 0,
        .remove_B = 0,
        .subsam_seed = 0,
        .subsam_frac = -1.,
        .library = NULL,
        .bed = NULL,
    };

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 'T'),
        { "threads", required_argument, NULL, '@' },
        { NULL, 0, NULL, 0 }
    };

    /* parse command-line options */
    strcpy(out_mode, "w");
    strcpy(out_un_mode, "w");
    while ((c = getopt_long(argc, argv,
                            "SbBcCt:h1Ho:O:q:f:F:ul:r:?T:R:L:s:@:m:x:U:",
                            lopts, NULL)) >= 0) {
        switch (c) {
        case 's':
            if ((settings.subsam_seed = strtol(optarg, &q, 10)) != 0) {
                srand(settings.subsam_seed);
                settings.subsam_seed = rand();
            }
            settings.subsam_frac = strtod(q, &q);
            break;
        case 'm': settings.min_qlen = atoi(optarg); break;
        case 'c': is_count = 1; break;
        case 'S': break;
        case 'b': out_format = "b"; break;
        case 'C': out_format = "c"; break;
        case 't': fn_list = strdup(optarg); break;
        case 'h': is_header = 1; break;
        case 'H': is_header_only = 1; break;
        case 'o': fn_out = strdup(optarg); break;
        case 'U': fn_un_out = strdup(optarg); break;
        case 'f': settings.flag_on |= strtol(optarg, 0, 0); break;
        case 'F': settings.flag_off |= strtol(optarg, 0, 0); break;
        case 'q': settings.min_mapQ = atoi(optarg); break;
        case 'u': compress_level = 0; break;
        case '1': compress_level = 1; break;
        case 'l': settings.library = strdup(optarg); break;
        case 'L':
            if ((settings.bed = bed_read(optarg)) == NULL) {
                print_error_errno("view", "Could not read file \"%s\"", optarg);
                ret = 1;
                goto view_end;
            }
            break;
        case 'r':
            if (add_read_group_single("view", &settings, optarg) != 0) {
                ret = 1;
                goto view_end;
            }
            break;
        case 'R':
            if (add_read_groups_file("view", &settings, optarg) != 0) {
                ret = 1;
                goto view_end;
            }
            break;
                /* REMOVED as htslib doesn't support this
        //case 'x': out_format = "x"; break;
        //case 'X': out_format = "X"; break;
                 */
        case '?': is_long_help = 1; break;
        case 'B': settings.remove_B = 1; break;
        case '@': n_threads = strtol(optarg, 0, 0); break;
        case 'x':
            {
                if (strlen(optarg) != 2) {
                    fprintf(stderr, "main_samview: Error parsing -x auxiliary tags should be exactly two characters long.\n");
                    return usage(stderr, EXIT_FAILURE, is_long_help);
                }
                settings.remove_aux = (char**)realloc(settings.remove_aux, sizeof(char*) * (++settings.remove_aux_len));
                settings.remove_aux[settings.remove_aux_len-1] = optarg;
            }
            break;

        default:
            if (parse_sam_global_opt(c, optarg, lopts, &ga) != 0)
                return usage(stderr, EXIT_FAILURE, is_long_help);
            break;
        }
    }
    if (compress_level >= 0 && !*out_format) out_format = "b";
    if (is_header_only) is_header = 1;
    // File format auto-detection first
    if (fn_out)    sam_open_mode(out_mode+1,    fn_out,    NULL);
    if (fn_un_out) sam_open_mode(out_un_mode+1, fn_un_out, NULL);
    // Overridden by manual -b, -C
    if (*out_format)
        out_mode[1] = out_un_mode[1] = *out_format;
    out_mode[2] = out_un_mode[2] = '\0';
    // out_(un_)mode now 1 or 2 bytes long, followed by nul.
    if (compress_level >= 0) {
        char tmp[2];
        tmp[0] = compress_level + '0'; tmp[1] = '\0';
        strcat(out_mode, tmp);
        strcat(out_un_mode, tmp);
    }
    if (argc == optind && isatty(STDIN_FILENO)) return usage(stdout, EXIT_SUCCESS, is_long_help); // potential memory leak...

    fn_in = (optind < argc)? argv[optind] : "-";
    // generate the fn_list if necessary
    if (fn_list == 0 && ga.reference) fn_list = samfaipath(ga.reference);
    // open file handlers
    if ((in = sam_open_format(fn_in, "r", &ga.in)) == 0) {
        print_error_errno("view", "failed to open \"%s\" for reading", fn_in);
        ret = 1;
        goto view_end;
    }

    if (fn_list) {
        if (hts_set_fai_filename(in, fn_list) != 0) {
            fprintf(stderr, "[main_samview] failed to use reference \"%s\".\n", fn_list);
            ret = 1;
            goto view_end;
        }
    }
    if ((header = sam_hdr_read(in)) == 0) {
        fprintf(stderr, "[main_samview] fail to read the header from \"%s\".\n", fn_in);
        ret = 1;
        goto view_end;
    }
    if (settings.rghash) { // FIXME: I do not know what "bam_header_t::n_text" is for...
        char *tmp;
        int l;
        tmp = drop_rg(header->text, settings.rghash, &l);
        free(header->text);
        header->text = tmp;
        header->l_text = l;
    }
    if (!is_count) {
        if ((out = sam_open_format(fn_out? fn_out : "-", out_mode, &ga.out)) == 0) {
            print_error_errno("view", "failed to open \"%s\" for writing", fn_out? fn_out : "standard output");
            ret = 1;
            goto view_end;
        }
        if (fn_list) {
            if (hts_set_fai_filename(out, fn_list) != 0) {
                fprintf(stderr, "[main_samview] failed to use reference \"%s\".\n", fn_list);
                ret = 1;
                goto view_end;
            }
        }
        if (*out_format || is_header ||
            out_mode[1] == 'b' || out_mode[1] == 'c' ||
            (ga.out.format != sam && ga.out.format != unknown_format))  {
            if (sam_hdr_write(out, header) != 0) {
                fprintf(stderr, "[main_samview] failed to write the SAM header\n");
                ret = 1;
                goto view_end;
            }
        }
        if (fn_un_out) {
            if ((un_out = sam_open_format(fn_un_out, out_un_mode, &ga.out)) == 0) {
                print_error_errno("view", "failed to open \"%s\" for writing", fn_un_out);
                ret = 1;
                goto view_end;
            }
            if (fn_list) {
                if (hts_set_fai_filename(un_out, fn_list) != 0) {
                    fprintf(stderr, "[main_samview] failed to use reference \"%s\".\n", fn_list);
                    ret = 1;
                    goto view_end;
                }
            }
            if (*out_format || is_header ||
                out_un_mode[1] == 'b' || out_un_mode[1] == 'c' ||
                (ga.out.format != sam && ga.out.format != unknown_format))  {
                if (sam_hdr_write(un_out, header) != 0) {
                    fprintf(stderr, "[main_samview] failed to write the SAM header\n");
                    ret = 1;
                    goto view_end;
                }
            }
        }
    }

    if (n_threads > 1) { if (out) hts_set_threads(out, n_threads); }
    if (is_header_only) goto view_end; // no need to print alignments

    if (optind + 1 >= argc) { // convert/print the entire file
        bam1_t *b = bam_init1();
        int r;
        while ((r = sam_read1(in, header, b)) >= 0) { // read one alignment from `in'
            if (!process_aln(header, b, &settings)) {
                if (!is_count) { if (check_sam_write1(out, header, b, fn_out, &ret) < 0) break; }
                count++;
            } else {
                if (un_out) { if (check_sam_write1(un_out, header, b, fn_un_out, &ret) < 0) break; }
            }
        }
        if (r < -1) {
            fprintf(stderr, "[main_samview] truncated file.\n");
            ret = 1;
        }
        bam_destroy1(b);
    } else { // retrieve alignments in specified regions
        int i;
        bam1_t *b;
        hts_idx_t *idx = sam_index_load(in, fn_in); // load index
        if (idx == 0) { // index is unavailable
            fprintf(stderr, "[main_samview] random alignment retrieval only works for indexed BAM or CRAM files.\n");
            ret = 1;
            goto view_end;
        }
        b = bam_init1();
        for (i = optind + 1; i < argc; ++i) {
            int result;
            hts_itr_t *iter = sam_itr_querys(idx, header, argv[i]); // parse a region in the format like `chr2:100-200'
            if (iter == NULL) { // region invalid or reference name not found
                int beg, end;
                if (hts_parse_reg(argv[i], &beg, &end))
                    fprintf(stderr, "[main_samview] region \"%s\" specifies an unknown reference name. Continue anyway.\n", argv[i]);
                else
                    fprintf(stderr, "[main_samview] region \"%s\" could not be parsed. Continue anyway.\n", argv[i]);
                continue;
            }
            // fetch alignments
            while ((result = sam_itr_next(in, iter, b)) >= 0) {
                if (!process_aln(header, b, &settings)) {
                    if (!is_count) { if (check_sam_write1(out, header, b, fn_out, &ret) < 0) break; }
                    count++;
                } else {
                    if (un_out) { if (check_sam_write1(un_out, header, b, fn_un_out, &ret) < 0) break; }
                }
            }
            hts_itr_destroy(iter);
            if (result < -1) {
                fprintf(stderr, "[main_samview] retrieval of region \"%s\" failed due to truncated file or corrupt BAM index file\n", argv[i]);
                ret = 1;
                break;
            }
        }
        bam_destroy1(b);
        hts_idx_destroy(idx); // destroy the BAM index
    }

view_end:
    if (is_count && ret == 0)
        printf("%" PRId64 "\n", count);

    // close files, free and return
    if (in) check_sam_close("view", in, fn_in, "standard input", &ret);
    if (out) check_sam_close("view", out, fn_out, "standard output", &ret);
    if (un_out) check_sam_close("view", un_out, fn_un_out, "file", &ret);

    free(fn_list); free(fn_out); free(settings.library);  free(fn_un_out);
    sam_global_args_free(&ga);
    if ( header ) bam_hdr_destroy(header);
    if (settings.bed) bed_destroy(settings.bed);
    if (settings.rghash) {
        khint_t k;
        for (k = 0; k < kh_end(settings.rghash); ++k)
            if (kh_exist(settings.rghash, k)) free((char*)kh_key(settings.rghash, k));
        kh_destroy(rg, settings.rghash);
    }
    if (settings.remove_aux_len) {
        free(settings.remove_aux);
    }
    return ret;
}

static int usage(FILE *fp, int exit_status, int is_long_help)
{
    fprintf(fp,
"\n"
"Usage: samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]\n"
"\n"
"Options:\n"
// output options
"  -b       output BAM\n"
"  -C       output CRAM (requires -T)\n"
"  -1       use fast BAM compression (implies -b)\n"
"  -u       uncompressed BAM output (implies -b)\n"
"  -h       include header in SAM output\n"
"  -H       print SAM header only (no alignments)\n"
"  -c       print only the count of matching records\n"
"  -o FILE  output file name [stdout]\n"
"  -U FILE  output reads not selected by filters to FILE [null]\n"
// extra input
"  -t FILE  FILE listing reference names and lengths (see long help) [null]\n"
// read filters
"  -L FILE  only include reads overlapping this BED FILE [null]\n"
"  -r STR   only include reads in read group STR [null]\n"
"  -R FILE  only include reads with read group listed in FILE [null]\n"
"  -q INT   only include reads with mapping quality >= INT [0]\n"
"  -l STR   only include reads in library STR [null]\n"
"  -m INT   only include reads with number of CIGAR operations consuming\n"
"           query sequence >= INT [0]\n"
"  -f INT   only include reads with all bits set in INT set in FLAG [0]\n"
"  -F INT   only include reads with none of the bits set in INT set in FLAG [0]\n"
// read processing
"  -x STR   read tag to strip (repeatable) [null]\n"
"  -B       collapse the backward CIGAR operation\n"
"  -s FLOAT integer part sets seed of random number generator [0];\n"
"           rest sets fraction of templates to subsample [no subsampling]\n"
// general options
"  -@, --threads INT\n"
"           number of BAM/CRAM compression threads [0]\n"
"  -?       print long help, including note about region specification\n"
"  -S       ignored (input format is auto-detected)\n");

    sam_global_opt_help(fp, "-.O.T");
    fprintf(fp, "\n");

    if (is_long_help)
        fprintf(fp,
"Notes:\n"
"\n"
"1. This command now auto-detects the input format (BAM/CRAM/SAM).\n"
"   Further control over the CRAM format can be specified by using the\n"
"   --output-fmt-option, e.g. to specify the number of sequences per slice\n"
"   and to use avoid reference based compression:\n"
"\n"
"\tsamtools view -C --output-fmt-option seqs_per_slice=5000 \\\n"
"\t   --output-fmt-option no_ref -o out.cram in.bam\n"
"\n"
"   Options can also be specified as a comma separated list within the\n"
"   --output-fmt value too.  For example this is equivalent to the above\n"
"\n"
"\tsamtools view --output-fmt cram,seqs_per_slice=5000,no_ref \\\n"
"\t   -o out.cram in.bam\n"
"\n"
"2. The file supplied with `-t' is SPACE/TAB delimited with the first\n"
"   two fields of each line consisting of the reference name and the\n"
"   corresponding sequence length. The `.fai' file generated by \n"
"   `samtools faidx' is suitable for use as this file. This may be an\n"
"   empty file if reads are unaligned.\n"
"\n"
"3. SAM->BAM conversion:  samtools view -bT ref.fa in.sam.gz\n"
"\n"
"4. BAM->SAM conversion:  samtools view -h in.bam\n"
"\n"
"5. A region should be presented in one of the following formats:\n"
"   `chr1', `chr2:1,000' and `chr3:1000-2,000'. When a region is\n"
"   specified, the input alignment file must be a sorted and indexed\n"
"   alignment (BAM/CRAM) file.\n"
"\n"
"6. Option `-u' is preferred over `-b' when the output is piped to\n"
"   another samtools command.\n"
"\n");

    return exit_status;
}

int main_import(int argc, char *argv[])
{
    int argc2, ret;
    char **argv2;
    if (argc != 4) {
        fprintf(stderr, "Usage: samtools import <in.ref_list> <in.sam> <out.bam>\n");
        return 1;
    }
    argc2 = 6;
    argv2 = calloc(6, sizeof(char*));
    argv2[0] = "import", argv2[1] = "-o", argv2[2] = argv[3], argv2[3] = "-bt", argv2[4] = argv[1], argv2[5] = argv[2];
    ret = main_samview(argc2, argv2);
    free(argv2);
    return ret;
}

int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };
static const char *copied_tags[] = { "RG", "BC", "QT", NULL };

static void bam2fq_usage(FILE *to, const char *command)
{
    int fq = strcasecmp("fastq", command) == 0 || strcasecmp("bam2fq", command) == 0;
    fprintf(to,
"Usage: samtools %s [options...] <in.bam>\n", command);
    fprintf(to,
"Options:\n"
"  -0 FILE   write paired reads flagged both or neither READ1 and READ2 to FILE\n"
"  -1 FILE   write paired reads flagged READ1 to FILE\n"
"  -2 FILE   write paired reads flagged READ2 to FILE\n"
"  -f INT    only include reads with all bits set in INT set in FLAG [0]\n"
"  -F INT    only include reads with none of the bits set in INT set in FLAG [0]\n"
"  -n        don't append /1 and /2 to the read name\n");
    if (fq) fprintf(to,
"  -O        output quality in the OQ tag if present\n");
    fprintf(to,
"  -s FILE   write singleton reads to FILE [assume single-end]\n"
"  -t        copy RG, BC and QT tags to the %s header line\n",
    fq ? "FASTQ" : "FASTA");
    if (fq) fprintf(to,
"  -v INT    default quality score if not given in file [1]\n");
    sam_global_opt_help(to, "-.--.");
}

typedef enum { READ_UNKNOWN = 0, READ_1 = 1, READ_2 = 2 } readpart;
typedef enum { FASTA, FASTQ } fastfile;
typedef struct bam2fq_opts {
    char *fnse;
    char *fnr[3];
    char *fn_input; // pointer to input filename in argv do not free
    bool has12, use_oq, copy_tags;
    int flag_on, flag_off;
    sam_global_args ga;
    fastfile filetype;
    int def_qual;
} bam2fq_opts_t;

typedef struct bam2fq_state {
    samFile *fp;
    FILE *fpse;
    FILE *fpr[3];
    bam_hdr_t *h;
    bool has12, use_oq, copy_tags;
    int flag_on, flag_off;
    fastfile filetype;
    int def_qual;
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

// Transform a bam1_t record into a string with the FASTQ representation of it
// @returns false for error, true for success
static bool bam1_to_fq(const bam1_t *b, kstring_t *linebuf, const bam2fq_state_t *state)
{
    int i;
    int32_t qlen = b->core.l_qseq;
    assert(qlen >= 0);
    uint8_t *seq;
    uint8_t *qual = bam_get_qual(b);
    const uint8_t *oq = NULL;
    if (state->use_oq) {
        oq = bam_aux_get(b, "OQ");
        if (oq) oq++; // skip tag type
    }
    bool has_qual = (qual[0] != 0xff || (state->use_oq && oq)); // test if there is quality

    linebuf->l = 0;
    // Write read name
    readpart readpart = which_readpart(b);
    kputc(state->filetype == FASTA? '>' : '@', linebuf);
    kputs(bam_get_qname(b), linebuf);
    // Add the /1 /2 if requested
    if (state->has12) {
        if (readpart == READ_1) kputs("/1", linebuf);
        else if (readpart == READ_2) kputs("/2", linebuf);
    }
    if (state->copy_tags) {
        for (i = 0; copied_tags[i]; ++i) {
            uint8_t *s;
            if ((s = bam_aux_get(b, copied_tags[i])) != 0) {
                kputc('\t', linebuf);
                kputsn(copied_tags[i], 2, linebuf);
                kputsn(":Z:", 3, linebuf);
                kputs(bam_aux2Z(s), linebuf);
            }
        }
    }
    kputc('\n', linebuf);

    seq = bam_get_seq(b);

    if (b->core.flag & BAM_FREVERSE) { // read is reverse complemented
        for (i = qlen-1; i > -1; --i) {
            char c = seq_nt16_str[seq_comp_table[bam_seqi(seq,i)]];
            kputc(c, linebuf);
        }
    } else {
        for (i = 0; i < qlen; ++i) {
            char c = seq_nt16_str[bam_seqi(seq,i)];
            kputc(c, linebuf);
        }
    }
    kputc('\n', linebuf);

    if (state->filetype == FASTQ) {
        // Write quality
        kputs("+\n", linebuf);
        if (has_qual) {
            if (state->use_oq && oq) {
                if (b->core.flag & BAM_FREVERSE) { // read is reverse complemented
                    for (i = qlen-1; i > -1; --i) {
                        kputc(oq[i], linebuf);
                    }
                } else {
                    kputs((char*)oq, linebuf);
                }
            } else {
                if (b->core.flag & BAM_FREVERSE) { // read is reverse complemented
                    for (i = qlen-1; i > -1; --i) {
                        kputc(33 + qual[i], linebuf);
                    }
                } else {
                    for (i = 0; i < qlen; ++i) {
                        kputc(33 + qual[i], linebuf);
                    }
                }
            }
        } else {
            for (i = 0; i < qlen; ++i) {
                kputc(33 + state->def_qual, linebuf);
            }
        }
        kputc('\n', linebuf);
    }
    return true;
}

// return true if valid
static bool parse_opts(int argc, char *argv[], bam2fq_opts_t** opts_out)
{
    // Parse args
    bam2fq_opts_t* opts = calloc(1, sizeof(bam2fq_opts_t));
    opts->has12 = true;
    opts->filetype = FASTQ;
    opts->def_qual = 1;

    int c;
    sam_global_args_init(&opts->ga);
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', 0),
        { NULL, 0, NULL, 0 }
    };
    while ((c = getopt_long(argc, argv, "0:1:2:f:F:nOs:tv:", lopts, NULL)) > 0) {
        switch (c) {
            case '0': opts->fnr[0] = optarg; break;
            case '1': opts->fnr[1] = optarg; break;
            case '2': opts->fnr[2] = optarg; break;
            case 'f': opts->flag_on |= strtol(optarg, 0, 0); break;
            case 'F': opts->flag_off |= strtol(optarg, 0, 0); break;
            case 'n': opts->has12 = false; break;
            case 'O': opts->use_oq = true; break;
            case 's': opts->fnse = optarg; break;
            case 't': opts->copy_tags = true; break;
            case 'v': opts->def_qual = atoi(optarg); break;
            case '?': bam2fq_usage(stderr, argv[0]); free(opts); return false;
            default:
                if (parse_sam_global_opt(c, optarg, lopts, &opts->ga) != 0) {
                    bam2fq_usage(stderr, argv[0]); free(opts); return false;
                }
                break;
        }
    }

    if (opts->fnr[1] || opts->fnr[2]) opts->has12 = false;

    if (opts->def_qual < 0 || 93 < opts->def_qual) {
        fprintf(stderr, "Invalid -v default quality %i, allowed range 0 to 93\n", opts->def_qual);
        bam2fq_usage(stderr, argv[0]);
        free(opts);
        return true;
    }

    const char* type_str = argv[0];
    if (strcasecmp("fastq", type_str) == 0 || strcasecmp("bam2fq", type_str) == 0) {
        opts->filetype = FASTQ;
    } else if (strcasecmp("fasta", type_str) == 0) {
        opts->filetype = FASTA;
    } else {
        print_error("bam2fq", "Unrecognised type call \"%s\", this should be impossible... but you managed it!", type_str);
        bam2fq_usage(stderr, argv[0]);
        free(opts);
        return false;
    }

    if ((argc - (optind)) == 0) {
        bam2fq_usage(stdout, argv[0]);
        free(opts);
        return false;
    }

    if ((argc - (optind)) != 1) {
        fprintf(stderr, "Too many arguments.\n");
        bam2fq_usage(stderr, argv[0]);
        free(opts);
        return false;
    }
    opts->fn_input = argv[optind];
    *opts_out = opts;
    return true;
}

static bool init_state(const bam2fq_opts_t* opts, bam2fq_state_t** state_out)
{
    bam2fq_state_t* state = calloc(1, sizeof(bam2fq_state_t));
    state->flag_on = opts->flag_on;
    state->flag_off = opts->flag_off;
    state->has12 = opts->has12;
    state->use_oq = opts->use_oq;
    state->copy_tags = opts->copy_tags;
    state->filetype = opts->filetype;
    state->def_qual = opts->def_qual;

    state->fp = sam_open(opts->fn_input, "r");
    if (state->fp == NULL) {
        print_error_errno("bam2fq","Cannot read file \"%s\"", opts->fn_input);
        free(state);
        return false;
    }
    uint32_t rf = SAM_QNAME | SAM_FLAG | SAM_SEQ | SAM_QUAL;
    if (opts->use_oq) rf |= SAM_AUX;
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
        state->fpse = fopen(opts->fnse,"w");
        if (state->fpse == NULL) {
            print_error_errno("bam2fq", "Cannot write to singleton file \"%s\"", opts->fnse);
            free(state);
            return false;
        }
    }

    int i;
    for (i = 0; i < 3; ++i) {
        if (opts->fnr[i]) {
            state->fpr[i] = fopen(opts->fnr[i], "w");
            if (state->fpr[i] == NULL) {
                print_error_errno("bam2fq", "Cannot write to r%d file \"%s\"", i, opts->fnr[i]);
                free(state);
                return false;
            }
        } else {
            state->fpr[i] = stdout;
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
    bam_hdr_destroy(state->h);
    check_sam_close("bam2fq", state->fp, opts->fn_input, "file", status);
    if (state->fpse && fclose(state->fpse)) { print_error_errno("bam2fq", "Error closing singleton file \"%s\"", opts->fnse); valid = false; }
    int i;
    for (i = 0; i < 3; ++i) {
        if (state->fpr[i] != stdout && fclose(state->fpr[i])) { print_error_errno("bam2fq", "Error closing r%d file \"%s\"", i, opts->fnr[i]); valid = false; }
    }
    free(state);
    return valid;
}

static inline bool filter_it_out(const bam1_t *b, const bam2fq_state_t *state)
{
    return (b->core.flag&(BAM_FSECONDARY|BAM_FSUPPLEMENTARY) // skip secondary and supplementary alignments
        || (b->core.flag&(state->flag_on)) != state->flag_on // or reads indicated by filter flags
        || (b->core.flag&(state->flag_off)) != 0);

}

static bool bam2fq_mainloop_singletontrack(bam2fq_state_t *state)
{
    bam1_t* b = bam_init1();
    char *current_qname = NULL;
    int64_t n_reads = 0, n_singletons = 0; // Statistics
    kstring_t linebuf[3] = {{0,0,NULL},{0,0,NULL},{0,0,NULL}};
    int score[3];
    int at_eof;
    if (b == NULL ) {
        perror("[bam2fq_mainloop_singletontrack] Malloc error for bam record buffer.");
        return false;
    }

    bool valid = true;
    while (true) {
        at_eof = sam_read1(state->fp, state->h, b) < 0;

        if (!at_eof && filter_it_out(b, state)) continue;
        if (!at_eof) ++n_reads;

        if (at_eof || !current_qname || (strcmp(current_qname, bam_get_qname(b)) != 0)) {
            if (current_qname) {
                if (score[1] > 0 && score[2] > 0) {
                    // print linebuf[1] to fpr[1], linebuf[2] to fpr[2]
                    if (fputs(linebuf[1].s, state->fpr[1]) == EOF) { valid = false; break; }
                    if (fputs(linebuf[2].s, state->fpr[2]) == EOF) { valid = false; break; }
                } else if (score[1] > 0 || score[2] > 0) {
                    // print whichever one exists to fpse
                    if (score[1] > 0) {
                        if (fputs(linebuf[1].s, state->fpse) == EOF) { valid = false; break; }
                    } else {
                        if (fputs(linebuf[2].s, state->fpse) == EOF) { valid = false; break; }
                    }
                    ++n_singletons;
                }
                if (score[0]) { // TODO: check this
                    // print linebuf[0] to fpr[0]
                    if (fputs(linebuf[0].s, state->fpr[0]) == EOF) { valid = false; break; }
                }
            }

            if (at_eof) break;

            free(current_qname);
            current_qname = strdup(bam_get_qname(b));
            score[0] = score[1] = score[2] = 0;
        }

        // Prefer a copy of the read that has base qualities
        int b_score = bam_get_qual(b)[0] != 0xff? 2 : 1;
        if (b_score > score[which_readpart(b)]) {
            if(!bam1_to_fq(b, &linebuf[which_readpart(b)], state)) {
                fprintf(stderr, "[%s] Error converting read to FASTA/Q\n", __func__);
                return false;
            }
            score[which_readpart(b)] = b_score;
        }
    }
    if (!valid)
    {
        perror("[bam2fq_mainloop_singletontrack] Error writing to FASTx files.");
    }
    bam_destroy1(b);
    free(current_qname);
    free(linebuf[0].s);
    free(linebuf[1].s);
    free(linebuf[2].s);
    fprintf(stderr, "[M::%s] discarded %" PRId64 " singletons\n", __func__, n_singletons);
    fprintf(stderr, "[M::%s] processed %" PRId64 " reads\n", __func__, n_reads);

    return valid;
}

static bool bam2fq_mainloop(bam2fq_state_t *state)
{
    // process a name collated BAM into fastq
    bam1_t* b = bam_init1();
    if (b == NULL) {
        perror(NULL);
        return false;
    }
    int64_t n_reads = 0; // Statistics
    kstring_t linebuf = { 0, 0, NULL }; // Buffer
    while (sam_read1(state->fp, state->h, b) >= 0) {
        if (b->core.flag&(BAM_FSECONDARY|BAM_FSUPPLEMENTARY) // skip secondary and supplementary alignments
            || (b->core.flag&(state->flag_on)) != state->flag_on             // or reads indicated by filter flags
            || (b->core.flag&(state->flag_off)) != 0) continue;
        ++n_reads;

        if (!bam1_to_fq(b, &linebuf, state)) return false;
        fputs(linebuf.s, state->fpr[which_readpart(b)]);
    }
    free(linebuf.s);
    bam_destroy1(b);

    fprintf(stderr, "[M::%s] processed %" PRId64 " reads\n", __func__, n_reads);
    return true;
}

int main_bam2fq(int argc, char *argv[])
{
    int status = EXIT_SUCCESS;
    bam2fq_opts_t* opts = NULL;
    bam2fq_state_t* state = NULL;

    bool valid = parse_opts(argc, argv, &opts);
    if (!valid || opts == NULL) return valid ? EXIT_SUCCESS : EXIT_FAILURE;

    if (!init_state(opts, &state)) return EXIT_FAILURE;

    if (state->fpse) {
        if (!bam2fq_mainloop_singletontrack(state)) status = EXIT_FAILURE;
    } else {
        if (!bam2fq_mainloop(state)) status = EXIT_FAILURE;
    }

    if (!destroy_state(opts, state, &status)) return EXIT_FAILURE;
    sam_global_args_free(&opts->ga);
    free(opts);

    return status;
}
