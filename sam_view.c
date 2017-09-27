/*  sam_view.c -- SAM<->BAM<->CRAM conversion.

    Copyright (C) 2009-2017 Genome Research Ltd.
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
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include <getopt.h>
#include <ctype.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/klist.h"
#include "htslib/thread_pool.h"
#include "htslib/bgzf.h"
#include "samtools.h"
#include "sam_opts.h"

#define DEFAULT_BARCODE_TAG "BC"
#define DEFAULT_QUALITY_TAG "QT"

KHASH_SET_INIT_STR(rg)
#define taglist_free(p)
KLIST_INIT(ktaglist, char*, taglist_free)

typedef khash_t(rg) *rghash_t;

// This structure contains the settings for a samview run
typedef struct samview_settings {
    rghash_t rghash;
    int min_mapQ;
    int flag_on;
    int flag_off;
    int flag_alloff;
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
    if (settings->flag_alloff && ((b->core.flag & settings->flag_alloff) == settings->flag_alloff))
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
    int is_long_help = 0;
    int64_t count = 0;
    samFile *in = 0, *out = 0, *un_out=0;
    FILE *fp_out = NULL;
    bam_hdr_t *header = NULL;
    char out_mode[5], out_un_mode[5], *out_format = "";
    char *fn_in = 0, *fn_out = 0, *fn_list = 0, *q, *fn_un_out = 0;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    htsThreadPool p = {NULL, 0};

    samview_settings_t settings = {
        .rghash = NULL,
        .min_mapQ = 0,
        .flag_on = 0,
        .flag_off = 0,
        .flag_alloff = 0,
        .min_qlen = 0,
        .remove_B = 0,
        .subsam_seed = 0,
        .subsam_frac = -1.,
        .library = NULL,
        .bed = NULL,
    };

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 'T', '@'),
        { NULL, 0, NULL, 0 }
    };

    /* parse command-line options */
    strcpy(out_mode, "w");
    strcpy(out_un_mode, "w");
    while ((c = getopt_long(argc, argv,
                            "SbBcCt:h1Ho:O:q:f:F:G:ul:r:?T:R:L:s:@:m:x:U:",
                            lopts, NULL)) >= 0) {
        switch (c) {
        case 's':
            if ((settings.subsam_seed = strtol(optarg, &q, 10)) != 0) {
                // Convert likely user input 0,1,2,... to pseudo-random
                // values with more entropy and more bits set
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
        case 'G': settings.flag_alloff |= strtol(optarg, 0, 0); break;
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
    else {
        if (fn_out) {
            fp_out = fopen(fn_out, "w");
            if (fp_out == NULL) {
                print_error_errno("view", "can't create \"%s\"", fn_out);
                ret = EXIT_FAILURE;
                goto view_end;
            }
        }
    }

    if (ga.nthreads > 1) {
        if (!(p.pool = hts_tpool_init(ga.nthreads))) {
            fprintf(stderr, "Error creating thread pool\n");
            ret = 1;
            goto view_end;
        }
        hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
        if (out) hts_set_opt(out, HTS_OPT_THREAD_POOL, &p);
    }
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
    if (is_count && ret == 0) {
        if (fprintf(fn_out? fp_out : stdout, "%" PRId64 "\n", count) < 0) {
            if (fn_out) print_error_errno("view", "writing to \"%s\" failed", fn_out);
            else print_error_errno("view", "writing to standard output failed");
            ret = EXIT_FAILURE;
        }
    }

    // close files, free and return
    if (in) check_sam_close("view", in, fn_in, "standard input", &ret);
    if (out) check_sam_close("view", out, fn_out, "standard output", &ret);
    if (un_out) check_sam_close("view", un_out, fn_un_out, "file", &ret);
    if (fp_out) fclose(fp_out);

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

    if (p.pool)
        hts_tpool_destroy(p.pool);

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
"  -f INT   only include reads with all  of the FLAGs in INT present [0]\n"       //   F&x == x
"  -F INT   only include reads with none of the FLAGS in INT present [0]\n"       //   F&x == 0
"  -G INT   only EXCLUDE reads with all  of the FLAGs in INT present [0]\n"       // !(F&x == x)
"  -s FLOAT subsample reads (given INT.FRAC option value, 0.FRAC is the\n"
"           fraction of templates/read pairs to keep; INT part sets seed)\n"
// read processing
"  -x STR   read tag to strip (repeatable) [null]\n"
"  -B       collapse the backward CIGAR operation\n"
// general options
"  -?       print long help, including note about region specification\n"
"  -S       ignored (input format is auto-detected)\n");

    sam_global_opt_help(fp, "-.O.T@");
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
"  -0 FILE              write paired reads flagged both or neither READ1 and READ2 to FILE\n"
"  -1 FILE              write paired reads flagged READ1 to FILE\n"
"  -2 FILE              write paired reads flagged READ2 to FILE\n"
"  -f INT               only include reads with all  of the FLAGs in INT present [0]\n"       //   F&x == x
"  -F INT               only include reads with none of the FLAGS in INT present [0]\n"       //   F&x == 0
"  -G INT               only EXCLUDE reads with all  of the FLAGs in INT present [0]\n"       // !(F&x == x)
"  -n                   don't append /1 and /2 to the read name\n"
"  -N                   always append /1 and /2 to the read name\n");
    if (fq) fprintf(to,
"  -O                   output quality in the OQ tag if present\n");
    fprintf(to,
"  -s FILE              write singleton reads to FILE [assume single-end]\n"
"  -t                   copy RG, BC and QT tags to the %s header line\n",
    fq ? "FASTQ" : "FASTA");
    fprintf(to,
"  -T TAGLIST           copy arbitrary tags to the %s header line\n",
    fq ? "FASTQ" : "FASTA");
    if (fq) fprintf(to,
"  -v INT               default quality score if not given in file [1]\n"
"  -i                   add Illumina Casava 1.8 format entry to header (eg 1:N:0:ATCACG)\n"
"  -c                   compression level [0..9] to use when creating gz or bgzf fastq files\n"
"  --i1 FILE            write first index reads to FILE\n"
"  --i2 FILE            write second index reads to FILE\n"
"  --barcode-tag TAG    Barcode tag [default: " DEFAULT_BARCODE_TAG "]\n"
"  --quality-tag TAG    Quality tag [default: " DEFAULT_QUALITY_TAG "]\n"
"  --index-format STR   How to parse barcode and quality tags\n\n");
    sam_global_opt_help(to, "-.--.@");
    fprintf(to,
"   \n"
"   The index-format string describes how to parse the barcode and quality tags, for example:\n"
"   i14i8       the first 14 characters are index 1, the next 8 characters are index 2\n"
"   n8i14       ignore the first 8 characters, and use the next 14 characters for index 1\n"
"   If the tag contains a separator, then the numeric part can be replaced with '*' to mean\n"
"   'read until the separator or end of tag', for example:\n"
"   n*i*        ignore the left part of the tag until the separator, then use the second part\n"
"               of the tag as index 1\n");
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
    bam_hdr_t *h;
    bool has12, use_oq, copy_tags, illumina_tag;
    int flag_on, flag_off, flag_alloff;
    fastfile filetype;
    int def_qual;
    klist_t(ktaglist) *taglist;
    char *index_sequence;
    char compression_level;
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

        if (action=='i' && *sub_tag && state->fpi[file_number]) {
            //if (file_number==0) state->index_sequence = strdup(sub_tag);    // we're going to need this later...
            state->index_sequence = strdup(sub_tag);    // we're going to need this later...
            if (!state->index_sequence) goto fail;
            if (!make_fq_line(rec, sub_tag, sub_qual, &linebuf, state)) goto fail;
            if (state->illumina_tag) {
                if (insert_index_sequence_into_linebuf(state->index_sequence, &linebuf, rec) < 0) {
                    goto fail;
                }
            }
            if (bgzf_write(state->fpi[file_number++], linebuf.s, linebuf.l) < 0)
                goto fail;
        }

    }

    free(sub_qual); free(sub_tag);
    free(linebuf.s);
    return true;

 fail:
    perror(__func__);
    free(sub_qual); free(sub_tag);
    free(linebuf.s);
    return true;
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

    if (state->use_oq) {
        oq = bam_aux_get(b, "OQ");
        if (oq) {
            oq++;
            qual = strdup(bam_aux2Z(oq));
            if (!qual) goto fail;
            if (b->core.flag & BAM_FREVERSE) { // read is reverse complemented
                reverse(qual);
            }
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
    while ((c = getopt_long(argc, argv, "0:1:2:f:F:G:niNOs:c:tT:v:@:", lopts, NULL)) > 0) {
        switch (c) {
            case 'b': opts->barcode_tag = strdup(optarg); break;
            case 'q': opts->quality_tag = strdup(optarg); break;
            case  1 : opts->index_file[0] = optarg; break;
            case  2 : opts->index_file[1] = optarg; break;
            case  3 : opts->index_format = strdup(optarg); break;
            case '0': opts->fnr[0] = optarg; break;
            case '1': opts->fnr[1] = optarg; break;
            case '2': opts->fnr[2] = optarg; break;
            case 'f': opts->flag_on |= strtol(optarg, 0, 0); break;
            case 'F': opts->flag_off |= strtol(optarg, 0, 0); break;
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

    if (nIndex==2 && !opts->index_file[1]) {
        fprintf(stderr, "index_format specifies two indexes, but only one index file given\n");
        bam2fq_usage(stderr, argv[0]);
        free_opts(opts);
        return false;
    }

    if (nIndex==1 && !opts->index_file[0]) {
        fprintf(stderr, "index_format specifies an index, but no index file given\n");
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

    if ((argc - (optind)) == 0) {
        fprintf(stderr, "No input file specified.\n");
        bam2fq_usage(stdout, argv[0]);
        free_opts(opts);
        return false;
    }

    if ((argc - (optind)) != 1) {
        fprintf(stderr, "Too many arguments.\n");
        bam2fq_usage(stderr, argv[0]);
        free_opts(opts);
        return false;
    }
    opts->fn_input = argv[optind];
    *opts_out = opts;
    return true;
}

static BGZF *open_fqfile(char *filename, int c)
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

    return bgzf_open(filename,mode);
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
    state->hstdout = bgzf_dopen(fileno(stdout), "wu");
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
    if (opts->ga.nthreads > 0)
        hts_set_threads(state->fp, opts->ga.nthreads);
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
        state->fpse = open_fqfile(opts->fnse, state->compression_level);
        if (state->fpse == NULL) {
            print_error_errno("bam2fq", "Cannot write to singleton file \"%s\"", opts->fnse);
            free(state);
            return false;
        }
    }

    int i;
    for (i = 0; i < 3; ++i) {
        if (opts->fnr[i]) {
            state->fpr[i] = open_fqfile(opts->fnr[i], state->compression_level);
            if (state->fpr[i] == NULL) {
                print_error_errno("bam2fq", "Cannot write to r%d file \"%s\"", i, opts->fnr[i]);
                free(state);
                return false;
            }
        } else {
            state->fpr[i] = state->hstdout;
        }
    }
    for (i = 0; i < 2; i++) {
        state->fpi[i] = NULL;
        if (opts->index_file[i]) {
            state->fpi[i] = open_fqfile(opts->index_file[i], state->compression_level);
            if (state->fpi[i] == NULL) {
                print_error_errno("bam2fq", "Cannot write to i%d file \"%s\"", i+1, opts->index_file[i]);
                free(state);
                return false;
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
    bam_hdr_destroy(state->h);
    check_sam_close("bam2fq", state->fp, opts->fn_input, "file", status);
    if (state->fpse && bgzf_close(state->fpse)) { print_error_errno("bam2fq", "Error closing singleton file \"%s\"", opts->fnse); valid = false; }
    int i;
    for (i = 0; i < 3; ++i) {
        if (state->fpr[i] == state->hstdout) {
            if (i==0 && bgzf_close(state->fpr[i])) { print_error_errno("bam2fq", "Error closing STDOUT"); valid = false; }
        } else {
            if (bgzf_close(state->fpr[i])) { print_error_errno("bam2fq", "Error closing r%d file \"%s\"", i, opts->fnr[i]); valid = false; }
        }
    }
    for (i = 0; i < 2; i++) {
        if (state->fpi[i] && bgzf_close(state->fpi[i])) {
            print_error_errno("bam2fq", "Error closing i%d file \"%s\"", i+1, opts->index_file[i]);
            valid = false;
        }
    }
    kl_destroy(ktaglist,state->taglist);
    free(state->index_sequence);
    free(state);
    return valid;
}

static inline bool filter_it_out(const bam1_t *b, const bam2fq_state_t *state)
{
    return (b->core.flag&(BAM_FSECONDARY|BAM_FSUPPLEMENTARY) // skip secondary and supplementary alignments
        || (b->core.flag&(state->flag_on)) != state->flag_on // or reads indicated by filter flags
        || (b->core.flag&(state->flag_off)) != 0
        || (b->core.flag&(state->flag_alloff) && (b->core.flag&(state->flag_alloff)) == state->flag_alloff));

}

static bool bam2fq_mainloop(bam2fq_state_t *state, bam2fq_opts_t* opts)
{
    int n;
    bam1_t *records[3];
    bam1_t* b = bam_init1();
    char *current_qname = NULL;
    int64_t n_reads = 0, n_singletons = 0; // Statistics
    kstring_t linebuf[3] = {{0,0,NULL},{0,0,NULL},{0,0,NULL}};
    int score[3];
    int at_eof;
    if (b == NULL ) {
        perror("[bam2fq_mainloop] Malloc error for bam record buffer.");
        return false;
    }

    bool valid = true;
    while (true) {
        int res = sam_read1(state->fp, state->h, b);
        if (res < -1) {
            fprintf(stderr, "[bam2fq_mainloop] Failed to read bam record.\n");
            return false;
        }
        at_eof = res < 0;

        if (!at_eof && filter_it_out(b, state)) continue;
        if (!at_eof) ++n_reads;

        if (at_eof || !current_qname || (strcmp(current_qname, bam_get_qname(b)) != 0)) {
            if (current_qname) {
                if (state->illumina_tag) {
                    for (n=0; valid && n<3; n++) {
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

            if (at_eof) break;

            free(current_qname);
            current_qname = strdup(bam_get_qname(b));
            if (!current_qname) { valid = false; break; }
            score[0] = score[1] = score[2] = 0;
        }

        // Prefer a copy of the read that has base qualities
        int b_score = bam_get_qual(b)[0] != 0xff? 2 : 1;
        if (b_score > score[which_readpart(b)]) {
            if (state->fpi[0]) if (!tags2fq(b, state, opts)) return false;
            records[which_readpart(b)] = b;
            if(!bam1_to_fq(b, &linebuf[which_readpart(b)], state)) {
                fprintf(stderr, "[%s] Error converting read to FASTA/Q\n", __func__);
                return false;
            }
            score[which_readpart(b)] = b_score;
        }
    }
    if (!valid)
    {
        perror("[bam2fq_mainloop] Error writing to FASTx files.");
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
