/*  sam_view.c -- SAM<->BAM<->CRAM conversion.

    Copyright (C) 2009-2014 Genome Research Ltd.
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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <assert.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "samtools.h"
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
        if (p && strcmp(p, settings->library) != 0) return 1;
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

static int usage(int is_long_help);

static int add_read_group_single(samview_settings_t *settings, char *name)
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
    print_error("Couldn't add \"%s\" to read group list: memory exhausted?", name);
    free(d);
    return -1;
}

static int add_read_groups_file(samview_settings_t *settings, char *fn)
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
        print_error_errno("failed to open \"%s\" for reading", fn);
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
        print_error_errno("failed to read \"%s\"", fn);
    }
    fclose(fp);
    return (ret != -1) ? 0 : -1;
}

static inline int check_sam_write1(samFile *fp, const bam_hdr_t *h, const bam1_t *b, const char *fname, int *retp)
{
    int r = sam_write1(fp, h, b);
    if (r >= 0) return r;

    if (fname) print_error_errno("writing to \"%s\" failed", fname);
    else print_error_errno("writing to standard output failed");

    *retp = EXIT_FAILURE;
    return r;
}

static void check_sam_close(samFile *fp, const char *fname, const char *null_fname, int *retp)
{
    int r = sam_close(fp);
    if (r >= 0) return;

    // TODO Need error infrastructure so we can print a message instead of r
    if (fname) print_error("error closing \"%s\": %d", fname, r);
    else print_error("error closing %s: %d", null_fname, r);

    *retp = EXIT_FAILURE;
}

int main_samview(int argc, char *argv[])
{
    int c, is_header = 0, is_header_only = 0, ret = 0, compress_level = -1, is_count = 0;
    int is_long_help = 0, n_threads = 0;
    int64_t count = 0;
    samFile *in = 0, *out = 0, *un_out=0;
    bam_hdr_t *header = NULL;
    char out_mode[5], *out_format = "", *fn_out = 0, *fn_list = 0, *fn_ref = 0, *q, *fn_un_out = 0;

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

    /* parse command-line options */
    /* TODO: convert this to getopt_long we're running out of letters */
    strcpy(out_mode, "w");
    while ((c = getopt(argc, argv, "SbBcCt:h1Ho:q:f:F:ul:r:?T:R:L:s:@:m:x:U:")) >= 0) {
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
                print_error_errno("Could not read file \"%s\"", optarg);
                ret = 1;
                goto view_end;
            }
            break;
        case 'r':
            if (add_read_group_single(&settings, optarg) != 0) {
                ret = 1;
                goto view_end;
            }
            break;
        case 'R':
            if (add_read_groups_file(&settings, optarg) != 0) {
                ret = 1;
                goto view_end;
            }
            break;
                /* REMOVED as htslib doesn't support this
        //case 'x': out_format = "x"; break;
        //case 'X': out_format = "X"; break;
                 */
        case '?': is_long_help = 1; break;
        case 'T': fn_ref = strdup(optarg); break;
        case 'B': settings.remove_B = 1; break;
        case '@': n_threads = strtol(optarg, 0, 0); break;
        case 'x':
            {
                if (strlen(optarg) != 2) {
                    fprintf(stderr, "main_samview: Error parsing -x auxiliary tags should be exactly two characters long.\n");
                    return usage(is_long_help);
                }
                settings.remove_aux = (char**)realloc(settings.remove_aux, sizeof(char*) * (++settings.remove_aux_len));
                settings.remove_aux[settings.remove_aux_len-1] = optarg;
            }
            break;
        default: return usage(is_long_help);
        }
    }
    if (compress_level >= 0) out_format = "b";
    if (is_header_only) is_header = 1;
    strcat(out_mode, out_format);
    if (compress_level >= 0) {
        char tmp[2];
        tmp[0] = compress_level + '0'; tmp[1] = '\0';
        strcat(out_mode, tmp);
    }
    if (argc == optind) return usage(is_long_help); // potential memory leak...

    // generate the fn_list if necessary
    if (fn_list == 0 && fn_ref) fn_list = samfaipath(fn_ref);
    // open file handlers
    if ((in = sam_open(argv[optind], "r")) == 0) {
        print_error_errno("failed to open \"%s\" for reading", argv[optind]);
        ret = 1;
        goto view_end;
    }
    if (fn_list) hts_set_fai_filename(in, fn_list);
    if ((header = sam_hdr_read(in)) == 0) {
        fprintf(stderr, "[main_samview] fail to read the header from \"%s\".\n", argv[optind]);
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
        if ((out = sam_open(fn_out? fn_out : "-", out_mode)) == 0) {
            print_error_errno("failed to open \"%s\" for writing", fn_out? fn_out : "standard output");
            ret = 1;
            goto view_end;
        }
        if (fn_list) hts_set_fai_filename(out, fn_list);
        if (*out_format || is_header)  {
            if (sam_hdr_write(out, header) != 0) {
                fprintf(stderr, "[main_samview] failed to write the SAM header\n");
                ret = 1;
                goto view_end;
            }
        }
        if (fn_un_out) {
            if ((un_out = sam_open(fn_un_out, out_mode)) == 0) {
                print_error_errno("failed to open \"%s\" for writing", fn_un_out);
                ret = 1;
                goto view_end;
            }
            if (*out_format || is_header) {
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

    if (argc == optind + 1) { // convert/print the entire file
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
        hts_idx_t *idx = sam_index_load(in, argv[optind]); // load index
        if (idx == 0) { // index is unavailable
            fprintf(stderr, "[main_samview] random alignment retrieval only works for indexed BAM or CRAM files.\n");
            ret = 1;
            goto view_end;
        }
        b = bam_init1();
        for (i = optind + 1; i < argc; ++i) {
            int result;
            hts_itr_t *iter = sam_itr_querys(idx, header, argv[i]); // parse a region in the format like `chr2:100-200'
            if (iter == NULL) { // reference name is not found
                fprintf(stderr, "[main_samview] region \"%s\" specifies an unknown reference name. Continue anyway.\n", argv[i]);
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
    if (in) check_sam_close(in, argv[optind], "standard input", &ret);
    if (out) check_sam_close(out, fn_out, "standard output", &ret);
    if (un_out) check_sam_close(un_out, fn_un_out, "file", &ret);

    free(fn_list); free(fn_ref); free(fn_out); free(settings.library);  free(fn_un_out);
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

static int usage(int is_long_help)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]\n\n");
    // output options
    fprintf(stderr, "Options: -b       output BAM\n");
    fprintf(stderr, "         -C       output CRAM (requires -T)\n");
    fprintf(stderr, "         -1       use fast BAM compression (implies -b)\n");
    fprintf(stderr, "         -u       uncompressed BAM output (implies -b)\n");
    fprintf(stderr, "         -h       include header in SAM output\n");
    fprintf(stderr, "         -H       print SAM header only (no alignments)\n");
    fprintf(stderr, "         -c       print only the count of matching records\n");
    fprintf(stderr, "         -o FILE  output file name [stdout]\n");
    fprintf(stderr, "         -U FILE  output reads not selected by filters to FILE [null]\n");
    // extra input
    fprintf(stderr, "         -t FILE  FILE listing reference names and lengths (see long help) [null]\n");
    fprintf(stderr, "         -T FILE  reference sequence FASTA FILE [null]\n");
    // read filters
    fprintf(stderr, "         -L FILE  only include reads overlapping this BED FILE [null]\n");
    fprintf(stderr, "         -r STR   only include reads in read group STR [null]\n");
    fprintf(stderr, "         -R FILE  only include reads with read group listed in FILE [null]\n");
    fprintf(stderr, "         -q INT   only include reads with mapping quality >= INT [0]\n");
    fprintf(stderr, "         -l STR   only include reads in library STR [null]\n");
    fprintf(stderr, "         -m INT   only include reads with number of CIGAR operations\n");
    fprintf(stderr, "                  consuming query sequence >= INT [0]\n");
    fprintf(stderr, "         -f INT   only include reads with all bits set in INT set in FLAG [0]\n");
    fprintf(stderr, "         -F INT   only include reads with none of the bits set in INT\n");
    fprintf(stderr, "                  set in FLAG [0]\n");
    // read processing
    fprintf(stderr, "         -x STR   read tag to strip (repeatable) [null]\n");
    fprintf(stderr, "         -B       collapse the backward CIGAR operation\n");
    fprintf(stderr, "         -s FLOAT integer part sets seed of random number generator [0];\n");
    fprintf(stderr, "                  rest sets fraction of templates to subsample [no subsampling]\n");
    // general options
    fprintf(stderr, "         -@ INT   number of BAM compression threads [0]\n");
    fprintf(stderr, "         -?       print long help, including note about region specification\n");
    fprintf(stderr, "         -S       ignored (input format is auto-detected)\n");
    fprintf(stderr, "\n");
    if (is_long_help)
        fprintf(stderr, "Notes:\n\
\n\
  1. This command now auto-detects the input format (BAM/CRAM/SAM).\n\
\n\
  2. The file supplied with `-t' is SPACE/TAB delimited with the first\n\
     two fields of each line consisting of the reference name and the\n\
     corresponding sequence length. The `.fai' file generated by \n\
     `samtools faidx' is suitable for use as this file. This may be an\n\
     empty file if reads are unaligned.\n\
\n\
  3. SAM->BAM conversion: `samtools view -bT ref.fa in.sam.gz'.\n\
\n\
  4. BAM->SAM conversion: `samtools view -h in.bam'.\n\
\n\
  5. A region should be presented in one of the following formats:\n\
     `chr1', `chr2:1,000' and `chr3:1000-2,000'. When a region is\n\
     specified, the input alignment file must be a sorted and indexed\n\
     alignment (BAM/CRAM) file.\n\
\n\
  6. Option `-u' is preferred over `-b' when the output is piped to\n\
     another samtools command.\n\
\n");
    return 1;
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

static void bam2fq_usage(FILE *to)
{
    fprintf(to, "\nUsage:   samtools bam2fq [-nO] [-s <outSE.fq>] <in.bam>\n\n");
    fprintf(to, "Options: -n        don't append /1 and /2 to the read name\n");
    fprintf(to, "         -O        output quality in the OQ tag if present\n");
    fprintf(to, "         -s FILE   write singleton reads to FILE [assume single-end]\n");
    fprintf(to, "\n");
}

int main_bam2fq(int argc, char *argv[])
{
    samFile *fp;
    bam_hdr_t *h;
    bam1_t *b;
    int8_t *buf;
    int status = EXIT_SUCCESS;
    size_t max_buf;
    FILE* fpse;
    // Parse args
    char* fnse = NULL;
    bool has12 = true, use_oq = false;
    int c;
    while ((c = getopt(argc, argv, "nOs:")) > 0) {
        switch (c) {
            case 'n': has12 = false; break;
            case 'O': use_oq = true; break;
            case 's': fnse = optarg; break;
            default: bam2fq_usage(stderr); return 1;
        }
    }

    if ((argc - (optind)) == 0) {
        bam2fq_usage(stdout);
        return 0;
    }

    if ((argc - (optind)) != 1) {
        fprintf(stderr, "Too many arguments.\n");
        bam2fq_usage(stderr);
        return 1;
    }

    fp = sam_open(argv[optind], "r");
    if (fp == NULL) {
        print_error_errno("Cannot read file \"%s\"", argv[optind]);
        return 1;
    }
    if (hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS,
                    SAM_QNAME | SAM_FLAG | SAM_SEQ | SAM_QUAL)) {
        fprintf(stderr, "Failed to set CRAM_OPT_REQUIRED_FIELDS value\n");
        return 1;
    }
    if (hts_set_opt(fp, CRAM_OPT_DECODE_MD, 0)) {
        fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
        return 1;
    }
    fpse = NULL;
    if (fnse) {
        fpse = fopen(fnse,"w");
        if (fpse == NULL) {
            print_error_errno("Cannot write to singleton file \"%s\"", fnse);
            return 1;
        }
    }

    h = sam_hdr_read(fp);
    b = bam_init1();
    buf = NULL;
    max_buf = 0;

    int64_t n_singletons = 0, n_reads = 0;
    char* previous = NULL;
    kstring_t linebuf = { 0, 0, NULL };
    kputsn("", 0, &linebuf);

    while (sam_read1(fp, h, b) >= 0) {
        if (b->core.flag&(BAM_FSECONDARY|BAM_FSUPPLEMENTARY)) continue; // skip secondary and supplementary alignments
        ++n_reads;

        int i;
        int32_t qlen = b->core.l_qseq;
        assert(qlen >= 0);
        uint8_t* seq;
        uint8_t* qual = bam_get_qual(b);
        const uint8_t *oq = NULL;
        if (use_oq) oq = bam_aux_get(b, "OQ");
        bool has_qual = (qual[0] != 0xff || (use_oq && oq)); // test if there is quality

        // If there was a previous readname
        if ( fpse && previous ) {
            if (!strcmp(bam_get_qname(b), previous ) ) {
                fputs(linebuf.s, stdout); // Write previous read
                free(previous);
                previous = NULL;
            } else { // Doesn't match it's a singleton
                ++n_singletons;
                fputs(linebuf.s, fpse);  // Write previous read to singletons
                free(previous);
                previous = strdup(bam_get_qname(b));
            }
        } else {
            fputs(linebuf.s, stdout); // Write pending read
            if (fpse) previous = strdup(bam_get_qname(b));
        }

        linebuf.l = 0;
        // Write read name
        kputc(!has_qual? '>' : '@', &linebuf);
        kputs(bam_get_qname(b), &linebuf);
        // Add the /1 /2 if requested
        if (has12) {
            if ((b->core.flag & BAM_FREAD1) && !(b->core.flag & BAM_FREAD2)) kputs("/1\n", &linebuf);
            else if ((b->core.flag & BAM_FREAD2) && !(b->core.flag & BAM_FREAD1)) kputs("/2\n", &linebuf);
            else kputc('\n', &linebuf);
        } else {
            kputc('\n', &linebuf);
        }

        if (max_buf < qlen + 1) {
            max_buf = qlen + 1;
            kroundup32(max_buf);
            buf = realloc(buf, max_buf);
            if (buf == NULL) {
                fprintf(stderr, "Out of memory");
                return 1;
            }
        }
        buf[qlen] = '\0';
        seq = bam_get_seq(b);
        for (i = 0; i < qlen; ++i)
            buf[i] = bam_seqi(seq, i);
        if (b->core.flag & BAM_FREVERSE) { // reverse complement
            for (i = 0; i < qlen>>1; ++i) {
                int8_t t = seq_comp_table[buf[qlen - 1 - i]];
                buf[qlen - 1 - i] = seq_comp_table[buf[i]];
                buf[i] = t;
            }
            if (qlen&1) buf[i] = seq_comp_table[buf[i]];
        }
        for (i = 0; i < qlen; ++i)
            buf[i] = seq_nt16_str[buf[i]];
        kputs((char*)buf, &linebuf);
        kputc('\n', &linebuf);

        if (has_qual) {
            // Write quality
            kputs("+\n", &linebuf);
            if (use_oq && oq) memcpy(buf, oq + 1, qlen);
            else {
                for (i = 0; i < qlen; ++i)
                    buf[i] = 33 + qual[i];
            }
            if (b->core.flag & BAM_FREVERSE) { // reverse
                for (i = 0; i < qlen>>1; ++i) {
                    int8_t t = buf[qlen - 1 - i];
                    buf[qlen - 1 - i] = buf[i];
                    buf[i] = t;
                }
            }
            kputs((char*)buf, &linebuf);
            kputc('\n', &linebuf);
        }
    }

    if (fpse) {
        if ( previous ) { // Nothing left to match it's a singleton
            ++n_singletons;
            fputs(linebuf.s, fpse);  // Write previous read to singletons
        } else {
            fputs(linebuf.s, stdout); // Write previous read
        }

        fprintf(stderr, "[M::%s] discarded %" PRId64 " singletons\n", __func__, n_singletons);
        fclose(fpse);
    } else {
        fputs(linebuf.s, stdout); // Write previous read
    }
    free(linebuf.s);
    free(previous);

    fprintf(stderr, "[M::%s] processed %" PRId64 " reads\n", __func__, n_reads);

    free(buf);
    bam_destroy1(b);
    bam_hdr_destroy(h);
    check_sam_close(fp, argv[optind], "file", &status);
    return status;
}
