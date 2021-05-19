/*  bam_samples -- print samples in a set of BAM files

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
#include <htslib/hts.h>
#include <htslib/kseq.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/hfile.h>
#include "htslib/khash.h"
#include <samtools.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

KHASH_MAP_INIT_STR(sm, int)

/** and chained struct containing the faidx and the fasta filename
    will be compared with the @SQ lines in the SAM header
 */
typedef struct FaidxPath {
    /** path to reference */
    char* filename;
    /** fasta index  */
    faidx_t * faidx;
    struct FaidxPath* next;
} FaidxPath;

/** program parameters */
typedef struct Params {
    /** output stream */
    FILE* out;
    /** tag in @RG line. default is "SM" */
    char tag[3];
    /** first faidx/path in chained list */
    FaidxPath* faidx;
    /** enable more than one sample in a bam */
    int enable_multiple;
    /** enable no sample in a bam */
    int enable_missing;
    /** show whether the bam is indexed */
    int test_index;
} Params;

/** print usage */
static void usage_samples(FILE *write_to) {
    fprintf(write_to,
            "Usage: samtools samples [options] <input> [...]\n"
            "       find dir1 dir2 -type f \\(-name \"*.bam\" -o -name \"*.cram\" \\) | samtools samples [options]\n"
            "\n"
            "Options:\n"
            "  -h              print help and exit\n"
            "  -H              print a header\n"
            "  -i              test if the file is indexed.\n"
            "  -T <tag>        set the TAG in the @RG line [SM].\n"
            "  -o <file>       output file [stdout].\n"
            "  -d              enable multiple samples in one bam\n"
            "  -m              enable missing sample in one bam. \".\" will be used as the default name.\n"
            "  -f <file.fa>    add an indexed fasta in the collection of reference. Can be used multiple times.\n"
            "  -F <file.txt>   read a file containing the path to the indexed fasta references. One path per line.\n"
            "\n"
            " Using -f or -F will add a column containing the path to the reference or \".\" if the reference was not found.\n"
            "\n"
    );
}


/** loads fasta fai file into FaidxPath, add it to params->faidx */
static int load_dictionary(struct Params* params, const char* filename) {
    FaidxPath* prev = params->faidx;
    FaidxPath* ptr = (struct FaidxPath*)malloc(sizeof(struct FaidxPath));
    if (ptr == NULL) {
        print_error("samples", "out of memory");
        return EXIT_FAILURE;
    }
    ptr->filename = strdup(filename);
    if (ptr->filename == NULL) {
        print_error("samples", "out of memory");
        return EXIT_FAILURE;
    }
    ptr->faidx = fai_load(filename);
    if (ptr->faidx == NULL) {
        print_error_errno("samples", "cannot load index from \"%s\".", filename);
        return EXIT_FAILURE;
    }
    /* insert in chained list */
    params->faidx = ptr;
    ptr->next = prev;
    return EXIT_SUCCESS;
}

/** load a faidx file and append it to params */
static int load_dictionaries(Params* params, const char* filename) {
    int ret;
    htsFile* in;
    int status = EXIT_SUCCESS;
    kstring_t  *ks = NULL;
    in = hts_open(filename, "r");

    ks = &in->line;
    if (in == NULL) {
        print_error_errno("samples", "cannot open \"%s\".", filename);
        return EXIT_FAILURE;
    }
    while ((ret = hts_getline(in, KS_SEP_LINE, ks)) >= 0) {
        if (load_dictionary(params, ks->s)!= EXIT_SUCCESS) {
            status = EXIT_FAILURE;
            break;
        }
    }
    hts_close(in);
    return status;
}

/** print the sample information, search for a reference */
static int print_sample(Params* params, sam_hdr_t *header, int has_index, const char* sample, const char* fname) {
    fputs(sample, params->out);
    fputc('\t', params->out);
    fputs(fname, params->out);
    if (params->test_index) {
        fprintf(params->out, "\t%d", has_index);
    }
    if (params->faidx!=NULL) {
        FaidxPath* ref = NULL;
        FaidxPath* curr = params->faidx;
        while(curr!=NULL) {
            /** check names and length are the same in the same order */
            if (faidx_nseq(curr->faidx) == header->n_targets) {
                int i;
                for (i = 0; i < faidx_nseq(curr->faidx); i++) {
                    /** check name if the same */
                    if (strcmp(faidx_iseq(curr->faidx, i), header->target_name[i])!=0) break;
                    /** check length if the same */
                    if (faidx_seq_len(curr->faidx, faidx_iseq(curr->faidx, i)) != header->target_len[i]) break;
                }
                if (i == faidx_nseq(curr->faidx)) {
                    ref = curr;
                    break;
                }
            }
            curr = curr->next;
        }
        fputc('\t', params->out);
        if (ref == NULL) {
            fputc('.', params->out);
        } else {
            fputs(curr->filename, params->out);
        }
    }
    fputc('\n', params->out);
    return 0;
}


/** open a sam file. Search for all samples in the @RG lines */
static int print_samples(Params* params, const char* fname) {
    samFile *in = 0;
    sam_hdr_t *header = NULL;
    int n_rg;
    int status = EXIT_SUCCESS;
    khash_t(sm) *sample_set = kh_init(sm);
    khint_t k;
    int count_samples = 0;
    int has_index = 0;


    /* open sam file */
    if ((in = sam_open_format(fname, "r", NULL)) == 0) {
        print_error_errno("samples", "failed to open \"%s\" for reading.", fname);
        status = EXIT_FAILURE;
        goto end_print;
    }
    /* load header */
    if ((header = sam_hdr_read(in)) == 0) {
        print_error("samples", "failed to read the header from \"%s\".", fname);
        status = EXIT_FAILURE;
        goto end_print;
    }

    /* try to load index if required */
    if (params->test_index) {
        hts_idx_t *bam_idx =  sam_index_load(in, fname);
        has_index = bam_idx!=NULL;
        if (bam_idx != NULL) hts_idx_destroy(bam_idx);
    }


    /* get the RG lines */
    n_rg = sam_hdr_count_lines(header, "RG");
    if (n_rg > 0) {
        int i, r, ret;
        char* sample;
        kstring_t sm_val = KS_INITIALIZE;
        /* loop over the RG lines and search for the params->tag */
        for (i = 0; i < n_rg; i++) {
            r = sam_hdr_find_tag_pos(header, "RG", i, params->tag, &sm_val);
            if (r < 0) continue;
            k = kh_get(sm, sample_set, sm_val.s);
            if (k != kh_end(sample_set)) continue;
            sample = strdup(sm_val.s);
            if (sample == NULL) {
                print_error("samples", "out of memory.");
                goto end_print;
            }
            kh_put(sm, sample_set, sample, &ret);
            ++count_samples;
        }
        ks_free(&sm_val);
    }
    if (count_samples == 0) {
        if (params->enable_missing) {
            print_sample(params, header, has_index, ".", fname);
        } else {
            print_error("samples", "no @RG:%s in \"%s\". Use option -m to enable missing.", params->tag, fname);
            status = EXIT_FAILURE;
            goto end_print;
        }
    } else if (count_samples > 1 && !params->enable_multiple) {
        print_error("samples", "multiple @RG:\"%s\" in \"%s\". Use option -d to enable multiple.", params->tag, fname);
        status = EXIT_FAILURE;
        goto end_print;
    } else {
        for (k = kh_begin(sample_set); k != kh_end(sample_set); ++k) {
            if (kh_exist(sample_set, k)) {
                char* sample = (char*)kh_key(sample_set, k);
                print_sample(params, header, has_index, sample, fname);
            }
        }
    }

end_print:
    for (k = kh_begin(sample_set); k != kh_end(sample_set); ++k) {
        if (kh_exist(sample_set, k)) {
            char* sample = (char*)kh_key(sample_set, k);
            if (kh_exist(sample_set, k)) free(sample);
        }
    }
    kh_destroy(sm, sample_set);
    if (header!=NULL) sam_hdr_destroy(header);
    if (in!=NULL) sam_close(in);

    return status;
}


int main_samples(int argc, char** argv) {
    int status = EXIT_SUCCESS;
    int print_header = 0;
    Params params;
    char* out_filename = NULL;
    strcpy(params.tag, "SM");
    params.faidx = NULL;
    params.enable_multiple = 0;
    params.enable_missing = 0;
    params.test_index =0;

    int opt;
    while ((opt = getopt(argc, argv,  "hHdmio:f:F:t:T:")) != -1) {
        switch (opt) {
        case 'H':
            print_header = 1;
            break;
        case 'o':
            out_filename = optarg;
            break;
        case 'd':
            params.enable_multiple = 1;
            break;
        case 'i':
            params.test_index = 1;
            break;
        case 'm':
            params.enable_missing = 1;
            break;
        case 'f':
            if (load_dictionary(&params, optarg) != EXIT_SUCCESS) {
                return EXIT_FAILURE;
            }
            break;
        case 'F':
            if (load_dictionaries(&params, optarg) != EXIT_SUCCESS) {
                return EXIT_FAILURE;
            }
            break;
        case 'T':
            if (strlen(optarg)!=2) {
                print_error("samples", "length ot a TAG must be 2 but got \"%s\".", optarg);
                return EXIT_FAILURE;
            }
            strcpy(params.tag, optarg);
            break;
        case 'h':
            usage_samples(stdout);
            return EXIT_SUCCESS;
        default:
            usage_samples(stderr);
            return EXIT_FAILURE;
        }
    }

    if(out_filename!=NULL) {
        params.out = fopen(out_filename, "w");
        if (params.out == NULL) {
            print_error_errno("samples", "cannot open \"%s\" for writing.", out_filename);
            return EXIT_FAILURE;
        }
    } else {
        params.out = stdout;
    }

    if (print_header) {
        fprintf(params.out, "#%s\tPATH", params.tag);
        if (params.test_index) fprintf(params.out, "\tINDEX");
        if (params.faidx != NULL) fprintf(params.out, "\tREFERENCE");
        fprintf(params.out, "\n");
    }

    /* input is stdin, each line contains the path to a bam file */
    if (argc == optind) {
        htsFile* fp = hts_open("-", "r");
        kstring_t  *ks = &fp->line;
        int ret;
        while ((ret = hts_getline(fp, KS_SEP_LINE, ks)) >= 0) {
            if (print_samples(&params, ks->s) != EXIT_SUCCESS) {
                status = EXIT_FAILURE;
                break;
            }
        }
        hts_close(fp);
    }
    else {
        /* loop over each bam file */
        int i;
        for (i = optind; i < argc; i++) {
            if (print_samples(&params, argv[i]) != EXIT_SUCCESS) {
                status = EXIT_FAILURE;
                break;
            }
        }
    }


    fflush(params.out);
    if (out_filename!=NULL)
        fclose(params.out);

    return status;
}
