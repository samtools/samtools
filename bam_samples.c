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
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <htslib/kseq.h>
#include <samtools.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

KHASH_MAP_INIT_STR(sm, int)

/** and chained struct containing the faidx and the fasta filename
    will be compared with the @SQ lines in the SAM header
 */
typedef struct FaidxPath {
    /** path to reference */
    char* filename;
    /** fasta index  */
    faidx_t* faidx;
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
    /** show whether the bam is indexed */
    int test_index;
} Params;

/** print usage */
static void usage_samples(FILE *write_to) {
    fprintf(write_to,
            "Usage: samtools samples [options] <input> [...]\n"
            "       samtools samples [options] -X f1.bam f2.bam f1.bam.bai f2.bai \n"
            "       find dir1 dir2 -type f \\(-name \"*.bam\" -o -name \"*.cram\" \\) | samtools samples [options]\n"
            "       find dir1 dir2 -type f \\(-name \"*.bam\" -o -name \"*.bai\" \\) | sort | paste - - | samtools samples -X [options]\n"
            "\n"
            "Options:\n"
            "  -?              print help and exit\n"
            "  -h              add the columns header before printing the results\n"
            "  -i              test if the file is indexed.\n"
            "  -T <tag>        provide the sample tag name from the @RG line [SM].\n"
            "  -o <file>       output file [stdout].\n"
            "  -f <file.fa>    load an indexed fasta file in the collection of references. Can be used multiple times.\n"
            "  -F <file.txt>   read a file containing the paths to indexed fasta files. One path per line.\n"
            "  -X              use a custom index file.\n"
            "\n"
            " Using -f or -F will add a column containing the path to the reference or \".\" if the reference was not found.\n"
            "\n"
    );
}


/** loads fasta fai file into FaidxPath, add it to params->faidx */
static int load_dictionary(struct Params* params, const char* filename) {
    FaidxPath* head = params->faidx;
    FaidxPath* ptr = (FaidxPath*)malloc(sizeof(FaidxPath));
    if (ptr == NULL) {
        print_error_errno("samples", "Out of memory");
        return EXIT_FAILURE;
    }
    ptr->filename = strdup(filename);
    if (ptr->filename == NULL) {
        free(ptr);
        print_error_errno("samples", "Out of memory");
        return EXIT_FAILURE;
    }
    ptr->faidx = fai_load(filename);
    if (ptr->faidx == NULL) {
        free(ptr->filename);
        free(ptr);
        print_error_errno("samples", "Cannot load index from \"%s\"", filename);
        return EXIT_FAILURE;
    }
    /* insert at the beginning of the linked list */
    params->faidx = ptr;
    ptr->next = head;
    return EXIT_SUCCESS;
}

/** load a faidx file and append it to params */
static int load_dictionaries(Params* params, const char* filename) {
    int ret;
    htsFile* in;
    int status = EXIT_SUCCESS;

    in = hts_open(filename, "r");
    if (in == NULL) {
        print_error_errno("samples", "Cannot open \"%s\"", filename);
        status = EXIT_FAILURE;
    } else {
        kstring_t ks = KS_INITIALIZE;
        while ((ret = hts_getline(in, KS_SEP_LINE, &ks)) >= 0) {
            if (load_dictionary(params, ks_str(&ks)) != EXIT_SUCCESS) {
                status = EXIT_FAILURE;
                break;
            }
        }
        ks_free(&ks);
        hts_close(in);
    }
    return status;
}

/** print the sample information, search for a reference */
static int print_sample(
        Params* params,
        sam_hdr_t *header,
        int has_index,
        const char* sample,
        const char* fname) {
    fputs(sample, params->out);
    fputc('\t', params->out);
    fputs(fname, params->out);
    if (params->test_index) {
        fprintf(params->out, "\t%c", has_index ? 'Y' : 'N');
    }
    if (params->faidx != NULL) {
        FaidxPath* ref = NULL;
        FaidxPath* curr = params->faidx;
        while (curr != NULL) {
            /** check names and length are the same in the same order */
            if (faidx_nseq(curr->faidx) == header->n_targets) {
                int i;
                for (i = 0; i < faidx_nseq(curr->faidx); i++) {
                    /** check name is the same */
                    if (strcmp(faidx_iseq(curr->faidx, i), header->target_name[i]) != 0) break;
                    /** check length is the same */
                    if (faidx_seq_len(curr->faidx, faidx_iseq(curr->faidx, i)) != header->target_len[i]) break;
                }
                /* the ref was found */
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
static int print_samples(Params* params, const char* fname, const char* baifname) {
    samFile *in = 0;
    sam_hdr_t *header = NULL;
    int n_rg;
    int status = EXIT_SUCCESS;
    khash_t(sm) *sample_set = NULL;
    khint_t k;
    int count_samples = 0;
    int has_index = 0;

    if ((sample_set = kh_init(sm)) == NULL) {
        print_error("samples", "Failed to initialise sample hash");
        status = EXIT_FAILURE;
        goto end_print;
    }

    if ((in = sam_open_format(fname, "r", NULL)) == 0) {
        print_error_errno("samples", "Failed to open \"%s\" for reading", fname);
        status = EXIT_FAILURE;
        goto end_print;
    }
    if ((header = sam_hdr_read(in)) == 0) {
        print_error("samples", "Failed to read the header from \"%s\"", fname);
        status = EXIT_FAILURE;
        goto end_print;
    }

    /* try to load index if required */
    if (params->test_index) {
        hts_idx_t *bam_idx;
        /* path to bam index was specified */
        if (baifname != NULL) {
            bam_idx = sam_index_load3(in, fname, baifname, HTS_IDX_SILENT_FAIL);
        }
        /* get default index */
        else {
            bam_idx = sam_index_load3(in, fname, NULL, HTS_IDX_SILENT_FAIL);
        }
        has_index = bam_idx != NULL;
        if (bam_idx != NULL) hts_idx_destroy(bam_idx);
        /* and we continue... we have tested the index file but we always test for the samples and the references */
    }

    /* get the RG lines */
    n_rg = sam_hdr_count_lines(header, "RG");
    if (n_rg > 0) {
        int i, r, ret;
        char* sample;
        kstring_t sm_val = KS_INITIALIZE;
        for (i = 0; i < n_rg; i++) {
            r = sam_hdr_find_tag_pos(header, "RG", i, params->tag, &sm_val);
            if (r < 0) continue;
            k = kh_get(sm, sample_set, ks_str(&sm_val));
            if (k != kh_end(sample_set)) continue;
            sample = strdup(ks_str(&sm_val));
            if (sample == NULL) {
                print_error_errno("samples", "Out of memory");
                status = EXIT_FAILURE;
                goto end_print;
            }
            kh_put(sm, sample_set, sample, &ret);
            if (ret < 0) {
                print_error("samples", "Failed to insert key '%s' into sample_set", sample);
                free(sample);
                status = EXIT_FAILURE;
                goto end_print;
            }
            ++count_samples;
        }
        ks_free(&sm_val);
    }
    if (count_samples == 0) {
        print_sample(params, header, has_index, ".", fname);
    } else {
        for (k = kh_begin(sample_set); k != kh_end(sample_set); ++k) {
            if (kh_exist(sample_set, k)) {
                char* sample = (char*)kh_key(sample_set, k);
                print_sample(params, header, has_index, sample, fname);
            }
        }
    }

end_print:
    if (sample_set != NULL) {
        for (k = kh_begin(sample_set); k != kh_end(sample_set); ++k) {
            if (kh_exist(sample_set, k)) {
                char* sample = (char*)kh_key(sample_set, k);
                free(sample);
            }
        }
        kh_destroy(sm, sample_set);
    }
    if (header != NULL) sam_hdr_destroy(header);
    if (in != NULL) sam_close(in);

    return status;
}


int main_samples(int argc, char** argv) {
    int status = EXIT_SUCCESS;
    int print_header = 0;
    int has_index_file = 0;
    Params params;
    char* out_filename = NULL;
    FaidxPath* fai;

    strcpy(params.tag, "SM");
    params.faidx = NULL;
    params.test_index =0;

    int opt;
    while ((opt = getopt_long(argc, argv,  "?hiXo:f:F:T:", NULL, NULL)) != -1) {
        switch (opt) {
        case 'h':
            print_header = 1;
            break;
        case 'o':
            out_filename = optarg;
            break;
        case 'i':
            params.test_index = 1;
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
            if (strlen(optarg) != 2) {
                print_error("samples", "Length of tag \"%s\" is not 2.", optarg);
                return EXIT_FAILURE;
            }
            strcpy(params.tag, optarg);
            break;
        case '?':
            usage_samples(stdout);
            return EXIT_SUCCESS;
        case 'X':
            has_index_file = 1;
            break;
        default:
            usage_samples(stderr);
            return EXIT_FAILURE;
        }
    }

    /* if no file was provided and input is the terminal, print the usage and exit */
    if (argc == optind && isatty(STDIN_FILENO)) {
       usage_samples(stderr);
       return EXIT_FAILURE;
    }

    if (out_filename != NULL) {
        params.out = fopen(out_filename, "w");
        if (params.out == NULL) {
            print_error_errno("samples", "Cannot open \"%s\" for writing", out_filename);
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

    /* no file was provided, input is stdin, each line contains the path to a bam file */
    if (argc == optind) {
        htsFile* fp = hts_open("-", "r");
        if (fp == NULL) {
            print_error_errno("samples", "Cannot read from stdin");
            status = EXIT_FAILURE;
        } else {
            kstring_t ks = KS_INITIALIZE;
            int ret;
            while ((ret = hts_getline(fp, KS_SEP_LINE, &ks)) >= 0) {
                char* bai_path = NULL;
                if (has_index_file) {
                    /* bam path and bam index file are separated by a tab */
                    char* tab = strchr(ks_str(&ks), '\t');
                    if (tab == NULL || *(tab+1) == '\0') {
                        print_error_errno("samples", "Expected path-to-bam(tab)path-to-index but got \"%s\"", ks_str(&ks));
                        status = EXIT_FAILURE;
                        break;
                    }
                    *tab=0;
                    bai_path = (tab + 1);
                }
                if (print_samples(&params, ks_str(&ks), bai_path) != EXIT_SUCCESS) {
                    status = EXIT_FAILURE;
                    break;
                }
            }
            ks_free(&ks);
            hts_close(fp);
        }
    }
    /* loop over each file in argc/argv bam index provided */
    else if (has_index_file) {
        /* calculate number of input BAM files */
        if ((argc - optind) % 2 != 0) {
            print_error("samples","Odd number of filenames detected! Each BAM file should have an index file");
            status = EXIT_FAILURE;
        } else {
            int i;
            int n = (argc - optind ) / 2;
            for (i = 0; i < n; i++) {
                if (print_samples(&params, argv[optind+i], argv[optind+i+n]) != EXIT_SUCCESS) {
                    status = EXIT_FAILURE;
                    break;
                }
            }
        }
    } else {
        int i;
        for (i = optind; i < argc; i++) {
            if (print_samples(&params, argv[i], NULL) != EXIT_SUCCESS) {
                status = EXIT_FAILURE;
                break;
            }
        }
    }

    fai = params.faidx;
    while (fai != NULL) {
        FaidxPath* next = fai -> next;
        free(fai->filename);
        fai_destroy(fai->faidx);
        free(fai);
        fai = next;
    }

    if (fflush(params.out) != 0) {
        print_error_errno("samples", "Cannot flush output");
        status = EXIT_FAILURE;
    }
    if (out_filename != NULL) {
        fclose(params.out);
    }

    return status;
}
