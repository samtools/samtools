/* bam_import -- Import of FASTQ files.
 *
 *   samtools import -1 a_1.fq -2 a_2.fq --i1 a_i1.fq --i2 a_i2.fq
 *   samtools import a_1.fq a_2.fq
 *   samtools import a_interleaved.fq
 *
 * Copyright (C) 2020-2021 Genome Research Ltd.
 *
 * Author: James Bonfield <jkb@sanger.ac.uk>
 */

/*
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
DEALINGS IN THE SOFTWARE.
*/

// TODO: Store other non-aux comments; in new sam tag?

#include <config.h>
#include <ctype.h>

#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#include "samtools.h"
#include "sam_opts.h"

static int usage(FILE *fp, int exit_status) {
    fprintf(fp, "Usage: samtools import [options] [file.fastq ...]\n");
    fprintf(fp, "\n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "  -s FILE      Read paired-ended data from single FILE\n");
    fprintf(fp, "  -0 FILE      Read single-ended data from FILE\n");
    fprintf(fp, "  -1 FILE      Read-1 from FILE\n");
    fprintf(fp, "  -2 FILE      Read-2 from FILE\n");
    fprintf(fp, "  --i1 FILE    Index-1 from FILE\n");
    fprintf(fp, "  --i2 FILE    Index-2 from FILE\n");
    fprintf(fp, "  -i           Parse CASAVA identifier\n");
    fprintf(fp, "  --barcode-tag TAG\n");
    fprintf(fp, "               Tag to use with barcode sequences [BC]\n");
    fprintf(fp, "  --quality-tag TAG\n");
    fprintf(fp, "               Tag to use with barcode qualities [QT]\n");
    fprintf(fp, "  -N, --name2  Use 2nd field as read name (SRA format)\n");
    fprintf(fp, "  -r STRING    Build up a complete @RG line\n");
    fprintf(fp, "  -R STRING    Add a simple RG line of \"@RG\\tID:STRING\"\n");
    fprintf(fp, "  -T TAGLIST   Parse tags in SAM format; list of '*' for all\n");
    fprintf(fp, "  -o FILE      Output to FILE instead of stdout\n");
    fprintf(fp, "  -u           Uncompressed output\n");
    fprintf(fp, "  --order TAG  Store Nth record count in TAG\n");
    fprintf(fp, "\n");
    sam_global_opt_help(fp, "-.O.-@--");

    fprintf(fp, "\nA single fastq file will be interpreted as -s, -0 or -1 depending on\n");
    fprintf(fp, "file contents, and a pair of fastq files as \"-1 FILE1 -2 FILE2\".\n");

    return exit_status;
}

// Order matters here as we want to read index elements before main
// sequences so on reading the seqs we can emit a fully annotated record.
enum fileno {
    FQ_I1, FQ_I2, // index seqs for R1 and R2
    FQ_R0,        // single file and unpaired data (singled-ended tech).
    FQ_R1, FQ_R2, // separate read1 and read2 files
    FQ_SINGLE,    // single file, but with read1 and/or read2 present.
    FQ_END
};

typedef struct {
    sam_global_args ga;
    int no_pg;
    char *fn[FQ_END], *fn_out;
    int idx_both;      // add index to READ2 too, not just READ1
    int casava;
    char *barcode_seq;
    char *barcode_qual;
    char *aux;
    char *rg;
    char *rg_line;
    char *order;
    int compress_level;
    htsThreadPool p;
    int name2;
} opts_t;

// Append a sequence and quality string from a BAM record to a BC:Z and
// QT:Z style aux tag string.
static int append_index(kstring_t *s, kstring_t *q, bam1_t *b) {
    char *sp, *qp;
    if (ks_resize(s, s->l + b->core.l_qseq+1 +1) < 0)
        return -1;
    if (ks_resize(q, q->l + b->core.l_qseq+1 +1) < 0)
        return -1;

    sp = s->s + s->l - (s->l > 0);
    qp = q->s + q->l - (q->l > 0);

    if (s->l)
        *sp++ = '-';

    if (q->l)
        *qp++ = ' ';

    int i;
    uint8_t *seq = bam_get_seq(b);
    uint8_t *qual = bam_get_qual(b);
    for (i = 0; i < b->core.l_qseq; i++) {
        *sp++ = seq_nt16_str[bam_seqi(seq, i)];
        *qp++ = qual[i] + '!';
    }
    *sp++ = 0;
    *qp++ = 0;

    s->l = sp - s->s;
    q->l = qp - q->s;

    return 0;
}

static int import_fastq(int argc, char **argv, opts_t *opts) {
    int i, n, ret = 0;
    samFile *fp_in[FQ_END] = {NULL};
    bam1_t *b = bam_init1();
    int ids[FQ_END];
    samFile *fp_out = NULL;
    sam_hdr_t *hdr_out = NULL;
    kstring_t index_str = {0,0};
    kstring_t read_str = {0,0};
    char *rg = opts->rg;
    kstring_t rg_line = {0,0};
    uint64_t read_num = 0;
    kstring_t idx_seq  = {0};
    kstring_t idx_qual = {0};

    // Any additional arguments are assumed to be r1 r2, as a
    // short cut. We support reading index tags out of those too (eg
    // Illumina CASAVA format), but if we do that we lack the barcode
    // quality string.
    //
    // We also consider a read name ending in /1 or /2 to be a single
    // file containing interleaved fastq records for both ends.
    // These will be labeled as fn[FQ_R1] but adjusted during reading.
    if (argc == 1)
        opts->fn[FQ_SINGLE] = argv[0];
    else
        for (i = 0; i < 4; i++)
            if (argc > i)
                opts->fn[FQ_R1+i] = argv[i];

    // Open all files
    for (i = n = 0; i < FQ_END; i++) {
        if (!opts->fn[i])
            continue;
        fp_in[i] = sam_open_format(opts->fn[i], "r", &opts->ga.in);
        if (!fp_in[i]) {
            perror(opts->fn[i]);
            ret = -1;
            goto err;
        }
        if (opts->p.pool)
            hts_set_thread_pool(fp_in[i], &opts->p);
        ids[n++] = i;

        if (opts->name2)
            hts_set_opt(fp_in[i], FASTQ_OPT_NAME2, 1);
        if (opts->casava)
            hts_set_opt(fp_in[i], FASTQ_OPT_CASAVA, 1);
        if (opts->barcode_seq) // for auto-CASAVA parsing
            hts_set_opt(fp_in[i], FASTQ_OPT_BARCODE, opts->barcode_seq);
        if (opts->aux)
            hts_set_opt(fp_in[i], FASTQ_OPT_AUX,
                        *opts->aux == '*' || *opts->aux == '\0'
                        ? NULL : opts->aux);

        switch (i) {
        case FQ_I1:
            kputs("--i1 I1.fastq ", &read_str);
            kputs("i*", &index_str);
            break;
        case FQ_I2:
            kputs("--i2 I2.fastq ", &read_str);
            kputs("i*", &index_str);
            break;

        case FQ_R0:
            kputs("-0 unpaired.fastq ", &read_str);
            break;

        case FQ_R1:
            kputs("-1 R1.fastq ", &read_str);
            break;

        case FQ_R2:
            kputs("-2 R2.fastq ", &read_str);
            break;

        case FQ_SINGLE:
            kputs("-N -o paired.fastq ", &read_str);
            break;

        default:
            ks_clear(&read_str); // not reversible
            kputs("", &read_str);
        }
    }
    if (n == 0) {
        bam_destroy1(b);
        return usage(stdout, EXIT_SUCCESS);
    }

    char out_mode[10] = {'w', 0, 0};
    if (opts->compress_level != -1)
        out_mode[1] = '0' + opts->compress_level;
    sam_open_mode(out_mode+strlen(out_mode), opts->fn_out, NULL);
    fp_out = sam_open_format(opts->fn_out, out_mode, &opts->ga.out);
    if (!fp_out) {
        perror(opts->fn_out);
        goto err;
    }
    autoflush_if_stdout(fp_out, opts->fn_out);
    if (opts->p.pool)
        hts_set_thread_pool(fp_out, &opts->p);

    // Create header
    if (ks_len(&read_str)) {
        char CO[2100];
        if (ks_len(&index_str))
            snprintf(CO, sizeof(CO), "@CO\tReverse with: samtools fastq %s "
                    "--index-format=\"%s\"\n",
                    ks_str(&read_str), ks_str(&index_str));
        else
            snprintf(CO, sizeof(CO), "@CO\tReverse with: samtools fastq %s\n",
                    ks_str(&read_str));

        hdr_out = sam_hdr_parse(strlen(CO), CO);
    } else {
        hdr_out = sam_hdr_init();
    }

    // Add a version line with the sort order to the output header
    if (sam_hdr_add_line(hdr_out, "HD", "VN", SAM_FORMAT_VERSION, "SO", "unsorted", "GO", "query", NULL) < 0) {
        fprintf(stderr, "Could not set SO and GO in the header.\n");
        goto err;
    }

    // Read group
    if (opts->rg_line) {
        if (*opts->rg_line != '@')
            ksprintf(&rg_line, "@RG\t%s", opts->rg_line);
        else
            kputs(opts->rg_line, &rg_line);
    } else if (opts->rg) {
        ksprintf(&rg_line, "@RG\tID:%s", opts->rg);
    }

    if (ks_len(&rg_line)) {
        if (sam_hdr_add_lines(hdr_out, ks_str(&rg_line), 0) < 0)
            goto err;
        rg = strstr(ks_str(&rg_line), "\tID:");
        if (!rg) {
            fprintf(stderr, "\"-r RG-LINE\" option contained no ID field\n");
            goto err;
        }
        rg += 4;

        i = 0;
        while (rg[i] != '\t' && rg[i] != '\0')
            i++;
        rg[i] = 0;
    }

    if ((ret = sam_hdr_write(fp_out, hdr_out)) < 0)
        goto err;


    // Interleave / combine from n files (ids[0..n-1]).
    int res;
    int eof = 0;
    do {
        idx_seq.l = idx_qual.l = 0;
        for (i = 0; i < n; i++) {
            if ((res = sam_read1(fp_in[ids[i]], NULL, b)) < 0) {
                if (res == -1) {
                    eof++;
                    continue;
                } else
                    break;
            }

            // index
            if (ids[i] == FQ_I1 || ids[i] == FQ_I2) {
                if (append_index(&idx_seq, &idx_qual, b) < 0) {
                    res = -1;
                    break;
                }
                continue;
            }

            // full read
            if (idx_seq.l) {
                if (opts->idx_both || ids[i] == FQ_SINGLE ||
                    ids[i] == FQ_R0 || ids[i] == FQ_R1) {
                    if (bam_aux_append(b, opts->barcode_seq, 'Z', idx_seq.l,
                                       (uint8_t *)idx_seq.s) ||
                        bam_aux_append(b, opts->barcode_qual, 'Z', idx_qual.l,
                                       (uint8_t *)idx_qual.s)) {
                        res = -1;
                        break;
                    }
                }
            }

            switch(ids[i]) {
            case FQ_R0:
                // unpaired; no flags to declare
                break;
            case FQ_SINGLE:
                // paired (but don't know if R1 or R2) or unpaired.
                // We rely on the /1 and /2 read suffix parsing in htslib
                // to distinguish the two cases, or CASAVA tags if
                // explicitly enabled.
                break;
            case FQ_R1:
                if ((b->core.flag & (BAM_FREAD1 | BAM_FREAD2)) == 0)
                    b->core.flag |= BAM_FREAD1;
                b->core.flag |= BAM_FPAIRED;
                if (i+1 < n && ids[i+1] == FQ_R2)
                    b->core.flag |= BAM_FMUNMAP;
                break;
            case FQ_R2:
                b->core.flag |= BAM_FPAIRED | BAM_FREAD2;
                if (i > 0 && ids[i-1] == FQ_R1)
                    b->core.flag |= BAM_FMUNMAP;
                break;
            }

            if (rg) {
                if (bam_aux_append(b, "RG", 'Z', strlen(rg)+1,
                                   (uint8_t *)rg) < 0) {
                    ret = -1;
                    goto err;
                }
            }

            if (opts->order) {
                if (bam_aux_update_int(b, opts->order, read_num++) < 0) {
                    ret = -1;
                    goto err;
                }
            }

            res = sam_write1(fp_out, hdr_out, b);
        }
    } while (res >= 0);

    if (res != -1) {
        print_error("import", "truncated file. Aborting");
        ret = res;
        goto err;
    }

    if (eof != n) {
        print_error("import", "input files with differing number of records");
        ret = -1;
        goto err;
    }

    // Close and return
    ret = 0;
err:
    bam_destroy1(b);
    sam_hdr_destroy(hdr_out);
    ks_free(&rg_line);
    ks_free(&index_str);
    ks_free(&read_str);
    if (fp_out) {
        release_autoflush(fp_out);
        if (sam_close(fp_out) < 0) {
            perror(opts->fn_out);
            ret |= -1;
        }
    }
    for (i = 0; i < FQ_END; i++) {
        if (fp_in[i] && sam_close(fp_in[i]) < 0) {
            perror(opts->fn[i]);
            ret |= -1;
        }
    }
    ks_free(&idx_seq);
    ks_free(&idx_qual);

    return ret;
}

int main_import(int argc, char *argv[]) {
    int c;
    opts_t opts = {
        .no_pg = 0,
        .ga = SAM_GLOBAL_ARGS_INIT,
        .fn = {NULL},
        .fn_out = "-",
        .casava = 0,
        .barcode_seq = "BC",
        .barcode_qual = "QT",
        .aux = NULL,
        .rg = NULL,
        .rg_line = NULL,
        .order = NULL,
        .compress_level = -1,
        .name2 = 0,
    };
    kstring_t rg = {0};

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, '-', '@'),
        {"no-PG", no_argument, NULL, 9},
        {"i1", required_argument, NULL, 1},
        {"i2", required_argument, NULL, 2},
        {"r1", required_argument, NULL, '1'},
        {"r2", required_argument, NULL, '2'},
        {"rg", required_argument, NULL, 'R'},
        {"rg-line", required_argument, NULL, 'r'},
        {"order", required_argument, NULL, 3},
        {"barcode-tag", required_argument, NULL, 4},
        {"quality-tag", required_argument, NULL, 5},
        {"name2", no_argument, NULL, 'N'},
        { NULL, 0, NULL, 0 }
    };

    while ((c = getopt_long(argc, argv, "1:2:s:0:bhiT:r:R:o:O:u@:N", lopts, NULL)) >= 0) {
        switch (c) {
        case 'b': opts.idx_both = 1; break;
        case '0': opts.fn[FQ_R0] = optarg; break;
        case '1': opts.fn[FQ_R1] = optarg; break;
        case '2': opts.fn[FQ_R2] = optarg; break;
        case  1:  opts.fn[FQ_I1] = optarg; break;
        case  2:  opts.fn[FQ_I2] = optarg; break;
        case 's': opts.fn[FQ_SINGLE] = optarg; break;
        case 'o': opts.fn_out = optarg; break;
        case 'i': opts.casava = 1; break;
        case  4:  opts.barcode_seq = optarg; break;
        case  5:  opts.barcode_qual = optarg; break;
        case 'T': opts.aux = optarg; break;
        case 'u': opts.compress_level = 0; break;
        case 'R': opts.rg = optarg; break;
        case 'r':
            if (*optarg != '@' && ks_len(&rg) == 0)
                kputs("@RG", &rg);
            if (ks_len(&rg))
                kputc_('\t', &rg);
            kputs(optarg, &rg);
            opts.rg_line = rg.s;
            break;

        case 'N': opts.name2 = 1; break;

        case 9: opts.no_pg = 1; break;
        case 3: opts.order = optarg; break;

        case 'h': return usage(stdout, EXIT_SUCCESS);
        case '?': return usage(stderr, EXIT_FAILURE);

        default:
            if (parse_sam_global_opt(c, optarg, lopts, &opts.ga) != 0)
                return usage(stderr, EXIT_FAILURE);
            break;
        }
    }

    if (opts.ga.nthreads > 0) {
        if (!(opts.p.pool = hts_tpool_init(opts.ga.nthreads))) {
            fprintf(stderr, "Failed to create thread pool\n");
            if (rg.s)
                free(rg.s);
            return -1;;
        }
    }

    int ret = import_fastq(argc-optind, argv+optind, &opts) ? 1 : 0;

    if (rg.s)
        free(rg.s);

    if (opts.p.pool)
        hts_tpool_destroy(opts.p.pool);

    return ret;
}
