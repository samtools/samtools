/*  bamtk.c -- main samtools command front-end.

    Copyright (C) 2008-2020 Genome Research Ltd.

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
DEALINGS IN THE SOFTWARE.  */

#include <config.h>

#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>

#include "htslib/hts.h"
#include "samtools.h"
#include "version.h"

int bam_taf2baf(int argc, char *argv[]);
int bam_mpileup(int argc, char *argv[]);
int bam_merge(int argc, char *argv[]);
int bam_index(int argc, char *argv[]);
int bam_sort(int argc, char *argv[]);
int bam_tview_main(int argc, char *argv[]);
int bam_mating(int argc, char *argv[]);
int bam_rmdup(int argc, char *argv[]);
int bam_flagstat(int argc, char *argv[]);
int bam_fillmd(int argc, char *argv[]);
int bam_idxstats(int argc, char *argv[]);
int bam_markdup(int argc, char *argv[]);
int main_samview(int argc, char *argv[]);
int main_reheader(int argc, char *argv[]);
int main_cut_target(int argc, char *argv[]);
int main_phase(int argc, char *argv[]);
int main_cat(int argc, char *argv[]);
int main_depth(int argc, char *argv[]);
int main_coverage(int argc, char *argv[]);
int main_bam2fq(int argc, char *argv[]);
int main_pad2unpad(int argc, char *argv[]);
int main_bedcov(int argc, char *argv[]);
int main_bamshuf(int argc, char *argv[]);
int main_stats(int argc, char *argv[]);
int main_flags(int argc, char *argv[]);
int main_split(int argc, char *argv[]);
int main_quickcheck(int argc, char *argv[]);
int main_addreplacerg(int argc, char *argv[]);
int faidx_main(int argc, char *argv[]);
int dict_main(int argc, char *argv[]);
int fqidx_main(int argc, char *argv[]);
int amplicon_clip_main(int argc, char *argv[]);
int main_ampliconstats(int argc, char *argv[]);

const char *samtools_version()
{
    return SAMTOOLS_VERSION;
}

static void usage(FILE *fp)
{
    /* Please improve the grouping */

    fprintf(fp,
"\n"
"Program: samtools (Tools for alignments in the SAM format)\n"
"Version: %s (using htslib %s)\n\n", samtools_version(), hts_version());
    fprintf(fp,
"Usage:   samtools <command> [options]\n"
"\n"
"Commands:\n"
"  -- Indexing\n"
"     dict           create a sequence dictionary file\n"
"     faidx          index/extract FASTA\n"
"     fqidx          index/extract FASTQ\n"
"     index          index alignment\n"
"\n"
"  -- Editing\n"
"     calmd          recalculate MD/NM tags and '=' bases\n"
"     fixmate        fix mate information\n"
"     reheader       replace BAM header\n"
"     targetcut      cut fosmid regions (for fosmid pool only)\n"
"     addreplacerg   adds or replaces RG tags\n"
"     markdup        mark duplicates\n"
"     ampliconclip   clip oligos from the end of reads\n"
"\n"
"  -- File operations\n"
"     collate        shuffle and group alignments by name\n"
"     cat            concatenate BAMs\n"
"     merge          merge sorted alignments\n"
"     mpileup        multi-way pileup\n"
"     sort           sort alignment file\n"
"     split          splits a file by read group\n"
"     quickcheck     quickly check if SAM/BAM/CRAM file appears intact\n"
"     fastq          converts a BAM to a FASTQ\n"
"     fasta          converts a BAM to a FASTA\n"
"\n"
"  -- Statistics\n"
"     bedcov         read depth per BED region\n"
"     coverage       alignment depth and percent coverage\n"
"     depth          compute the depth\n"
"     flagstat       simple stats\n"
"     idxstats       BAM index stats\n"
"     phase          phase heterozygotes\n"
"     stats          generate stats (former bamcheck)\n"
"     ampliconstats  generate amplicon specific stats\n"
"\n"
"  -- Viewing\n"
"     flags          explain BAM flags\n"
"     tview          text alignment viewer\n"
"     view           SAM<->BAM<->CRAM conversion\n"
"     depad          convert padded BAM to unpadded BAM\n"
"\n");
}

// This is a tricky one, but on Windows the filename wildcard expansion is done by
// the application and not by the shell, as traditionally it never had a "shell".
// Even now, DOS and Powershell do not do this expansion (but bash does).
//
// This means that Mingw/Msys implements code before main() that takes e.g. "*" and
// expands it up to a list of matching filenames.  This in turn breaks things like
// specifying "*" as a region (all the unmapped reads).  We take a hard line here -
// filename expansion is the task of the shell, not our application!
#ifdef _WIN32
int _CRT_glob = 0;
#endif

static void bam_import_err(void) {
    fprintf(stderr, "[main] \"samtools import\" has been removed. "
            "Please use \"samtools view\" instead.\n");
}

int main(int argc, char *argv[])
{
#ifdef _WIN32
    setmode(fileno(stdout), O_BINARY);
    setmode(fileno(stdin),  O_BINARY);
#endif
    if (argc < 2) { usage(stderr); return 1; }

    if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0) {
        if (argc == 2) { usage(stdout); return 0; }

        // Otherwise change "samtools help COMMAND [...]" to "samtools COMMAND";
        // main_xyz() functions by convention display the subcommand's usage
        // when invoked without any arguments.
        argv++;
        argc = 2;
    }

    int ret = 0;
    if (strcmp(argv[1], "view") == 0)           ret = main_samview(argc-1, argv+1);
    else if (strcmp(argv[1], "import") == 0)    { bam_import_err(); return 1; }
    else if (strcmp(argv[1], "mpileup") == 0)   ret = bam_mpileup(argc-1, argv+1);
    else if (strcmp(argv[1], "merge") == 0)     ret = bam_merge(argc-1, argv+1);
    else if (strcmp(argv[1], "sort") == 0)      ret = bam_sort(argc-1, argv+1);
    else if (strcmp(argv[1], "index") == 0)     ret = bam_index(argc-1, argv+1);
    else if (strcmp(argv[1], "idxstat") == 0 ||
             strcmp(argv[1], "idxstats") == 0)  ret = bam_idxstats(argc-1, argv+1);
    else if (strcmp(argv[1], "faidx") == 0)     ret = faidx_main(argc-1, argv+1);
    else if (strcmp(argv[1], "fqidx") == 0)     ret = fqidx_main(argc-1, argv+1);
    else if (strcmp(argv[1], "dict") == 0)      ret = dict_main(argc-1, argv+1);
    else if (strcmp(argv[1], "fixmate") == 0)   ret = bam_mating(argc-1, argv+1);
    else if (strcmp(argv[1], "rmdup") == 0)     ret = bam_rmdup(argc-1, argv+1);
    else if (strcmp(argv[1], "markdup") == 0)   ret = bam_markdup(argc-1, argv+1);
    else if (strcmp(argv[1], "ampliconclip") == 0) ret = amplicon_clip_main(argc-1, argv+1);
    else if (strcmp(argv[1], "flagstat") == 0 ||
             strcmp(argv[1], "flagstats") == 0) ret = bam_flagstat(argc-1, argv+1);
    else if (strcmp(argv[1], "calmd") == 0)     ret = bam_fillmd(argc-1, argv+1);
    else if (strcmp(argv[1], "fillmd") == 0)    ret = bam_fillmd(argc-1, argv+1);
    else if (strcmp(argv[1], "reheader") == 0)  ret = main_reheader(argc-1, argv+1);
    else if (strcmp(argv[1], "cat") == 0)       ret = main_cat(argc-1, argv+1);
    else if (strcmp(argv[1], "targetcut") == 0) ret = main_cut_target(argc-1, argv+1);
    else if (strcmp(argv[1], "phase") == 0)     ret = main_phase(argc-1, argv+1);
    else if (strcmp(argv[1], "depth") == 0)     ret = main_depth(argc-1, argv+1);
    else if (strcmp(argv[1], "coverage") == 0)  ret = main_coverage(argc-1, argv+1);
    else if (strcmp(argv[1], "bam2fq") == 0 ||
             strcmp(argv[1], "fastq") == 0 ||
             strcmp(argv[1], "fasta") == 0)     ret = main_bam2fq(argc-1, argv+1);
    else if (strcmp(argv[1], "pad2unpad") == 0) ret = main_pad2unpad(argc-1, argv+1);
    else if (strcmp(argv[1], "depad") == 0)     ret = main_pad2unpad(argc-1, argv+1);
    else if (strcmp(argv[1], "bedcov") == 0)    ret = main_bedcov(argc-1, argv+1);
    else if (strcmp(argv[1], "bamshuf") == 0)   ret = main_bamshuf(argc-1, argv+1);
    else if (strcmp(argv[1], "collate") == 0)   ret = main_bamshuf(argc-1, argv+1);
    else if (strcmp(argv[1], "stat") == 0 ||
             strcmp(argv[1], "stats") == 0)     ret = main_stats(argc-1, argv+1);
    else if (strcmp(argv[1], "flag") == 0 ||
             strcmp(argv[1], "flags") == 0)     ret = main_flags(argc-1, argv+1);
    else if (strcmp(argv[1], "split") == 0)     ret = main_split(argc-1, argv+1);
    else if (strcmp(argv[1], "quickcheck") == 0)  ret = main_quickcheck(argc-1, argv+1);
    else if (strcmp(argv[1], "addreplacerg") == 0) ret = main_addreplacerg(argc-1, argv+1);
    else if (strcmp(argv[1], "pileup") == 0) {
        fprintf(stderr, "[main] The `pileup' command has been removed. Please use `mpileup' instead.\n");
        return 1;
    }
    else if (strcmp(argv[1], "tview") == 0)   ret = bam_tview_main(argc-1, argv+1);
    else if (strcmp(argv[1], "ampliconstats") == 0)     ret = main_ampliconstats(argc-1, argv+1);
    else if (strcmp(argv[1], "--version") == 0) {
        printf(
"samtools %s\n"
"Using htslib %s\n"
"Copyright (C) 2020 Genome Research Ltd.\n",
               samtools_version(), hts_version());
    }
    else if (strcmp(argv[1], "--version-only") == 0) {
        printf("%s+htslib-%s\n", samtools_version(), hts_version());
    }
    else {
        fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
        return 1;
    }
    return ret;
}
