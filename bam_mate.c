/*  bam_mate.c -- fix mate pairing information and clean up flags.

    Copyright (C) 2009, 2011-2017 Genome Research Ltd.
    Portions copyright (C) 2011 Broad Institute.
    Portions copyright (C) 2012 Peter Cock, The James Hutton Institute.

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

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "htslib/thread_pool.h"
#include "sam_opts.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "samtools.h"

/*
 * This function calculates ct tag for two bams, it assumes they are from the same template and
 * writes the tag to the first read in position terms.
 */
static void bam_template_cigar(bam1_t *b1, bam1_t *b2, kstring_t *str)
{
    bam1_t *swap;
    int i, end;
    uint32_t *cigar;
    str->l = 0;
    if (b1->core.tid != b2->core.tid || b1->core.tid < 0 || b1->core.pos < 0 || b2->core.pos < 0 || b1->core.flag&BAM_FUNMAP || b2->core.flag&BAM_FUNMAP) return; // coordinateless or not on the same chr; skip
    if (b1->core.pos > b2->core.pos) swap = b1, b1 = b2, b2 = swap; // make sure b1 has a smaller coordinate
    kputc((b1->core.flag & BAM_FREAD1)? '1' : '2', str); // segment index
    kputc((b1->core.flag & BAM_FREVERSE)? 'R' : 'F', str); // strand
    for (i = 0, cigar = bam_get_cigar(b1); i < b1->core.n_cigar; ++i) {
        kputw(bam_cigar_oplen(cigar[i]), str);
        kputc(bam_cigar_opchr(cigar[i]), str);
    }
    end = bam_endpos(b1);
    kputw(b2->core.pos - end, str);
    kputc('T', str);
    kputc((b2->core.flag & BAM_FREAD1)? '1' : '2', str); // segment index
    kputc((b2->core.flag & BAM_FREVERSE)? 'R' : 'F', str); // strand
    for (i = 0, cigar = bam_get_cigar(b2); i < b2->core.n_cigar; ++i) {
        kputw(bam_cigar_oplen(cigar[i]), str);
        kputc(bam_cigar_opchr(cigar[i]), str);
    }

    uint8_t* data;
    if ((data = bam_aux_get(b1,"ct")) != NULL) bam_aux_del(b1, data);
    if ((data = bam_aux_get(b2,"ct")) != NULL) bam_aux_del(b2, data);

    bam_aux_append(b1, "ct", 'Z', str->l+1, (uint8_t*)str->s);
}

/*
 * What This Program is Supposed To Do:
 * Fill in mate coordinates, ISIZE and mate related flags from a name-sorted
 * alignment.
 *
 * How We Handle Input
 *
 * Secondary and supplementary Reads:
 * -write to output unchanged
 * All Reads:
 * -if pos == 0 (1 based), tid == -1 set UNMAPPED flag
 * single Reads:
 * -if pos == 0 (1 based), tid == -1, or UNMAPPED then set UNMAPPED, pos = 0,
 *  tid = -1
 * -clear bad flags (PAIRED, MREVERSE, PROPER_PAIR)
 * -set mpos = 0 (1 based), mtid = -1 and isize = 0
 * -write to output
 * Paired Reads:
 * -if read is unmapped and mate is not, set pos and tid to equal that of mate
 * -sync mate flags (MREVERSE, MUNMAPPED), mpos, mtid
 * -recalculate ISIZE if possible, otherwise set it to 0
 * -optionally clear PROPER_PAIR flag from reads where mapping or orientation
 *  indicate this is not possible (Illumina orientation only)
 * -calculate ct and apply to lowest positioned read
 * -write to output
 * Limitations
 * -Does not handle tandem reads
 * -Should mark supplementary reads the same as primary.
 * Notes
 * -CT definition appears to be something else in spec, this was in here before
 *  I started tampering with it, anyone know what is going on here? To work
 *  around this I have demoted the CT this tool generates to ct.
 */

static void sync_unmapped_pos_inner(bam1_t* src, bam1_t* dest) {
    if ((dest->core.flag & BAM_FUNMAP) && !(src->core.flag & BAM_FUNMAP)) {
        // Set unmapped read's RNAME and POS to those of its mapped mate
        // (recommended best practice, ensures if coord sort will be together)
        dest->core.tid = src->core.tid;
        dest->core.pos = src->core.pos;
    }
}

static void sync_mate_inner(bam1_t* src, bam1_t* dest)
{
    // sync mate pos information
    dest->core.mtid = src->core.tid; dest->core.mpos = src->core.pos;
    // sync flag info
    if (src->core.flag&BAM_FREVERSE)
        dest->core.flag |= BAM_FMREVERSE;
    else
        dest->core.flag &= ~BAM_FMREVERSE;
    if (src->core.flag & BAM_FUNMAP) {
        dest->core.flag |= BAM_FMUNMAP;
    }
}

// Is it plausible that these reads are properly paired?
// Can't really give definitive answer without checking isize
static bool plausibly_properly_paired(bam1_t* a, bam1_t* b)
{
    if ((a->core.flag & BAM_FUNMAP) || (b->core.flag & BAM_FUNMAP)) return false;
    assert(a->core.tid >= 0); // This should never happen if FUNMAP is set correctly

    if (a->core.tid != b->core.tid) return false;

    bam1_t* first = a;
    bam1_t* second = b;
    int32_t a_pos = a->core.flag&BAM_FREVERSE ? bam_endpos(a) : a->core.pos;
    int32_t b_pos = b->core.flag&BAM_FREVERSE ? bam_endpos(b) : b->core.pos;
    if (a_pos > b_pos) {
        first = b;
        second = a;
    } else {
        first = a;
        second = b;
    }

    if (!(first->core.flag&BAM_FREVERSE) && (second->core.flag&BAM_FREVERSE))
        return true;
    else
        return false;
}

// Returns 0 on success, -1 on failure.
static int bam_format_cigar(const bam1_t* b, kstring_t* str)
{
    // An empty cigar is a special case return "*" rather than ""
    if (b->core.n_cigar == 0) {
        return (kputc('*', str) == EOF) ? -1 : 0;
    }

    const uint32_t *cigar = bam_get_cigar(b);
    uint32_t i;

    for (i = 0; i < b->core.n_cigar; ++i) {
        if (kputw(bam_cigar_oplen(cigar[i]), str) == EOF) return -1;
        if (kputc(bam_cigar_opchr(cigar[i]), str) == EOF) return -1;
    }

    return 0;
}

// Returns 0 on success, -1 on failure.
static int sync_mq_mc(bam1_t* src, bam1_t* dest)
{
    if ( (src->core.flag & BAM_FUNMAP) == 0 ) { // If mapped
        // Copy Mate Mapping Quality
        uint32_t mq = src->core.qual;
        uint8_t* data;
        if ((data = bam_aux_get(dest,"MQ")) != NULL) {
            bam_aux_del(dest, data);
        }

        bam_aux_append(dest, "MQ", 'i', sizeof(uint32_t), (uint8_t*)&mq);
    }
    // Copy mate cigar if either read is mapped
    if ( (src->core.flag & BAM_FUNMAP) == 0 || (dest->core.flag & BAM_FUNMAP) == 0 ) {
        uint8_t* data_mc;
        if ((data_mc = bam_aux_get(dest,"MC")) != NULL) {
            bam_aux_del(dest, data_mc);
        }

        // Convert cigar to string
        kstring_t mc = { 0, 0, NULL };
        if (bam_format_cigar(src, &mc) < 0) return -1;

        bam_aux_append(dest, "MC", 'Z', ks_len(&mc)+1, (uint8_t*)ks_str(&mc));
        free(mc.s);
    }
    return 0;
}

// Copy flags.
// Returns 0 on success, -1 on failure.
static int sync_mate(bam1_t* a, bam1_t* b)
{
    sync_unmapped_pos_inner(a,b);
    sync_unmapped_pos_inner(b,a);
    sync_mate_inner(a,b);
    sync_mate_inner(b,a);
    if (sync_mq_mc(a,b) < 0) return -1;
    if (sync_mq_mc(b,a) < 0) return -1;
    return 0;
}


static uint32_t calc_mate_score(bam1_t *b)
{
    uint32_t score = 0;
    uint8_t  *qual = bam_get_qual(b);
    int i;

    for (i = 0; i < b->core.l_qseq; i++) {
        if (qual[i] >= 15) score += qual[i];
    }

    return score;
}


static int add_mate_score(bam1_t *src, bam1_t *dest)
{
    uint8_t *data_ms;
    uint32_t mate_score = calc_mate_score(src);

    if ((data_ms = bam_aux_get(dest, "ms")) != NULL) {
        bam_aux_del(dest, data_ms);
    }

    if (bam_aux_append(dest, "ms", 'i', sizeof(uint32_t), (uint8_t*)&mate_score) == -1) {
        return -1;
    }

    return 0;
}

// currently, this function ONLY works if each read has one hit
static int bam_mating_core(samFile *in, samFile *out, int remove_reads, int proper_pair_check, int add_ct, int do_mate_scoring)
{
    bam_hdr_t *header;
    bam1_t *b[2] = { NULL, NULL };
    int curr, has_prev, pre_end = 0, cur_end = 0;
    kstring_t str;

    str.l = str.m = 0; str.s = 0;
    header = sam_hdr_read(in);
    if (header == NULL) {
        fprintf(stderr, "[bam_mating_core] ERROR: Couldn't read header\n");
        return 1;
    }
    // Accept unknown, unsorted, or queryname sort order, but error on coordinate sorted.
    if ((header->l_text > 3) && (strncmp(header->text, "@HD", 3) == 0)) {
        char *p, *q;
        p = strstr(header->text, "\tSO:coordinate");
        q = strchr(header->text, '\n');
        // Looking for SO:coordinate within the @HD line only
        // (e.g. must ignore in a @CO comment line later in header)
        if ((p != 0) && (p < q)) {
            fprintf(stderr, "[bam_mating_core] ERROR: Coordinate sorted, require grouped/sorted by queryname.\n");
            goto fail;
        }
    }
    if (sam_hdr_write(out, header) < 0) goto write_fail;

    b[0] = bam_init1();
    b[1] = bam_init1();
    curr = 0; has_prev = 0;
    while (sam_read1(in, header, b[curr]) >= 0) {
        bam1_t *cur = b[curr], *pre = b[1-curr];
        if (cur->core.flag & BAM_FSECONDARY)
        {
            if ( !remove_reads ) {
                if (sam_write1(out, header, cur) < 0) goto write_fail;
            }
            continue; // skip secondary alignments
        }
        if (cur->core.flag & BAM_FSUPPLEMENTARY)
        {
            if (sam_write1(out, header, cur) < 0) goto write_fail;
            continue; // pass supplementary alignments through unchanged (TODO:make them match read they came from)
        }
        if (cur->core.tid < 0 || cur->core.pos < 0) // If unmapped set the flag
        {
            cur->core.flag |= BAM_FUNMAP;
        }
        if ((cur->core.flag&BAM_FUNMAP) == 0) // If mapped calculate end
        {
            cur_end = bam_endpos(cur);

            // Check cur_end isn't past the end of the contig we're on, if it is set the UNMAP'd flag
            if (cur_end > (int)header->target_len[cur->core.tid]) cur->core.flag |= BAM_FUNMAP;
        }
        if (has_prev) { // do we have a pair of reads to examine?
            if (strcmp(bam_get_qname(cur), bam_get_qname(pre)) == 0) { // identical pair name
                pre->core.flag |= BAM_FPAIRED;
                cur->core.flag |= BAM_FPAIRED;
                if (sync_mate(pre, cur)) goto fail;

                if (pre->core.tid == cur->core.tid && !(cur->core.flag&(BAM_FUNMAP|BAM_FMUNMAP))
                    && !(pre->core.flag&(BAM_FUNMAP|BAM_FMUNMAP))) // if safe set TLEN/ISIZE
                {
                    uint32_t cur5, pre5;
                    cur5 = (cur->core.flag&BAM_FREVERSE)? cur_end : cur->core.pos;
                    pre5 = (pre->core.flag&BAM_FREVERSE)? pre_end : pre->core.pos;
                    cur->core.isize = pre5 - cur5; pre->core.isize = cur5 - pre5;
                } else cur->core.isize = pre->core.isize = 0;
                if (add_ct) bam_template_cigar(pre, cur, &str);
                // TODO: Add code to properly check if read is in a proper pair based on ISIZE distribution
                if (proper_pair_check && !plausibly_properly_paired(pre,cur)) {
                    pre->core.flag &= ~BAM_FPROPER_PAIR;
                    cur->core.flag &= ~BAM_FPROPER_PAIR;
                }

                if (do_mate_scoring) {
                    if ((add_mate_score(pre, cur) == -1) || (add_mate_score(cur, pre) == -1)) {
                        fprintf(stderr, "[bam_mating_core] ERROR: unable to add mate score.\n");
                        goto fail;
                    }
                }

                // Write out result
                if ( !remove_reads ) {
                    if (sam_write1(out, header, pre) < 0) goto write_fail;
                    if (sam_write1(out, header, cur) < 0) goto write_fail;
                } else {
                    // If we have to remove reads make sure we do it in a way that doesn't create orphans with bad flags
                    if(pre->core.flag&BAM_FUNMAP) cur->core.flag &= ~(BAM_FPAIRED|BAM_FMREVERSE|BAM_FPROPER_PAIR);
                    if(cur->core.flag&BAM_FUNMAP) pre->core.flag &= ~(BAM_FPAIRED|BAM_FMREVERSE|BAM_FPROPER_PAIR);
                    if(!(pre->core.flag&BAM_FUNMAP)) {
                        if (sam_write1(out, header, pre) < 0) goto write_fail;
                    }
                    if(!(cur->core.flag&BAM_FUNMAP)) {
                        if (sam_write1(out, header, cur) < 0) goto write_fail;
                    }
                }
                has_prev = 0;
            } else { // unpaired?  clear bad info and write it out
                if (pre->core.tid < 0 || pre->core.pos < 0 || pre->core.flag&BAM_FUNMAP) { // If unmapped
                    pre->core.flag |= BAM_FUNMAP;
                    pre->core.tid = -1;
                    pre->core.pos = -1;
                }
                pre->core.mtid = -1; pre->core.mpos = -1; pre->core.isize = 0;
                pre->core.flag &= ~(BAM_FPAIRED|BAM_FMREVERSE|BAM_FPROPER_PAIR);
                if ( !remove_reads || !(pre->core.flag&BAM_FUNMAP) ) {
                    if (sam_write1(out, header, pre) < 0) goto write_fail;
                }
            }
        } else has_prev = 1;
        curr = 1 - curr;
        pre_end = cur_end;
    }
    if (has_prev && !remove_reads) { // If we still have a BAM in the buffer it must be unpaired
        bam1_t *pre = b[1-curr];
        if (pre->core.tid < 0 || pre->core.pos < 0 || pre->core.flag&BAM_FUNMAP) { // If unmapped
            pre->core.flag |= BAM_FUNMAP;
            pre->core.tid = -1;
            pre->core.pos = -1;
        }
        pre->core.mtid = -1; pre->core.mpos = -1; pre->core.isize = 0;
        pre->core.flag &= ~(BAM_FPAIRED|BAM_FMREVERSE|BAM_FPROPER_PAIR);

        if (sam_write1(out, header, pre) < 0) goto write_fail;
    }
    bam_hdr_destroy(header);
    bam_destroy1(b[0]);
    bam_destroy1(b[1]);
    free(str.s);
    return 0;

 write_fail:
    print_error_errno("fixmate", "Couldn't write to output file");
 fail:
    bam_hdr_destroy(header);
    bam_destroy1(b[0]);
    bam_destroy1(b[1]);
    return 1;
}

void usage(FILE* where)
{
    fprintf(where,
"Usage: samtools fixmate <in.nameSrt.bam> <out.nameSrt.bam>\n"
"Options:\n"
"  -r           Remove unmapped reads and secondary alignments\n"
"  -p           Disable FR proper pair check\n"
"  -c           Add template cigar ct tag\n"
"  -m           Add mate score tag\n");

    sam_global_opt_help(where, "-.O..@");

    fprintf(where,
"\n"
"As elsewhere in samtools, use '-' as the filename for stdin/stdout. The input\n"
"file must be grouped by read name (e.g. sorted by name). Coordinated sorted\n"
"input is not accepted.\n");
}

int bam_mating(int argc, char *argv[])
{
    htsThreadPool p = {NULL, 0};
    samFile *in = NULL, *out = NULL;
    int c, remove_reads = 0, proper_pair_check = 1, add_ct = 0, res = 1, mate_score = 0;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    char wmode[3] = {'w', 'b', 0};
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 0, '@'),
        { NULL, 0, NULL, 0 }
    };

    // parse args
    if (argc == 1) { usage(stdout); return 0; }
    while ((c = getopt_long(argc, argv, "rpcmO:@:", lopts, NULL)) >= 0) {
        switch (c) {
            case 'r': remove_reads = 1; break;
            case 'p': proper_pair_check = 0; break;
            case 'c': add_ct = 1; break;
            case 'm': mate_score = 1; break;
            default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
                      /* else fall-through */
            case '?': usage(stderr); goto fail;
        }
    }
    if (optind+1 >= argc) { usage(stderr); goto fail; }

    // init
    if ((in = sam_open_format(argv[optind], "rb", &ga.in)) == NULL) {
        print_error_errno("fixmate", "cannot open input file");
        goto fail;
    }
    sam_open_mode(wmode+1, argv[optind+1], NULL);
    if ((out = sam_open_format(argv[optind+1], wmode, &ga.out)) == NULL) {
        print_error_errno("fixmate", "cannot open output file");
        goto fail;
    }

    if (ga.nthreads > 0) {
        if (!(p.pool = hts_tpool_init(ga.nthreads))) {
            fprintf(stderr, "Error creating thread pool\n");
            goto fail;
        }
        hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
        hts_set_opt(out, HTS_OPT_THREAD_POOL, &p);
    }

    // run
    res = bam_mating_core(in, out, remove_reads, proper_pair_check, add_ct, mate_score);

    // cleanup
    sam_close(in);
    if (sam_close(out) < 0) {
        fprintf(stderr, "[bam_mating] error while closing output file\n");
        res = 1;
    }

    if (p.pool) hts_tpool_destroy(p.pool);
    sam_global_args_free(&ga);
    return res;

 fail:
    if (in) sam_close(in);
    if (out) sam_close(out);
    if (p.pool) hts_tpool_destroy(p.pool);
    sam_global_args_free(&ga);
    return 1;
}


