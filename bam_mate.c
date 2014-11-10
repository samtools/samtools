/*  bam_mate.c -- fix mate pairing information and clean up flags.

    Copyright (C) 2009, 2011-2014 Genome Research Ltd.
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

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "htslib/kstring.h"
#include "htslib/sam.h"

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

static void sync_mq(bam1_t* src, bam1_t* dest)
{
    if ( (src->core.flag & BAM_FUNMAP) == 0 ) { // If mapped
        uint32_t mq = src->core.qual;
        uint8_t* data;
        if ((data = bam_aux_get(dest,"MQ")) != NULL) {
            bam_aux_del(dest, data);
        }

        bam_aux_append(dest, "MQ", 'i', sizeof(uint32_t), (uint8_t*)&mq);
    }
}

// copy flags
static void sync_mate(bam1_t* a, bam1_t* b)
{
    sync_unmapped_pos_inner(a,b);
    sync_unmapped_pos_inner(b,a);
    sync_mate_inner(a,b);
    sync_mate_inner(b,a);
    sync_mq(a,b);
    sync_mq(b,a);
}

// currently, this function ONLY works if each read has one hit
static void bam_mating_core(samFile* in, samFile* out, int remove_reads, int proper_pair_check, int add_ct)
{
    bam_hdr_t *header;
    bam1_t *b[2];
    int curr, has_prev, pre_end = 0, cur_end = 0;
    kstring_t str;

    str.l = str.m = 0; str.s = 0;
    header = sam_hdr_read(in);
    // Accept unknown, unsorted, or queryname sort order, but error on coordinate sorted.
    if ((header->l_text > 3) && (strncmp(header->text, "@HD", 3) == 0)) {
        char *p, *q;
        p = strstr(header->text, "\tSO:coordinate");
        q = strchr(header->text, '\n');
        // Looking for SO:coordinate within the @HD line only
        // (e.g. must ignore in a @CO comment line later in header)
        if ((p != 0) && (p < q)) {
            fprintf(stderr, "[bam_mating_core] ERROR: Coordinate sorted, require grouped/sorted by queryname.\n");
            exit(1);
        }
    }
    sam_hdr_write(out, header);

    b[0] = bam_init1();
    b[1] = bam_init1();
    curr = 0; has_prev = 0;
    while (sam_read1(in, header, b[curr]) >= 0) {
        bam1_t *cur = b[curr], *pre = b[1-curr];
        if (cur->core.flag & BAM_FSECONDARY)
        {
            if ( !remove_reads ) sam_write1(out, header, cur);
            continue; // skip secondary alignments
        }
        if (cur->core.flag & BAM_FSUPPLEMENTARY)
        {
            sam_write1(out, header, cur);
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
                sync_mate(pre, cur);

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

                // Write out result
                if ( !remove_reads ) {
                    sam_write1(out, header, pre);
                    sam_write1(out, header, cur);
                } else {
                    // If we have to remove reads make sure we do it in a way that doesn't create orphans with bad flags
                    if(pre->core.flag&BAM_FUNMAP) cur->core.flag &= ~(BAM_FPAIRED|BAM_FMREVERSE|BAM_FPROPER_PAIR);
                    if(cur->core.flag&BAM_FUNMAP) pre->core.flag &= ~(BAM_FPAIRED|BAM_FMREVERSE|BAM_FPROPER_PAIR);
                    if(!(pre->core.flag&BAM_FUNMAP)) sam_write1(out, header, pre);
                    if(!(cur->core.flag&BAM_FUNMAP)) sam_write1(out, header, cur);
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
                if ( !remove_reads || !(pre->core.flag&BAM_FUNMAP) ) sam_write1(out, header, pre);
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

        sam_write1(out, header, pre);
    }
    bam_hdr_destroy(header);
    bam_destroy1(b[0]);
    bam_destroy1(b[1]);
    free(str.s);
}

void usage(FILE* where)
{
    fprintf(where,"Usage: samtools fixmate <in.nameSrt.bam> <out.nameSrt.bam>\n\n");
    fprintf(where,"Options:\n");
    fprintf(stderr,"  -r         Remove unmapped reads and secondary alignments\n");
    fprintf(stderr,"  -p         Disable FR proper pair check\n");
    fprintf(stderr,"  -c         Add template cigar ct tag\n");
    fprintf(stderr,"  -O FORMAT  Write output as FORMAT ('sam'/'bam'/'cram')\n");
    fprintf(stderr,"As elsewhere in samtools, use '-' as the filename for stdin/stdout. The input\n");
    fprintf(stderr,"file must be grouped by read name (e.g. sorted by name). Coordinated sorted\n");
    fprintf(stderr,"input is not accepted.\n");
}

int bam_mating(int argc, char *argv[])
{
    samFile *in, *out;
    int c, remove_reads = 0, proper_pair_check = 1, add_ct = 0;
    char* fmtout = NULL;
    char modeout[12];
    // parse args
    if (argc == 1) { usage(stdout); return 0; }
    while ((c = getopt(argc, argv, "rpcO:")) >= 0) {
        switch (c) {
            case 'r': remove_reads = 1; break;
            case 'p': proper_pair_check = 0; break;
            case 'c': add_ct = 1; break;
            case 'O': fmtout = optarg; break;
            default: usage(stderr); return 1;
        }
    }
    if (optind+1 >= argc) { usage(stderr); return 1; }
    strcpy(modeout, "w");
    if (sam_open_mode(&modeout[1], argv[optind+1], fmtout) < 0) {
        if (fmtout) fprintf(stderr, "[bam_mating] cannot parse output format \"%s\"\n", fmtout);
        else fprintf(stderr, "[bam_mating] cannot determine output format\n");
        return 1;
    }

    // init
    if ((in = sam_open(argv[optind], "r")) == NULL) {
        fprintf(stderr, "[bam_mating] cannot open input file\n");
        return 1;
    }
    if ((out = sam_open(argv[optind+1], modeout)) == NULL) {
        fprintf(stderr, "[bam_mating] cannot open output file\n");
        return 1;
    }

    // run
    bam_mating_core(in, out, remove_reads, proper_pair_check, add_ct);
    // cleanup
    sam_close(in); sam_close(out);
    return 0;
}


