// Copyright (c) 2009-2013 Genome Research Limited.
//
// This file is part of samtools.
//
// samtools is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// this program. If not, see L<http://www.gnu.org/licenses/>.

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "htslib/kstring.h"
#include "bam.h"

/*
 * This function calculates CT tag for two bams, it assumes they are from the same template and
 * writes the tag to the first read in position terms.
 */
static void bam_template_cigar(bam1_t *b1, bam1_t *b2, kstring_t *str)
{
	bam1_t *swap;
	int i, end;
	uint32_t *cigar;
	str->l = 0;
	if (b1->core.tid != b2->core.tid || b1->core.tid < 0 || b1->core.pos == 0 || b2->core.pos == 0 || b1->core.flag&BAM_FUNMAP || b2->core.flag&BAM_FUNMAP) return; // coordinateless or not on the same chr; skip
	if (b1->core.pos > b2->core.pos) swap = b1, b1 = b2, b2 = swap; // make sure b1 has a smaller coordinate
	kputc((b1->core.flag & BAM_FREAD1)? '1' : '2', str); // segment index
	kputc((b1->core.flag & BAM_FREVERSE)? 'R' : 'F', str); // strand
	for (i = 0, cigar = bam1_cigar(b1); i < b1->core.n_cigar; ++i) {
		kputw(bam_cigar_oplen(cigar[i]), str);
		kputc(bam_cigar_opchr(cigar[i]), str);
	}
	end = bam_calend(&b1->core, cigar);
	kputw(b2->core.pos - end, str);
	kputc('T', str);
	kputc((b2->core.flag & BAM_FREAD1)? '1' : '2', str); // segment index
	kputc((b2->core.flag & BAM_FREVERSE)? 'R' : 'F', str); // strand
	for (i = 0, cigar = bam1_cigar(b2); i < b2->core.n_cigar; ++i) {
		kputw(bam_cigar_oplen(cigar[i]), str);
		kputc(bam_cigar_opchr(cigar[i]), str);
	}
	bam_aux_append(b1, "CT", 'Z', str->l+1, (uint8_t*)str->s); 
}

/*
 * What This Program is Supposed To Do:
 * Fill in mate coordinates, ISIZE and mate related flags from a name-sorted alignment.
 *
 * How We Handle Input
 *
 * Secondary Reads:
 * -write to output unchanged
 * All Reads:
 * -if pos == 0, tid == -1 set UNMAPPED flag
 * single Reads:
 * -if pos == 0, tid == -1, or UNMAPPED then set UNMAPPED, pos = 0, tid = -1
 * -clear flags (PAIRED, MREVERSE, PROPER_PAIR)
 * -set mpos = 0, mtid = -1 and isize = 0
 * -write to output
 * Paired Reads:
 * -if read is unmapped and mate is not, set pos and tid to equal that of mate
 * -sync mate flags (MREVERSE, MUNMAPPED), mpos, mtid
 * -recalculate ISIZE if possible, otherwise set it to 0
 * -optionally clear PROPER_PAIR flag from reads where mapping or orientation indicate this is not possible
 * -calculate CT and apply to lowest positioned read
 * -write to output
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
	if (a->core.pos > b->core.pos) {
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

// copy flags
static void sync_mate(bam1_t* a, bam1_t* b)
{
	sync_unmapped_pos_inner(a,b);
	sync_unmapped_pos_inner(b,a);
	sync_mate_inner(a,b);
	sync_mate_inner(b,a);
}

// currently, this function ONLY works if each read has one hit
static void bam_mating_core(bamFile in, bamFile out, int remove_reads, int proper_pair_check)
{
	bam_header_t *header;
	bam1_t *b[2];
	int curr, has_prev, pre_end = 0, cur_end;
	kstring_t str;

	str.l = str.m = 0; str.s = 0;
	header = bam_header_read(in);
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
	bam_header_write(out, header);

	b[0] = bam_init1();
	b[1] = bam_init1();
	curr = 0; has_prev = 0;
	while (bam_read1(in, b[curr]) >= 0) {
		bam1_t *cur = b[curr], *pre = b[1-curr];
		if (cur->core.flag & BAM_FSECONDARY)
		{
			if ( !remove_reads ) bam_write1(out, cur);
			continue; // skip secondary alignments
		}
		if (cur->core.tid < 0 || cur->core.pos == 0) // If unmapped set the flag
		{
			cur->core.flag |= BAM_FUNMAP;
		}
		if ((cur->core.flag&BAM_FUNMAP) != 0) // If mapped calculate end
		{
			cur_end = bam_calend(&cur->core, bam1_cigar(cur));

			// Check cur_end isn't past the end of the contig we're on, if it is set the UNMAP'd flag
			if (cur_end > (int)header->target_len[cur->core.tid]) cur->core.flag |= BAM_FUNMAP;
		}
		if (has_prev) { // do we have a pair of reads to examine?
			if (strcmp(bam1_qname(cur), bam1_qname(pre)) == 0) { // identical pair name
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
				bam_template_cigar(pre, cur, &str);
				// TODO: Add code to properly check if read is in a proper pair based on ISIZE distribution
				if (proper_pair_check && !plausibly_properly_paired(pre,cur)) {
					pre->core.flag &= ~BAM_FPROPER_PAIR;
					cur->core.flag &= ~BAM_FPROPER_PAIR;
				}
				
				// Write out result
				if ( !remove_reads ) {
					bam_write1(out, pre);
					bam_write1(out, cur);
				} else {
					// If we have to remove reads make sure we do it in a way that doesn't create orphans with bad flags
					if(pre->core.flag&BAM_FUNMAP) cur->core.flag &= ~(BAM_FPAIRED|BAM_FMREVERSE|BAM_FPROPER_PAIR);
					if(cur->core.flag&BAM_FUNMAP) pre->core.flag &= ~(BAM_FPAIRED|BAM_FMREVERSE|BAM_FPROPER_PAIR);
					if(!(pre->core.flag&BAM_FUNMAP)) bam_write1(out, pre);
					if(!(cur->core.flag&BAM_FUNMAP)) bam_write1(out, cur);
				}
				has_prev = 0;
			} else { // unpaired?  clear bad info and write it out
				if (pre->core.tid < 0 || pre->core.pos == 0 || pre->core.flag&BAM_FUNMAP) { // If unmapped
					pre->core.flag |= BAM_FUNMAP;
					pre->core.tid = -1;
					pre->core.pos = 0;
				}
				pre->core.mtid = -1; pre->core.mpos = -1; pre->core.isize = 0;
				pre->core.flag &= ~(BAM_FPAIRED|BAM_FMREVERSE|BAM_FPROPER_PAIR);
				if ( !remove_reads || !(pre->core.flag&BAM_FUNMAP) ) bam_write1(out, pre);
			}
		} else has_prev = 1;
		curr = 1 - curr;
		pre_end = cur_end;
	}
	if (has_prev && !remove_reads) { // If we still have a BAM in the buffer it must be unpaired
		bam1_t *pre = b[1-curr];
		if (pre->core.tid < 0 || pre->core.pos == 0 || pre->core.flag&BAM_FUNMAP) { // If unmapped
			pre->core.flag |= BAM_FUNMAP;
			pre->core.tid = -1;
			pre->core.pos = 0;
		}
		pre->core.mtid = -1; pre->core.mpos = -1; pre->core.isize = 0;
		pre->core.flag &= ~(BAM_FPAIRED|BAM_FMREVERSE|BAM_FPROPER_PAIR);

		bam_write1(out, pre);
	}
	bam_header_destroy(header);
	bam_destroy1(b[0]);
	bam_destroy1(b[1]);
	free(str.s);
}

void usage()
{
	fprintf(stderr,"Usage: samtools fixmate <in.nameSrt.bam> <out.nameSrt.bam>\n\n");
	fprintf(stderr,"Options:\n");
	fprintf(stderr,"       -r    remove unmapped reads and secondary alignments\n");
	fprintf(stderr,"       -p    disable FR proper pair check\n\n");
	fprintf(stderr,"As elsewhere in samtools, use '-' as the filename for stdin/stdout. The input\n");
	fprintf(stderr,"file must be grouped by read name (e.g. sorted by name). Coordinated sorted\n");
	fprintf(stderr,"input is not accepted.\n");
	exit(1);
}

int bam_mating(int argc, char *argv[])
{
	bamFile in, out;
	int c, remove_reads = 0, proper_pair_check = 1;
	while ((c = getopt(argc, argv, "rp")) >= 0) {
		switch (c) {
			case 'r': remove_reads = 1; break;
			case 'p': proper_pair_check = 0; break;
		}
	}
	if (optind+1 >= argc) usage();
	in = (strcmp(argv[optind], "-") == 0)? bam_dopen(fileno(stdin), "r") : bam_open(argv[optind], "r");
	out = (strcmp(argv[optind+1], "-") == 0)? bam_dopen(fileno(stdout), "w") : bam_open(argv[optind+1], "w");
	bam_mating_core(in, out, remove_reads, proper_pair_check);
	bam_close(in); bam_close(out);
	return 0;
}


