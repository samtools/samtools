#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "kstring.h"
#include "bam.h"

void bam_template_cigar(bam1_t *b1, bam1_t *b2, kstring_t *str)
{
	bam1_t *swap;
	int i, end;
	uint32_t *cigar;
	str->l = 0;
	if (b1->core.tid != b2->core.tid || b1->core.tid < 0) return; // coordinateless or not on the same chr; skip
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

// currently, this function ONLY works if each read has one hit
void bam_mating_core(bamFile in, bamFile out, int remove_reads)
{
	bam_header_t *header;
	bam1_t *b[2];
	int curr, has_prev, pre_end = 0, cur_end;
	kstring_t str;

	str.l = str.m = 0; str.s = 0;
	header = bam_header_read(in);
	bam_header_write(out, header);

	b[0] = bam_init1();
	b[1] = bam_init1();
	curr = 0; has_prev = 0;
	while (bam_read1(in, b[curr]) >= 0) {
		bam1_t *cur = b[curr], *pre = b[1-curr];
		if (cur->core.tid < 0) 
        {
            if ( !remove_reads ) bam_write1(out, cur);
            continue;
        }
		cur_end = bam_calend(&cur->core, bam1_cigar(cur));
		if (cur_end > (int)header->target_len[cur->core.tid]) cur->core.flag |= BAM_FUNMAP;
		if (cur->core.flag & BAM_FSECONDARY) 
        {
            if ( !remove_reads ) bam_write1(out, cur);
            continue; // skip secondary alignments
        }
		if (has_prev) {
			if (strcmp(bam1_qname(cur), bam1_qname(pre)) == 0) { // identical pair name
				cur->core.mtid = pre->core.tid; cur->core.mpos = pre->core.pos;
				pre->core.mtid = cur->core.tid; pre->core.mpos = cur->core.pos;
				if (pre->core.tid == cur->core.tid && !(cur->core.flag&(BAM_FUNMAP|BAM_FMUNMAP))
					&& !(pre->core.flag&(BAM_FUNMAP|BAM_FMUNMAP))) // set TLEN/ISIZE
				{
					uint32_t cur5, pre5;
					cur5 = (cur->core.flag&BAM_FREVERSE)? cur_end : cur->core.pos;
					pre5 = (pre->core.flag&BAM_FREVERSE)? pre_end : pre->core.pos;
					cur->core.isize = pre5 - cur5; pre->core.isize = cur5 - pre5;
				} else cur->core.isize = pre->core.isize = 0;
				if (pre->core.flag&BAM_FREVERSE) cur->core.flag |= BAM_FMREVERSE;
				else cur->core.flag &= ~BAM_FMREVERSE;
				if (cur->core.flag&BAM_FREVERSE) pre->core.flag |= BAM_FMREVERSE;
				else pre->core.flag &= ~BAM_FMREVERSE;
				if (cur->core.flag & BAM_FUNMAP) { pre->core.flag |= BAM_FMUNMAP; pre->core.flag &= ~BAM_FPROPER_PAIR; }
				if (pre->core.flag & BAM_FUNMAP) { cur->core.flag |= BAM_FMUNMAP; cur->core.flag &= ~BAM_FPROPER_PAIR; }
				bam_template_cigar(pre, cur, &str);
				bam_write1(out, pre);
				bam_write1(out, cur);
				has_prev = 0;
			} else { // unpaired or singleton
				pre->core.mtid = -1; pre->core.mpos = -1; pre->core.isize = 0;
				if (pre->core.flag & BAM_FPAIRED) {
					pre->core.flag |= BAM_FMUNMAP;
					pre->core.flag &= ~BAM_FMREVERSE & ~BAM_FPROPER_PAIR;
				}
				bam_write1(out, pre);
			}
		} else has_prev = 1;
		curr = 1 - curr;
		pre_end = cur_end;
	}
	if (has_prev) bam_write1(out, b[1-curr]);
	bam_header_destroy(header);
	bam_destroy1(b[0]);
	bam_destroy1(b[1]);
	free(str.s);
}

void usage()
{
    fprintf(stderr,"Usage: samtools fixmate <in.nameSrt.bam> <out.nameSrt.bam>\n");
    fprintf(stderr,"Options:\n");
    fprintf(stderr,"       -r    remove unmapped reads and secondary alignments\n");
    exit(1);
}

int bam_mating(int argc, char *argv[])
{
	bamFile in, out;
    int c, remove_reads=0;
    while ((c = getopt(argc, argv, "r")) >= 0) {
        switch (c) {
            case 'r': remove_reads=1; break;
        }
    }
    if (optind+1 >= argc) usage();
	in = (strcmp(argv[optind], "-") == 0)? bam_dopen(fileno(stdin), "r") : bam_open(argv[optind], "r");
    out = (strcmp(argv[optind+1], "-") == 0)? bam_dopen(fileno(stdout), "w") : bam_open(argv[optind+1], "w");
	bam_mating_core(in, out, remove_reads);
	bam_close(in); bam_close(out);
	return 0;
}


