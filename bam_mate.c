#include <stdlib.h>
#include <string.h>
#include "bam.h"

// currently, this function ONLY works if each read has one hit
void bam_mating_core(bamFile in, bamFile out)
{
	bam_header_t *header;
	bam1_t *b[2];
	int curr, has_prev;

	header = bam_header_read(in);
	bam_header_write(out, header);

	b[0] = bam_init1();
	b[1] = bam_init1();
	curr = 0; has_prev = 0;
	while (bam_read1(in, b[curr]) >= 0) {
		bam1_t *cur = b[curr], *pre = b[1-curr];
		if (has_prev) {
			if (strcmp(bam1_qname(cur), bam1_qname(pre)) == 0) { // identical pair name
				cur->core.mtid = pre->core.tid; cur->core.mpos = pre->core.pos;
				pre->core.mtid = cur->core.tid; pre->core.mpos = cur->core.pos;
				if (pre->core.tid == cur->core.tid && !(cur->core.flag&(BAM_FUNMAP|BAM_FMUNMAP))
					&& !(pre->core.flag&(BAM_FUNMAP|BAM_FMUNMAP)))
				{
					uint32_t cur5, pre5;
					cur5 = (cur->core.flag&BAM_FREVERSE)? bam_calend(&cur->core, bam1_cigar(cur)) : cur->core.pos;
					pre5 = (pre->core.flag&BAM_FREVERSE)? bam_calend(&pre->core, bam1_cigar(pre)) : pre->core.pos;
					cur->core.isize = pre5 - cur5; pre->core.isize = cur5 - pre5;
				} else cur->core.isize = pre->core.isize = 0;
				if (pre->core.flag&BAM_FREVERSE) cur->core.flag |= BAM_FMREVERSE;
				else cur->core.flag &= ~BAM_FMREVERSE;
				if (cur->core.flag&BAM_FREVERSE) pre->core.flag |= BAM_FMREVERSE;
				else pre->core.flag &= ~BAM_FMREVERSE;
				if (cur->core.flag & BAM_FUNMAP) { pre->core.flag |= BAM_FMUNMAP; pre->core.flag &= ~BAM_FPROPER_PAIR; }
				if (pre->core.flag & BAM_FUNMAP) { cur->core.flag |= BAM_FMUNMAP; cur->core.flag &= ~BAM_FPROPER_PAIR; }
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
	}
	if (has_prev) bam_write1(out, b[1-curr]);
	bam_header_destroy(header);
	bam_destroy1(b[0]);
	bam_destroy1(b[1]);
}

int bam_mating(int argc, char *argv[])
{
	bamFile in, out;
	if (argc < 3) {
		fprintf(stderr, "samtools fixmate <in.nameSrt.bam> <out.nameSrt.bam>\n");
		return 1;
	}
	in = (strcmp(argv[1], "-") == 0)? bam_dopen(fileno(stdin), "r") : bam_open(argv[1], "r");
    out = (strcmp(argv[2], "-") == 0)? bam_dopen(fileno(stdout), "w") : bam_open(argv[2], "w");
	bam_mating_core(in, out);
	bam_close(in); bam_close(out);
	return 0;
}
