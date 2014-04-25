#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>

#include "test.h"

void xfreopen(const char *path, const char *mode, FILE *stream)
{
	if (freopen(path, mode, stream) == NULL) {
		fprintf(stderr, __FILE__": error reopening %s: %s\n",
				path, strerror(errno));
		exit(2);
	}
}

void dump_hdr(const bam_hdr_t* hdr)
{
	printf("n_targets: %d\n", hdr->n_targets);
	printf("ignore_sam_err: %d\n", hdr->ignore_sam_err);
	printf("l_text: %u\n", hdr->l_text);
	printf("idx\ttarget_len\ttarget_name:\n");
	int32_t target;
	for (target = 0; target < hdr->n_targets; ++target) {
		printf("%d\t%u\t\"%s\"\n", target, hdr->target_len[target], hdr->target_name[target]);
	}
	printf("text: \"%s\"\n", hdr->text);
}

void dump_read(const bam1_t* read)
{
	printf("->core.tid: %d\n"
		   "->core.pos: %d\n"
		   "->core.bin: %u\n"
		   "->core.qual: %u\n"
		   "->core.l_qname: %u\n"
		   "->core.flag: %u\n"
		   "->core.n_cigar: %u\n"
		   "->core.l_qseq: %d\n"
		   "->core.mtid: %d\n"
		   "->core.mpos: %d\n"
		   "->core.isize: %d\n",
		   read->core.tid,
		   read->core.pos,
		   read->core.bin,
		   read->core.qual,
		   read->core.l_qname,
		   read->core.flag,
		   read->core.n_cigar,
		   read->core.l_qseq,
		   read->core.mtid,
		   read->core.mpos,
		   read->core.isize
	);
	if (read->data) {
		printf("->data:");
		int i;
		for (i = 0; i < read->l_data; ++i) {
			printf("%x ", read->data[i]);
		}
		printf("\n");
	}
	if (read->core.l_qname) {
		printf("qname: %s\n",bam_get_qname(read));
	}
	if (read->core.l_qseq) {
		printf("qseq:");
		int i;
		for (i = 0; i < read->core.l_qseq; ++i) {
			printf("%c",seq_nt16_str[seq_nt16_table[bam_seqi(bam_get_seq(read),i)]]);
		}
		printf("\n");
		printf("qual:");
		for (i = 0; i < read->core.l_qseq; ++i) {
			printf("%c",bam_get_qual(read)[i]);
		}
		printf("\n");
	}

	if (bam_get_l_aux(read)) {
		int i = 0;
		uint8_t* aux = bam_get_aux(read);
		
		while (i < bam_get_l_aux(read)) {
			printf("%.2s:%c:",aux+i,*(aux+i+2));
			i += 2;
			switch (*(aux+i)) {
				case 'Z':
					while (*(aux+1+i) != '\0') { putc(*(aux+1+i), stdout); ++i; }
					break;
			}
			putc('\n',stdout);
			++i;++i;
		}
	}
	printf("\n");
}
