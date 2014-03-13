#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>

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
