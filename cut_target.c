#include <unistd.h>
#include <stdlib.h>
#include "bam.h"

int min_baseQ = 20;

int main(int argc, char *argv[])
{
	bamFile fp;
	int c, tid, pos, n, lasttid = -1, lastpos = -1;
	const bam_pileup1_t *p;
	bam_plp_t plp;
	while ((c = getopt(argc, argv, "Q:")) >= 0) {
		switch (c) {
			case 'Q': min_baseQ = atoi(optarg); break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: cut_target <in.bam>\n");
		return 1;
	}
	fp = bam_open(argv[optind], "r");
	plp = bam_plp_init(bam_read1, fp);
	while ((p = bam_plp_next(plp, &tid, &pos, &n)) >= 0) {
		if (tid != lasttid) {
		}
	}
	bam_plp_destroy(plp);
	bam_close(fp);
	return 0;
}
