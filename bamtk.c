#include <stdio.h>
#include <unistd.h>
#include "bam.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1.1-18"
#endif

int bam_taf2baf(int argc, char *argv[]);
int bam_pileup(int argc, char *argv[]);
int bam_merge(int argc, char *argv[]);
int bam_index(int argc, char *argv[]);
int bam_sort(int argc, char *argv[]);
int bam_tview_main(int argc, char *argv[]);
int bam_mating(int argc, char *argv[]);
int bam_rmdup(int argc, char *argv[]);

int faidx_main(int argc, char *argv[]);
int glf_view_main(int argc, char *argv[]);

static int view_aux(const bam1_t *b, void *data)
{
	bam_view1((bam_header_t*)data, b);
	return 0;
}
static int view_auxb(const bam1_t *b, void *data)
{
	bam_write1((bamFile)data, b);
	return 0;
}

int bam_view(int argc, char *argv[])
{
	bamFile fp, fpout = 0;
	bam_header_t *header;
	bam1_t *b;
	int ret, c, is_bam = 0;
	while ((c = getopt(argc, argv, "b")) >= 0) {
		switch (c) {
		case 'b': is_bam = 1; break;
		default: fprintf(stderr, "Unrecognized option: -%c\n", c); return 1;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: samtools view [-b] <in.bam> [<region> [...]]\n");
		return 1;
	}
	fp = strcmp(argv[optind], "-")? bam_open(argv[optind], "r") : bam_dopen(fileno(stdin), "r");
	assert(fp);
	header = bam_header_read(fp);
	if (is_bam) {
		assert(fpout = bam_dopen(fileno(stdout), "w"));
		bam_header_write(fpout, header);
	}
	if (optind + 1 == argc) {
		b = (bam1_t*)calloc(1, sizeof(bam1_t));
		while ((ret = bam_read1(fp, b)) >= 0) bam_view1(header, b);
		if (ret < -1) fprintf(stderr, "[bam_view] truncated file? Continue anyway. (%d)\n", ret);
		free(b->data); free(b);
	} else {
		int i;
		bam_index_t *idx;
		idx = bam_index_load(argv[optind]);
		for (i = optind + 1; i < argc; ++i) {
			int tid, beg, end;
			bam_parse_region(header, argv[i], &tid, &beg, &end);
			if (is_bam) bam_fetch(fp, idx, tid, beg, end, fpout, view_auxb);
			else bam_fetch(fp, idx, tid, beg, end, header, view_aux);
		}
		bam_index_destroy(idx);
	}
	bam_header_destroy(header);
	bam_close(fp);
	if (is_bam) bam_close(fpout);
	return 0;
}

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: samtools (Tools for alignments in the SAM format)\n");
	fprintf(stderr, "Version: %s\n\n", PACKAGE_VERSION);
	fprintf(stderr, "Usage:   samtools <command> [options]\n\n");
	fprintf(stderr, "Command: import      import from the text format\n");
	fprintf(stderr, "         view        export to the text format\n");
	fprintf(stderr, "         sort        sort alignment file\n");
	fprintf(stderr, "         merge       merge multiple sorted alignment files\n");
	fprintf(stderr, "         pileup      generate pileup output\n");
	fprintf(stderr, "         faidx       index/extract FASTA\n");
#ifndef _NO_CURSES
	fprintf(stderr, "         tview       text alignment viewer\n");
#endif
	fprintf(stderr, "         index       index alignment\n");
	fprintf(stderr, "         fixmate     fix mate information\n");
	fprintf(stderr, "         rmdup       remove PCR duplicates\n");
	fprintf(stderr, "         glfview     print GLFv2 file\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();
	if (strcmp(argv[1], "view") == 0) return bam_view(argc-1, argv+1);
	else if (strcmp(argv[1], "import") == 0) return bam_taf2baf(argc-1, argv+1);
	else if (strcmp(argv[1], "pileup") == 0) return bam_pileup(argc-1, argv+1);
	else if (strcmp(argv[1], "merge") == 0) return bam_merge(argc-1, argv+1);
	else if (strcmp(argv[1], "sort") == 0) return bam_sort(argc-1, argv+1);
	else if (strcmp(argv[1], "index") == 0) return bam_index(argc-1, argv+1);
	else if (strcmp(argv[1], "faidx") == 0) return faidx_main(argc-1, argv+1);
	else if (strcmp(argv[1], "fixmate") == 0) return bam_mating(argc-1, argv+1);
	else if (strcmp(argv[1], "rmdup") == 0) return bam_rmdup(argc-1, argv+1);
	else if (strcmp(argv[1], "glfview") == 0) return glf_view_main(argc-1, argv+1);
#ifndef _NO_CURSES
	else if (strcmp(argv[1], "tview") == 0) return bam_tview_main(argc-1, argv+1);
#endif
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;	
}
