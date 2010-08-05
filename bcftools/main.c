#include <string.h>
#include <stdlib.h>
#include "bcf.h"

int bcfview(int argc, char *argv[]);
int bcf_main_index(int argc, char *argv[]);

int main(int argc, char *argv[])
{
	if (argc == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bcftools <command> <arguments>\n\n");
		fprintf(stderr, "Command: view      print, extract, convert and call SNPs from BCF\n");
		fprintf(stderr, "         index     index BCF\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if (strcmp(argv[1], "view") == 0) return bcfview(argc-1, argv+1);
	if (strcmp(argv[1], "index") == 0) return bcf_main_index(argc-1, argv+1);
	if (strcmp(argv[1], "test") == 0) return vcf_test(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] Unrecognized command.\n");
		return 1;
	}
	return 0;
}
