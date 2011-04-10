/* This program demonstrates how to generate pileup from multiple BAMs
 * simutaneously, to achieve random access and to use the BED interface.
 * To compile this program separately, you may:
 *
 *   gcc -g -O2 -Wall -o bam2depth -D_MAIN_BAM2DEPTH bam2depth.c -L. -lbam -lz
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "bam.h"

typedef struct {     // auxiliary data structure
	bamFile fp;      // the file handler
	bam_iter_t iter; // NULL if a region not specified
	int min_mapQ;    // mapQ filter
} aux_t;

void *bed_read(const char *fn); // read a BED or position list file
void bed_destroy(void *_h);     // destroy the BED data structure
int bed_overlap(const void *_h, const char *chr, int beg, int end); // test if chr:beg-end overlaps

// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret = aux->iter? bam_iter_read(aux->fp, aux->iter, b) : bam_read1(aux->fp, b);
	if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
	return ret;
}

#ifdef _MAIN_BAM2DEPTH
int main(int argc, char *argv[])
#else
int main_depth(int argc, char *argv[])
#endif
{
	int i, n, tid, beg, end, pos, *n_plp, baseQ = 0, mapQ = 0;
	const bam_pileup1_t **plp;
	char *reg = 0; // specified region
	void *bed = 0; // BED data structure
	bam_header_t *h = 0; // BAM header of the 1st input
	aux_t **data;
	bam_mplp_t mplp;

	// parse the command line
	while ((n = getopt(argc, argv, "r:b:q:Q:")) >= 0) {
		switch (n) {
			case 'r': reg = strdup(optarg); break;   // parsing a region requires a BAM header
			case 'b': bed = bed_read(optarg); break; // BED or position list file can be parsed now
			case 'q': baseQ = atoi(optarg); break;   // base quality threshold
			case 'Q': mapQ = atoi(optarg); break;    // mapping quality threshold
		}
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: bam2depth [-r reg] [-q baseQthres] [-Q mapQthres] [-b in.bed] <in1.bam> [...]\n");
		return 1;
	}

	// initialize the auxiliary data structures
	n = argc - optind; // the number of BAMs on the command line
	data = calloc(n, sizeof(void*)); // data[i] for the i-th input
	beg = 0; end = 1<<30; tid = -1;  // set the default region
	for (i = 0; i < n; ++i) {
		bam_header_t *htmp;
		data[i] = calloc(1, sizeof(aux_t));
		data[i]->fp = bam_open(argv[optind+i], "r"); // open BAM
		data[i]->min_mapQ = mapQ;                    // set the mapQ filter
		htmp = bam_header_read(data[i]->fp);         // read the BAM header
		if (i == 0) {
			h = htmp; // keep the header of the 1st BAM
			if (reg) bam_parse_region(h, reg, &tid, &beg, &end); // also parse the region
		} else bam_header_destroy(htmp); // if not the 1st BAM, trash the header
		if (tid >= 0) { // if a region is specified and parsed successfully
			bam_index_t *idx = bam_index_load(argv[optind+i]);  // load the index
			data[i]->iter = bam_iter_query(idx, tid, beg, end); // set the iterator
			bam_index_destroy(idx); // the index is not needed any more; phase out of the memory
		}
	}

	// the core multi-pileup loop
	mplp = bam_mplp_init(n, read_bam, (void**)data); // initialization
	n_plp = calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
	plp = calloc(n, sizeof(void*)); // plp[i] points to the array of covering reads (internal in mplp)
	while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) { // come to the next covered position
		if (pos < beg || pos >= end) continue; // out of range; skip
		if (bed && bed_overlap(bed, h->target_name[tid], pos, pos + 1) == 0) continue; // not in BED; skip
		fputs(h->target_name[tid], stdout); printf("\t%d", pos+1); // a customized printf() would be faster
		for (i = 0; i < n; ++i) { // base level filters have to go here
			int j, m = 0;
			for (j = 0; j < n_plp[i]; ++j) {
				const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
				if (p->is_del || p->is_refskip) ++m; // having dels or refskips at tid:pos
				else if (bam1_qual(p->b)[p->qpos] < baseQ) ++m; // low base quality
			}
			printf("\t%d", n_plp[i] - m); // this the depth to output
		}
		putchar('\n');
	}
	free(n_plp); free(plp);
	bam_mplp_destroy(mplp);

	bam_header_destroy(h);
	for (i = 0; i < n; ++i) {
		bam_close(data[i]->fp);
		if (data[i]->iter) bam_iter_destroy(data[i]->iter);
		free(data[i]);
	}
	free(data); free(reg);
	if (bed) bed_destroy(bed);
	return 0;
}
