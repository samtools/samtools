#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include "sam_header.h"
#include "sam.h"
#include "faidx.h"
#include "kstring.h"
#include "khash.h"
KHASH_SET_INIT_STR(rg)

// When counting records instead of printing them,
// data passed to the bam_fetch callback is encapsulated in this struct.
typedef struct {
	bam_header_t *header;
	int64_t *count;  // int does overflow for very big BAMs
} count_func_data_t;

typedef khash_t(rg) *rghash_t;

// FIXME: we'd better use no global variables...
static rghash_t g_rghash = 0;
static int g_min_mapQ = 0, g_flag_on = 0, g_flag_off = 0, g_qual_scale = 0, g_min_qlen = 0;
static uint32_t g_subsam_seed = 0;
static double g_subsam_frac = -1.;
static char *g_library, *g_rg;
static void *g_bed;

void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);

static int process_aln(const bam_header_t *h, bam1_t *b)
{
	if (g_qual_scale > 1) {
		int i;
		uint8_t *qual = bam1_qual(b);
		for (i = 0; i < b->core.l_qseq; ++i) {
			int c = qual[i] * g_qual_scale;
			qual[i] = c < 93? c : 93;
		}
	}
	if (g_min_qlen > 0) {
		int k, qlen = 0;
		uint32_t *cigar = bam1_cigar(b);
		for (k = 0; k < b->core.n_cigar; ++k)
			if ((bam_cigar_type(bam_cigar_op(cigar[k]))&1) || bam_cigar_op(cigar[k]) == BAM_CHARD_CLIP)
				qlen += bam_cigar_oplen(cigar[k]);
		if (qlen < g_min_qlen) return 1;
	}
	if (b->core.qual < g_min_mapQ || ((b->core.flag & g_flag_on) != g_flag_on) || (b->core.flag & g_flag_off))
		return 1;
	if (g_bed && b->core.tid >= 0 && !bed_overlap(g_bed, h->target_name[b->core.tid], b->core.pos, bam_calend(&b->core, bam1_cigar(b))))
		return 1;
	if (g_subsam_frac > 0.) {
		uint32_t k = __ac_X31_hash_string(bam1_qname(b)) + g_subsam_seed;
		if ((double)(k&0xffffff) / 0x1000000 >= g_subsam_frac) return 1;
	}
	if (g_rg || g_rghash) {
		uint8_t *s = bam_aux_get(b, "RG");
		if (s) {
			if (g_rg) return (strcmp(g_rg, (char*)(s + 1)) == 0)? 0 : 1;
			if (g_rghash) {
				khint_t k = kh_get(rg, g_rghash, (char*)(s + 1));
				return (k != kh_end(g_rghash))? 0 : 1;
			}
		}
	}
	if (g_library) {
		const char *p = bam_get_library((bam_header_t*)h, b);
		return (p && strcmp(p, g_library) == 0)? 0 : 1;
	}
	return 0;
}

static char *drop_rg(char *hdtxt, rghash_t h, int *len)
{
	char *p = hdtxt, *q, *r, *s;
	kstring_t str;
	memset(&str, 0, sizeof(kstring_t));
	while (1) {
		int toprint = 0;
		q = strchr(p, '\n');
		if (q == 0) q = p + strlen(p);
		if (q - p < 3) break; // the line is too short; then stop
		if (strncmp(p, "@RG\t", 4) == 0) {
			int c;
			khint_t k;
			if ((r = strstr(p, "\tID:")) != 0) {
				r += 4;
				for (s = r; *s != '\0' && *s != '\n' && *s != '\t'; ++s);
				c = *s; *s = '\0';
				k = kh_get(rg, h, r);
				*s = c;
				if (k != kh_end(h)) toprint = 1;
			}
		} else toprint = 1;
		if (toprint) {
			kputsn(p, q - p, &str); kputc('\n', &str);
		}
		p = q + 1;
	}
	*len = str.l;
	return str.s;
}

// callback function for bam_fetch() that prints nonskipped records
static int view_func(const bam1_t *b, void *data)
{
	if (!process_aln(((samfile_t*)data)->header, (bam1_t*)b))
		samwrite((samfile_t*)data, b);
	return 0;
}

// callback function for bam_fetch() that counts nonskipped records
static int count_func(const bam1_t *b, void *data)
{
	if (!process_aln(((count_func_data_t*)data)->header, (bam1_t*)b)) {
		(*((count_func_data_t*)data)->count)++;
	}
	return 0;
}

static int usage(int is_long_help);

int main_samview(int argc, char *argv[])
{
	int c, is_header = 0, is_header_only = 0, is_bamin = 1, ret = 0, compress_level = -1, is_bamout = 0, is_count = 0;
	int of_type = BAM_OFDEC, is_long_help = 0, n_threads = 0;
	int64_t count = 0;
	samfile_t *in = 0, *out = 0;
	char in_mode[5], out_mode[5], *fn_out = 0, *fn_list = 0, *fn_ref = 0, *fn_rg = 0, *q;

	/* parse command-line options */
	strcpy(in_mode, "r"); strcpy(out_mode, "w");
	while ((c = getopt(argc, argv, "SbBct:h1Ho:q:f:F:ul:r:xX?T:R:L:s:Q:@:m:")) >= 0) {
		switch (c) {
		case 's':
			if ((g_subsam_seed = strtol(optarg, &q, 10)) != 0) {
				srand(g_subsam_seed);
				g_subsam_seed = rand();
			}
			g_subsam_frac = strtod(q, &q);
			break;
		case 'm': g_min_qlen = atoi(optarg); break;
		case 'c': is_count = 1; break;
		case 'S': is_bamin = 0; break;
		case 'b': is_bamout = 1; break;
		case 't': fn_list = strdup(optarg); is_bamin = 0; break;
		case 'h': is_header = 1; break;
		case 'H': is_header_only = 1; break;
		case 'o': fn_out = strdup(optarg); break;
		case 'f': g_flag_on = strtol(optarg, 0, 0); break;
		case 'F': g_flag_off = strtol(optarg, 0, 0); break;
		case 'q': g_min_mapQ = atoi(optarg); break;
		case 'u': compress_level = 0; break;
		case '1': compress_level = 1; break;
		case 'l': g_library = strdup(optarg); break;
		case 'L': g_bed = bed_read(optarg); break;
		case 'r': g_rg = strdup(optarg); break;
		case 'R': fn_rg = strdup(optarg); break;
		case 'x': of_type = BAM_OFHEX; break;
		case 'X': of_type = BAM_OFSTR; break;
		case '?': is_long_help = 1; break;
		case 'T': fn_ref = strdup(optarg); is_bamin = 0; break;
		case 'B': bam_no_B = 1; break;
		case 'Q': g_qual_scale = atoi(optarg); break;
		case '@': n_threads = strtol(optarg, 0, 0); break;
		default: return usage(is_long_help);
		}
	}
	if (compress_level >= 0) is_bamout = 1;
	if (is_header_only) is_header = 1;
	if (is_bamout) strcat(out_mode, "b");
	else {
		if (of_type == BAM_OFHEX) strcat(out_mode, "x");
		else if (of_type == BAM_OFSTR) strcat(out_mode, "X");
	}
	if (is_bamin) strcat(in_mode, "b");
	if (is_header) strcat(out_mode, "h");
	if (compress_level >= 0) {
		char tmp[2];
		tmp[0] = compress_level + '0'; tmp[1] = '\0';
		strcat(out_mode, tmp);
	}
	if (argc == optind) return usage(is_long_help); // potential memory leak...

	// read the list of read groups
	if (fn_rg) {
		FILE *fp_rg;
		char buf[1024];
		int ret;
		g_rghash = kh_init(rg);
		fp_rg = fopen(fn_rg, "r");
		while (!feof(fp_rg) && fscanf(fp_rg, "%s", buf) > 0) // this is not a good style, but bear me...
			kh_put(rg, g_rghash, strdup(buf), &ret); // we'd better check duplicates...
		fclose(fp_rg);
	}

	// generate the fn_list if necessary
	if (fn_list == 0 && fn_ref) fn_list = samfaipath(fn_ref);
	// open file handlers
	if ((in = samopen(argv[optind], in_mode, fn_list)) == 0) {
		fprintf(stderr, "[main_samview] fail to open \"%s\" for reading.\n", argv[optind]);
		ret = 1;
		goto view_end;
	}
	if (in->header == 0) {
		fprintf(stderr, "[main_samview] fail to read the header from \"%s\".\n", argv[optind]);
		ret = 1;
		goto view_end;
	}
	if (g_rghash) { // FIXME: I do not know what "bam_header_t::n_text" is for...
		char *tmp;
		int l;
		tmp = drop_rg(in->header->text, g_rghash, &l);
		free(in->header->text);
		in->header->text = tmp;
		in->header->l_text = l;
	}
	if (!is_count && (out = samopen(fn_out? fn_out : "-", out_mode, in->header)) == 0) {
		fprintf(stderr, "[main_samview] fail to open \"%s\" for writing.\n", fn_out? fn_out : "standard output");
		ret = 1;
		goto view_end;
	}
	if (n_threads > 1) samthreads(out, n_threads, 256); 
	if (is_header_only) goto view_end; // no need to print alignments

	if (argc == optind + 1) { // convert/print the entire file
		bam1_t *b = bam_init1();
		int r;
		while ((r = samread(in, b)) >= 0) { // read one alignment from `in'
			if (!process_aln(in->header, b)) {
				if (!is_count) samwrite(out, b); // write the alignment to `out'
				count++;
			}
		}
		if (r < -1) {
			fprintf(stderr, "[main_samview] truncated file.\n");
			ret = 1;
		}
		bam_destroy1(b);
	} else { // retrieve alignments in specified regions
		int i;
		bam_index_t *idx = 0;
		if (is_bamin) idx = bam_index_load(argv[optind]); // load BAM index
		if (idx == 0) { // index is unavailable
			fprintf(stderr, "[main_samview] random alignment retrieval only works for indexed BAM files.\n");
			ret = 1;
			goto view_end;
		}
		for (i = optind + 1; i < argc; ++i) {
			int tid, beg, end, result;
			bam_parse_region(in->header, argv[i], &tid, &beg, &end); // parse a region in the format like `chr2:100-200'
			if (tid < 0) { // reference name is not found
				fprintf(stderr, "[main_samview] region \"%s\" specifies an unknown reference name. Continue anyway.\n", argv[i]);
				continue;
			}
			// fetch alignments
			if (is_count) {
				count_func_data_t count_data = { in->header, &count };
				result = bam_fetch(in->x.bam, idx, tid, beg, end, &count_data, count_func);
			} else
				result = bam_fetch(in->x.bam, idx, tid, beg, end, out, view_func);
			if (result < 0) {
				fprintf(stderr, "[main_samview] retrieval of region \"%s\" failed due to truncated file or corrupt BAM index file\n", argv[i]);
				ret = 1;
				break;
			}
		}
		bam_index_destroy(idx); // destroy the BAM index
	}

view_end:
	if (is_count && ret == 0) 
		printf("%" PRId64 "\n", count);

	// close files, free and return
	free(fn_list); free(fn_ref); free(fn_out); free(g_library); free(g_rg); free(fn_rg);
	if (g_bed) bed_destroy(g_bed);
	if (g_rghash) {
		khint_t k;
		for (k = 0; k < kh_end(g_rghash); ++k)
			if (kh_exist(g_rghash, k)) free((char*)kh_key(g_rghash, k));
		kh_destroy(rg, g_rghash);
	}
	samclose(in);
	if (!is_count)
		samclose(out);
	return ret;
}

static int usage(int is_long_help)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   samtools view [options] <in.bam>|<in.sam> [region1 [...]]\n\n");
	fprintf(stderr, "Options: -b       output BAM\n");
	fprintf(stderr, "         -h       print header for the SAM output\n");
	fprintf(stderr, "         -H       print header only (no alignments)\n");
	fprintf(stderr, "         -S       input is SAM\n");
	fprintf(stderr, "         -u       uncompressed BAM output (force -b)\n");
	fprintf(stderr, "         -1       fast compression (force -b)\n");
	fprintf(stderr, "         -x       output FLAG in HEX (samtools-C specific)\n");
	fprintf(stderr, "         -X       output FLAG in string (samtools-C specific)\n");
	fprintf(stderr, "         -c       print only the count of matching records\n");
	fprintf(stderr, "         -B       collapse the backward CIGAR operation\n");
	fprintf(stderr, "         -@ INT   number of BAM compression threads [0]\n");
	fprintf(stderr, "         -L FILE  output alignments overlapping the input BED FILE [null]\n");
	fprintf(stderr, "         -t FILE  list of reference names and lengths (force -S) [null]\n");
	fprintf(stderr, "         -T FILE  reference sequence file (force -S) [null]\n");
	fprintf(stderr, "         -o FILE  output file name [stdout]\n");
	fprintf(stderr, "         -R FILE  list of read groups to be outputted [null]\n");
	fprintf(stderr, "         -f INT   required flag, 0 for unset [0]\n");
	fprintf(stderr, "         -F INT   filtering flag, 0 for unset [0]\n");
	fprintf(stderr, "         -q INT   minimum mapping quality [0]\n");
	fprintf(stderr, "         -l STR   only output reads in library STR [null]\n");
	fprintf(stderr, "         -r STR   only output reads in read group STR [null]\n");
	fprintf(stderr, "         -s FLOAT fraction of templates to subsample; integer part as seed [-1]\n");
	fprintf(stderr, "         -?       longer help\n");
	fprintf(stderr, "\n");
	if (is_long_help)
		fprintf(stderr, "Notes:\n\
\n\
  1. By default, this command assumes the file on the command line is in\n\
     the BAM format and it prints the alignments in SAM. If `-t' is\n\
     applied, the input file is assumed to be in the SAM format. The\n\
     file supplied with `-t' is SPACE/TAB delimited with the first two\n\
     fields of each line consisting of the reference name and the\n\
     corresponding sequence length. The `.fai' file generated by `faidx'\n\
     can be used here. This file may be empty if reads are unaligned.\n\
\n\
  2. SAM->BAM conversion: `samtools view -bT ref.fa in.sam.gz'.\n\
\n\
  3. BAM->SAM conversion: `samtools view in.bam'.\n\
\n\
  4. A region should be presented in one of the following formats:\n\
     `chr1', `chr2:1,000' and `chr3:1000-2,000'. When a region is\n\
     specified, the input alignment file must be an indexed BAM file.\n\
\n\
  5. Option `-u' is preferred over `-b' when the output is piped to\n\
     another samtools command.\n\
\n\
  6. In a string FLAG, each character represents one bit with\n\
     p=0x1 (paired), P=0x2 (properly paired), u=0x4 (unmapped),\n\
     U=0x8 (mate unmapped), r=0x10 (reverse), R=0x20 (mate reverse)\n\
     1=0x40 (first), 2=0x80 (second), s=0x100 (not primary), \n\
     f=0x200 (failure) and d=0x400 (duplicate). Note that `-x' and\n\
     `-X' are samtools-C specific. Picard and older samtools do not\n\
     support HEX or string flags.\n\
\n");
	return 1;
}

int main_import(int argc, char *argv[])
{
	int argc2, ret;
	char **argv2;
	if (argc != 4) {
		fprintf(stderr, "Usage: bamtk import <in.ref_list> <in.sam> <out.bam>\n");
		return 1;
	}
	argc2 = 6;
	argv2 = calloc(6, sizeof(char*));
	argv2[0] = "import", argv2[1] = "-o", argv2[2] = argv[3], argv2[3] = "-bt", argv2[4] = argv[1], argv2[5] = argv[2];
	ret = main_samview(argc2, argv2);
	free(argv2);
	return ret;
}

int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };

int main_bam2fq(int argc, char *argv[])
{
	bamFile fp;
	bam_header_t *h;
	bam1_t *b;
	int8_t *buf;
	int max_buf, c, no12 = 0;
	while ((c = getopt(argc, argv, "n")) > 0)
		if (c == 'n') no12 = 1;
	if (argc == 1) {
		fprintf(stderr, "Usage: samtools bam2fq <in.bam>\n");
		return 1;
	}
	fp = strcmp(argv[optind], "-")? bam_open(argv[optind], "r") : bam_dopen(fileno(stdin), "r");
	if (fp == 0) return 1;
	h = bam_header_read(fp);
	b = bam_init1();
	buf = 0;
	max_buf = 0;
	while (bam_read1(fp, b) >= 0) {
		int i, qlen = b->core.l_qseq;
		uint8_t *seq;
		putchar('@'); fputs(bam1_qname(b), stdout);
		if (no12) putchar('\n');
		else {
			if ((b->core.flag & 0x40) && !(b->core.flag & 0x80)) puts("/1");
			else if ((b->core.flag & 0x80) && !(b->core.flag & 0x40)) puts("/2");
			else putchar('\n');
		}
		if (max_buf < qlen + 1) {
			max_buf = qlen + 1;
			kroundup32(max_buf);
			buf = realloc(buf, max_buf);
		}
		buf[qlen] = 0;
		seq = bam1_seq(b);
		for (i = 0; i < qlen; ++i)
			buf[i] = bam1_seqi(seq, i);
		if (b->core.flag & 16) { // reverse complement
			for (i = 0; i < qlen>>1; ++i) {
				int8_t t = seq_comp_table[buf[qlen - 1 - i]];
				buf[qlen - 1 - i] = seq_comp_table[buf[i]];
				buf[i] = t;
			}
			if (qlen&1) buf[i] = seq_comp_table[buf[i]];
		}
		for (i = 0; i < qlen; ++i)
			buf[i] = bam_nt16_rev_table[buf[i]];
		puts((char*)buf);
		puts("+");
		seq = bam1_qual(b);
		for (i = 0; i < qlen; ++i)
			buf[i] = 33 + seq[i];
		if (b->core.flag & 16) { // reverse
			for (i = 0; i < qlen>>1; ++i) {
				int8_t t = buf[qlen - 1 - i];
				buf[qlen - 1 - i] = buf[i];
				buf[i] = t;
			}
		}
		puts((char*)buf);
	}
	free(buf);
	bam_destroy1(b);
	bam_header_destroy(h);
	bam_close(fp);
	return 0;
}
