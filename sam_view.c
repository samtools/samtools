#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
KHASH_SET_INIT_STR(rg)

typedef khash_t(rg) *rghash_t;

// This structure contains the settings for a samview run
typedef struct samview_settings {
	rghash_t rghash;
	int min_mapQ;
	int flag_on;
	int flag_off;
	int qual_scale;
	int min_qlen;
	int remove_B;
	uint32_t subsam_seed;
	double subsam_frac;
	char* library;
	char* rg;
	void* bed;
	size_t remove_aux_len;
	char** remove_aux;
} samview_settings_t;


// TODO Add declarations of these to a viable htslib or samtools header
extern const char *bam_get_library(bam_hdr_t *header, const bam1_t *b);
extern int bam_remove_B(bam1_t *b);
extern char *samfaipath(const char *fn_ref);
void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);

// Returns 0 to indicate read should be output 1 otherwise
static int process_aln(const bam_hdr_t *h, bam1_t *b, samview_settings_t* settings)
{
	if (settings->remove_B) bam_remove_B(b);
	if (settings->qual_scale > 1) {
		int i;
		uint8_t *qual = bam_get_qual(b);
		for (i = 0; i < b->core.l_qseq; ++i) {
			int c = qual[i] * settings->qual_scale;
			qual[i] = c < 93? c : 93;
		}
	}
	if (settings->min_qlen > 0) {
		int k, qlen = 0;
		uint32_t *cigar = bam_get_cigar(b);
		for (k = 0; k < b->core.n_cigar; ++k)
			if ((bam_cigar_type(bam_cigar_op(cigar[k]))&1) || bam_cigar_op(cigar[k]) == BAM_CHARD_CLIP)
				qlen += bam_cigar_oplen(cigar[k]);
		if (qlen < settings->min_qlen) return 1;
	}
	if (b->core.qual < settings->min_mapQ || ((b->core.flag & settings->flag_on) != settings->flag_on) || (b->core.flag & settings->flag_off))
		return 1;
	if (settings->bed && b->core.tid >= 0 && !bed_overlap(settings->bed, h->target_name[b->core.tid], b->core.pos, bam_endpos(b)))
		return 1;
	if (settings->subsam_frac > 0.) {
		uint32_t k = __ac_X31_hash_string(bam_get_qname(b)) + settings->subsam_seed;
		if ((double)(k&0xffffff) / 0x1000000 >= settings->subsam_frac) return 1;
	}
	if (settings->rg || settings->rghash) {
		uint8_t *s = bam_aux_get(b, "RG");
		if (s) {
			if (settings->rg) return (strcmp(settings->rg, (char*)(s + 1)) == 0)? 0 : 1;
			if (settings->rghash) {
				khint_t k = kh_get(rg, settings->rghash, (char*)(s + 1));
				return (k != kh_end(settings->rghash))? 0 : 1;
			}
		}
	}
	if (settings->library) {
		const char *p = bam_get_library((bam_hdr_t*)h, b);
		return (p && strcmp(p, settings->library) == 0)? 0 : 1;
	}
	if (settings->remove_aux_len) {
		size_t i;
		for (i = 0; i < settings->remove_aux_len; ++i) {
			uint8_t *s = bam_aux_get(b, settings->remove_aux[i]);
			if (s) {
				bam_aux_del(b, s);
			}
		}
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

static int usage(int is_long_help);

int main_samview(int argc, char *argv[])
{
	int c, is_header = 0, is_header_only = 0, ret = 0, compress_level = -1, is_count = 0;
	int is_long_help = 0, n_threads = 0;
	int64_t count = 0;
	samFile *in = 0, *out = 0;
	bam_hdr_t *header;
	char out_mode[5], *out_format = "", *fn_out = 0, *fn_list = 0, *fn_ref = 0, *fn_rg = 0, *q;
	
	samview_settings_t settings = {
		.rghash = NULL,
		.min_mapQ = 0,
		.flag_on = 0,
		.flag_off = 0,
		.qual_scale = 0,
		.min_qlen = 0,
		.remove_B = 0,
		.subsam_seed = 0,
		.subsam_frac = -1.,
		.library = NULL,
		.rg = NULL,
		.bed = NULL,
	};

	/* parse command-line options */
	/* TODO: convert this to getopt_long we're running out of letters */
	strcpy(out_mode, "w");
	while ((c = getopt(argc, argv, "SbBcCt:h1Ho:q:f:F:ul:r:?T:R:L:s:Q:@:m:x:")) >= 0) {
		switch (c) {
		case 's':
			if ((settings.subsam_seed = strtol(optarg, &q, 10)) != 0) {
				srand(settings.subsam_seed);
				settings.subsam_seed = rand();
			}
			settings.subsam_frac = strtod(q, &q);
			break;
		case 'm': settings.min_qlen = atoi(optarg); break;
		case 'c': is_count = 1; break;
		case 'S': break;
		case 'b': out_format = "b"; break;
		case 'C': out_format = "c"; break;
		case 't': fn_list = strdup(optarg); break;
		case 'h': is_header = 1; break;
		case 'H': is_header_only = 1; break;
		case 'o': fn_out = strdup(optarg); break;
		case 'f': settings.flag_on |= strtol(optarg, 0, 0); break;
		case 'F': settings.flag_off |= strtol(optarg, 0, 0); break;
		case 'q': settings.min_mapQ = atoi(optarg); break;
		case 'u': compress_level = 0; break;
		case '1': compress_level = 1; break;
		case 'l': settings.library = strdup(optarg); break;
		case 'L': settings.bed = bed_read(optarg); break;
		case 'r': settings.rg = strdup(optarg); break;
		case 'R': fn_rg = strdup(optarg); break;
				/* REMOVED as htslib doesn't support this
		//case 'x': out_format = "x"; break;
		//case 'X': out_format = "X"; break;
				 */
		case '?': is_long_help = 1; break;
		case 'T': fn_ref = strdup(optarg); break;
		case 'B': settings.remove_B = 1; break;
		case 'Q': settings.qual_scale = atoi(optarg); break;
		case '@': n_threads = strtol(optarg, 0, 0); break;
		case 'x':
			{
				if (strlen(optarg) != 2) {
					fprintf(stderr, "main_samview: Error parsing -x auxiliary tags should be exactly two characters long.\n");
					return usage(is_long_help);
				}
				settings.remove_aux = (char**)realloc(settings.remove_aux, sizeof(char*) * (++settings.remove_aux_len));
				settings.remove_aux[settings.remove_aux_len-1] = optarg;
			}
			break;
		default: return usage(is_long_help);
		}
	}
	if (compress_level >= 0) out_format = "b";
	if (is_header_only) is_header = 1;
	strcat(out_mode, out_format);
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
		settings.rghash = kh_init(rg);
		fp_rg = fopen(fn_rg, "r");
		while (!feof(fp_rg) && fscanf(fp_rg, "%s", buf) > 0) // this is not a good style, but bear me...
			kh_put(rg, settings.rghash, strdup(buf), &ret); // we'd better check duplicates...
		fclose(fp_rg);
	}

	// generate the fn_list if necessary
	if (fn_list == 0 && fn_ref) fn_list = samfaipath(fn_ref);
	// open file handlers
	if ((in = sam_open(argv[optind], "r")) == 0) {
		fprintf(stderr, "[main_samview] fail to open \"%s\" for reading.\n", argv[optind]);
		ret = 1;
		goto view_end;
	}
	if (fn_list) hts_set_fai_filename(in, fn_list);
	if ((header = sam_hdr_read(in)) == 0) {
		fprintf(stderr, "[main_samview] fail to read the header from \"%s\".\n", argv[optind]);
		ret = 1;
		goto view_end;
	}
	if (settings.rghash) { // FIXME: I do not know what "bam_header_t::n_text" is for...
		char *tmp;
		int l;
		tmp = drop_rg(header->text, settings.rghash, &l);
		free(header->text);
		header->text = tmp;
		header->l_text = l;
	}
	if (!is_count) {
		if ((out = sam_open(fn_out? fn_out : "-", out_mode)) == 0) {
			fprintf(stderr, "[main_samview] fail to open \"%s\" for writing.\n", fn_out? fn_out : "standard output");
			ret = 1;
			goto view_end;
		}
		if (*out_format || is_header) sam_hdr_write(out, header);
	}
#if 0
	// TODO Add function for setting I/O threading to htslib API
	if (n_threads > 1) { samthreads(out, n_threads, 256); }
#endif
	if (is_header_only) goto view_end; // no need to print alignments

	if (argc == optind + 1) { // convert/print the entire file
		bam1_t *b = bam_init1();
		int r;
		while ((r = sam_read1(in, header, b)) >= 0) { // read one alignment from `in'
			if (!process_aln(header, b, &settings)) {
				if (!is_count) sam_write1(out, header, b); // write the alignment to `out'
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
		bam1_t *b;
		hts_idx_t *idx = sam_index_load(in, argv[optind]); // load index
		if (idx == 0) { // index is unavailable
			fprintf(stderr, "[main_samview] random alignment retrieval only works for indexed BAM or CRAM files.\n");
			ret = 1;
			goto view_end;
		}
		b = bam_init1();
		for (i = optind + 1; i < argc; ++i) {
			int result;
			hts_itr_t *iter = sam_itr_querys(idx, header, argv[i]); // parse a region in the format like `chr2:100-200'
			if (iter == NULL) { // reference name is not found
				fprintf(stderr, "[main_samview] region \"%s\" specifies an unknown reference name. Continue anyway.\n", argv[i]);
				continue;
			}
			// fetch alignments
			while ((result = sam_itr_next(in, iter, b)) >= 0) {
				if (!process_aln(header, b, &settings)) {
					if (!is_count) sam_write1(out, header, b); // write the alignment to `out'
					count++;
				}
			}
			hts_itr_destroy(iter);
			if (result < -1) {
				fprintf(stderr, "[main_samview] retrieval of region \"%s\" failed due to truncated file or corrupt BAM index file\n", argv[i]);
				ret = 1;
				break;
			}
		}
		bam_destroy1(b);
		hts_idx_destroy(idx); // destroy the BAM index
	}

view_end:
	if (is_count && ret == 0) 
		printf("%" PRId64 "\n", count);

	// close files, free and return
	free(fn_list); free(fn_ref); free(fn_out); free(settings.library); free(settings.rg); free(fn_rg);
	if (settings.bed) bed_destroy(settings.bed);
	if (settings.rghash) {
		khint_t k;
		for (k = 0; k < kh_end(settings.rghash); ++k)
			if (kh_exist(settings.rghash, k)) free((char*)kh_key(settings.rghash, k));
		kh_destroy(rg, settings.rghash);
	}
	if (settings.remove_aux_len) {
		free(settings.remove_aux);
	}
	sam_close(in);
	if (!is_count)
		sam_close(out);
	return ret;
}

static int usage(int is_long_help)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   samtools view [options] <in.bam>|<in.sam> [region1 [...]]\n\n");
	fprintf(stderr, "Options: -b       output BAM\n");
	fprintf(stderr, "         -C       output CRAM\n");
	fprintf(stderr, "         -h       print header for the SAM output\n");
	fprintf(stderr, "         -H       print header only (no alignments)\n");
	fprintf(stderr, "         -S       ignored (input format is auto-detected)\n");
	fprintf(stderr, "         -u       uncompressed BAM output (force -b)\n");
	fprintf(stderr, "         -1       fast compression (force -b)\n");
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
	fprintf(stderr, "         -x STR   read tag to strip (repeatable) [null]\n");
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
		fprintf(stderr, "Usage: samtools import <in.ref_list> <in.sam> <out.bam>\n");
		return 1;
	}
	argc2 = 6;
	argv2 = calloc(6, sizeof(char*));
	argv2[0] = "import", argv2[1] = "-o", argv2[2] = argv[3], argv2[3] = "-bt", argv2[4] = argv[1], argv2[5] = argv[2];
	ret = main_samview(argc2, argv2);
	free(argv2);
	return ret;
}

int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

int main_bam2fq(int argc, char *argv[])
{
	samFile *fp;
	bam_hdr_t *h;
	bam1_t *b;
	int8_t *buf;
	int max_buf, c, no12 = 0;
	while ((c = getopt(argc, argv, "n")) > 0)
		if (c == 'n') no12 = 1;
	if (argc == 1) {
		fprintf(stderr, "Usage: samtools bam2fq <in.bam>\n");
		return 1;
	}
	fp = sam_open(argv[optind], "r");
	if (fp == 0) return 1;
	h = sam_hdr_read(fp);
	b = bam_init1();
	buf = 0;
	max_buf = 0;
	while (sam_read1(fp, h, b) >= 0) {
		int i, qlen = b->core.l_qseq;
		uint8_t *seq;
		putchar('@'); fputs(bam_get_qname(b), stdout);
		if (no12) putchar('\n');
		else {
			if ((b->core.flag & BAM_FREAD1) && !(b->core.flag & BAM_FREAD2)) puts("/1");
			else if ((b->core.flag & BAM_FREAD2) && !(b->core.flag & BAM_FREAD1)) puts("/2");
			else putchar('\n');
		}
		if (max_buf < qlen + 1) {
			max_buf = qlen + 1;
			kroundup32(max_buf);
			buf = realloc(buf, max_buf);
		}
		buf[qlen] = 0;
		seq = bam_get_seq(b);
		for (i = 0; i < qlen; ++i)
			buf[i] = bam_seqi(seq, i);
		if (b->core.flag & BAM_FREVERSE) { // reverse complement
			for (i = 0; i < qlen>>1; ++i) {
				int8_t t = seq_comp_table[buf[qlen - 1 - i]];
				buf[qlen - 1 - i] = seq_comp_table[buf[i]];
				buf[i] = t;
			}
			if (qlen&1) buf[i] = seq_comp_table[buf[i]];
		}
		for (i = 0; i < qlen; ++i)
			buf[i] = seq_nt16_str[buf[i]];
		puts((char*)buf);
		puts("+");
		seq = bam_get_qual(b);
		for (i = 0; i < qlen; ++i)
			buf[i] = 33 + seq[i];
		if (b->core.flag & BAM_FREVERSE) { // reverse
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
	bam_hdr_destroy(h);
	sam_close(fp);
	return 0;
}
