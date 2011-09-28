/* The MIT License

   Copyright (c) 2011  Heng Li <lh3@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "kstring.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

#define N_TMPSTR 5
#define LINE_LEN 60

// append a CIGAR operation plus length
#define write_cigar(_c, _n, _m, _v) do { \
		if (_n == _m) { \
			_m = _m? _m<<1 : 4; \
			_c = realloc(_c, _m * sizeof(unsigned)); \
		} \
		_c[_n++] = (_v); \
	} while (0)

// a fatal error
static void fatal(const char *msg)
{
	fprintf(stderr, "E %s\n", msg);
	exit(1);
}
// remove pads
static void remove_pads(const kstring_t *src, kstring_t *dst)
{
	int i, j;
	dst->l = 0;
	kputsn(src->s, src->l, dst);
	for (i = j = 0; i < dst->l; ++i)
		if (dst->s[i] != '*') dst->s[j++] = dst->s[i];
	dst->s[j] = 0;
	dst->l = j;
}

int main(int argc, char *argv[])
{
	gzFile fp;
	kstream_t *ks;
	kstring_t s, t[N_TMPSTR];
	int dret, i, k, af_n, af_max, af_i, c, is_padded = 0, write_cns = 0, *p2u = 0;
	long m_cigar = 0, n_cigar = 0;
	unsigned *af, *cigar = 0;

	while ((c = getopt(argc, argv, "pc")) >= 0) {
		switch (c) {
			case 'p': is_padded = 1; break;
			case 'c': write_cns = 1; break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "\nUsage:   ace2sam [-pc] <in.ace>\n\n");
		fprintf(stderr, "Options: -p     output padded SAM\n");
		fprintf(stderr, "         -c     write the contig sequence in SAM\n\n");
		fprintf(stderr, "Notes: 1. Fields must appear in the following order: (CO->[BQ]->(AF)->(RD->QA))\n");
		fprintf(stderr, "       2. The order of reads in AF and in RD must be identical\n");
		fprintf(stderr, "       3. Except in BQ, words and numbers must be separated by a single SPACE or TAB\n");
		fprintf(stderr, "       4. This program writes the headerless SAM to stdout and header to stderr\n\n");
		return 1;
	}

	s.l = s.m = 0; s.s = 0;
	af_n = af_max = af_i = 0; af = 0;
	for (i = 0; i < N_TMPSTR; ++i) t[i].l = t[i].m = 0, t[i].s = 0;
	fp = strcmp(argv[1], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	while (ks_getuntil(ks, 0, &s, &dret) >= 0) {
		if (strcmp(s.s, "CO") == 0) { // contig sequence
			kstring_t *cns;
			t[0].l = t[1].l = t[2].l = t[3].l = t[4].l = 0; // 0: name; 1: padded ctg; 2: unpadded ctg/padded read; 3: unpadded read; 4: SAM line
			af_n = af_i = 0; // reset the af array
			ks_getuntil(ks, 0, &s, &dret); kputs(s.s, &t[0]); // contig name
			ks_getuntil(ks, '\n', &s, &dret); // read the whole line
			while (ks_getuntil(ks, '\n', &s, &dret) >= 0 && s.l > 0) kputsn(s.s, s.l, &t[1]); // read the padded consensus sequence
			remove_pads(&t[1], &t[2]); // construct the unpadded sequence
			// compute the array for mapping padded positions to unpadded positions
			p2u = realloc(p2u, t[1].m * sizeof(int));
			for (i = k = 0; i < t[1].l; ++i) {
				p2u[i] = k;
				if (t[1].s[i] != '*') ++k;
			}
			// write out the SAM header and contig sequences
			fprintf(stderr, "H @SQ\tSN:%s\tLN:%ld\n", t[0].s, t[is_padded?1:2].l); // The SAM header line
			cns = &t[is_padded?1:2];
			fprintf(stderr, "S >%s\n", t[0].s);
			for (i = 0; i < cns->l; i += LINE_LEN) {
				fputs("S ", stderr);
				for (k = 0; k < LINE_LEN && i + k < cns->l; ++k)
					fputc(cns->s[i + k], stderr);
				fputc('\n', stderr);
			}

#define __padded2cigar(sp) do { \
		int i, l_M = 0, l_D = 0; \
		for (i = 0; i < sp.l; ++i) { \
			if (sp.s[i] == '*') { \
				if (l_M) write_cigar(cigar, n_cigar, m_cigar, l_M<<4); \
				++l_D; l_M = 0; \
			} else { \
				if (l_D) write_cigar(cigar, n_cigar, m_cigar, l_D<<4 | 2); \
				++l_M; l_D = 0; \
			} \
		} \
		if (l_M) write_cigar(cigar, n_cigar, m_cigar, l_M<<4); \
		else write_cigar(cigar, n_cigar, m_cigar, l_D<<4 | 2); \
	} while (0)

			if (write_cns) { // write the consensus SAM line (dummy read)
				n_cigar = 0;
				if (is_padded) __padded2cigar(t[1]);
				else write_cigar(cigar, n_cigar, m_cigar, t[2].l<<4);
				kputsn(t[0].s, t[0].l, &t[4]); kputs("\t516\t", &t[4]); kputsn(t[0].s, t[0].l, &t[4]); kputs("\t1\t60\t", &t[4]);
				for (i = 0; i < n_cigar; ++i) {
					kputw(cigar[i]>>4, &t[4]); kputc("MIDNSHP=X"[cigar[i]&0xf], &t[4]);
				}
				kputs("\t*\t0\t0\t", &t[4]); kputsn(t[2].s, t[2].l, &t[4]); kputs("\t*", &t[4]);
			}
		} else if (strcmp(s.s, "BQ") == 0) { // contig quality
			if (t[0].l == 0) fatal("come to 'BQ' before reading 'CO'");
			if (dret != '\n') ks_getuntil(ks, '\n', &s, &dret); // read the entire "BQ" line
			if (write_cns) t[4].s[--t[4].l] = 0; // remove the trailing "*"
			for (i = 0; i < t[2].l; ++i) { // read the consensus quality
				int q;
				if (ks_getuntil(ks, 0, &s, &dret) < 0) fprintf(stderr, "E truncated contig quality\n");
				if (s.l) {
					q = atoi(s.s) + 33;
					if (q > 126) q = 126;
					if (write_cns) kputc(q, &t[4]);
				} else --i;
			}
			if (dret != '\n') ks_getuntil(ks, '\n', &s, &dret);
			ks_getuntil(ks, '\n', &s, &dret); // skip the empty line
			if (write_cns) puts(t[4].s); t[4].l = 0;
		} else if (strcmp(s.s, "AF") == 0) { // padded read position
			int reversed, neg, pos;
			if (t[0].l == 0) fatal("come to 'AF' before reading 'CO'");
			if (write_cns) {
				if (t[4].l) puts(t[4].s);
				t[4].l = 0;
			}
			ks_getuntil(ks, 0, &s, &dret); // read name
			ks_getuntil(ks, 0, &s, &dret); reversed = s.s[0] == 'C'? 1 : 0; // strand
			ks_getuntil(ks, 0, &s, &dret); pos = atoi(s.s); neg = pos < 0? 1 : 0; pos = pos < 0? -pos : pos; // position
			if (af_n == af_max) { // double the af array
				af_max = af_max? af_max<<1 : 4;
				af = realloc(af, af_max * sizeof(unsigned));
			}
			af[af_n++] = pos << 2 | neg << 1 | reversed; // keep the placement information
		} else if (strcmp(s.s, "RD") == 0) { // read sequence
			if (af_i >= af_n) fatal("more 'RD' records than 'AF'");
			t[2].l = t[3].l = t[4].l = 0;
			ks_getuntil(ks, 0, &t[4], &dret); // QNAME
			if (dret != '\n') ks_getuntil(ks, '\n', &s, &dret); // read the entire RD line
			while (ks_getuntil(ks, '\n', &s, &dret) >= 0 && s.l > 0) kputs(s.s, &t[2]); // read the read sequence
		} else if (strcmp(s.s, "QA") == 0) { // clipping
			if (af_i >= af_n) fatal("more 'QA' records than 'AF'");
			int beg, end, pos, op;
			ks_getuntil(ks, 0, &s, &dret); ks_getuntil(ks, 0, &s, &dret); // skip quality clipping
			ks_getuntil(ks, 0, &s, &dret); beg = atoi(s.s) - 1; // align clipping start
			ks_getuntil(ks, 0, &s, &dret); end = atoi(s.s); // clipping end
			if (dret != '\n') ks_getuntil(ks, '\n', &s, &dret);
			// compute 1-based POS
			pos = af[af_i]>>2; // retrieve the position information
			if (af[af_i]>>1&1) pos = -pos;
			pos += beg; // now pos is the true padded position
			// generate CIGAR
			remove_pads(&t[2], &t[3]); // backup the unpadded read sequence
			n_cigar = 0;
			if (beg) write_cigar(cigar, n_cigar, m_cigar, beg<<4|4);
			if (is_padded) {
				__padded2cigar(t[2]);
				if (beg && n_cigar > 1) cigar[1] -= beg<<4; // fix the left-hand CIGAR 
				if (end < t[2].l && n_cigar) cigar[n_cigar-1] -= (t[2].l - end)<<4; // fix the right-hand CIGAR
			} else {
				// generate flattened CIGAR string
				for (i = beg, k = pos - 1; i < end; ++i, ++k)
					t[2].s[i] = t[2].s[i] != '*'? (t[1].s[k] != '*'? 0 : 1) : (t[1].s[k] != '*'? 2 : 6);
				// generate the proper CIGAR
				for (i = beg + 1, k = 1, op = t[2].s[beg]; i < end; ++i) {
					if (op != t[2].s[i]) {
						write_cigar(cigar, n_cigar, m_cigar, k<<4|op);
						op = t[2].s[i]; k = 1;
					} else ++k;
				}
				write_cigar(cigar, n_cigar, m_cigar, k<<4|op);
				// remove unnecessary "P" and possibly merge adjacent operations
				for (i = 2; i < n_cigar; ++i) {
					if ((cigar[i]&0xf) != 1 && (cigar[i-1]&0xf) == 6 && (cigar[i-2]&0xf) != 1) {
						cigar[i-1] = 0;
						if ((cigar[i]&0xf) == (cigar[i-2]&0xf)) // merge operations
							cigar[i] += cigar[i-2], cigar[i-2] = 0;
					}
				}
				for (i = k = 0; i < n_cigar; ++i) // squeeze out dumb operations
					if (cigar[i]) cigar[k++] = cigar[i];
				n_cigar = k;
			}
			if (end < t[2].l) write_cigar(cigar, n_cigar, m_cigar, (t[2].l - end)<<4|4);
			// write the SAM line for the read
			kputc('\t', &t[4]); // QNAME has already been written
			kputw((af[af_i]&1)? 16 : 0, &t[4]); kputc('\t', &t[4]); // FLAG
			kputsn(t[0].s, t[0].l, &t[4]); kputc('\t', &t[4]); // RNAME
			kputw(is_padded? pos : p2u[pos-1]+1, &t[4]); // POS
			kputs("\t60\t", &t[4]); // MAPQ
			for (i = 0; i < n_cigar; ++i) { // CIGAR
				kputw(cigar[i]>>4, &t[4]); kputc("MIDNSHP=X"[cigar[i]&0xf], &t[4]);
			}
			kputs("\t*\t0\t0\t", &t[4]); // empty MRNM, MPOS and TLEN
			kputsn(t[3].s, t[3].l, &t[4]); // unpadded SEQ
			kputs("\t*", &t[4]); // QUAL
			puts(t[4].s); // print to stdout
			++af_i;
		} else if (dret != '\n') ks_getuntil(ks, '\n', &s, &dret);
	}
	ks_destroy(ks);
	gzclose(fp);
	free(af); free(s.s); free(cigar); free(p2u);
	for (i = 0; i < N_TMPSTR; ++i) free(t[i].s);
	return 0;
}
