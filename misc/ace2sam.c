#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "kstring.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

#define N_TMPSTR 5

#define write_cigar(_c, _n, _m, _v) do { \
		if (_n == _m) { \
			_m = _m? _m<<1 : 4; \
			_c = realloc(_c, _m * sizeof(unsigned)); \
		} \
		_c[_n++] = (_v); \
	} while (0)

int main(int argc, char *argv[])
{
	gzFile fp;
	kstream_t *ks;
	kstring_t s, t[N_TMPSTR];
	int dret, i, af_n, af_max, af_i;
	long n_ctgs = 0, n_reads = 0, m_cigar = 0, n_cigar = 0;
	unsigned *af, *cigar = 0;

	if (argc == 1) {
		fprintf(stderr, "\nUsage: ace2sam <in.ace>\n\n");
		fprintf(stderr, "  ace2sam in.ace > a2s.out 2> a2s.err\n");
		fprintf(stderr, "  (grep ^H a2s.err | sed s,^..,,; cat a2s.out) > a2s.sam\n");
		fprintf(stderr, "  grep ^S a2s.err | sed s,^..,, > a2s.fa\n\n");
		return 1;
	}
	s.l = s.m = 0; s.s = 0;
	af_n = af_max = af_i = 0; af = 0;
	for (i = 0; i < N_TMPSTR; ++i) t[i].l = t[i].m = 0, t[i].s = 0;
	fp = strcmp(argv[1], "-")? gzopen(argv[1], "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	while (ks_getuntil(ks, 0, &s, &dret) >= 0) {
		if (strcmp(s.s, "AS") == 0) {
			ks_getuntil(ks, 0, &s, &dret); n_ctgs = atol(s.s);
			ks_getuntil(ks, 0, &s, &dret); n_reads = atol(s.s);
			if (dret != '\n') ks_getuntil(ks, '\n', &s, &dret);
		} else if (strcmp(s.s, "CO") == 0) {
			t[0].l = t[1].l = t[2].l = t[3].l = t[4].l = 0; // 0: name; 1: padded; 2: unpadded; 3: CIGAR/flat-cigar; 4: SAM line
			af_n = af_i = 0;
			ks_getuntil(ks, 0, &s, &dret); kputs(s.s, &t[0]);
			ks_getuntil(ks, 0, &s, &dret);
			ks_getuntil(ks, 0, &s, &dret); n_reads = atoi(s.s);
			ks_getuntil(ks, '\n', &s, &dret);
			fprintf(stderr, "S >%s\n", t[0].s);
			while (ks_getuntil(ks, '\n', &s, &dret) >= 0 && s.l > 0) {
				kputs(s.s, &t[1]);
				fputs("S ", stderr); fputs(s.s, stderr); fputc('\n', stderr);
			}

#define __padded2cigar(sp, su) do { \
		int i, l_M = 0, l_D = 0; \
		kputsn(sp.s, sp.l, &su); su.l = 0; \
		for (i = 0; i < sp.l; ++i) { \
			if (sp.s[i] == '*') { \
				if (l_M) write_cigar(cigar, n_cigar, m_cigar, l_M<<4); \
				++l_D; l_M = 0; \
			} else { \
				if (l_D) write_cigar(cigar, n_cigar, m_cigar, l_D<<4 | 2); \
				++l_M; l_D = 0; \
				su.s[su.l++] = sp.s[i]; \
			} \
		} \
		su.s[su.l] = 0; \
		if (l_M) write_cigar(cigar, n_cigar, m_cigar, l_M<<4); \
		else write_cigar(cigar, n_cigar, m_cigar, l_D<<4 | 2); \
	} while (0)

			n_cigar = 0;
			__padded2cigar(t[1], t[2]);
			for (i = 0; i < n_cigar; ++i) {
				kputw(cigar[i]>>4, &t[3]); kputc("MIDNSHP=X"[cigar[i]&0xf], &t[3]);
			}
			kputsn(t[0].s, t[0].l, &t[4]); kputs("\t516\t", &t[4]); kputsn(t[0].s, t[0].l, &t[4]); kputs("\t1\t60\t", &t[4]);
			kputsn(t[3].s, t[3].l, &t[4]); kputs("\t*\t0\t0\t", &t[4]); kputsn(t[2].s, t[2].l, &t[4]); kputs("\t*", &t[4]);
			fprintf(stderr, "H @SQ\tSN:%s\tLN:%ld\n", t[0].s, t[1].l);
		} else if (strcmp(s.s, "BQ") == 0) {
			if (dret != '\n') ks_getuntil(ks, '\n', &s, &dret);
			t[4].s[--t[4].l] = 0;
			for (i = 0; i < t[2].l; ++i) {
				int q;
				if (ks_getuntil(ks, 0, &s, &dret) < 0) fprintf(stderr, "E truncated contig quality\n");
				q = atoi(s.s) + 33;
				if (q > 126) q = 126;
				kputc(q, &t[4]);
			}
			if (dret != '\n') ks_getuntil(ks, '\n', &s, &dret);
			ks_getuntil(ks, '\n', &s, &dret);
			puts(t[4].s); t[4].l = 0;
		} else if (strcmp(s.s, "AF") == 0) { // read coordinate
			int reversed, neg, pos;
			if (t[4].l) puts(t[4].s);
			t[4].l = 0;
			ks_getuntil(ks, 0, &s, &dret);
			ks_getuntil(ks, 0, &s, &dret); reversed = s.s[0] == 'C'? 1 : 0;
			ks_getuntil(ks, 0, &s, &dret); pos = atoi(s.s); neg = pos < 0? 1 : 0; pos = pos < 0? -pos : pos;
			if (af_n == af_max) {
				af_max = af_max? af_max<<1 : 4;
				af = realloc(af, af_max * sizeof(unsigned));
			}
			af[af_n++] = pos << 2 | neg << 1 | reversed;
		} else if (strcmp(s.s, "RD") == 0) { // read sequence
			t[2].l = t[3].l = t[4].l = 0;
			ks_getuntil(ks, 0, &t[4], &dret); kputc('\t', &t[4]); // read name
			kputw((af[af_i]&1)? 16 : 0, &t[4]); kputc('\t', &t[4]);
			kputsn(t[0].s, t[0].l, &t[4]); kputc('\t', &t[4]); // the SAM line stops at RNAME
			if (dret != '\n') ks_getuntil(ks, '\n', &s, &dret);
			while (ks_getuntil(ks, '\n', &s, &dret) >= 0 && s.l > 0) kputs(s.s, &t[2]);
		} else if (strcmp(s.s, "QA") == 0) { // clipping
			int beg, end, pos;
			ks_getuntil(ks, 0, &s, &dret); ks_getuntil(ks, 0, &s, &dret); // skip quality clipping
			ks_getuntil(ks, 0, &s, &dret); beg = atoi(s.s) - 1; // align clipping start
			ks_getuntil(ks, 0, &s, &dret); end = atoi(s.s);
			pos = af[af_i]>>2;
			if (af[af_i]>>1&1) pos = -pos;
			pos += beg;
			kputw(pos, &t[4]); kputs("\t60\t", &t[4]);
			n_cigar = 0; // start to generate CIGAR
			if (beg) write_cigar(cigar, n_cigar, m_cigar, beg<<4|4);
			__padded2cigar(t[2], t[3]);
			if (end < t[2].l) write_cigar(cigar, n_cigar, m_cigar, (t[2].l-end)<<4|4);
			if (n_cigar > 1 && (cigar[0]&0xf) == 4) cigar[1] -= beg<<4;
			if (n_cigar > 1 && (cigar[n_cigar-1]&0xf) == 4) cigar[n_cigar-2] -= cigar[n_cigar-1]>>4<<4;
			for (i = 0; i < n_cigar; ++i) {
				kputw(cigar[i]>>4, &t[4]); kputc("MIDNSHP=X"[cigar[i]&0xf], &t[4]);
			}
			kputs("\t*\t0\t0\t", &t[4]);
			kputsn(t[3].s, t[3].l, &t[4]); kputs("\t*", &t[4]);
			puts(t[4].s);
			++af_i;
		} else if (dret != '\n') ks_getuntil(ks, '\n', &s, &dret);
	}
	ks_destroy(ks);
	gzclose(fp);
	free(af); free(s.s); free(cigar);
	for (i = 0; i < N_TMPSTR; ++i) free(t[i].s);
	return 0;
}
