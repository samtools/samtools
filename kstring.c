#include <stdarg.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include "kstring.h"

int ksprintf(kstring_t *s, const char *fmt, ...)
{
	va_list ap;
	int l;
	va_start(ap, fmt);
	l = vsnprintf(s->s + s->l, s->m - s->l, fmt, ap); // This line does not work with glibc 2.0. See `man snprintf'.
	va_end(ap);
	if (l + 1 > s->m - s->l) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
		va_start(ap, fmt);
		l = vsnprintf(s->s + s->l, s->m - s->l, fmt, ap);
	}
	va_end(ap);
	s->l += l;
	return l;
}

// s MUST BE a null terminated string; l = strlen(s)
int ksplit_core(char *s, int delimiter, int *_max, int **_offsets)
{
	int i, n, max, last_char, last_start, *offsets, l;
	n = 0; max = *_max; offsets = *_offsets;
	l = strlen(s);
	
#define __ksplit_aux do {												\
		if (_offsets) {													\
			s[i] = 0;													\
			if (n == max) {												\
				max = max? max<<1 : 2;									\
				offsets = (int*)realloc(offsets, sizeof(int) * max);	\
			}															\
			offsets[n++] = last_start;									\
		} else ++n;														\
	} while (0)

	for (i = 0, last_char = last_start = 0; i <= l; ++i) {
		if (delimiter == 0) {
			if (isspace(s[i]) || s[i] == 0) {
				if (isgraph(last_char)) __ksplit_aux; // the end of a field
			} else {
				if (isspace(last_char) || last_char == 0) last_start = i;
			}
		} else {
			if (s[i] == delimiter || s[i] == 0) {
				if (last_char != 0 && last_char != delimiter) __ksplit_aux; // the end of a field
			} else {
				if (last_char == delimiter || last_char == 0) last_start = i;
			}
		}
		last_char = s[i];
	}
	*_max = max; *_offsets = offsets;
	return n;
}

/**********************
 * Boyer-Moore search *
 **********************/

// reference: http://www-igm.univ-mlv.fr/~lecroq/string/node14.html
int *ksBM_prep(const uint8_t *pat, int m)
{
	int i, *suff, *prep, *bmGs, *bmBc;
	prep = calloc(m + 256, 1);
	bmGs = prep; bmBc = prep + m;
	{ // preBmBc()
		for (i = 0; i < 256; ++i) bmBc[i] = m;
		for (i = 0; i < m - 1; ++i) bmBc[pat[i]] = m - i - 1;
	}
	suff = calloc(m, sizeof(int));
	{ // suffixes()
		int f = 0, g;
		suff[m - 1] = m;
		g = m - 1;
		for (i = m - 2; i >= 0; --i) {
			if (i > g && suff[i + m - 1 - f] < i - g)
				suff[i] = suff[i + m - 1 - f];
			else {
				if (i < g) g = i;
				f = i;
				while (g >= 0 && pat[g] == pat[g + m - 1 - f]) --g;
				suff[i] = f - g;
			}
		}
	}
	{ // preBmGs()
		int j = 0;
		for (i = 0; i < m; ++i) bmGs[i] = m;
		for (i = m - 1; i >= 0; --i)
			if (suff[i] == i + 1)
				for (; j < m - 1 - i; ++j)
					if (bmGs[j] == m)
						bmGs[j] = m - 1 - i;
		for (i = 0; i <= m - 2; ++i)
			bmGs[m - 1 - suff[i]] = m - 1 - i;
	}
	free(suff);
	return prep;
}

int *ksBM_search(const uint8_t *str, int n, const uint8_t *pat, int m, int *_prep, int *n_matches)
{
	int i, j, *prep, *bmGs, *bmBc;
	int *matches = 0, mm = 0, nm = 0;
	prep = _prep? _prep : ksBM_prep(pat, m);
	bmGs = prep; bmBc = prep + m;
	j = 0;
	while (j <= n - m) {
		for (i = m - 1; i >= 0 && pat[i] == str[i+j]; --i);
		if (i < 0) {
			if (nm == mm) {
				mm = mm? mm<<1 : 1;
				matches = realloc(matches, mm * sizeof(int));
			}
			matches[nm++] = j;
			j += bmGs[0];
		} else {
			int max = bmBc[str[i+j]] - m + 1 + i;
			if (max < bmGs[i]) max = bmGs[i];
			j += max;
		}
	}
	*n_matches = nm;
	if (_prep == 0) free(prep);
	return matches;
}

#ifdef KSTRING_MAIN
#include <stdio.h>
int main()
{
	kstring_t *s;
	int *fields, n, i;
	s = (kstring_t*)calloc(1, sizeof(kstring_t));
	// test ksprintf()
	ksprintf(s, " abcdefg:    %d ", 100);
	printf("'%s'\n", s->s);
	// test ksplit()
	fields = ksplit(s, 0, &n);
	for (i = 0; i < n; ++i)
		printf("field[%d] = '%s'\n", i, s->s + fields[i]);
	free(s);

	{
		static char *str = "abcdefgcdg";
		static char *pat = "cd";
		int n, *matches;
		matches = ksBM_search(str, strlen(str), pat, strlen(pat), 0, &n);
		printf("%d: \n", n);
		for (i = 0; i < n; ++i)
			printf("- %d\n", matches[i]);
		free(matches);
	}
	return 0;
}
#endif
