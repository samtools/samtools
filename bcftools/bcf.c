#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include "kstring.h"
#include "bcf.h"

bcf_t *bcf_open(const char *fn, const char *mode)
{
	bcf_t *b;
	b = calloc(1, sizeof(bcf_t));
	if (strchr(mode, 'w')) {
		b->fp = strcmp(fn, "-")? bgzf_open(fn, mode) : bgzf_fdopen(fileno(stdout), mode);
	} else {
		b->fp = strcmp(fn, "-")? bgzf_open(fn, mode) : bgzf_fdopen(fileno(stdin), mode);
	}
	return b;
}

int bcf_close(bcf_t *b)
{
	int ret;
	if (b == 0) return 0;
	ret = bgzf_close(b->fp);
	free(b);
	return ret;
}

int bcf_hdr_write(bcf_t *b, const bcf_hdr_t *h)
{
	if (b == 0 || h == 0) return -1;
	bgzf_write(b->fp, "BCF\4", 4);
	bgzf_write(b->fp, &h->l_nm, 4);
	bgzf_write(b->fp, h->name, h->l_nm);
	bgzf_write(b->fp, &h->l_smpl, 4);
	bgzf_write(b->fp, h->sname, h->l_smpl);
	bgzf_write(b->fp, &h->l_txt, 4);
	bgzf_write(b->fp, h->txt, h->l_txt);
	bgzf_flush(b->fp);
	return 16 + h->l_nm + h->l_smpl + h->l_txt;
}

bcf_hdr_t *bcf_hdr_read(bcf_t *b)
{
	uint8_t magic[4];
	bcf_hdr_t *h;
	if (b == 0) return 0;
	h = calloc(1, sizeof(bcf_hdr_t));
	bgzf_read(b->fp, magic, 4);
	bgzf_read(b->fp, &h->l_nm, 4);
	h->name = malloc(h->l_nm);
	bgzf_read(b->fp, h->name, h->l_nm);
	bgzf_read(b->fp, &h->l_smpl, 4);
	h->sname = malloc(h->l_smpl);
	bgzf_read(b->fp, h->sname, h->l_smpl);
	bgzf_read(b->fp, &h->l_txt, 4);
	h->txt = malloc(h->l_txt);
	bgzf_read(b->fp, h->txt, h->l_txt);
	bcf_hdr_sync(h);
	return h;
}

void bcf_hdr_destroy(bcf_hdr_t *h)
{
	if (h == 0) return;
	free(h->name); free(h->sname); free(h->txt); free(h->ns); free(h->sns);
	free(h);
}

static inline char **cnt_null(int l, char *str, int *_n)
{
	int n = 0;
	char *p, **list;
	*_n = 0;
	if (l == 0 || str == 0) return 0;
	for (p = str; p != str + l; ++p)
		if (*p == 0) ++n;
	*_n = n;
	list = calloc(n, sizeof(void*));
	list[0] = str;
	for (p = str, n = 1; p < str + l - 1; ++p)
		if (*p == 0) list[n++] = p + 1;
	return list;
}

int bcf_hdr_sync(bcf_hdr_t *b)
{
	if (b == 0) return -1;
	if (b->ns) free(b->ns);
	if (b->sns) free(b->sns);
	if (b->l_nm) b->ns = cnt_null(b->l_nm, b->name, &b->n_ref);
	else b->ns = 0, b->n_ref = 0;
	b->sns = cnt_null(b->l_smpl, b->sname, &b->n_smpl);
	return 0;
}

int bcf_sync(bcf1_t *b)
{
	char *p, *tmp[5];
	int i, n, n_smpl = b->n_smpl;
	ks_tokaux_t aux;
	// set ref, alt, flt, info, fmt
	b->ref = b->alt = b->flt = b->info = b->fmt = 0;
	for (p = b->str, n = 0; p < b->str + b->l_str; ++p) {
		if (*p == 0 && p+1 != b->str + b->l_str) {
			if (n == 5) {
				++n;
				break;
			} else tmp[n++] = p + 1;
		}
	}
	if (n != 5) {
		fprintf(stderr, "[%s] incorrect number of fields (%d != 5) at %d:%d\n", __func__, n, b->tid, b->pos);
		return -1;
	}
	b->ref = tmp[0]; b->alt = tmp[1]; b->flt = tmp[2]; b->info = tmp[3]; b->fmt = tmp[4];
	// set n_alleles
	if (*b->alt == 0) b->n_alleles = 1;
	else {
		for (p = b->alt, n = 1; *p; ++p)
			if (*p == ',') ++n;
		b->n_alleles = n + 1;
	}
	// set n_gi and gi[i].fmt
	for (p = b->fmt, n = 1; *p; ++p)
		if (*p == ':') ++n;
	if (n > b->m_gi) {
		int old_m = b->m_gi;
		b->m_gi = n;
		kroundup32(b->m_gi);
		b->gi = realloc(b->gi, b->m_gi * sizeof(bcf_ginfo_t));
		memset(b->gi + old_m, 0, (b->m_gi - old_m) * sizeof(bcf_ginfo_t));
	}
	b->n_gi = n;
	for (p = kstrtok(b->fmt, ":", &aux), n = 0; p; p = kstrtok(0, 0, &aux))
		b->gi[n++].fmt = bcf_str2int(p, aux.p - p);
	// set gi[i].len
	for (i = 0; i < b->n_gi; ++i) {
		if (b->gi[i].fmt == bcf_str2int("PL", 2)) {
			b->gi[i].len = b->n_alleles * (b->n_alleles + 1) / 2;
		} else if (b->gi[i].fmt == bcf_str2int("DP", 2) || b->gi[i].fmt == bcf_str2int("HQ", 2) || b->gi[i].fmt == bcf_str2int("DV", 2)) {
			b->gi[i].len = 2;
		} else if (b->gi[i].fmt == bcf_str2int("GQ", 2) || b->gi[i].fmt == bcf_str2int("GT", 2)) {
			b->gi[i].len = 1;
		} else if (b->gi[i].fmt == bcf_str2int("SP", 2)) {
			b->gi[i].len = 4;
		} else if (b->gi[i].fmt == bcf_str2int("GL", 2)) {
			b->gi[i].len = b->n_alleles * (b->n_alleles + 1) / 2 * 4;
		}
		b->gi[i].data = realloc(b->gi[i].data, n_smpl * b->gi[i].len);
	}
	return 0;
}

int bcf_write(bcf_t *bp, const bcf_hdr_t *h, const bcf1_t *b)
{
	int i, l = 0;
	if (b == 0) return -1;
	bgzf_write(bp->fp, &b->tid, 4);
	bgzf_write(bp->fp, &b->pos, 4);
	bgzf_write(bp->fp, &b->qual, 4);
	bgzf_write(bp->fp, &b->l_str, 4);
	bgzf_write(bp->fp, b->str, b->l_str);
	l = 12 + b->l_str;
	for (i = 0; i < b->n_gi; ++i) {
		bgzf_write(bp->fp, b->gi[i].data, b->gi[i].len * h->n_smpl);
		l += b->gi[i].len * h->n_smpl;
	}
	return l;
}

int bcf_read(bcf_t *bp, const bcf_hdr_t *h, bcf1_t *b)
{
	int i, l = 0;
	if (b == 0) return -1;
	if (bgzf_read(bp->fp, &b->tid, 4) == 0) return -1;
	b->n_smpl = h->n_smpl;
	bgzf_read(bp->fp, &b->pos, 4);
	bgzf_read(bp->fp, &b->qual, 4);
	bgzf_read(bp->fp, &b->l_str, 4);
	if (b->l_str > b->m_str) {
		b->m_str = b->l_str;
		kroundup32(b->m_str);
		b->str = realloc(b->str, b->m_str);
	}
	bgzf_read(bp->fp, b->str, b->l_str);
	l = 12 + b->l_str;
	if (bcf_sync(b) < 0) return -2;
	for (i = 0; i < b->n_gi; ++i) {
		bgzf_read(bp->fp, b->gi[i].data, b->gi[i].len * h->n_smpl);
		l += b->gi[i].len * h->n_smpl;
	}
	return l;
}

int bcf_destroy(bcf1_t *b)
{
	int i;
	if (b == 0) return -1;
	free(b->str);
	for (i = 0; i < b->m_gi; ++i)
		free(b->gi[i].data);
	free(b->gi);
	free(b);
	return 0;
}

static inline void fmt_str(const char *p, kstring_t *s)
{
	if (*p == 0) kputc('.', s);
	else kputs(p, s);
}

void bcf_fmt_core(const bcf_hdr_t *h, bcf1_t *b, kstring_t *s)
{
	int i, j, x;
	s->l = 0;
	if (h->n_ref) kputs(h->ns[b->tid], s);
	else kputw(b->tid, s);
	kputc('\t', s);
	kputw(b->pos + 1, s); kputc('\t', s);
	fmt_str(b->str, s); kputc('\t', s);
	fmt_str(b->ref, s); kputc('\t', s);
	fmt_str(b->alt, s); kputc('\t', s);
	ksprintf(s, "%.3g", b->qual); kputc('\t', s);
	fmt_str(b->flt, s); kputc('\t', s);
	fmt_str(b->info, s);
	if (b->fmt[0]) {
		kputc('\t', s);
		fmt_str(b->fmt, s);
	}
	x = b->n_alleles * (b->n_alleles + 1) / 2;
	if (b->n_gi == 0) return;
    int iPL = -1;
    if ( b->n_alleles > 2 ) {
        for (i=0; i<b->n_gi; i++) {
            if ( b->gi[i].fmt == bcf_str2int("PL", 2) ) {
                iPL = i;
                break;
            }
        }
    }
	for (j = 0; j < h->n_smpl; ++j) {
        int ploidy = b->ploidy ? b->ploidy[j] : 2;
		kputc('\t', s);
		for (i = 0; i < b->n_gi; ++i) {
			if (i) kputc(':', s);
			if (b->gi[i].fmt == bcf_str2int("PL", 2)) {
				uint8_t *d = (uint8_t*)b->gi[i].data + j * x;
				int k;
                if ( ploidy==1 )
                    for (k=0; k<b->n_alleles; k++)
                    {
                        if (k>0) kputc(',', s);
                        kputw(d[(k+1)*(k+2)/2-1], s);
                    }
                else
                    for (k = 0; k < x; ++k) {
                        if (k > 0) kputc(',', s);
                        kputw(d[k], s);
                    }
			} else if (b->gi[i].fmt == bcf_str2int("DP", 2) || b->gi[i].fmt == bcf_str2int("DV", 2)) {
				kputw(((uint16_t*)b->gi[i].data)[j], s);
			} else if (b->gi[i].fmt == bcf_str2int("GQ", 2)) {
				kputw(((uint8_t*)b->gi[i].data)[j], s);
			} else if (b->gi[i].fmt == bcf_str2int("SP", 2)) {
				kputw(((int32_t*)b->gi[i].data)[j], s);
			} else if (b->gi[i].fmt == bcf_str2int("GT", 2)) {
                int y = ((uint8_t*)b->gi[i].data)[j];
                if ( ploidy==1 )
                {
                    if ( y>>7&1 )
                        kputc('.', s);
                    else 
                        kputc('0' + (y>>3&7), s);
                }
                else
                {
                    if ( y>>7&1 )
                        kputsn("./.", 3, s);
                    else { 
                        kputc('0' + (y>>3&7), s);
                        kputc("/|"[y>>6&1], s);
                        kputc('0' + (y&7), s);
                    }
                }
			} else if (b->gi[i].fmt == bcf_str2int("GL", 2)) {
				float *d = (float*)b->gi[i].data + j * x;
				int k;
				//printf("- %lx\n", d);
				for (k = 0; k < x; ++k) {
					if (k > 0) kputc(',', s);
					ksprintf(s, "%.2f", d[k]);
				}
			} else kputc('.', s); // custom fields
		}
	}
}

char *bcf_fmt(const bcf_hdr_t *h, bcf1_t *b)
{
	kstring_t s;
	s.l = s.m = 0; s.s = 0;
	bcf_fmt_core(h, b, &s);
	return s.s;
}

int bcf_append_info(bcf1_t *b, const char *info, int l)
{
	int shift = b->fmt - b->str;
	int l_fmt = b->l_str - shift;
	char *ori = b->str;
	if (b->l_str + l > b->m_str) { // enlarge if necessary
		b->m_str = b->l_str + l;
		kroundup32(b->m_str);
		b->str = realloc(b->str, b->m_str);
	}
	memmove(b->str + shift + l, b->str + shift, l_fmt); // move the FORMAT field
	memcpy(b->str + shift - 1, info, l); // append to the INFO field
	b->str[shift + l - 1] = '\0';
	b->fmt = b->str + shift + l;
	b->l_str += l;
	if (ori != b->str) bcf_sync(b); // synchronize when realloc changes the pointer
	return 0;
}

int remove_tag(char *str, const char *tag, char delim)
{
    char *tmp = str, *p;
    int len_diff = 0, ori_len = strlen(str);
    while ( *tmp && (p = strstr(tmp,tag)) )
    {
        if ( p>str )
        {
            if ( *(p-1)!=delim ) { tmp=p+1; continue; } // shared substring
            p--;
        }
        char *q=p+1;
        while ( *q && *q!=delim ) q++;
        if ( p==str && *q ) q++;        // the tag is first, don't move the delim char
        len_diff += q-p;
        if ( ! *q ) { *p = 0; break; }  // the tag was last, no delim follows
        else
            memmove(p,q,ori_len-(int)(p-str)-(int)(q-p));  // *q==delim
    }
    if ( len_diff==ori_len )
        str[0]='.', str[1]=0, len_diff--;

    return len_diff;
}


void rm_info(kstring_t *s, const char *key)
{
    char *p = s->s; 
    int n = 0;
    while ( n<4 )
    {
        if ( !*p ) n++;
        p++;
    }
    char *q = p+1; 
    while ( *q && q-s->s<s->l ) q++;

    int nrm = remove_tag(p, key, ';');
    if ( nrm )
        memmove(q-nrm, q, s->s+s->l-q+1);
    s->l -= nrm;
}

int bcf_cpy(bcf1_t *r, const bcf1_t *b)
{
	char *t1 = r->str;
	bcf_ginfo_t *t2 = r->gi;
	int i, t3 = r->m_str, t4 = r->m_gi;
	*r = *b;
	r->str = t1; r->gi = t2; r->m_str = t3; r->m_gi = t4;
	if (r->m_str < b->m_str) {
		r->m_str = b->m_str;
		r->str = realloc(r->str, r->m_str);
	}
	memcpy(r->str, b->str, r->m_str);
	bcf_sync(r); // calling bcf_sync() is simple but inefficient
	for (i = 0; i < r->n_gi; ++i)
		memcpy(r->gi[i].data, b->gi[i].data, r->n_smpl * r->gi[i].len);
	return 0;
}

int bcf_is_indel(const bcf1_t *b)
{
	char *p;
	if (strlen(b->ref) > 1) return 1;
	for (p = b->alt; *p; ++p)
		if (*p != ',' && p[1] != ',' && p[1] != '\0')
			return 1;
	return 0;
}
