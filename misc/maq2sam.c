/*  maq2sam.c -- convert MAQ output to SAM.

    Copyright (C) 2008, 2009 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <string.h>
#include <zlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <assert.h>

#define PACKAGE_VERSION "r439"

//#define MAQ_LONGREADS

#ifdef MAQ_LONGREADS
#  define MAX_READLEN 128
#else
#  define MAX_READLEN 64
#endif

#define MAX_NAMELEN 36
#define MAQMAP_FORMAT_OLD 0
#define MAQMAP_FORMAT_NEW -1

#define PAIRFLAG_FF      0x01
#define PAIRFLAG_FR      0x02
#define PAIRFLAG_RF      0x04
#define PAIRFLAG_RR      0x08
#define PAIRFLAG_PAIRED  0x10
#define PAIRFLAG_DIFFCHR 0x20
#define PAIRFLAG_NOMATCH 0x40
#define PAIRFLAG_SW      0x80

typedef struct
{
    uint8_t seq[MAX_READLEN]; /* the last base is the single-end mapping quality. */
    uint8_t size, map_qual, info1, info2, c[2], flag, alt_qual;
    uint32_t seqid, pos;
    int dist;
    char name[MAX_NAMELEN];
} maqmap1_t;

typedef struct
{
    int format, n_ref;
    char **ref_name;
    uint64_t n_mapped_reads;
    maqmap1_t *mapped_reads;
} maqmap_t;

maqmap_t *maq_new_maqmap(void)
{
    maqmap_t *mm = (maqmap_t*)calloc(1, sizeof(maqmap_t));
    mm->format = MAQMAP_FORMAT_NEW;
    return mm;
}
void maq_delete_maqmap(maqmap_t *mm)
{
    int i;
    if (mm == 0) return;
    for (i = 0; i < mm->n_ref; ++i)
        free(mm->ref_name[i]);
    free(mm->ref_name);
    free(mm->mapped_reads);
    free(mm);
}
maqmap_t *maqmap_read_header(gzFile fp)
{
    maqmap_t *mm;
    int k, len;
    mm = maq_new_maqmap();
    gzread(fp, &mm->format, sizeof(int));
    if (mm->format != MAQMAP_FORMAT_NEW) {
        if (mm->format > 0) {
            fprintf(stderr, "** Obsolete map format is detected. Please use 'mapass2maq' command to convert the format.\n");
            exit(3);
        }
        assert(mm->format == MAQMAP_FORMAT_NEW);
    }
    gzread(fp, &mm->n_ref, sizeof(int));
    mm->ref_name = (char**)calloc(mm->n_ref, sizeof(char*));
    for (k = 0; k != mm->n_ref; ++k) {
        gzread(fp, &len, sizeof(int));
        mm->ref_name[k] = (char*)malloc(len * sizeof(char));
        gzread(fp, mm->ref_name[k], len);
    }
    /* read number of mapped reads */
    gzread(fp, &mm->n_mapped_reads, sizeof(uint64_t));
    return mm;
}

void maq2tam_core(gzFile fp, const char *rg)
{
    maqmap_t *mm;
    maqmap1_t mm1, *m1;
    int ret;
    m1 = &mm1;
    mm = maqmap_read_header(fp);
    while ((ret = gzread(fp, m1, sizeof(maqmap1_t))) == sizeof(maqmap1_t)) {
        int j, flag = 0, se_mapq = m1->seq[MAX_READLEN-1];
        if (m1->flag) flag |= 1;
        if ((m1->flag&PAIRFLAG_PAIRED) || ((m1->flag&PAIRFLAG_SW) && m1->flag != 192)) flag |= 2;
        if (m1->flag == 192) flag |= 4;
        if (m1->flag == 64) flag |= 8;
        if (m1->pos&1) flag |= 0x10;
        if ((flag&1) && m1->dist != 0) {
            int c;
            if (m1->dist > 0) {
                if (m1->flag&(PAIRFLAG_FF|PAIRFLAG_RF)) c = 0;
                else if (m1->flag&(PAIRFLAG_FR|PAIRFLAG_RR)) c = 1;
                else c = m1->pos&1;
            } else {
                if (m1->flag&(PAIRFLAG_FF|PAIRFLAG_FR)) c = 0;
                else if (m1->flag&(PAIRFLAG_RF|PAIRFLAG_RR)) c = 1;
                else c = m1->pos&1;
            }
            if (c) flag |= 0x20;
        }
        if (m1->flag) {
            int l = strlen(m1->name);
            if (m1->name[l-2] == '/') {
                flag |= (m1->name[l-1] == '1')? 0x40 : 0x80;
                m1->name[l-2] = '\0';
            }
        }
        printf("%s\t%d\t", m1->name, flag);
        printf("%s\t%d\t", mm->ref_name[m1->seqid], (m1->pos>>1)+1);
        if (m1->flag == 130) {
            int c = (int8_t)m1->seq[MAX_READLEN-1];
            printf("%d\t", m1->alt_qual);
            if (c == 0) printf("%dM\t", m1->size);
            else {
                if (c > 0) printf("%dM%dI%dM\t", m1->map_qual, c, m1->size - m1->map_qual - c);
                else printf("%dM%dD%dM\t", m1->map_qual, -c, m1->size - m1->map_qual);
            }
            se_mapq = 0; // zero SE mapQ for reads aligned by SW
        } else {
            if (flag&4) printf("0\t*\t");
            else printf("%d\t%dM\t", m1->map_qual, m1->size);
        }
        printf("*\t0\t%d\t", m1->dist);
        for (j = 0; j != m1->size; ++j) {
            if (m1->seq[j] == 0) putchar('N');
            else putchar("ACGT"[m1->seq[j]>>6&3]);
        }
        putchar('\t');
        for (j = 0; j != m1->size; ++j)
            putchar((m1->seq[j]&0x3f) + 33);
        putchar('\t');
        if (rg) printf("RG:Z:%s\t", rg);
        if (flag&4) { // unmapped
            printf("MF:i:%d\n", m1->flag);
        } else {
            printf("MF:i:%d\t", m1->flag);
            if (m1->flag) printf("AM:i:%d\tSM:i:%d\t", m1->alt_qual, se_mapq);
            printf("NM:i:%d\tUQ:i:%d\tH0:i:%d\tH1:i:%d\n", m1->info1&0xf, m1->info2, m1->c[0], m1->c[1]);
        }
    }
    if (ret > 0)
        fprintf(stderr, "Truncated! Continue anyway.\n");
    maq_delete_maqmap(mm);
}

int main(int argc, char *argv[])
{
    gzFile fp;
    if (argc == 1) {
        fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
        fprintf(stderr, "Usage: maq2sam <in.map> [<readGroup>]\n");
        return 1;
    }
    fp = strcmp(argv[1], "-")? gzopen(argv[1], "r") : gzdopen(fileno(stdin), "r");
    maq2tam_core(fp, argc > 2? argv[2] : 0);
    gzclose(fp);
    return 0;
}
