/*  bam_tview.h -- tview subcommand.

    Copyright (C) 2008, 2013 Genome Research Ltd.
    Portions copyright (C) 2013 Pierre Lindenbaum, Institut du Thorax, INSERM U1087, Université de Nantes.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notices and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef BAM_TVIEW_H
#define BAM_TVIEW_H

#include <ctype.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdarg.h>
#include <htslib/sam.h>
#include "bam2bcf.h"
#include <htslib/khash.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>
#include "bam_lpileup.h"


KHASH_MAP_INIT_STR(kh_rg, const char *)

/* Holds state of Tview */
typedef struct AbstractTview {
    int mrow, mcol;

    hts_idx_t* idx;
    bam_lplbuf_t* lplbuf;
    bam_hdr_t* header;
    samFile* fp;
    int curr_tid, left_pos;
    faidx_t* fai;
    bcf_callaux_t* bca;

    int ccol, last_pos, row_shift, base_for, color_for, is_dot, l_ref, ins;
    int no_skip, show_name, inverse;
    char *ref;
    /* maps @RG ID => SM (sample), in practice only used to determine whether a particular RG is in the list of allowed ones */
    khash_t(kh_rg) *rg_hash;
    /* callbacks */
    void (*my_destroy)(struct AbstractTview* );
    void (*my_mvprintw)(struct AbstractTview* ,int,int,const char*,...);
    void (*my_mvaddch)(struct AbstractTview*,int,int,int);
    void (*my_attron)(struct AbstractTview*,int);
    void (*my_attroff)(struct AbstractTview*,int);
    void (*my_clear)(struct AbstractTview*);
    int (*my_colorpair)(struct AbstractTview*,int);
    int (*my_drawaln)(struct AbstractTview*,int,int);
    int (*my_loop)(struct AbstractTview*);
    int (*my_underline)(struct AbstractTview*);
} tview_t;


char bam_aux_getCEi(bam1_t *b, int i);
char bam_aux_getCSi(bam1_t *b, int i);
char bam_aux_getCQi(bam1_t *b, int i);

#define TV_MIN_ALNROW 2
#define TV_MAX_GOTO  40
#define TV_LOW_MAPQ  10

#define TV_COLOR_MAPQ   0
#define TV_COLOR_BASEQ  1
#define TV_COLOR_NUCL   2
#define TV_COLOR_COL    3
#define TV_COLOR_COLQ   4

#define TV_BASE_NUCL 0
#define TV_BASE_COLOR_SPACE 1

int tv_pl_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);
int base_tv_init(tview_t*,const char *fn, const char *fn_fa,
                 const char *samples, const htsFormat *fmt);
void base_tv_destroy(tview_t*);
int base_draw_aln(tview_t *tv, int tid, int pos);

typedef struct Tixel
    {
    int ch;
    int attributes;
    }tixel_t;

#endif

