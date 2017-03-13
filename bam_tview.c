/*  bam_tview.c -- tview subcommand.

    Copyright (C) 2008-2016 Genome Research Ltd.
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

#include <config.h>

#include <regex.h>
#include <assert.h>
#include "bam_tview.h"
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include "samtools.h"
#include "sam_opts.h"

khash_t(kh_rg)* get_rg_sample(const char* header, const char* sample)
{
    khash_t(kh_rg)* rg_hash = kh_init(kh_rg);
    // given sample id return all the RD ID's
    const char rg_regex[] = "^@RG.*\tID:([!-)+-<>-~][ !-~]*)(\t.*$|$)";

    regex_t rg_id;
    regmatch_t* matches = (regmatch_t*)calloc(2, sizeof(regmatch_t));
    if (matches == NULL) { perror("out of memory"); exit(-1); }
    regcomp(&rg_id, rg_regex, REG_EXTENDED|REG_NEWLINE);
    char* text = strdup(header);
    char* end = text + strlen(header);
    char* tofree = text;
    while (end > text && regexec(&rg_id, text, 2, matches, 0) == 0) { //    foreach rg id in  header
        int ret;
        text[matches[1].rm_eo] = '\0';
        kh_put(kh_rg, rg_hash, strdup(text+matches[1].rm_so), &ret); // Add the RG to the list
        text += matches[0].rm_eo + 1; // Move search pointer forward
    }
    free(tofree);
    return rg_hash;
}

int base_tv_init(tview_t* tv, const char *fn, const char *fn_fa,
                 const char *samples, const htsFormat *fmt)
{
    assert(tv!=NULL);
    assert(fn!=NULL);
    tv->mrow = 24; tv->mcol = 80;
    tv->color_for = TV_COLOR_MAPQ;
    tv->is_dot = 1;

    tv->fp = sam_open_format(fn, "r", fmt);
    if(tv->fp == NULL)
    {
        print_error_errno("tview", "can't open \"%s\"", fn);
        exit(EXIT_FAILURE);
    }
    // TODO bgzf_set_cache_size(tv->fp->fp.bgzf, 8 * 1024 *1024);
    assert(tv->fp);

    tv->header = sam_hdr_read(tv->fp);
    if(tv->header == NULL)
    {
        print_error("tview", "cannot read \"%s\"", fn);
        exit(EXIT_FAILURE);
    }
    tv->idx = sam_index_load(tv->fp, fn);
    if (tv->idx == NULL)
    {
        print_error("tview", "cannot read index for \"%s\"", fn);
        exit(EXIT_FAILURE);
    }
    tv->lplbuf = bam_lplbuf_init(tv_pl_func, tv);
    if (fn_fa) tv->fai = fai_load(fn_fa);
    tv->bca = bcf_call_init(0.83, 13);
    tv->ins = 1;

    // If the user has asked for specific samples find out create a list of readgroups make up these samples
    if ( samples )
    {
        tv->rg_hash = get_rg_sample(tv->header->text, samples); // Init the list of rg's
    }

    return 0;
}


void base_tv_destroy(tview_t* tv)
{
    bam_lplbuf_destroy(tv->lplbuf);
    bcf_call_destroy(tv->bca);
    hts_idx_destroy(tv->idx);
    if (tv->fai) fai_destroy(tv->fai);
    free(tv->ref);
    bam_hdr_destroy(tv->header);
    sam_close(tv->fp);
}


int tv_pl_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
    tview_t *tv = (tview_t*)data;
    int i, j, c, rb, attr, max_ins = 0;
    uint32_t call = 0;
    if (pos < tv->left_pos || tv->ccol > tv->mcol) return 0; // out of screen
    // print reference
    rb = (tv->ref && pos - tv->left_pos < tv->l_ref)? tv->ref[pos - tv->left_pos] : 'N';
    for (i = tv->last_pos + 1; i < pos; ++i) {
        if (i%10 == 0 && tv->mcol - tv->ccol >= 10) tv->my_mvprintw(tv,0, tv->ccol, "%-d", i+1);
        c = tv->ref? tv->ref[i - tv->left_pos] : 'N';
        tv->my_mvaddch(tv,1, tv->ccol++, c);
    }
    if (pos%10 == 0 && tv->mcol - tv->ccol >= 10) tv->my_mvprintw(tv,0, tv->ccol, "%-d", pos+1);
    { // call consensus
        bcf_callret1_t bcr;
        memset(&bcr, 0, sizeof bcr);
        int qsum[4], a1, a2, tmp;
        double p[3], prior = 30;
        bcf_call_glfgen(n, pl, seq_nt16_table[rb], tv->bca, &bcr);
        for (i = 0; i < 4; ++i) qsum[i] = ((int)bcr.qsum[i])<<2 | i;
        for (i = 1; i < 4; ++i) // insertion sort
            for (j = i; j > 0 && qsum[j] > qsum[j-1]; --j)
                tmp = qsum[j], qsum[j] = qsum[j-1], qsum[j-1] = tmp;
        a1 = qsum[0]&3; a2 = qsum[1]&3;
        p[0] = bcr.p[a1*5+a1]; p[1] = bcr.p[a1*5+a2] + prior; p[2] = bcr.p[a2*5+a2];
        if ("ACGT"[a1] != toupper(rb)) p[0] += prior + 3;
        if ("ACGT"[a2] != toupper(rb)) p[2] += prior + 3;
        if (p[0] < p[1] && p[0] < p[2]) call = (1<<a1)<<16 | (int)((p[1]<p[2]?p[1]:p[2]) - p[0] + .499);
        else if (p[2] < p[1] && p[2] < p[0]) call = (1<<a2)<<16 | (int)((p[0]<p[1]?p[0]:p[1]) - p[2] + .499);
        else call = (1<<a1|1<<a2)<<16 | (int)((p[0]<p[2]?p[0]:p[2]) - p[1] + .499);
    }
    attr = tv->my_underline(tv);
    c = ",ACMGRSVTWYHKDBN"[call>>16&0xf];
    i = (call&0xffff)/10+1;
    if (i > 4) i = 4;
    attr |= tv->my_colorpair(tv,i);
    if (c == toupper(rb)) c = '.';
    tv->my_attron(tv,attr);
    tv->my_mvaddch(tv,2, tv->ccol, c);
    tv->my_attroff(tv,attr);
    if(tv->ins) {
        // calculate maximum insert
        for (i = 0; i < n; ++i) {
            const bam_pileup1_t *p = pl + i;
            if (p->indel > 0 && max_ins < p->indel) max_ins = p->indel;
        }
    }
    // core loop
    for (j = 0; j <= max_ins; ++j) {
        for (i = 0; i < n; ++i) {
            const bam_pileup1_t *p = pl + i;
            int row = TV_MIN_ALNROW + p->level - tv->row_shift;
            if (j == 0) {
                if (!p->is_del) {
                    if (tv->base_for == TV_BASE_COLOR_SPACE &&
                            (c = bam_aux_getCSi(p->b, p->qpos))) {
                        // assume that if we found one color, we will be able to get the color error
                        if (tv->is_dot && '-' == bam_aux_getCEi(p->b, p->qpos)) c = bam_is_rev(p->b)? ',' : '.';
                    } else {
                        if (tv->show_name) {
                            char *name = bam_get_qname(p->b);
                            c = (p->qpos + 1 >= p->b->core.l_qname)? ' ' : name[p->qpos];
                        } else {
                            c = seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)];
                            if (tv->is_dot && toupper(c) == toupper(rb)) c = bam_is_rev(p->b)? ',' : '.';
                        }
                    }
                } else c = p->is_refskip? (bam_is_rev(p->b)? '<' : '>') : '*';
            } else { // padding
                if (j > p->indel) c = '*';
                else { // insertion
                    if (tv->base_for ==  TV_BASE_NUCL) {
                        if (tv->show_name) {
                            char *name = bam_get_qname(p->b);
                            c = (p->qpos + j + 1 >= p->b->core.l_qname)? ' ' : name[p->qpos + j];
                        } else {
                            c = seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos + j)];
                            if (j == 0 && tv->is_dot && toupper(c) == toupper(rb)) c = bam_is_rev(p->b)? ',' : '.';
                        }
                    } else {
                        c = bam_aux_getCSi(p->b, p->qpos + j);
                        if (tv->is_dot && '-' == bam_aux_getCEi(p->b, p->qpos + j)) c = bam_is_rev(p->b)? ',' : '.';
                    }
                }
            }
            if (row > TV_MIN_ALNROW && row < tv->mrow) {
                int x;
                attr = 0;
                if (((p->b->core.flag&BAM_FPAIRED) && !(p->b->core.flag&BAM_FPROPER_PAIR))
                        || (p->b->core.flag & BAM_FSECONDARY)) attr |= tv->my_underline(tv);
                if (tv->color_for == TV_COLOR_BASEQ) {
                    x = bam_get_qual(p->b)[p->qpos]/10 + 1;
                    if (x > 4) x = 4;
                    attr |= tv->my_colorpair(tv,x);
                } else if (tv->color_for == TV_COLOR_MAPQ) {
                    x = p->b->core.qual/10 + 1;
                    if (x > 4) x = 4;
                    attr |= tv->my_colorpair(tv,x);
                } else if (tv->color_for == TV_COLOR_NUCL) {
                    x = seq_nt16_int[bam_seqi(bam_get_seq(p->b), p->qpos)] + 5;
                    attr |= tv->my_colorpair(tv,x);
                } else if(tv->color_for == TV_COLOR_COL) {
                    x = 0;
                    switch(bam_aux_getCSi(p->b, p->qpos)) {
                        case '0': x = 0; break;
                        case '1': x = 1; break;
                        case '2': x = 2; break;
                        case '3': x = 3; break;
                        case '4': x = 4; break;
                        default: x = seq_nt16_int[bam_seqi(bam_get_seq(p->b), p->qpos)]; break;
                    }
                    x+=5;
                    attr |= tv->my_colorpair(tv,x);
                } else if(tv->color_for == TV_COLOR_COLQ) {
                    x = bam_aux_getCQi(p->b, p->qpos);
                    if(0 == x) x = bam_get_qual(p->b)[p->qpos];
                    x = x/10 + 1;
                    if (x > 4) x = 4;
                    attr |= tv->my_colorpair(tv,x);
                }
                tv->my_attron(tv,attr);
                tv->my_mvaddch(tv,row, tv->ccol, bam_is_rev(p->b)? tolower(c) : toupper(c));
                tv->my_attroff(tv,attr);
            }
        }
        c = j? '*' : rb;
        if (c == '*') {
            attr = tv->my_colorpair(tv,8);
            tv->my_attron(tv,attr);
            tv->my_mvaddch(tv,1, tv->ccol++, c);
            tv->my_attroff(tv,attr);
        } else tv->my_mvaddch(tv,1, tv->ccol++, c);
    }
    tv->last_pos = pos;
    return 0;
}




static int tv_push_aln(const bam1_t *b, tview_t *tv)
{
    /* If we are restricted to specific readgroups check RG is in the list */
    if ( tv->rg_hash )
    {
        const uint8_t *rg = bam_aux_get(b, "RG");
        if ( !rg ) return 0; // If we don't have an RG tag exclude read
        khiter_t k = kh_get(kh_rg, tv->rg_hash, (const char*)(rg + 1));
        if ( k == kh_end(tv->rg_hash) ) return 0; // if RG tag is not in list of allowed tags exclude read
    }
    if (tv->no_skip) {
        uint32_t *cigar = bam_get_cigar(b); // this is cheating...
        int i;
        for (i = 0; i <b->core.n_cigar; ++i) {
            if ((cigar[i]&0xf) == BAM_CREF_SKIP)
                cigar[i] = cigar[i]>>4<<4 | BAM_CDEL;
        }
    }
    bam_lplbuf_push(b, tv->lplbuf);
    return 0;
}

int base_draw_aln(tview_t *tv, int tid, int pos)
{
    assert(tv!=NULL);
    // reset
    tv->my_clear(tv);
    tv->curr_tid = tid; tv->left_pos = pos;
    tv->last_pos = tv->left_pos - 1;
    tv->ccol = 0;
    // print ref and consensus
    if (tv->fai) {
        char *str;
        if (tv->ref) free(tv->ref);
        assert(tv->curr_tid>=0);

        str = (char*)calloc(strlen(tv->header->target_name[tv->curr_tid]) + 30, 1);
        assert(str!=NULL);
        sprintf(str, "%s:%d-%d", tv->header->target_name[tv->curr_tid], tv->left_pos + 1, tv->left_pos + tv->mcol);
        tv->ref = fai_fetch(tv->fai, str, &tv->l_ref);
        free(str);
        if ( !tv->ref )
        {
            fprintf(stderr,"Could not read the reference sequence. Is it seekable (plain text or compressed + .gzi indexed with bgzip)?\n");
            exit(1);
        }
    }
    // draw aln
    bam_lplbuf_reset(tv->lplbuf);
    hts_itr_t *iter = sam_itr_queryi(tv->idx, tv->curr_tid, tv->left_pos, tv->left_pos + tv->mcol);
    bam1_t *b = bam_init1();
    while (sam_itr_next(tv->fp, iter, b) >= 0) tv_push_aln(b, tv);
    bam_destroy1(b);
    hts_itr_destroy(iter);
    bam_lplbuf_push(0, tv->lplbuf);

    while (tv->ccol < tv->mcol) {
        int pos = tv->last_pos + 1;
        if (pos%10 == 0 && tv->mcol - tv->ccol >= 10) tv->my_mvprintw(tv,0, tv->ccol, "%-d", pos+1);
        tv->my_mvaddch(tv,1, tv->ccol++, (tv->ref && pos < tv->l_ref)? tv->ref[pos - tv->left_pos] : 'N');
        ++tv->last_pos;
    }
    return 0;
}




static void error(const char *format, ...)
{
    if ( !format )
    {
        fprintf(stderr,
"Usage: samtools tview [options] <aln.bam> [ref.fasta]\n"
"Options:\n"
"   -d display      output as (H)tml or (C)urses or (T)ext \n"
"   -p chr:pos      go directly to this position\n"
"   -s STR          display only reads from this sample or group\n");
        sam_global_opt_help(stderr, "-.--.-");
    }
    else
    {
        va_list ap;
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
    }
    exit(-1);
}

enum dipsay_mode {display_ncurses,display_html,display_text};
extern tview_t* curses_tv_init(const char *fn, const char *fn_fa,
                               const char *samples, const htsFormat *fmt);
extern tview_t* html_tv_init(const char *fn, const char *fn_fa,
                             const char *samples, const htsFormat *fmt);
extern tview_t* text_tv_init(const char *fn, const char *fn_fa,
                             const char *samples, const htsFormat *fmt);

int bam_tview_main(int argc, char *argv[])
{
    int view_mode=display_ncurses;
    tview_t* tv=NULL;
    char *samples=NULL, *position=NULL, *ref;
    int c;

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', 0, '-'),
        { NULL, 0, NULL, 0 }
    };

    while ((c = getopt_long(argc, argv, "s:p:d:", lopts, NULL)) >= 0) {
        switch (c) {
            case 's': samples=optarg; break;
            case 'p': position=optarg; break;
            case 'd':
            {
                switch(optarg[0])
                {
                    case 'H': case 'h': view_mode=display_html;break;
                    case 'T': case 't': view_mode=display_text;break;
                    case 'C': case 'c': view_mode=display_ncurses;break;
                    default: view_mode=display_ncurses;break;
                }
                break;
            }
            default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
                      /* else fall-through */
            case '?': error(NULL);
        }
    }
    if (argc==optind) error(NULL);

    ref = (optind+1>=argc)? ga.reference : argv[optind+1];

    switch(view_mode)
    {
        case display_ncurses:
            tv = curses_tv_init(argv[optind], ref, samples, &ga.in);
            break;

        case display_text:
            tv = text_tv_init(argv[optind], ref, samples, &ga.in);
            break;

        case display_html:
            tv = html_tv_init(argv[optind], ref, samples, &ga.in);
            break;
    }
    if (tv==NULL)
    {
        error("cannot create view");
        return EXIT_FAILURE;
    }

    if ( position )
    {
        int tid, beg, end;
        char *name_lim = (char *) hts_parse_reg(position, &beg, &end);
        if (name_lim) *name_lim = '\0';
        else beg = 0; // region parsing failed, but possibly a seq named "foo:a"
        tid = bam_name2id(tv->header, position);
        if (tid >= 0) { tv->curr_tid = tid; tv->left_pos = beg; }
    }
    else if ( tv->fai )
    {
        // find the first sequence present in both BAM and the reference file
        int i;
        for (i=0; i<tv->header->n_targets; i++)
        {
            if ( faidx_has_seq(tv->fai, tv->header->target_name[i]) ) break;
        }
        if ( i==tv->header->n_targets )
        {
            fprintf(stderr,"None of the BAM sequence names present in the fasta file\n");
            exit(EXIT_FAILURE);
        }
        tv->curr_tid = i;
    }
    tv->my_drawaln(tv, tv->curr_tid, tv->left_pos);
    tv->my_loop(tv);
    tv->my_destroy(tv);

    return EXIT_SUCCESS;
}
