#include <assert.h>
#include "bam_tview.h"

int base_tv_init(tview_t* tv,const char *fn, const char *fn_fa, const char *samples)
	{
	assert(tv!=NULL);
	assert(fn!=NULL);
	tv->mrow = 24; tv->mcol = 80;
	tv->color_for = TV_COLOR_MAPQ;
	tv->is_dot = 1;
	
	tv->fp = bam_open(fn, "r");
	if(tv->fp==0)
		{
		fprintf(stderr,"bam_open %s. %s\n", fn,fn_fa);
            	exit(EXIT_FAILURE);
		}
	bgzf_set_cache_size(tv->fp, 8 * 1024 *1024);
	assert(tv->fp);
	
	tv->header = bam_header_read(tv->fp);
	if(tv->header==0)
		{
		fprintf(stderr,"Cannot read '%s'.\n", fn);
            	exit(EXIT_FAILURE);
		}
	tv->idx = bam_index_load(fn);
	if (tv->idx == 0)
		{
		fprintf(stderr,"Cannot read index for '%s'.\n", fn);
		exit(EXIT_FAILURE);
		}
	tv->lplbuf = bam_lplbuf_init(tv_pl_func, tv);
	if (fn_fa) tv->fai = fai_load(fn_fa);
	tv->bca = bcf_call_init(0.83, 13);
	tv->ins = 1;

    if ( samples ) 
    {
        if ( !tv->header->dict ) tv->header->dict = sam_header_parse2(tv->header->text);
        void *iter = tv->header->dict;
        const char *key, *val;
        int n = 0;
        tv->rg_hash = kh_init(kh_rg);
        while ( (iter = sam_header2key_val(iter, "RG","ID","SM", &key, &val)) )
        {
            if ( !strcmp(samples,key) || (val && !strcmp(samples,val)) )
            {
                khiter_t k = kh_get(kh_rg, tv->rg_hash, key);
                if ( k != kh_end(tv->rg_hash) ) continue;
                int ret;
                k = kh_put(kh_rg, tv->rg_hash, key, &ret);
                kh_value(tv->rg_hash, k) = val;
                n++;
            }
        }
        if ( !n )
        {
            fprintf(stderr,"The sample or read group \"%s\" not present.\n", samples);
            exit(EXIT_FAILURE);
        }
    }

	return 0;
	}


void base_tv_destroy(tview_t* tv)
	{
	bam_lplbuf_destroy(tv->lplbuf);
	bcf_call_destroy(tv->bca);
	bam_index_destroy(tv->idx);
	if (tv->fai) fai_destroy(tv->fai);
	free(tv->ref);
	bam_header_destroy(tv->header);
	bam_close(tv->fp);
	}


int tv_pl_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
	extern unsigned char bam_nt16_table[256];
	tview_t *tv = (tview_t*)data;
	int i, j, c, rb, attr, max_ins = 0;
	uint32_t call = 0;
	if (pos < tv->left_pos || tv->ccol > tv->mcol) return 0; // out of screen
	// print referece
	rb = (tv->ref && pos - tv->left_pos < tv->l_ref)? tv->ref[pos - tv->left_pos] : 'N';
	for (i = tv->last_pos + 1; i < pos; ++i) {
		if (i%10 == 0 && tv->mcol - tv->ccol >= 10) tv->my_mvprintw(tv,0, tv->ccol, "%-d", i+1);
		c = tv->ref? tv->ref[i - tv->left_pos] : 'N';
		tv->my_mvaddch(tv,1, tv->ccol++, c);
	}
	if (pos%10 == 0 && tv->mcol - tv->ccol >= 10) tv->my_mvprintw(tv,0, tv->ccol, "%-d", pos+1);
	{ // call consensus
		bcf_callret1_t bcr;
		int qsum[4], a1, a2, tmp;
		double p[3], prior = 30;
		bcf_call_glfgen(n, pl, bam_nt16_table[rb], tv->bca, &bcr);
		for (i = 0; i < 4; ++i) qsum[i] = bcr.qsum[i]<<2 | i;
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
						if (tv->is_dot && '-' == bam_aux_getCEi(p->b, p->qpos)) c = bam1_strand(p->b)? ',' : '.';
					} else {
						if (tv->show_name) {
							char *name = bam1_qname(p->b);
							c = (p->qpos + 1 >= p->b->core.l_qname)? ' ' : name[p->qpos];
						} else {
							c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
							if (tv->is_dot && toupper(c) == toupper(rb)) c = bam1_strand(p->b)? ',' : '.';
						}
					}
				} else c = p->is_refskip? (bam1_strand(p->b)? '<' : '>') : '*';
			} else { // padding
				if (j > p->indel) c = '*';
				else { // insertion
					if (tv->base_for ==  TV_BASE_NUCL) {
						if (tv->show_name) {
							char *name = bam1_qname(p->b);
							c = (p->qpos + j + 1 >= p->b->core.l_qname)? ' ' : name[p->qpos + j];
						} else {
							c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos + j)];
							if (j == 0 && tv->is_dot && toupper(c) == toupper(rb)) c = bam1_strand(p->b)? ',' : '.';
						}
					} else {
						c = bam_aux_getCSi(p->b, p->qpos + j);
						if (tv->is_dot && '-' == bam_aux_getCEi(p->b, p->qpos + j)) c = bam1_strand(p->b)? ',' : '.';
					}
				}
			}
			if (row > TV_MIN_ALNROW && row < tv->mrow) {
				int x;
				attr = 0;
				if (((p->b->core.flag&BAM_FPAIRED) && !(p->b->core.flag&BAM_FPROPER_PAIR))
						|| (p->b->core.flag & BAM_FSECONDARY)) attr |= tv->my_underline(tv);
				if (tv->color_for == TV_COLOR_BASEQ) {
					x = bam1_qual(p->b)[p->qpos]/10 + 1;
					if (x > 4) x = 4;
					attr |= tv->my_colorpair(tv,x);
				} else if (tv->color_for == TV_COLOR_MAPQ) {
					x = p->b->core.qual/10 + 1;
					if (x > 4) x = 4;
					attr |= tv->my_colorpair(tv,x);
				} else if (tv->color_for == TV_COLOR_NUCL) {
					x = bam_nt16_nt4_table[bam1_seqi(bam1_seq(p->b), p->qpos)] + 5;
					attr |= tv->my_colorpair(tv,x);
				} else if(tv->color_for == TV_COLOR_COL) {
					x = 0;
					switch(bam_aux_getCSi(p->b, p->qpos)) {
						case '0': x = 0; break;
						case '1': x = 1; break;
						case '2': x = 2; break;
						case '3': x = 3; break;
						case '4': x = 4; break;
						default: x = bam_nt16_nt4_table[bam1_seqi(bam1_seq(p->b), p->qpos)]; break;
					}
					x+=5;
					attr |= tv->my_colorpair(tv,x);
				} else if(tv->color_for == TV_COLOR_COLQ) {
					x = bam_aux_getCQi(p->b, p->qpos);
					if(0 == x) x = bam1_qual(p->b)[p->qpos];
					x = x/10 + 1;
					if (x > 4) x = 4;
					attr |= tv->my_colorpair(tv,x);
				}
				tv->my_attron(tv,attr);
				tv->my_mvaddch(tv,row, tv->ccol, bam1_strand(p->b)? tolower(c) : toupper(c));
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




int tv_fetch_func(const bam1_t *b, void *data)
{
	tview_t *tv = (tview_t*)data;
    if ( tv->rg_hash )
    {
        const uint8_t *rg = bam_aux_get(b, "RG");
        if ( !rg ) return 0; 
        khiter_t k = kh_get(kh_rg, tv->rg_hash, (const char*)(rg + 1));
        if ( k == kh_end(tv->rg_hash) ) return 0;
    }
	if (tv->no_skip) {
		uint32_t *cigar = bam1_cigar(b); // this is cheating...
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
	}
	// draw aln
	bam_lplbuf_reset(tv->lplbuf);
	bam_fetch(tv->fp, tv->idx, tv->curr_tid, tv->left_pos, tv->left_pos + tv->mcol, tv, tv_fetch_func);
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
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: bamtk tview [options] <aln.bam> [ref.fasta]\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "   -d display      output as (H)tml or (C)urses or (T)ext \n");
        fprintf(stderr, "   -p chr:pos      go directly to this position\n");
        fprintf(stderr, "   -s STR          display only reads from this sample or group\n");
        fprintf(stderr, "\n\n");
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
extern tview_t* curses_tv_init(const char *fn, const char *fn_fa, const char *samples);
extern tview_t* html_tv_init(const char *fn, const char *fn_fa, const char *samples);
extern tview_t* text_tv_init(const char *fn, const char *fn_fa, const char *samples);

int bam_tview_main(int argc, char *argv[])
	{
	int view_mode=display_ncurses;
	tview_t* tv=NULL;
    char *samples=NULL, *position=NULL;
    int c;
    while ((c = getopt(argc, argv, "s:p:d:")) >= 0) {
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
            default: error(NULL);
        }
    }
	if (argc==optind) error(NULL);
	
	switch(view_mode)
		{
		case display_ncurses:
			{
			tv = curses_tv_init(argv[optind], (optind+1>=argc)? 0 : argv[optind+1], samples);
			break;
			}
		case display_text:
			{
			tv = text_tv_init(argv[optind], (optind+1>=argc)? 0 : argv[optind+1], samples);
			break;
			}
		case display_html:
			{
			tv = html_tv_init(argv[optind], (optind+1>=argc)? 0 : argv[optind+1], samples);
			break;
			}
		}
	if(tv==NULL)
		{
		error("cannot create view");
		return EXIT_FAILURE;
		}
	
	if ( position )
	   	 {
		int _tid = -1, _beg, _end;
		bam_parse_region(tv->header, position, &_tid, &_beg, &_end);
		if (_tid >= 0) { tv->curr_tid = _tid; tv->left_pos = _beg; }
	    	}
	tv->my_drawaln(tv, tv->curr_tid, tv->left_pos);
	tv->my_loop(tv);
	tv->my_destroy(tv);
	
	return EXIT_SUCCESS;
	}
