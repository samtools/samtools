#undef _HAVE_CURSES

#if _CURSES_LIB == 0
#elif _CURSES_LIB == 1
#include <curses.h>
#ifndef NCURSES_VERSION
#warning "_CURSES_LIB=1 but NCURSES_VERSION not defined; tview is NOT compiled"
#else
#define _HAVE_CURSES
#endif
#elif _CURSES_LIB == 2
#include <xcurses.h>
#define _HAVE_CURSES
#else
#warning "_CURSES_LIB is not 0, 1 or 2; tview is NOT compiled"
#endif

#ifdef _HAVE_CURSES
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "bam.h"
#include "faidx.h"
#include "bam_maqcns.h"

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

typedef struct {
	int mrow, mcol;
	WINDOW *wgoto, *whelp;

	bam_index_t *idx;
	bam_lplbuf_t *lplbuf;
	bam_header_t *header;
	bamFile fp;
	int curr_tid, left_pos;
	faidx_t *fai;
	bam_maqcns_t *bmc;

	int ccol, last_pos, row_shift, base_for, color_for, is_dot, l_ref, ins, no_skip, show_name;
	char *ref;
} tview_t;

int tv_pl_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
	tview_t *tv = (tview_t*)data;
	int i, j, c, rb, attr, max_ins = 0;
	uint32_t call = 0;
	if (pos < tv->left_pos || tv->ccol > tv->mcol) return 0; // out of screen
	// print referece
	rb = (tv->ref && pos - tv->left_pos < tv->l_ref)? tv->ref[pos - tv->left_pos] : 'N';
	for (i = tv->last_pos + 1; i < pos; ++i) {
		if (i%10 == 0 && tv->mcol - tv->ccol >= 10) mvprintw(0, tv->ccol, "%-d", i+1);
		c = tv->ref? tv->ref[i - tv->left_pos] : 'N';
		mvaddch(1, tv->ccol++, c);
	}
	if (pos%10 == 0 && tv->mcol - tv->ccol >= 10) mvprintw(0, tv->ccol, "%-d", pos+1);
	// print consensus
	call = bam_maqcns_call(n, pl, tv->bmc);
	attr = A_UNDERLINE;
	c = ",ACMGRSVTWYHKDBN"[call>>28&0xf];
	i = (call>>8&0xff)/10+1;
	if (i > 4) i = 4;
	attr |= COLOR_PAIR(i);
	if (c == toupper(rb)) c = '.';
	attron(attr);
	mvaddch(2, tv->ccol, c);
	attroff(attr);
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
						c = bam_aux_getCSi(p->b, p->qpos);
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
						|| (p->b->core.flag & BAM_FSECONDARY)) attr |= A_UNDERLINE;
				if (tv->color_for == TV_COLOR_BASEQ) {
					x = bam1_qual(p->b)[p->qpos]/10 + 1;
					if (x > 4) x = 4;
					attr |= COLOR_PAIR(x);
				} else if (tv->color_for == TV_COLOR_MAPQ) {
					x = p->b->core.qual/10 + 1;
					if (x > 4) x = 4;
					attr |= COLOR_PAIR(x);
				} else if (tv->color_for == TV_COLOR_NUCL) {
					x = bam_nt16_nt4_table[bam1_seqi(bam1_seq(p->b), p->qpos)] + 5;
					attr |= COLOR_PAIR(x);
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
					attr |= COLOR_PAIR(x);
				} else if(tv->color_for == TV_COLOR_COLQ) {
					x = bam_aux_getCQi(p->b, p->qpos);
					if(0 == x) x = bam1_qual(p->b)[p->qpos];
					x = x/10 + 1;
					if (x > 4) x = 4;
					attr |= COLOR_PAIR(x);
				}
				attron(attr);
				mvaddch(row, tv->ccol, bam1_strand(p->b)? tolower(c) : toupper(c));
				attroff(attr);
			}
		}
		c = j? '*' : rb;
		if (c == '*') {
			attr = COLOR_PAIR(8);
			attron(attr);
			mvaddch(1, tv->ccol++, c);
			attroff(attr);
		} else mvaddch(1, tv->ccol++, c);
	}
	tv->last_pos = pos;
	return 0;
}

tview_t *tv_init(const char *fn, const char *fn_fa)
{
	tview_t *tv = (tview_t*)calloc(1, sizeof(tview_t));
	tv->is_dot = 1;
	tv->fp = bam_open(fn, "r");
	bgzf_set_cache_size(tv->fp, 8 * 1024 *1024);
	assert(tv->fp);
	tv->header = bam_header_read(tv->fp);
	tv->idx = bam_index_load(fn);
	if (tv->idx == 0) exit(1);
	tv->lplbuf = bam_lplbuf_init(tv_pl_func, tv);
	if (fn_fa) tv->fai = fai_load(fn_fa);
	tv->bmc = bam_maqcns_init();
	tv->ins = 1;
	bam_maqcns_prepare(tv->bmc);

	initscr();
	keypad(stdscr, TRUE);
	clear();
	noecho();
	cbreak();
	tv->mrow = 24; tv->mcol = 80;
	getmaxyx(stdscr, tv->mrow, tv->mcol);
	tv->wgoto = newwin(3, TV_MAX_GOTO + 10, 10, 5);
	tv->whelp = newwin(29, 40, 5, 5);
	tv->color_for = TV_COLOR_MAPQ;
	start_color();
	init_pair(1, COLOR_BLUE, COLOR_BLACK);
	init_pair(2, COLOR_GREEN, COLOR_BLACK);
	init_pair(3, COLOR_YELLOW, COLOR_BLACK);
	init_pair(4, COLOR_WHITE, COLOR_BLACK);
	init_pair(5, COLOR_GREEN, COLOR_BLACK);
	init_pair(6, COLOR_CYAN, COLOR_BLACK);
	init_pair(7, COLOR_YELLOW, COLOR_BLACK);
	init_pair(8, COLOR_RED, COLOR_BLACK);
	init_pair(9, COLOR_BLUE, COLOR_BLACK);
	return tv;
}

void tv_destroy(tview_t *tv)
{
	delwin(tv->wgoto); delwin(tv->whelp);
	endwin();

	bam_lplbuf_destroy(tv->lplbuf);
	bam_maqcns_destroy(tv->bmc);
	bam_index_destroy(tv->idx);
	if (tv->fai) fai_destroy(tv->fai);
	free(tv->ref);
	bam_header_destroy(tv->header);
	bam_close(tv->fp);
	free(tv);
}

int tv_fetch_func(const bam1_t *b, void *data)
{
	tview_t *tv = (tview_t*)data;
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

int tv_draw_aln(tview_t *tv, int tid, int pos)
{
	// reset
	clear();
	tv->curr_tid = tid; tv->left_pos = pos;
	tv->last_pos = tv->left_pos - 1;
	tv->ccol = 0;
	// print ref and consensus
	if (tv->fai) {
		char *str;
		if (tv->ref) free(tv->ref);
		str = (char*)calloc(strlen(tv->header->target_name[tv->curr_tid]) + 30, 1);
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
		if (pos%10 == 0 && tv->mcol - tv->ccol >= 10) mvprintw(0, tv->ccol, "%-d", pos+1);
		mvaddch(1, tv->ccol++, (tv->ref && pos < tv->l_ref)? tv->ref[pos - tv->left_pos] : 'N');
		++tv->last_pos;
	}
	return 0;
}

static void tv_win_goto(tview_t *tv, int *tid, int *pos)
{
	char str[256], *p;
	int i, l = 0;
	wborder(tv->wgoto, '|', '|', '-', '-', '+', '+', '+', '+');
	mvwprintw(tv->wgoto, 1, 2, "Goto: ");
	for (;;) {
		int c = wgetch(tv->wgoto);
		wrefresh(tv->wgoto);
		if (c == KEY_BACKSPACE || c == '\010' || c == '\177') {
			--l;
		} else if (c == KEY_ENTER || c == '\012' || c == '\015') {
			int _tid = -1, _beg, _end;
			if (str[0] == '=') {
				_beg = strtol(str+1, &p, 10) - 1;
				if (_beg > 0) {
					*pos = _beg;
					return;
				}
			} else {
				bam_parse_region(tv->header, str, &_tid, &_beg, &_end);
				if (_tid >= 0) {
					*tid = _tid; *pos = _beg;
					return;
				}
			}
		} else if (isgraph(c)) {
			if (l < TV_MAX_GOTO) str[l++] = c;
		} else if (c == '\027') l = 0;
		else if (c == '\033') return;
		str[l] = '\0';
		for (i = 0; i < TV_MAX_GOTO; ++i) mvwaddch(tv->wgoto, 1, 8 + i, ' ');
		mvwprintw(tv->wgoto, 1, 8, "%s", str);
	}
}

static void tv_win_help(tview_t *tv) {
	int r = 1;
	WINDOW *win = tv->whelp;
	wborder(win, '|', '|', '-', '-', '+', '+', '+', '+');
	mvwprintw(win, r++, 2, "        -=-    Help    -=- ");
	r++;
	mvwprintw(win, r++, 2, "?          This window");
	mvwprintw(win, r++, 2, "Arrows     Small scroll movement");
	mvwprintw(win, r++, 2, "h,j,k,l    Small scroll movement");
	mvwprintw(win, r++, 2, "H,J,K,L    Large scroll movement");
	mvwprintw(win, r++, 2, "ctrl-H     Scroll 1k left");
	mvwprintw(win, r++, 2, "ctrl-L     Scroll 1k right");
	mvwprintw(win, r++, 2, "space      Scroll one screen");
	mvwprintw(win, r++, 2, "backspace  Scroll back one screen");
	mvwprintw(win, r++, 2, "g          Go to specific location");
	mvwprintw(win, r++, 2, "m          Color for mapping qual");
	mvwprintw(win, r++, 2, "n          Color for nucleotide");
	mvwprintw(win, r++, 2, "b          Color for base quality");
	mvwprintw(win, r++, 2, "c          Color for cs color");
	mvwprintw(win, r++, 2, "z          Color for cs qual");
	mvwprintw(win, r++, 2, ".          Toggle on/off dot view");
	mvwprintw(win, r++, 2, "s          Toggle on/off ref skip");
	mvwprintw(win, r++, 2, "r          Toggle on/off rd name");
	mvwprintw(win, r++, 2, "N          Turn on nt view");
	mvwprintw(win, r++, 2, "C          Turn on cs view");
	mvwprintw(win, r++, 2, "i          Toggle on/off ins");
	mvwprintw(win, r++, 2, "q          Exit");
	r++;
	mvwprintw(win, r++, 2, "Underline:      Secondary or orphan");
	mvwprintw(win, r++, 2, "Blue:    0-9    Green: 10-19");
	mvwprintw(win, r++, 2, "Yellow: 20-29   White: >=30");
	wrefresh(win);
	wgetch(win);
}

void tv_loop(tview_t *tv)
{
	int tid, pos;
	tid = tv->curr_tid; pos = tv->left_pos;
	while (1) {
		int c = getch();
		switch (c) {
			case '?': tv_win_help(tv); break;
			case '\033':
			case 'q': goto end_loop;
			case '/': 
			case 'g': tv_win_goto(tv, &tid, &pos); break;
			case 'm': tv->color_for = TV_COLOR_MAPQ; break;
			case 'b': tv->color_for = TV_COLOR_BASEQ; break;
			case 'n': tv->color_for = TV_COLOR_NUCL; break;
			case 'c': tv->color_for = TV_COLOR_COL; break;
			case 'z': tv->color_for = TV_COLOR_COLQ; break;
			case 's': tv->no_skip = !tv->no_skip; break;
			case 'r': tv->show_name = !tv->show_name; break;
			case KEY_LEFT:
			case 'h': --pos; break;
			case KEY_RIGHT:
			case 'l': ++pos; break;
			case KEY_SLEFT:
			case 'H': pos -= 20; break;
			case KEY_SRIGHT:
			case 'L': pos += 20; break;
			case '.': tv->is_dot = !tv->is_dot; break;
			case 'N': tv->base_for = TV_BASE_NUCL; break;
			case 'C': tv->base_for = TV_BASE_COLOR_SPACE; break;
			case 'i': tv->ins = !tv->ins; break;
			case '\010': pos -= 1000; break;
			case '\014': pos += 1000; break;
			case ' ': pos += tv->mcol; break;
			case KEY_UP:
			case 'j': --tv->row_shift; break;
			case KEY_DOWN:
			case 'k': ++tv->row_shift; break;
			case KEY_BACKSPACE:
			case '\177': pos -= tv->mcol; break;
			case KEY_RESIZE: getmaxyx(stdscr, tv->mrow, tv->mcol); break;
			default: continue;
		}
		if (pos < 0) pos = 0;
		if (tv->row_shift < 0) tv->row_shift = 0;
		tv_draw_aln(tv, tid, pos);
	}
end_loop:
	return;
}

int bam_tview_main(int argc, char *argv[])
{
	tview_t *tv;
	if (argc == 1) {
		fprintf(stderr, "Usage: bamtk tview <aln.bam> [ref.fasta]\n");
		return 1;
	}
	tv = tv_init(argv[1], (argc == 2)? 0 : argv[2]);
	tv_draw_aln(tv, 0, 0);
	tv_loop(tv);
	tv_destroy(tv);
	return 0;
}
#else // #ifdef _HAVE_CURSES
#include <stdio.h>
#warning "No curses library is available; tview is disabled."
int bam_tview_main(int argc, char *argv[])
{
	fprintf(stderr, "[bam_tview_main] The ncurses library is unavailable; tview is not compiled.\n");
	return 1;
}
#endif // #ifdef _HAVE_CURSES
