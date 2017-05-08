/*  bam_tview_curses.c -- curses tview implementation.

    Copyright (C) 2008-2013 Genome Research Ltd.
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

#include "bam_tview.h"

#ifdef HAVE_CURSES

#if defined HAVE_NCURSESW_CURSES_H
#include <ncursesw/curses.h>
#elif defined HAVE_NCURSESW_H
#include <ncursesw.h>
#elif defined HAVE_NCURSES_CURSES_H
#include <ncurses/curses.h>
#elif defined HAVE_NCURSES_H
#include <ncurses.h>
#elif defined HAVE_CURSES_H
#include <curses.h>
#endif

typedef struct CursesTview {
    tview_t view;
    WINDOW *wgoto, *whelp;
    } curses_tview_t;

#define FROM_TV(ptr) ((curses_tview_t*)ptr)

static void curses_destroy(tview_t* base)
    {
    curses_tview_t* tv=(curses_tview_t*)base;


    delwin(tv->wgoto); delwin(tv->whelp);
    endwin();

    base_tv_destroy(base);

    free(tv);
    }

/*
 void (*my_mvprintw)(struct AbstractTview* ,int,int,const char*,...);
    void (*my_)(struct AbstractTview*,int,int,int);
    void (*my_attron)(struct AbstractTview*,int);
    void (*my_attroff)(struct AbstractTview*,int);
    void (*my_clear)(struct AbstractTview*);
    int (*my_colorpair)(struct AbstractTview*,int);
*/

static void curses_mvprintw(struct AbstractTview* tv,int y ,int x,const char* fmt,...)
    {
    unsigned int size=tv->mcol+2;
    char* str=malloc(size);
    if(str==0) exit(EXIT_FAILURE);
    va_list argptr;
    va_start(argptr, fmt);
    vsnprintf(str,size, fmt, argptr);
    va_end(argptr);
    mvprintw(y,x,str);
    free(str);
    }

static void curses_mvaddch(struct AbstractTview* tv,int y,int x,int ch)
    {
    mvaddch(y,x,ch);
    }

static void curses_attron(struct AbstractTview* tv,int flag)
    {
    attron(flag);
    }
static void curses_attroff(struct AbstractTview* tv,int flag)
    {
    attroff(flag);
    }
static void curses_clear(struct AbstractTview* tv)
    {
    clear();
    }

static int curses_init_colors(int inverse)
{
    if (inverse) {
        init_pair(1, COLOR_WHITE, COLOR_BLUE);
        init_pair(2, COLOR_BLACK, COLOR_GREEN);
        init_pair(3, COLOR_BLACK, COLOR_YELLOW);
        init_pair(4, COLOR_BLACK, COLOR_WHITE);
        init_pair(5, COLOR_BLACK, COLOR_GREEN);
        init_pair(6, COLOR_BLACK, COLOR_CYAN);
        init_pair(7, COLOR_WHITE, COLOR_MAGENTA);
        init_pair(8, COLOR_WHITE, COLOR_RED);
        init_pair(9, COLOR_WHITE, COLOR_BLUE);
    } else {
        init_pair(1, COLOR_BLUE, COLOR_BLACK);
        init_pair(2, COLOR_GREEN, COLOR_BLACK);
        init_pair(3, COLOR_YELLOW, COLOR_BLACK);
        init_pair(4, COLOR_WHITE, COLOR_BLACK);
        init_pair(5, COLOR_GREEN, COLOR_BLACK);
        init_pair(6, COLOR_CYAN, COLOR_BLACK);
        init_pair(7, COLOR_MAGENTA, COLOR_BLACK);
        init_pair(8, COLOR_RED, COLOR_BLACK);
        init_pair(9, COLOR_BLUE, COLOR_BLACK);
    }

    return 0;
}

static int curses_colorpair(struct AbstractTview* tv,int flag)
    {
    return COLOR_PAIR(flag);
    }

static int curses_drawaln(struct AbstractTview* tv, int tid, int pos)
    {
    return base_draw_aln(tv,  tid, pos);
    }



static void tv_win_goto(curses_tview_t *tv, int *tid, int *pos)
    {
    char str[256], *p;
    int i, l = 0;
    tview_t *base=(tview_t*)tv;
    wborder(tv->wgoto, '|', '|', '-', '-', '+', '+', '+', '+');
    mvwprintw(tv->wgoto, 1, 2, "Goto: ");
    for (;;) {
        int invalid = 0;
        int c = wgetch(tv->wgoto);
        wrefresh(tv->wgoto);
        if (c == KEY_BACKSPACE || c == '\010' || c == '\177') {
            if(l > 0) --l;
        } else if (c == KEY_ENTER || c == '\012' || c == '\015') {
            int _tid = -1, _beg, _end;
            if (str[0] == '=') {
                _beg = strtol(str+1, &p, 10) - 1;
                if (_beg > 0) {
                    *pos = _beg;
                    return;
                }
            } else {
                char *name_lim = (char *) hts_parse_reg(str, &_beg, &_end);
                if (name_lim) {
                    char name_terminator = *name_lim;
                    *name_lim = '\0';
                    _tid = bam_name2id(base->header, str);
                    *name_lim = name_terminator;
                }
                else {
                    // Unparsable region, but possibly a sequence named "foo:a"
                    _tid = bam_name2id(base->header, str);
                    _beg = 0;
                }

                if (_tid >= 0) {
                    *tid = _tid; *pos = _beg;
                    return;
                }
            }

            // If we get here, the region string is invalid
            invalid = 1;
        } else if (isgraph(c)) {
            if (l < TV_MAX_GOTO) str[l++] = c;
        } else if (c == '\027') l = 0;
        else if (c == '\033') return;
        str[l] = '\0';
        for (i = 0; i < TV_MAX_GOTO; ++i) mvwaddch(tv->wgoto, 1, 8 + i, ' ');
        if (invalid) mvwprintw(tv->wgoto, 1, TV_MAX_GOTO - 1, "[Invalid]");
        mvwprintw(tv->wgoto, 1, 8, "%s", str);
    }
}




static void tv_win_help(curses_tview_t *tv) {
    int r = 1;
    tview_t* base=(tview_t*)base;
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
    mvwprintw(win, r++, 2, "v          Inverse video");
    mvwprintw(win, r++, 2, "q          Exit");
    r++;
    mvwprintw(win, r++, 2, "Underline:      Secondary or orphan");
    mvwprintw(win, r++, 2, "Blue:    0-9    Green: 10-19");
    mvwprintw(win, r++, 2, "Yellow: 20-29   White: >=30");
    wrefresh(win);
    wgetch(win);
}

static int curses_underline(tview_t* tv)
    {
    return A_UNDERLINE;
    }

static int curses_loop(tview_t* tv)
    {
    int tid, pos;
    curses_tview_t *CTV=(curses_tview_t *)tv;
    tid = tv->curr_tid; pos = tv->left_pos;
    while (1) {
        int c = getch();
        switch (c) {
            case '?': tv_win_help(CTV); break;
            case '\033':
            case 'q': goto end_loop;
            case '/':
            case 'g': tv_win_goto(CTV, &tid, &pos); break;
            case 'm': tv->color_for = TV_COLOR_MAPQ; break;
            case 'b': tv->color_for = TV_COLOR_BASEQ; break;
            case 'n': tv->color_for = TV_COLOR_NUCL; break;
            case 'c': tv->color_for = TV_COLOR_COL; break;
            case 'z': tv->color_for = TV_COLOR_COLQ; break;
            case 'v': curses_init_colors(tv->inverse = !tv->inverse); break;
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
            case 'J': tv->row_shift -= 20; break;
            case KEY_DOWN:
            case 'k': ++tv->row_shift; break;
            case 'K': tv->row_shift += 20; break;
            case KEY_BACKSPACE:
            case '\177': pos -= tv->mcol; break;
            case KEY_RESIZE: getmaxyx(stdscr, tv->mrow, tv->mcol); break;
            default: continue;
        }
        if (pos < 0) pos = 0;
        if (tv->row_shift < 0) tv->row_shift = 0;
        tv->my_drawaln(tv, tid, pos);
    }
end_loop:
    return 0;
}




tview_t* curses_tv_init(const char *fn, const char *fn_fa, const char *samples,
                        const htsFormat *fmt)
    {
    curses_tview_t *tv = (curses_tview_t*)calloc(1, sizeof(curses_tview_t));
    tview_t* base=(tview_t*)tv;
    if(tv==0)
        {
        fprintf(stderr,"Calloc failed\n");
        return 0;
        }

    base_tv_init(base,fn,fn_fa,samples,fmt);
    /* initialize callbacks */
#define SET_CALLBACK(fun) base->my_##fun=curses_##fun;
    SET_CALLBACK(destroy);
    SET_CALLBACK(mvprintw);
    SET_CALLBACK(mvaddch);
    SET_CALLBACK(attron);
    SET_CALLBACK(attroff);
    SET_CALLBACK(clear);
    SET_CALLBACK(colorpair);
    SET_CALLBACK(drawaln);
    SET_CALLBACK(loop);
    SET_CALLBACK(underline);
#undef SET_CALLBACK

    initscr();
    keypad(stdscr, TRUE);
    clear();
    noecho();
    cbreak();

    getmaxyx(stdscr, base->mrow, base->mcol);
    tv->wgoto = newwin(3, TV_MAX_GOTO + 10, 10, 5);
    tv->whelp = newwin(30, 40, 5, 5);

    start_color();
    curses_init_colors(0);
    return base;
    }

#else // !HAVE_CURSES

#warning "No curses library is available; tview with curses is disabled."

extern tview_t* text_tv_init(const char *fn, const char *fn_fa, const char *samples,
                             const htsFormat *fmt);

tview_t* curses_tv_init(const char *fn, const char *fn_fa, const char *samples,
                        const htsFormat *fmt)
    {
    return text_tv_init(fn,fn_fa,samples,fmt);
    }

#endif
