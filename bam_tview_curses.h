#ifndef BAM_TVIEW_CURSES_H
#define BAM_TVIEW_CURSES_H

#include <curses.h>
#include "bam_tview.h"

typedef struct CursesTview {
	tview_t view;
	WINDOW *wgoto, *whelp;
	} curses_tview_t;


curses_tview_t* curses_tv_init(const char *fn, const char *fn_fa,const char *samples);

#endif

