#include "bam_tview.h"

typedef struct HtmlTview {
	tview_t view;
	int row_count;
	tixel_t** screen;
	FILE* out;
	int attributes;/* color... */
	} html_tview_t;

#define FROM_TV(ptr) ((html_tview_t*)ptr)

static void html_destroy(tview_t* base)
	{
	int i;
	html_tview_t* tv=(html_tview_t*)base;
	if(tv->screen!=NULL)
		{
		for(i=0;i< tv->row_count;++i) free(tv->screen[i]);
		free(tv->screen);
		}
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

static void html_mvprintw(struct AbstractTview* tv,int y ,int x,const char* fmt,...)
	{
	int i,nchars=0;
	unsigned int size=tv->mcol+2;
	char* str=malloc(size);
	if(str==0) exit(EXIT_FAILURE);
	va_list argptr;
  	va_start(argptr, fmt);
	nchars=vsnprintf(str,size, fmt, argptr);
	va_end(argptr);
	
	for(i=0;i< nchars;++i)
		{
		tv->my_mvaddch(tv,y,x+i,str[i]);
		}
	free(str);
	}

static void html_mvaddch(struct AbstractTview* tv,int y,int x,int ch)
	{
	tixel_t* row=NULL;
	html_tview_t* ptr=FROM_TV(tv);
	if( x >= tv->mcol ) return; //out of screen
	while(ptr->row_count<=y)
		{
		int x;
		row=(tixel_t*)calloc(tv->mcol,sizeof(tixel_t));
		if(row==0)  exit(EXIT_FAILURE);
		for(x=0;x<tv->mcol;++x) {row[x].ch=' ';row[x].attributes=0;}
		ptr->screen=(tixel_t**)realloc(ptr->screen,sizeof(tixel_t*)*(ptr->row_count+1));
		ptr->screen[ptr->row_count++]=row;
		}
	row=ptr->screen[y];
	row[x].ch=ch;
	row[x].attributes=ptr->attributes;
	}
    	
static void html_attron(struct AbstractTview* tv,int flag)
    {
    html_tview_t* ptr=FROM_TV(tv);
    ptr->attributes |= 1 << flag;
    }
   
static void html_attroff(struct AbstractTview* tv,int flag)
    {
    html_tview_t* ptr=FROM_TV(tv);
    ptr->attributes &= ~(1 << flag);
    }
    
static void html_clear(struct AbstractTview* tv)
    {
    html_tview_t* ptr=FROM_TV(tv);
    if(ptr->screen!=NULL)
	{
	int i;
	for(i=0;i< ptr->row_count;++i) free(ptr->screen[i]);
	free(ptr->screen);
	ptr->screen=NULL;
	}
    ptr->row_count=0;
    ptr->attributes=0;
    }
    
static int html_colorpair(struct AbstractTview* tv,int flag)
    {
    return flag;
    }

static int html_drawaln(struct AbstractTview* tv, int tid, int pos)
    {
    int y,x;
    html_tview_t* ptr=FROM_TV(tv);
    html_clear(tv);
    base_draw_aln(tv,  tid, pos);
    fputs("<html><head>",ptr->out);
    
    //style
    fputs("<style type='text/css'>",ptr->out);
    #define CSS(id,col) fprintf(ptr->out,".tviewc%d {color:%s;}\n.tviewcu%d {color:%s;text-decoration:underline;}\n",id,col,id,col);
        CSS(0, "black");
    	CSS(1, "blue");
	CSS(2, "green");
	CSS(3, "yellow");
	CSS(4, "black");
	CSS(5, "green");
	CSS(6, "cyan");
	CSS(7, "yellow");
	CSS(8, "red");
	CSS(9, "blue");
    #undef CSS
    fputs("</style>",ptr->out);
    
    fputs("</head><body>",ptr->out);
    fputs("<pre class='tview'>\n",ptr->out);
    for(y=0;y< ptr->row_count;++y)
    	{
    	
    	for(x=0;x< tv->mcol;++x)
	    	{
	    	
		
		if(x== 0 || ptr->screen[y][x].attributes!=ptr->screen[y][x-1].attributes)
	    		{
	    		int css=0;
	    		fputs("<span class='",ptr->out);
	    		while(css<10)
	    			{
	    			if((ptr->screen[y][x].attributes & (1 << css))!=0)
	    				{
	    				break;
	    				}
	    			++css;
	    			}
	    		if(css>=10) css=0;
	    		fprintf(ptr->out,"tviewc%d",css);
	    		fputs("'>",ptr->out);
	    		}
		
		int ch=ptr->screen[y][x].ch;
		switch(ch)
			{
			case '<': fputs("&lt;",ptr->out);break;
			case '>': fputs("&gt;",ptr->out);break;
			case '&': fputs("&amp;",ptr->out);break;
			default: fputc(ch,ptr->out); break;
			}
	    	
	    	
	    	if(x+1 == tv->mcol  || ptr->screen[y][x].attributes!=ptr->screen[y][x+1].attributes)
	    		{
	    		fputs("</span>",ptr->out);
	    		}
	    	}
    	if(y+1 < ptr->row_count) fputs("<br/>",ptr->out);
    	}
    fputs("</pre></body></html>",ptr->out);
    return 0;
    }

static int text_drawaln(struct AbstractTview* tv, int tid, int pos)
    {
    int y,x;
    html_tview_t* ptr=FROM_TV(tv);
    html_clear(tv);
    base_draw_aln(tv,  tid, pos);
    for(y=0;y< ptr->row_count;++y)
    	{
    	for(x=0;x< tv->mcol;++x)
	    	{
	    	int ch=ptr->screen[y][x].ch;

	    	fputc(ch,ptr->out);
	    	}
    	fputc('\n',ptr->out);
    	}
    return 0;
    }


static int html_loop(tview_t* tv)
	{
	//tv->my_drawaln(tv, tv->curr_tid, tv->left_pos);
	return 0;	
	}

static int html_underline(tview_t* tv)
	{
	return 11;	
	}

/*
static void init_pair(html_tview_t *tv,int id_ge_1, const char* pen, const char* paper)
	{
	
	}
*/

tview_t* html_tv_init(const char *fn, const char *fn_fa, const char *samples)
	{
	char* colstr=getenv("COLUMNS");
	html_tview_t *tv = (html_tview_t*)calloc(1, sizeof(html_tview_t));
	tview_t* base=(tview_t*)tv;
	if(tv==0)
		{
		fprintf(stderr,"Calloc failed\n");
		return 0;
		}
	tv->row_count=0;
	tv->screen=NULL;
	tv->out=stdout;
	tv->attributes=0;
	base_tv_init(base,fn,fn_fa,samples);
	/* initialize callbacks */
#define SET_CALLBACK(fun) base->my_##fun=html_##fun;
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

	
	if(colstr!=0)
		{
		base->mcol=atoi(colstr);
		if(base->mcol<10) base->mcol=80;
		}
	base->mrow=99999;
	
/*
	init_pair(tv,1, "blue", "white");
	init_pair(tv,2, "green", "white");
	init_pair(tv,3, "yellow", "white");
	init_pair(tv,4, "white", "white");
	init_pair(tv,5, "green", "white");
	init_pair(tv,6, "cyan", "white");
	init_pair(tv,7, "yellow", "white");
	init_pair(tv,8, "red", "white");
	init_pair(tv,9, "blue", "white");
	*/
	return base;
	}


tview_t* text_tv_init(const char *fn, const char *fn_fa, const char *samples)
	{
	tview_t* tv=html_tv_init(fn,fn_fa,samples);
	tv->my_drawaln=text_drawaln;
	return tv;
	}

