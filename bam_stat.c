#include <unistd.h>
#include <errno.h>
#include <assert.h>
#include "bam.h"

typedef struct {
	long long n_reads[2], n_mapped[2], n_pair_all[2], n_pair_map[2], n_pair_good[2];
	long long n_sgltn[2], n_read1[2], n_read2[2];
	long long n_dup[2];
	long long n_diffchr[2], n_diffhigh[2];
} bam_flagstat_t;

typedef struct bamStatApp
	{
	FILE* out;
	char* filename;
	bamFile fpin;
	bamFile fpout;
	bam_flagstat_t *stats;
	int number_printed;
	void (*my_init)( struct bamStatApp* );
	void (*my_print)( struct bamStatApp*);
	void (*my_finish)( struct bamStatApp* );
	} BamStatApp;


static bam_flagstat_t *bam_flagstat_core(BamStatApp* app)
	{
	int w;
	bam_flagstat_t *s;
	bam1_t *b;
	bam1_core_t *c;
	int ret;
	s = (bam_flagstat_t*)calloc(1, sizeof(bam_flagstat_t));
	if(s==NULL) return NULL;
	b = bam_init1();
	c = &b->core;
	while ((ret = bam_read1(app->fpin, b)) >= 0)
		{
		if( app->fpout!=NULL) bam_write1(app->fpout, b);
		
		w = ((c)->flag & BAM_FQCFAIL)? 1 : 0;
		++(s)->n_reads[w];												\
		if ((c)->flag & BAM_FPAIRED) {
			++(s)->n_pair_all[w];
			if ((c)->flag & BAM_FPROPER_PAIR) ++(s)->n_pair_good[w];
			if ((c)->flag & BAM_FREAD1) ++(s)->n_read1[w];
			if ((c)->flag & BAM_FREAD2) ++(s)->n_read2[w];
			if (((c)->flag & BAM_FMUNMAP) && !((c)->flag & BAM_FUNMAP)) ++(s)->n_sgltn[w];
			if (!((c)->flag & BAM_FUNMAP) && !((c)->flag & BAM_FMUNMAP))
				{
				++(s)->n_pair_map[w];
				if ((c)->mtid != (c)->tid)
					{
					++(s)->n_diffchr[w];
					if ((c)->qual >= 5) ++(s)->n_diffhigh[w];
					}														
				}															
			}																
		if (!((c)->flag & BAM_FUNMAP)) ++(s)->n_mapped[w];				
		if ((c)->flag & BAM_FDUP) ++(s)->n_dup[w];						
			
		}
	bam_destroy1(b);
	if (ret != -1)
		fprintf(stderr, "[bam_flagstat_core] Truncated file? Continue anyway.\n");
	return s;
	}

static void print_bam_flagstat_t(BamStatApp* app)
	{
	FILE* out =app->out;
	const bam_flagstat_t *s=app->stats;
	if(app->number_printed>0) fprintf(out,"\n\n");
	fprintf(out,"File: %s\n",app->filename);
	fprintf(out,"%lld + %lld in total (QC-passed reads + QC-failed reads)\n", s->n_reads[0], s->n_reads[1]);
	fprintf(out,"%lld + %lld duplicates\n", s->n_dup[0], s->n_dup[1]);
	fprintf(out,"%lld + %lld mapped (%.2f%%:%.2f%%)\n", s->n_mapped[0], s->n_mapped[1], (float)s->n_mapped[0] / s->n_reads[0] * 100.0, (float)s->n_mapped[1] / s->n_reads[1] * 100.0);
	fprintf(out,"%lld + %lld paired in sequencing\n", s->n_pair_all[0], s->n_pair_all[1]);
	fprintf(out,"%lld + %lld read1\n", s->n_read1[0], s->n_read1[1]);
	fprintf(out,"%lld + %lld read2\n", s->n_read2[0], s->n_read2[1]);
	fprintf(out,"%lld + %lld properly paired (%.2f%%:%.2f%%)\n", s->n_pair_good[0], s->n_pair_good[1], (float)s->n_pair_good[0] / s->n_pair_all[0] * 100.0, (float)s->n_pair_good[1] / s->n_pair_all[1] * 100.0);
	fprintf(out,"%lld + %lld with itself and mate mapped\n", s->n_pair_map[0], s->n_pair_map[1]);
	fprintf(out,"%lld + %lld singletons (%.2f%%:%.2f%%)\n", s->n_sgltn[0], s->n_sgltn[1], (float)s->n_sgltn[0] / s->n_pair_all[0] * 100.0, (float)s->n_sgltn[1] / s->n_pair_all[1] * 100.0);
	fprintf(out,"%lld + %lld with mate mapped to a different chr\n", s->n_diffchr[0], s->n_diffchr[1]);
	fprintf(out,"%lld + %lld with mate mapped to a different chr (mapQ>=5)\n", s->n_diffhigh[0], s->n_diffhigh[1]);
	
	}

static void scan(BamStatApp* app)
	{
	bam_header_t *header;
	
	header = bam_header_read(app->fpin);
	if(app->fpout!=NULL) bam_header_write(app->fpout, header);
	
	app->stats = bam_flagstat_core(app);
	app->my_print(app);
	free(app->stats);
	bam_header_destroy(header);
	}

static void do_nothing( BamStatApp* app)
	{
	}

static void json_init( BamStatApp* app)
	{
	fputc('[',app->out);
	}

static void json_finish( BamStatApp* app)
	{
	fputc(']',app->out);
	fputc('\n',app->out);
	}

static void json_print(BamStatApp* app)
	{
	FILE* out =app->out;
	const bam_flagstat_t *s=app->stats;
	if(app->number_printed>0) fputc(',',out);
	
	fprintf(out,"{\"file\":\"%s\" ",app->filename);
	fprintf(out,",\"total\":{\"pass\":%lld,\"fail\":%lld}", s->n_reads[0], s->n_reads[1]);
	fprintf(out,",\"duplicates\":{\"pass\":%lld,\"fail\":%lld}", s->n_dup[0], s->n_dup[1]);


	fprintf(out,",\"mapped\":{\"percent_pass\":%.2f,\"pass\":%lld,\"percent_fail\":%.2f,\"fail\":%lld}",
		(s->n_reads[0]==0?0.0f:(float)s->n_mapped[0] / s->n_reads[0] * 100.0),
		s->n_mapped[0],
		(s->n_reads[1]==0?0.0f:(float)s->n_mapped[1] / s->n_reads[1] * 100.0),
		s->n_mapped[1]
		);
		
	fprintf(out,",\"paired\":{\"pass\":%lld,\"fail\":%lld}", s->n_pair_all[0], s->n_pair_all[1]);
	fprintf(out,",\"read1\":{\"pass\":%lld,\"fail\":%lld}", s->n_read1[0], s->n_read1[1]);
	fprintf(out,",\"read2\":{\"pass\":%lld,\"fail\":%lld}", s->n_read2[0], s->n_read2[1]);

	fprintf(out,",\"properly-paired\":{\"percent_pass\":%.2f,\"pass\":%lld,\"percent_fail\":%.2f,\"fail\":%lld}",
		(s->n_pair_all[0]==0?0.0f:(float)s->n_pair_good[0] / s->n_pair_all[0] * 100.0),
		s->n_pair_good[0],
		(s->n_pair_all[1]==0?0.0f:(float)s->n_pair_good[1] / s->n_pair_all[1] * 100.0),
		s->n_pair_good[1]
		);
	fprintf(out,",\"pair_map\":{\"pass\":%lld,\"fail\":%lld}", s->n_pair_map[0], s->n_pair_map[1]);
	
	fprintf(out,",\"singleton\":{\"percent_pass\":%.2f,\"pass\":%lld,\"percent_fail\":%.2f,\"fail\":%lld}",
		(s->n_pair_all[0]==0?0.0f:(float)s->n_sgltn[0] / s->n_pair_all[0] * 100.0),
		s->n_sgltn[0],
		(s->n_pair_all[1]==0?0.0f:(float)s->n_sgltn[1] / s->n_pair_all[1] * 100.0),
		s->n_sgltn[1]
		);
	fprintf(out,",\"diffchr\":{\"pass\":%lld,\"fail\":%lld}", s->n_diffchr[0], s->n_diffchr[1]);
	fprintf(out,",\"diffchrhigh\":{\"pass\":%lld,\"fail\":%lld}", s->n_diffhigh[0], s->n_diffhigh[1]);

	fprintf(out,"}");
	}


static void xml_init( BamStatApp* app)
	{
	fputs("<flagstat>",app->out);
	}

static void xml_print(BamStatApp* app)
	{
	FILE* out =app->out;
	const bam_flagstat_t *s=app->stats;

	fprintf(out,"<input file=\"%s\">",app->filename);
	fprintf(out,"<property key=\"total\"><pass>%lld</pass><fail>%lld</fail></property>", s->n_reads[0], s->n_reads[1]);
	fprintf(out,"<property key=\"duplicates\"><pass>%lld</pass><fail>%lld</fail></property>", s->n_dup[0], s->n_dup[1]);
	fprintf(out,"<property key=\"mapped\"><pass percent=\"%.2f\">%lld</pass><fail percent=\"%.2f\">%lld</fail></property>",
		(s->n_reads[0]==0?0.0f:(float)s->n_mapped[0] / s->n_reads[0] * 100.0),
		s->n_mapped[0],
		(s->n_reads[1]==0?0.0f:(float)s->n_mapped[1] / s->n_reads[1] * 100.0),
		s->n_mapped[1]
		);
	fprintf(out,"<property key=\"paired\"><pass>%lld</pass><fail>%lld</fail></property>", s->n_pair_all[0], s->n_pair_all[1]);
	fprintf(out,"<property key=\"read1\"><pass>%lld</pass><fail>%lld</fail></property>", s->n_read1[0], s->n_read1[1]);
	fprintf(out,"<property key=\"read2\"><pass>%lld</pass><fail>%lld</fail></property>", s->n_read2[0], s->n_read2[1]);

	fprintf(out,"<property key=\"properly-paired\"><pass percent=\"%.2f\">%lld</pass><fail percent=\"%.2f\">%lld</fail></property>",
		(s->n_pair_all[0]==0?0.0f:(float)s->n_pair_good[0] / s->n_pair_all[0] * 100.0),
		s->n_pair_good[0],
		(s->n_pair_all[1]==0?0.0f:(float)s->n_pair_good[1] / s->n_pair_all[1] * 100.0),
		s->n_pair_good[1]
		);
	fprintf(out,"<property key=\"pair_map\"><pass>%lld</pass><fail>%lld</fail></property>", s->n_pair_map[0], s->n_pair_map[1]);
	
	fprintf(out,"<property key=\"singleton\"><pass percent=\"%.2f\">%lld</pass><fail percent=\"%.2f\">%lld</fail></property>",
		(s->n_pair_all[0]==0?0.0f:(float)s->n_sgltn[0] / s->n_pair_all[0] * 100.0),
		s->n_sgltn[0],
		(s->n_pair_all[1]==0?0.0f:(float)s->n_sgltn[1] / s->n_pair_all[1] * 100.0),
		s->n_sgltn[1]
		);
	fprintf(out,"<property key=\"diffchr\"><pass>%lld</pass><fail>%lld</fail></property>", s->n_diffchr[0], s->n_diffchr[1]);
	fprintf(out,"<property key=\"diffchrhigh\"><pass>%lld</pass><fail>%lld</fail></property>", s->n_diffhigh[0], s->n_diffhigh[1]);

	fprintf(out,"</input>\n");
	}

static void xml_finish( BamStatApp* app)
	{
	fputs("</flagstat>\n",app->out);
	}

int bam_flagstat(int argc, char *argv[])
	{
	char* filenameout=NULL;
	int c;
	BamStatApp app;
	memset(&app, 0, sizeof(BamStatApp));
	app.my_init=do_nothing;
	app.my_print=print_bam_flagstat_t;
	app.my_finish=do_nothing;
	app.out=stdout;
	int streaming=0;
	while ((c = getopt(argc, argv, "o:sf:")) >= 0) {
		switch (c)
			{
			case 'o': filenameout=optarg; break;
			case 's': streaming=1;break;
			case 'f':
				{
				switch(optarg[0])
					{
					case 'j':case 'J':
						{
						app.my_init=json_init;
						app.my_print=json_print;
						app.my_finish=json_finish;
						break;
						}
					case 'x':case 'X':
						{
						app.my_init=xml_init;
						app.my_print=xml_print;
						app.my_finish=xml_finish;
						break;
						}
					default:
						{
						app.my_init=do_nothing;
						app.my_print=print_bam_flagstat_t;
						app.my_finish=do_nothing;
						break;
						}
					}
				break;
				}
			case ':': fputs("argument missing\n",stderr); return EXIT_FAILURE;
			case '?': fputs("unknown argument.\n",stderr); return EXIT_FAILURE;
			}
		}
	

	
	if (argc == optind)
		{
		fprintf(stderr, "Usage:\n");
		fprintf(stderr, "    samtools flagstat <in.bam>\n");
		fprintf(stderr, "    samtools flagstat -o report.txt <in1.bam> <in2.bam> ... \n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, " -o (filename) save report to file:\n");
		fprintf(stderr, " -s write input to stdout (-o required)\n");
		fprintf(stderr, " -f (format) (j)son (x)ml. default: text\n");
		return 1;
		}
	
	/* open report file */
	if(filenameout!=NULL)
		{
		app.out=fopen(filenameout,"w");
		if(app.out==NULL)
			{
			fprintf(stderr,"Cannot open \"%s\" : %s\n",filenameout,strerror(errno));
			return EXIT_FAILURE;
			}
		}	
		
	
	if(optind+1==argc)
		{
		/* streaming bam to stdout */
		if(streaming==1)
			{
			/** report need to be defined to be saved somewhere */
			if(filenameout==NULL)
				{
				fputs("streaming but output filename undefined.\n",stderr);
				return EXIT_FAILURE;
				}
			/* open bam to stdout */
			app.fpout= bam_dopen(fileno(stdout), "wb");
			if( app.fpout== 0)
				{
				fprintf(stderr,"cannot write BAM to stdout %s.\n",strerror(errno));
				return EXIT_FAILURE;
				}
			}
		app.filename=argv[optind];
		/* open stdin */
		if(strcmp(argv[optind], "-")==0)
			{
			app.fpin = bam_dopen(fileno(stdin), "r");
			}
		else /* open file */
			{
			app.fpin = bam_open(argv[optind], "r");
			}
		/* check fileopen */
		if( app.fpin ==0)
			{
			fprintf(stderr,"cannot read \"%s\" %s.\n",
				argv[optind],
				strerror(errno))
				;
			return EXIT_FAILURE;
			}
		/* scan file */
		app.my_init(&app);
		scan(	&app );
		app.my_finish(&app);
		/* close BAM input */
		bam_close(app.fpin);
		/* close BAM output */
		if(app.fpout!=NULL) bam_close(app.fpout);
		}
	else
		{
		if(streaming==1)
			{
			
			fputs("multiple files : Cannot use streaming.\n",stderr);
			return EXIT_FAILURE;
			}
		app.my_init(&app);
		while(optind< argc)
			{
			/* open BAM file */
			app.fpin = bam_open(argv[optind], "r");
			if( app.fpin ==0)
				{
				fprintf(stderr,"cannot read \"%s\" %s.\n",
					argv[optind],
					strerror(errno))
					;
				return EXIT_FAILURE;
				}
			app.filename= argv[optind];
			/* scan bam input */
			scan(	&app );
			/* close bam input */
			bam_close(app.fpin);
			++optind;
			app.number_printed++;
			}
		app.my_finish(&app);
		}
	/* close report */
	if(filenameout!=NULL)
		{
		fclose(app.out);
		}
	return 0;
	}
