/*  bam2tsc.c -- depth subcommand.


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

#include <config.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <unistd.h>
#include "htslib/sam.h"
#include "samtools.h"
#include "bedidx.h"
#include "sam_opts.h"
#include "htslib/thread_pool.h"
#include "htslib/khash.h"

// defined in bam_markdup.c
extern hts_pos_t unclipped_start(bam1_t *b);

// From bam_plcmd.c
int read_file_list(const char *file_list, int *n, char **argv[]);



typedef struct base_info {
	// current record
    bam1_t *b;
    //
    hts_pos_t ref1;
    hts_pos_t read0;
    char op;
    char base;
    char qual;
    char ref;
} BaseInfo,*BaseInfoPtr;

typedef struct tsv_opt {
    FILE *out;
    htsFile *in;
	sam_hdr_t *header;
	char* query;
	//current reference interval loaded
	int ref_tid;
	int ref_start;
	int ref_end;
	char* ref_seq;
} TsvOpt,*TsvOptPtr;

#define REF_BUFFER_SIZE 1000000
static char get_reference_base_at(TsvOptPtr opt, int tid, int ref1) {
	if(ref1<1) return 'N';
	if(!(tid==opt->ref_tid && ref1>=opt->ref_start && ref1<= opt->ref_end)) {
		int ref_len = opt->header[tid]->len; 
		if (ref1 > ref_len) return 'N';
		free(opt->ref_seq);
		opt->ref_start = (ref1 <= REF_BUFFER_SIZE ? 1: ref1 - REF_BUFFER_SIZE);
		opt->ref_end = (ref1 + REF_BUFFER_SIZE >= ref_len ? ref_len : ref1 + REF_BUFFER_SIZE);
		}
	return return opt->ref_seq[tid - opt->ref_start];
	}

static int sam2tsv_base(TsvOptPtr ctx, BaseInfoPtr aln) {
	char* p = ctx->query;
	while(*p!=0) {
		if(p!=ctx->query) fputc('\t',ctx->out);
		switch(*p) {
			#define CASE_OPCODE(OPCODE,DESC,FUN) case OPCODE: FUN; break
			#include "bam2tsv.h"
			#undef CASE_OPCODE
			}
		p++;
		}
	return fputc('\n',ctx->out);
	}

#define READ_BASE_AT(i) read_bases==NULL?'N':read_bases[i]
#define REF_BASE_AT(i) get_reference_base_at(opt,b->core.tid,i);
#define READ_QUAL_AT(i) read_quals==NULL?'.':read_quals[i]

static int sam2tsv_aln(TsvOptPtr opt,bam1_t* b) {
// return value
int ret = 0, i,j;
// read sequence
uint8_t *read_bases = NULL;
uint8_t *read_quals = NULL;
BaseInfo aln;
// skip unmapped records
if((b->core.flag & BAM_FMUNMAP)) return 0;
// seq and qual
if (b->core.l_qseq) { 
        read_bases = bam_get_seq(b);
        read_quals = bam_get_qual(b);
	}
uint32_t *cig = bam_get_cigar(b);
// one based reference positon
int ref1 = unclipped_start(b);
int n_cigar = b->core.n_cigar;
int read0 = 0;
aln.b = b;

//loop over each cigar element
for(i=0;i< n_cigar;i++) {
	aln.op   = bam_cigar_op(cig[j]);
	int oplen = bam_cigar_oplen(cig[j]);
	switch(aln.op) {
			case BAM_PADDING: break;
			case BAM_HCLIP:
				for(j=0;j< oplen;j++) {
					aln.base = 'N';
					aln.qual = '.';
					aln.ref = REF_BASE_AT(ref1);
					aln.ref1 = ref1;
					aln.read0 = -1;
					if (sam2tsv_base(opt,&aln)!=0) {
						ret = -1;
						break;
						}
					ref1++;
					}
				break;
			case BAM_INS:
		    	for(j=0;j< oplen;j++) {
					aln.base = READ_BASE_AT(read0);
					aln.qual = READ_QUAL_AT(read0);
					aln.ref = 'N';
					aln.ref1 = -1;
					aln.read0 = read0;
					if (sam2tsv_base(opt,&aln)!=0) {
						ret = -1;
						break;
						}
					read0++;
					}
				break;
			case BAM_CREF_SKIP:
				if(opt_>skip_N) {
		    		ref1 += oplen;
		    		break;
		    		}
		    	// NO break here, continue
			case BAM_CDEL:
					for(j=0;j< oplen;j++) {
						aln.base = '.';
						aln.qual = '.';
						aln.ref = REF_BASE_AT(ref1);
						aln.ref1 = ref1;
						aln.read0 = -1;
						if (sam2tsv_base(opt,&aln)!=0) {
							ret = -1;
							break;
							}
						ref1++;
						}
				break;
			case BAM_SCLIP:
		    case BAM_CMATCH:
		    case BAM_CEQUAL:
		    case BAM_CDIFF:
				for(j=0;j< oplen;j++) {
					aln.base = READ_BASE_AT(read0);
					aln.qual = READ_QUAL_AT(read0);
					aln.ref = REF_BASE_AT(ref1);
					aln.ref1 = ref1;
					aln.read0 = read0;
					if (sam2tsv_base(opt,&aln)!=0) {
						ret = -1;
						break;
						}
					ref1++;
					read0++;
				}
				break; 

			default: {
					print_error("tsv", "Unsupported cigar op '%c'", aln.op);
					ret =  -1;
					break;
					}
			}

		}

return ret;
}

static int sam2tsv_core(TsvOptPtr param) {
    int ret=0;
    bam1_t *b = bam_init1();
	while ((r = sam_read1(param->in, param->header, b)) >= 0) {
		if ( sam2tsv_aln(param,b) != 0) {
			ret = -1;
			break;
			}
		}
	bam_destroy1(b);
    return ret;
	}

static void usage_exit(FILE *fp, int exit_status)
{
    fprintf(fp, "Usage: samtools tsv [options] (in.bam|stdin)\n");
    fprintf(fp, "\nOptions:\n");
    fprintf(fp, "  -o FILE      Write output to FILE [stdout]\n");
    sam_global_opt_help(fp, "-.---@-.");
    exit(exit_status);
}

int main_bam2tsv(int argc, char *argv[])
{
    int nfiles, i, c, ret = EXIT_SUCCESS;
    TsvOpt param;
    char* out_fname = NULL;
	htsThreadPool p = {NULL, 0};
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', '-', '@'),
        {NULL, 0, NULL, 0}
    };

	param.out = stdout;


    while ((c = getopt_long(argc, argv, "o:",
                            lopts, NULL)) >= 0) {
        switch (c) {
			case 'o': out_fname = optarg; break;
			case 'l':
        		{
        		#define CASE_OPCODE(OPCODE,DESC,FUN) printf("\t%c\t%s\n",OPCODE,DESC)
				#include "bam2tsv.h"
				#undef CASE_OPCODE
        		return EXIT_SUCCESS;
        		}
        	
        	default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
        case '?':
            usage_exit(stderr, EXIT_FAILURE);
        }
    }

    param.in = sam_open_format(argv[optind], "r", &ga.in);

    if (!param.in) {
        print_error_errno("tsv", "failed to open \"%s\" for input", argv[optind]);
        return 1;
    }

	if (out_fname!=NULL) {
		params.out = fopen(out_fname,"w");
		if(!params.out) {
			print_error_errno("tsv", "failed to open \"%s\" for writing", out_fname);
        	return EXIT_FAILURE;
			}
   		}

	//print header
	char* p = params.query;
	while(*p!=0) {
		if(p!=ctx->query) fputc('\t',params.out);
		switch(*p) {
			#define CASE_OPCODE(OPCODE,DESC,FUN) case OPCODE: fputs(DESC,params.out); break
			#include "bam2tsv.h"
			#undef CASE_OPCODE
			}
		p++;
		}
	fputc('\n',params.out);

    if (ga.nthreads > 0)  {
        if (!(p.pool = hts_tpool_init(ga.nthreads))) {
            fprintf(stderr, "[tsv] error creating thread pool\n");
            return 1;
        }

        hts_set_opt(param.in,  HTS_OPT_THREAD_POOL, &p);
    }



    ret = sam2tsv_core(&param);

    sam_close(param.in);

	if (out_fname!=NULL) {
		fflush(params.out);
		fclose(params.out);
		free(out_fname);
		} 

    if (p.pool) hts_tpool_destroy(p.pool);
    sam_global_args_free(&ga);

    return ret;
}

