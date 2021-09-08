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
#include <inttypes.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "samtools.h"
#include "bedidx.h"
#include "sam_opts.h"
#include "htslib/thread_pool.h"
#include "htslib/khash.h"

#define WHERE do{fprintf(stderr,"[%s:%d]",__FILE__,__LINE__);} while(0)
#define DEBUG(...) do{WHERE;fprintf(stderr, __VA_ARGS__);fputc('\n',stderr);} while(0)

typedef struct base_info {
     // current record
    bam1_t *b;
    /* do not use hts_pos_t because can be negative due to clipping */
    int64_t ref1;
    hts_pos_t read0;
    int64_t unclipped_read0;
    char op_chr;
    char op_int;
    char base;
    char qual;
    char ref;
} BaseInfo,*BaseInfoPtr;

typedef struct tsv_opt {
    FILE *out;
    htsFile *in;
    faidx_t *fai;
     sam_hdr_t *header;
    int skip_N;
     char* query;
     //current reference interval loaded
     int ref_tid;
     hts_pos_t ref_start;
     hts_pos_t ref_end;
     char* ref_seq;
} TsvOpt,*TsvOptPtr;


/** code from markdup.c but use a signed position */
static int64_t unclipped_start(bam1_t *b) {
    uint32_t *cigar = bam_get_cigar(b);
    int64_t clipped = 0;
    uint32_t i;

    for (i = 0; i < b->core.n_cigar; i++) {
        char c = bam_cigar_opchr(cigar[i]);

        if (c == 'S' || c == 'H') { // clips
            clipped += (int64_t)bam_cigar_oplen(cigar[i]);
        } else {
            break;
        }
    }

    return (int64_t)b->core.pos - clipped + 1;
}




#define REF_BUFFER_SIZE 1000000
static int get_reference_base_at(TsvOptPtr opt, int tid, int ref1,char* base) {
     if (ref1<1 || opt->fai==NULL) return 'N';
     if (!(tid==opt->ref_tid && ref1>=opt->ref_start && ref1<= opt->ref_end)) {
        int len;
          hts_pos_t ref_len = sam_hdr_tid2len(opt->header,tid);
          if (ref1 > ref_len) {
		*base = 'N';
		return 0;
		}
          free(opt->ref_seq);
          opt->ref_start = (ref1 <= REF_BUFFER_SIZE ? 1: ref1 - REF_BUFFER_SIZE);
          opt->ref_end = (ref1 + REF_BUFFER_SIZE >= ref_len ? ref_len : ref1 + REF_BUFFER_SIZE);
        opt->ref_seq = faidx_fetch_seq(opt->fai, sam_hdr_tid2name(opt->header,tid),opt->ref_start-1,opt->ref_end,&len);
        if (opt->ref_seq==NULL) {
            fprintf(stderr,"[tsv] warning. Cannot fetch reference.\n");
            *base = 'N';
	    return -1;
            }
          }
     *base = opt->ref_seq[ref1 - opt->ref_start];
     return 0;
     }

static int sam2tsv_base(TsvOptPtr ctx, BaseInfoPtr aln) {
     DEBUG("ici");
     char* p = ctx->query;
     while(*p!=0) {
          if (p!=ctx->query) fputc('\t',ctx->out);
          switch(*p) {
               #define CASE_OPCODE(OPCODE,LABEL,DESC,FUN) case OPCODE: FUN; break
               #include "bam2tsv.h"
               #undef CASE_OPCODE
               }
          p++;
          }
     return fputc('\n',ctx->out) == EOF ? -1: 0;
     }

#define READ_BASE_AT(i) read_bases==NULL?'N':seq_nt16_str[bam_seqi(read_bases, i)]
#define READ_QUAL_AT(i) (read_quals==NULL?'*':(read_quals[0] == 0xff?'B':read_quals[i]+33))
#define REF_BASE_AT(i) if (get_reference_base_at(opt,b->core.tid,i,&aln.ref)!=0) \
	{\
	print_error("tsv","Cannot fetch base."); \
	ret=-1; \
	break; \
	}

static int sam2tsv_aln(TsvOptPtr opt,bam1_t* b) {
// return value
int ret = 0, i,j;
// read sequence
uint8_t *read_bases = NULL;
uint8_t *read_quals = NULL;
BaseInfo aln;

// skip unmapped records
if ((b->core.flag & BAM_FMUNMAP)) return 0;
// seq and qual
if (b->core.l_qseq) { 
        read_bases = bam_get_seq(b);
        read_quals = bam_get_qual(b);
     }
uint32_t *cig = bam_get_cigar(b);
// one based reference positon
int64_t ref1 = unclipped_start(b);
int n_cigar = b->core.n_cigar;
int64_t read0 = 0;
int64_t unclipped_read0 = 0;
aln.b = b;

//loop over each cigar element
for(i=0;i< n_cigar;i++) {
     aln.op_chr   = bam_cigar_opchr(cig[i]);
     aln.op_int = bam_cigar_op(cig[i]);
     int oplen = bam_cigar_oplen(cig[i]);
     switch(aln.op_int) {
               case BAM_CPAD: break;
               case BAM_CHARD_CLIP:
                    for(j=0;j< oplen;j++) {
                         aln.base = 'N';
                         aln.qual = '.';
                         REF_BASE_AT(ref1)
                         aln.ref1 = ref1;
                         aln.read0 = -1;
			 aln.unclipped_read0 = unclipped_read0;
                         if (sam2tsv_base(opt,&aln)!=0) {
                              ret = -1;
                              break;
                              }
                         ref1++;
                         unclipped_read0++;
                        }
                    break;
               case BAM_CINS:
                   for(j=0;j< oplen;j++) {
                         aln.base = READ_BASE_AT(read0);
                         aln.qual = READ_QUAL_AT(read0);
                         aln.ref = '-';
                         aln.ref1 = -1;
                         aln.read0 = read0;
			 aln.unclipped_read0 = unclipped_read0;

                         if (sam2tsv_base(opt,&aln)!=0) {
                              ret = -1;
                              break;
                              }
                         read0++;
                         unclipped_read0++;
                         }
                    break;
               case BAM_CREF_SKIP:
                    if (opt->skip_N) {
                             ref1 += oplen;
                             break;
                             }
                   // NO break here, continue
               case BAM_CDEL:
                         for(j=0;j< oplen;j++) {
                              aln.base = '-';
                              aln.qual = '-';
                              REF_BASE_AT(ref1)
                              aln.ref1 = ref1;
                              aln.read0 = -1;
			      aln.unclipped_read0 = -1;
                              if (sam2tsv_base(opt,&aln)!=0) {
                                   ret = -1;
                                   break;
                                   }
                              ref1++;
                              }
                    break;
               case BAM_CSOFT_CLIP:
              case BAM_CMATCH:
              case BAM_CEQUAL:
              case BAM_CDIFF:
                    for(j=0;j< oplen;j++) {
                         aln.base = READ_BASE_AT(read0);
                         aln.qual = READ_QUAL_AT(read0);
                         REF_BASE_AT(ref1)
                         aln.ref1 = ref1;
                         aln.read0 = read0;
                         aln.unclipped_read0 = unclipped_read0;
                         if (sam2tsv_base(opt,&aln)!=0) {
                              ret = -1;
                              break;
                              }
                         ref1++;
                         read0++;
                        unclipped_read0++;
                        }
                    break;

               default: {
                         print_error("tsv", "Unsupported cigar operator '%c'.", aln.op_chr);
                         ret =  -1;
                         break;
                         }
               }

          }

return ret;
}

static int sam2tsv_core(TsvOptPtr param) {
    int ret=0, r;
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
    sam_global_opt_help(fp, "-.--T@-.");
    exit(exit_status);
}

int main_bam2tsv(int argc, char *argv[])
{
    int c, ret = EXIT_SUCCESS;
    char* p = NULL;
    TsvOpt param;
    char* out_fname = NULL;
     htsThreadPool pool = {NULL, 0};
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 0, 0,'T', '@'),
        {NULL, 0, NULL, 0}
    };

    memset((void*)&param,0,sizeof(TsvOpt));
     param.out = stdout;
    param.query = strdup("NQFTOBbqRpu");

    while ((c = getopt_long(argc, argv, "Nlo:q:T:",
                            lopts, NULL)) >= 0) {
        switch (c) {
               case 'o': out_fname = optarg; break;
               case 'N': param.skip_N = 1; break;
            case 'q':
                free(param.query);
                param.query = strdup(optarg);
                break;
               case 'l':
                  {
                  #define CASE_OPCODE(OPCODE,LABEL,DESC,FUN) printf("\t%c\t%s\t%s\n",OPCODE,LABEL,DESC)
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

    if (optind==argc  && isatty(STDIN_FILENO)) {
         usage_exit(stderr, EXIT_FAILURE);
          return -1;
          }

    // check query
     p = param.query;
     while(*p!=0) {
          switch(*p) {
               #define CASE_OPCODE(OPCODE,LABEL,DESC,FUN) case OPCODE:break
               #include "bam2tsv.h"
               #undef CASE_OPCODE
            default: print_error("tsv","In query \"%s\" unknown opcode \"%c\".",param.query,*p); return EXIT_FAILURE; break;
               }
          p++;
          }


    /* load reference index */
    if (ga.reference==NULL) {
        print_error("tsv", "undefined reference.\n");
        ret = EXIT_FAILURE;
        goto cleanup;
         }

    char* fn_fai = fai_path(ga.reference);
    if (fn_fai==NULL) {
         print_error("tsv", "Cannot get fasta index.");
        ret = EXIT_FAILURE;
        goto cleanup;
         }
    
    param.fai = fai_load3(ga.reference, fn_fai, NULL, FAI_CREATE);
     if (param.fai==NULL) {
         print_error("tsv", "Cannot load fasta index.");
        ret = EXIT_FAILURE;
        goto cleanup;
         }

     /* open the SAM */
    param.in = sam_open_format(optind==argc?"-":argv[optind], "r", &ga.in);
    if (!param.in) {
        print_error_errno("tsv", "failed to open \"%s\" for input.", argv[optind]);
        ret = EXIT_FAILURE;
        goto cleanup;
         }

    if (fn_fai && hts_set_fai_filename(param.in, fn_fai) != 0) {
        fprintf(stderr, "[tsv] failed to load reference file \"%s\".\n", fn_fai);
        ret = EXIT_FAILURE;
        goto cleanup;
     }

    param.header= sam_hdr_read(param.in);
    if (param.header==NULL) {
          print_error("tsv", "cannot read SAM header.\n");
        ret = EXIT_FAILURE;
        goto cleanup;
     }


     if (out_fname!=NULL) {
          param.out = fopen(out_fname,"w");
          if (!param.out) {
               print_error_errno("tsv", "failed to open \"%s\" for writing", out_fname);
               ret = EXIT_FAILURE;
               goto cleanup;
               }
             }

     //print header
     p = param.query;
     while(*p!=0) {
          fputc(p==param.query?'#':'\t',param.out);
          switch(*p) {
               #define CASE_OPCODE(OPCODE,LABEL,DESC,FUN) case OPCODE: fputs(LABEL,param.out); break
               #include "bam2tsv.h"
               #undef CASE_OPCODE
               }
          p++;
          }
     fputc('\n',param.out);

    if (ga.nthreads > 0)  {
        if (!(pool.pool = hts_tpool_init(ga.nthreads))) {
            fprintf(stderr, "[tsv] error creating thread pool\n");
            return 1;
        }

        hts_set_opt(param.in,  HTS_OPT_THREAD_POOL, &pool);
    }


    ret = sam2tsv_core(&param);

     cleanup:
    if (param.in != NULL) {
         sam_close(param.in);
          }
     if (out_fname!=NULL) {
          fflush(param.out);
          fclose(param.out);
          free(out_fname);
          } 
    if (param.header!=NULL) {
         sam_hdr_destroy(param.header);
         }
    if (param.fai) fai_destroy(param.fai);
    if (pool.pool) hts_tpool_destroy(pool.pool);
    sam_global_args_free(&ga);
    free(param.query);
    return ret;
}

