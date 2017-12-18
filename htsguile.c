#include "htslib/kstring.h"
#include "htsguile.h"

#ifdef HAVE_GUILE
#define UNUSED
/*


 ./samtools view -g '(define read-filter (lambda (R) (integer? (string-contains  (hts-read-seq R) "GCCTTGGCTCTTGTTCC") )))' ~/jeter.bam
 
 
 */
 
static SCM make_mod = SCM_EOL;


static SCM hts_read_name(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=(HtsGuileCtxPtr)scm_to_pointer(scm_ctx);
	return  scm_from_locale_string(bam_get_qname(ptr->b));
	}

static SCM hts_read_flags(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=(HtsGuileCtxPtr)scm_to_pointer(scm_ctx);
	return  scm_from_signed_integer(ptr->b->core.flag);
	}
static SCM hts_read_tid(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=(HtsGuileCtxPtr)scm_to_pointer(scm_ctx);
	return  scm_from_signed_integer(ptr->b->core.tid);
	}

static SCM hts_read_mapq(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=(HtsGuileCtxPtr)scm_to_pointer(scm_ctx);
	return  scm_from_signed_integer((int)ptr->b->core.qual);
	}

static SCM hts_read_contig(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=(HtsGuileCtxPtr)scm_to_pointer(scm_ctx);
    	if ( ptr->header !=NULL &&
    	     ptr->b->core.tid >=0 &&
    	     ptr->b->core.tid < ptr->header->n_targets
    	     )
    		{
    		scm_from_locale_string(ptr->header->target_name[ptr->b->core.tid]);
    		}
    	return SCM_UNDEFINED;
    	}

static SCM hts_read_pos(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=(HtsGuileCtxPtr)scm_to_pointer(scm_ctx);
	return  scm_from_signed_integer(ptr->b->core.pos + 1);
	}


static SCM hts_read_cigar_string(SCM scm_ctx) {
    SCM cigar_str;
    HtsGuileCtxPtr ptr=(HtsGuileCtxPtr)scm_to_pointer(scm_ctx);
    if (ptr->b->core.n_cigar>0) { // cigar
         int i;
          kstring_t str;
         memset(&str, 0, sizeof(kstring_t));
        uint32_t *cigar = bam_get_cigar(ptr->b);
        for (i = 0; i < ptr->b->core.n_cigar; ++i) {
            kputw(bam_cigar_oplen(cigar[i]), &str);
            kputc(bam_cigar_opchr(cigar[i]), &str);
        	}
        cigar_str =  scm_from_locale_string(str.s);
        free(str.s);
        return cigar_str;
        }
    else
    	{
    	return SCM_UNDEFINED;
    	}
     }

static SCM hts_read_seq_length(SCM scm_ctx) {
	HtsGuileCtxPtr ptr=(HtsGuileCtxPtr)scm_to_pointer(scm_ctx);
	return  scm_from_signed_integer(ptr->b->core.l_qseq);
	}

static SCM hts_read_seq(SCM scm_ctx) {
	  HtsGuileCtxPtr ptr=(HtsGuileCtxPtr)scm_to_pointer(scm_ctx);
	  if (ptr->b->core.l_qseq) { // seq and qual
	  	kstring_t str;
	  	SCM seq;
	  	int i;
         	memset(&str, 0, sizeof(kstring_t));
		uint8_t *s = bam_get_seq(ptr->b);
		for (i = 0; i < ptr->b->core.l_qseq; ++i) {
			kputc("=ACMGRSVTWYHKDBN"[bam_seqi(s, i)], &str);
			}

		seq =  scm_from_locale_string(str.s);		
        	free(str.s);
        	return seq;
	   	}	
	    else
	    	{
	    	return SCM_UNDEFINED;
	    	}

	}

static SCM hts_read_is_paired(SCM scm_ctx)
	{
	HtsGuileCtxPtr ptr=(HtsGuileCtxPtr)scm_to_pointer(scm_ctx);
	return scm_from_bool(ptr->b->core.flag & BAM_FPAIRED );
	}
/*

    if ( flag&BAM_FPAIRED ) ksprintf(&str,"%s%s", str.l?",":"","PAIRED");
    if ( flag&BAM_FPROPER_PAIR ) ksprintf(&str,"%s%s", str.l?",":"","PROPER_PAIR");
    if ( flag&BAM_FUNMAP ) ksprintf(&str,"%s%s", str.l?",":"","UNMAP");
    if ( flag&BAM_FMUNMAP ) ksprintf(&str,"%s%s", str.l?",":"","MUNMAP");
    if ( flag&BAM_FREVERSE ) ksprintf(&str,"%s%s", str.l?",":"","REVERSE");
    if ( flag&BAM_FMREVERSE ) ksprintf(&str,"%s%s", str.l?",":"","MREVERSE");
    if ( flag&BAM_FREAD1 ) ksprintf(&str,"%s%s", str.l?",":"","READ1");
    if ( flag&BAM_FREAD2 ) ksprintf(&str,"%s%s", str.l?",":"","READ2");
    if ( flag&BAM_FSECONDARY ) ksprintf(&str,"%s%s", str.l?",":"","SECONDARY");
    if ( flag&BAM_FQCFAIL ) ksprintf(&str,"%s%s", str.l?",":"","QCFAIL");
    if ( flag&BAM_FDUP ) ksprintf(&str,"%s%s", str.l?",":"","DUP");
    if ( flag&BAM_FSUPPLEMENTARY ) ksprintf(&str,"%s%s", str.l?",":"","SUPPLEMENTARY");*/

static void hts_guile_define_module(void *data UNUSED)
	{

	scm_c_define_gsubr ("hts-read-length", 1, 0, 0, hts_read_seq_length);
	scm_c_define_gsubr ("hts-read-name", 1, 0, 0, hts_read_name);
	scm_c_define_gsubr ("hts-read-pos", 1, 0, 0, hts_read_pos);
	scm_c_define_gsubr ("hts-read-seq", 1, 0, 0, hts_read_seq);

	scm_c_export(
		"hts-read-length",
		"hts-read-name",
		"hts-read-pos",
		"hts-read-seq",
		NULL);

	}

void* hts_guile_init() {
	if(scm_is_eq(make_mod,SCM_EOL)) {
		make_mod = scm_c_define_module ("hts", hts_guile_define_module, NULL);
		}
	return NULL;
	}



#endif

