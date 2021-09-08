/**
 This file keep the columns definitions of bam2tsv.c and the way to extract them.
 This file will be included several times by bam2tsv.c in a 'switch' statement , and , depending
 of the macro definition, it will be expanded to different outputs.

 CASE_OPCODE ( char code, column header, description, lambda )
*/
CASE_OPCODE('N',"ReadName","Read Bame",fputs(bam_get_qname(aln->b), ctx->out); if (aln->b->core.flag & BAM_FPAIRED) fputs((aln->b->core.flag & BAM_FREAD1 ?"/1":"/2"),ctx->out));
CASE_OPCODE('Q',"ReadQual","Read Mapping Quality",fprintf(ctx->out,"%"PRId32,aln->b->core.qual));
CASE_OPCODE('F',"ReadFlag","Read Flag",fprintf(ctx->out,"%"PRId32,aln->b->core.flag));
CASE_OPCODE('T',"CHROM","Contig name",fputs(sam_hdr_tid2name(ctx->header,aln->b->core.tid),ctx->out));
CASE_OPCODE('R',"REF","Ref base.",fputc(aln->ref,ctx->out));
CASE_OPCODE('p',"POS","1-based Ref position.",if(aln->op_int == BAM_CINS) fputc('.',ctx->out);  else fprintf(ctx->out,"%"PRId64,aln->ref1));
CASE_OPCODE('O',"CigarOperator","Cigar Operator",fputc(aln->op_chr,ctx->out));
CASE_OPCODE('B',"ReadPos0","0-based Read Position", if(aln->op_int == BAM_CREF_SKIP || aln->op_int == BAM_CDEL || aln->op_int == BAM_CHARD_CLIP)  fputc('.',ctx->out); else fprintf(ctx->out,"%"PRId64,aln->read0));
CASE_OPCODE('u',"UReadPos0","0-based Unclipped Read Position", if(aln->op_int == BAM_CREF_SKIP || aln->op_int == BAM_CDEL )  fputc('.',ctx->out); else fprintf(ctx->out,"%"PRId64,aln->unclipped_read0));
CASE_OPCODE('b',"ALT","Read Base", fputc(aln->base,ctx->out));
CASE_OPCODE('q',"ReadQual","Quality base", fputc(aln->qual,ctx->out));

