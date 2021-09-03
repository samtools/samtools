CASE_OPCODE('N',"ReadName","Read Bame",fputs(bam_get_qname(aln->b), ctx->out));
CASE_OPCODE('Q',"ReadQual","Read Qual",fprintf(ctx->out,"%d",aln->b->core.qual));
CASE_OPCODE('F',"ReadFlag","Read Flag",fprintf(ctx->out,"%d",aln->b->core.flag));
CASE_OPCODE('R',"CHROM","Contig name",fputs(sam_hdr_tid2name(ctx->header,aln->b->core.tid),ctx->out));
CASE_OPCODE('O',"Op","Cigar Operator",fputc(aln->op_chr,ctx->out));
CASE_OPCODE('b',"ReadBase","Read Base", fputc(aln->base,ctx->out));
CASE_OPCODE('q',"QualBase","Quality base", fputc(aln->qual,ctx->out));
