CASE_OPCODE('N',"ReadName","Read Bame",fputs(bam_get_qname(aln->b), ctx->out));
CASE_OPCODE('Q',"ReadQual","Read Qual",fprintf(ctx->out,"%d",aln->b->core.qual));
CASE_OPCODE('F',"ReadFlag","Read Flag",fprintf(ctx->out,"%d",aln->b->core.flag));
CASE_OPCODE('T',"CHROM","Contig name",fputs(sam_hdr_tid2name(ctx->header,aln->b->core.tid),ctx->out));
CASE_OPCODE('R',"REF","Ref base.",fputc(aln->ref,ctx->out));
CASE_OPCODE('p',"POS","1-based Ref position.",if(aln->op_int == BAM_CINS) fputc('.',ctx->out);  else fprintf(ctx->out,"%d",aln->ref1));
CASE_OPCODE('O',"CigarOperator","Cigar Operator",fputc(aln->op_chr,ctx->out));
CASE_OPCODE('B',"ReadPos0","0-based Read Position", if(aln->op_int == BAM_CREF_SKIP || aln->op_int == BAM_CDEL || aln->op_int == BAM_CHARD_CLIP)  fputc('.',ctx->out); else fprintf(ctx->out,"%d",aln->read0));
CASE_OPCODE('b',"ReadBase","Read Base", fputc(aln->base,ctx->out));
CASE_OPCODE('q',"ReadQual","Quality base", fputc(aln->qual,ctx->out));

