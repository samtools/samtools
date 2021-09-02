CASE_OPCODE('N',"Read Bame",fputs(bam_get_qname(aln->b), ctx->out));
CASE_OPCODE('Q',"Read Qual",fprintf(ctx->out,"%d",aln->b->core.qual));
CASE_OPCODE('F',"Read Flag",fprintf(ctx->out,"%d",aln->b->core.flag));
CASE_OPCODE('R',"Contig name",fprintf(ctx->out,"%s",sam_hdr_tid2name(ctx->header,aln->b->core.tid)));
CASE_OPCODE('O',"Cigar Operator",fprintf(ctx->out,"%c",aln->op));
