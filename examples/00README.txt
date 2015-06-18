File ex1.fa contains two sequences cut from the human genome
build36. They were extracted with command:

    samtools faidx human_b36.fa 2:2043966-2045540 20:67967-69550

Sequence names were changed manually for simplicity. File ex1.sam.gz
contains MAQ alignments extracted with:

   (samtools view NA18507_maq.bam 2:2044001-2045500;
    samtools view NA18507_maq.bam 20:68001-69500)

and processed with `samtools fixmate' to make it self-consistent as a
standalone alignment.

To try samtools, you may run the following commands.

Index the reference FASTA.
    samtools faidx ex1.fa

Convert the (headerless) SAM file to BAM.  Note if we had used
"samtools view -h" above to create the ex1.sam.gz then we could omit the
"-t ex1.fa.fai" option here.
    samtools view -S -b -t ex1.fa.fai -o ex1.bam ex1.sam.gz

Build an index for the BAM file:
    samtools index ex1.bam

View a portion of the BAM file:
    samtools view ex1.bam seq2:450-550

Visually inspect the alignments at the same location:
    samtools tview -p seq2:450 ex1.bam ex1.fa

View the data in pileup format:
    samtools mpileup -f ex1.fa ex1.bam

Generate an uncompressed VCF file of variants:
    samtools mpileup -vu -f ex1.fa ex1.bam > ex1.vcf

Generate a compressed VCF file of variants:
    samtools mpileup -g -f ex1.fa ex1.bam > ex1.bcf
