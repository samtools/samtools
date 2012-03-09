File ex1.fa contains two sequences cut from the human genome
build36. They were exatracted with command:

  samtools faidx human_b36.fa 2:2043966-2045540 20:67967-69550

Sequence names were changed manually for simplicity. File ex1.sam.gz
contains MAQ alignments exatracted with:

  (samtools view NA18507_maq.bam 2:2044001-2045500;
   samtools view NA18507_maq.bam 20:68001-69500)

and processed with `samtools fixmate' to make it self-consistent as a
standalone alignment.

To try samtools, you may run the following commands:

  samtools faidx ex1.fa                 # index the reference FASTA
  samtools import ex1.fa.fai ex1.sam.gz ex1.bam   # SAM->BAM
  samtools index ex1.bam                # index BAM
  samtools tview ex1.bam ex1.fa         # view alignment
  samtools pileup -cf ex1.fa ex1.bam    # pileup and consensus
  samtools pileup -cf ex1.fa -t ex1.fa.fai ex1.sam.gz

