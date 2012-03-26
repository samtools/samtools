#include <stdio.h>
#include "radix.h"
#include "sam.h"

typedef struct
{
  int printAll,doMedian,maxCoverage;
} Options;

/**
 * Check if read is properly mapped
 * @return true if read mapped, false otherwise
 */
static inline int is_mapped(const bam1_core_t *core)
{
  return !(core->flag&BAM_FUNMAP);
}

/**
 * Print usage instructions
 */
static int print_usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:   samtools qa [options] <in.bam>\n");
  fprintf(stderr, "Options: -a            Don't print alternate assemblies to the output file (for human genome)\n");
  fprintf(stderr, "         -m            Also compute median coverage\n");
  fprintf(stderr, "         -c [INT]      Maximum coverage to consider in histogram [30]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Note: Input file should be sorted\n\n");
  return 1;
}

static void compute_print_cov(FILE* outputFile, Options userOpt, int* data, char* name,const uint32_t chrSize, int64_t* coverageHist,const int currentTid)
{
  int32_t covVal = 0;
  int64_t covSum = 0;
  int32_t i;

  //Go through chromosome and count avarage covarage.
  for (i=0; i<chrSize; ++i){
    covVal += data[i];
    //This will be sorted later.
    //If -m was not defined, this is useless, but cheaper than an 'if'
    data[i] = covVal;
    covSum += covVal;
    //Add value to histogram
    if (covVal > userOpt.maxCoverage) {
      ++coverageHist[userOpt.maxCoverage];
    } else {
      ++coverageHist[covVal];
    }

  }
  if (userOpt.doMedian)
    //Sort entireChr
    radix_sort(data, chrSize);

  //Printout avarage coverage over this chrom
  if (userOpt.printAll == 1) {
    if (userOpt.doMedian)
      fprintf(outputFile, "%s\t%d\t%3.5f\t%d\n", name, chrSize, (double)covSum / chrSize, data[chrSize/2]);
    else
      fprintf(outputFile, "%s\t%d\t%3.5f\n", name, chrSize, (double)covSum / chrSize);
  } else if (currentTid < 24) {
    //Don't print alternate assemblies to the file
    //This is human genome specific
    if (userOpt.doMedian)
      fprintf(outputFile, "%s\t%d\t%3.5f\t%d\n", name, chrSize, (double)covSum / chrSize, data[chrSize/2]);
    else
      fprintf(outputFile, "%s\t%d\t%3.5f\n", name, chrSize, (double)covSum / chrSize);
  }
}

/**
 * Main of app
 */
int main_qa(int argc, char *argv[])
{
  samfile_t *fp;
  FILE *outputFile;
  Options userOpt;
  userOpt.printAll = 1;
  userOpt.doMedian = 0;
  userOpt.maxCoverage = 30;
  int arg;
  //Get args
  while ((arg = getopt(argc, argv, "amc:")) >= 0) {
    switch (arg) {
    case 'a': userOpt.printAll = 0; break;
    case 'm': userOpt.doMedian = 1; break;
    case 'c': userOpt.maxCoverage = atoi(optarg); break;
    }
  }

  if (argc-optind != 1) {
    print_usage();
    return 1;
  }

  //Note that file is supposed to have been ordered beforehand!
  if ((fp = samopen(argv[optind], "rb", 0)) == 0) {
    fprintf(stderr, "qaCompute: Fail to open BAM file %s\n", argv[1]);
    return 1;
  }
  outputFile = stdout;


    //Initialize bam entity
    bam1_t *b = bam_init1();

    //All var declarations
    int64_t totalGenomeLength = 0;
    int32_t unmappedReads = 0;
    int32_t zeroQualityReads = 0;
    int32_t totalNumberOfReads = 0;
    int32_t totalProperPaires = 0;
    uint32_t chrSize = 0;

    int32_t duplicates = 0;

    int *entireChr = NULL;
    //Keep header for further reference
    bam_header_t* head = fp->header;

    int32_t currentTid = -1;

    //Create "map" vector for histogram
    int64_t* coverageHist = (int64_t*)malloc((userOpt.maxCoverage+1)*sizeof(int64_t));
    memset( coverageHist, 0, (userOpt.maxCoverage+1)*sizeof(int64_t));

    //Write file table header
    if (userOpt.doMedian == 1)
      fprintf(outputFile, "Chromosome\tSeq_len\tAvg_Cov\tMedian_Cov\n");
    else
      fprintf(outputFile, "Chromosome\tSeq_lem\tAvg_Cov\n");

    while (samread(fp, b) >= 0) {

      //uint32_t* cigar = bam1_cigar(b);

      //Get bam core.
      const bam1_core_t *core = &b->core;

      if (core == NULL) {
    //There is something wrong with the read/file
    //Leak everything and exit!
    return -1;
      }

      //BAM block has been read
      if (!is_mapped(core))
    ++unmappedReads;
      else {

    if (core->tid != currentTid) {

      //Count coverage!
      if (currentTid != -1) {
        compute_print_cov(outputFile, userOpt, entireChr, head->target_name[currentTid], chrSize, coverageHist, currentTid);
      }

      //Get length of next section
          chrSize = head->target_len[core->tid];
          totalGenomeLength += chrSize;

      //Done with current section.
      //Allocate memory
      entireChr = (int*)realloc(entireChr, (chrSize+1)*sizeof(int));

      if (entireChr == NULL) {
        return -1;
      }
      memset(entireChr, 0, (chrSize+1)*sizeof(int));

      currentTid = core->tid;

    }

    //If read has quality == 0, we won't count it as mapped
    if (core->qual != 0) {
     if (core->flag&BAM_FPROPER_PAIR) {
        //Is part of a proper pair
        ++totalProperPaires;
      }

     if (core->flag&BAM_FDUP) {
       //This is a duplicate. Don't count it!.
       ++duplicates;
     } else {
       //All entries in SAM file are represented on the forward strand! (See specs of SAM format for details)
       ++entireChr[core->pos];

       if (core->pos+core->l_qseq >= chrSize)
         --entireChr[chrSize-1];
       else
         --entireChr[core->pos+core->l_qseq];
     }

    } else {
      //Count is as unmapped?
      ++zeroQualityReads;
    }
      }

      ++totalNumberOfReads;

    }

    //Compute coverage for the last "chromosome"
    compute_print_cov(outputFile, userOpt, entireChr, head->target_name[currentTid], chrSize, coverageHist, currentTid);

    bam_destroy1(b);
    free(entireChr);

    //Print header for next table in output file
    fprintf(outputFile,"\nCov*X\tPercentage\tNr. of bases\n");

    //Compute procentages of genome cover!
    int i = 0;
    for (; i <= userOpt.maxCoverage; ++i) {
      if (i == 0) {
        //Non-covered!
      } else {
        int64_t coverage = 0;
        //All that has been covered i, had been covered i+1, i+2 and so on times. Thus, do this addition
        int x = i;
        for (; x <= userOpt.maxCoverage; ++x)
            coverage += coverageHist[x];
        fprintf(outputFile,"%d\t%3.5f\t%ld\n",i, (double)(100*coverage)/totalGenomeLength, (long)coverageHist[i]);
      }
    }

    fprintf(outputFile,"\nOther\n");

    //Printout procentage of mapped/unmapped reads
    double procentageOfUnmapped = (100*unmappedReads)/totalNumberOfReads;
    double procentageOfZeroQuality = (100*zeroQualityReads)/totalNumberOfReads;
    fprintf(outputFile,"Total number of reads: %d\n", totalNumberOfReads);
    fprintf(outputFile,"Total number of duplicates found and ignored: %d\n", duplicates);
    fprintf(outputFile,"Percentage of unmapped reads: %3.5f\n", procentageOfUnmapped);
    fprintf(outputFile,"Percentage of zero quality mappings: %3.5f\n", procentageOfZeroQuality);
    int32_t nrOfPaires = totalNumberOfReads/2;
    double procOfProperPaires = (double)(100*(double)totalProperPaires/2)/nrOfPaires;
    fprintf(outputFile,"Number of proper paired reads: %d\n", totalProperPaires);
    fprintf(outputFile,"Percentage of proper pairs: %3.5f\n", procOfProperPaires);

    free(coverageHist);


  samclose(fp);
  fclose(outputFile);
  return 0;
}
