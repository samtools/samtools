#include "sam.h"
#include <time.h>

void bam_sample_core(samfile_t *in, samfile_t *out, float r)
{
       bam1_t *b;
       b = bam_init1();
       while (samread(in, b) >= 0)
       {
               // generate a random float
               if ((float)drand48() <= r)
                       samwrite(out, b);
       }
}


int bam_sample(int argc, char *argv[])
{
       int c = 0;
       int sample_seed = 0;
       float sample_rate = 1.0;
       samfile_t *in, *out;
       while (( c = getopt(argc, argv, "r:s:")) >= 0) {
               switch (c) {
                       case 'r': sample_rate = atof(optarg); break;
                       case 's': sample_seed = atoi(optarg); break;
                       default: fprintf(stderr, "Unrecognizd option '-%c'.\n", c); return 1;
               }
       }

       if (optind + 2 > argc) {
               fprintf(stderr, "\n");
               fprintf(stderr, "Usage:  samtools sample [-r sample rate] [-s seed] <input.srt.bam> <output.bam>\n\n");
               fprintf(stderr, "Option: -r    sample rate (default: 1.0)\n");
               fprintf(stderr, "        -s    seed for random number generator\n");
               return 1;
       }
       // some checks
       if (sample_rate < 0) sample_rate = 0.0;
       if (sample_rate > 1) sample_rate = 1.0;
       if (sample_seed < 0) sample_seed *= -1;

       if (sample_seed) {
               // a seed has been specified, initialize the seed
               srand48(sample_seed);
       } else {
               // otherwise initialize random
               srand48((unsigned)time(0)+(unsigned)getpid());
       }
       in = samopen(argv[optind], "rb", 0);
       out = samopen(argv[optind + 1], "wb", in->header);
       if (in == 0 || out == 0) {
               fprintf(stderr, "[bam_sample] fail to read/write input files\n");
               return 1;
       }
       bam_sample_core(in, out, sample_rate);
       samclose(in);
       samclose(out);
       return 0;
}
