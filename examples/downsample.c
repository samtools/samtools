#include <stdio.h>
#include "sam.h"

struct fetchdata {
    float sample;
    long int total, filtered;
};

// callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data)
{
        struct fetchdata *d = (struct fetchdata *) data;
        d->total++;

        // Validate b against criteria
        uint32_t flag = b->core.flag;

        // If read is not mapped or flaged as duplicate skip
        if(flag & BAM_FUNMAP || flag & BAM_FDUP) {
            return 0;
        }

        // Use only primary alignments
        if(flag & BAM_FSECONDARY) {
            return 0; // Don't use secondary alignments
        }

        d->filtered++;

        if( ( (float)rand() / (float) RAND_MAX ) > d->sample ) return 0;

        return 1;
}

int main(int argc, char *argv[])
{
        /* FIXME - stdin > stdout */
	if (argc != 3) {
		fprintf(stderr, "Usage: downsample <in.bam> <out.bam>\n");
		return 1;
	}

	float sample = 0.50f;

	samfile_t *in = samopen(argv[1], "rb", 0);
	if(in == 0) {
		fprintf(stderr, "Failed to open input BAM file %s\n", argv[1]);
		return 1;
	}

	samfile_t *out;
	/* Provide bam_header_t from source file as header for output */
	out = samopen(argv[2], "wbh", in->header);
	if(out == 0) {
		fprintf(stderr, "Failed to open output file %s\n", argv[1]);
		return 1;
	}

        /* NB: rand is seeded with constant - result should be reproducible */
        srand(42);

        struct fetchdata data = { sample, 0, 0 };

        bam1_t *b = bam_init1();
        unsigned long count = 0;
        while(samread(in, b) > 0) {
            if(fetch_func(b, &data) != 0) continue;

            // Write to out
            samwrite(out, b);
            count++;
        }
        bam_destroy1(b);

        fprintf(stderr, "Sampled %ld reads out of %ld total, %ld filtered [mapped, primary, non-duplicate]\n",count, data.total, data.filtered);

	samclose(in);
	samclose(out);
	return 0;
}
