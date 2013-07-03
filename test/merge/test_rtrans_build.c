#include "../../bam_sort.c"

void dump_rtrans(int* rtrans, int n, int n_targets) {
	printf("->n_targets:(%d)\n", n_targets);
	int i, j;
	for (i = 0; i < n_targets; ++i) {
		for (j = 0; j < n; ++j)
			fprintf(stderr, "%x",i,rtrans[i*n_targets+n]);
		fprintf(stderr, "\n");
	}
}

void setup_test_1(trans_tbl_t* tbl) {
}

bool check_test_1(trans_tbl_t* tbl, int* rtrans) {
	// Check input is unchanged
	
	// Check output
	
	return true;
}


int main(int argc, char**argv)
{
	const int NUM_TESTS = 1;
	int verbose = 0;
	int success = 0;
	int failure = 0;
	int getopt_char;
	while ((getopt_char = getopt(argc, argv, "v")) != -1) {
		switch (getopt_char) {
			case 'v':
				verbose = 1;
				break;
			default:
				break;
		}
	}
	const long GIMIC_SEED = 0x1234abcd330e;
	srand48(GIMIC_SEED);

	if (verbose) printf("BEGIN test 1\n");
	// setup
	trans_tbl_t tbl_1;
	int* rtrans_1 = NULL;
	setup_test_1(&tbl_1);
	// test
	if (verbose) {
		printf("rtrans\n");
		dump_rtrans(rtrans_1, 1, 1);
	}
	if (verbose) printf("RUN test 1\n");
	rtrans_1 = rtrans_build(1, 1, NULL);
	if (verbose) printf("END RUN test 1\n");
	if (verbose) {
		printf("rtrans\n");
		dump_rtrans(rtrans_1, 1, 1);
	}
	if (check_test_1(&tbl_1, rtrans_1)) { ++success; } else { ++failure; }
	// teardown
	trans_tbl_destroy(&tbl_1);
	free(rtrans_1);
	if (verbose) printf("END test 1\n");
	
	if (success == NUM_TESTS) {
		return 0;
	} else {
		fprintf(stderr, "%d failures %d successes\n", failure, success);
		return 1;
	}
}
