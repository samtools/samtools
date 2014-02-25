#include "../../bam_sort.c"

void setup_test_1(char** input) {
	*input = strdup(
	"@HD\n"
	"@SQ\n"
	"@RG\n"
	"@PG\n"
	"@CO\n");
}

bool check_test_1(char* input) {
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
				++verbose;
				break;
			default:
				break;
		}
	}

	if (verbose) printf("BEGIN test 1\n");
	// setup
	char* input;
	setup_test_1(&input);
	// test
	if (verbose > 1) {
		printf("input:\n%s",input);
	}
	if (verbose) printf("RUN test 1\n");
	pretty_header(&input, strlen(input));
	if (verbose) printf("END RUN test 1\n");
	if (verbose > 1) {
		printf("input:\n%s",input);
	}
	if (check_test_1(input)) { ++success; } else { ++failure; }
	// teardown
	free(input);
	if (verbose) printf("END test 1\n");
	
	if (success == NUM_TESTS) {
		return 0;
	} else {
		fprintf(stderr, "%d failures %d successes\n", failure, success);
		return 1;
	}
}
