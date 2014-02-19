#include "../../bam_tview.c"
#include <stdbool.h>

const char header_1[] =
"@HD	VN:1.4	SO:undefined\n"
"@SQ	SN:dummy\n"
"@RG	ID:blah	SM:po\n";

void setup_test_1(char** header)
{
	*header = strndup(header_1, sizeof(header_1));
}

khash_t(kh_rg)* run_test_1(char* header)
{
	khash_t(kh_rg)* test_result = get_rg_sample(header,"po");
	return test_result;
}

bool check_test_1(khash_t(kh_rg)* test_result, char* header)
{
	if (strcmp(header_1, header)) return false;
	// test blah is in there
	if (kh_get(kh_rg, test_result, "blah") == kh_end(test_result))
	{
		return false;
	}
	return true;
}

void teardown_1(khash_t(kh_rg)* test_result, char* header)
{
	free(header);
}

int main(int argc, char** argv)
{
	const int NUM_TESTS = 1;
	int success = 0;
	int failure = 0;

	char* test_header_1;
	setup_test_1(&test_header_1);
	khash_t(kh_rg)* test_result_1 = run_test_1(test_header_1);
	if (!check_test_1(test_result_1, test_header_1))
		failure++;
	else
		success++;
	teardown_1(test_result_1, test_header_1);

	if (success == NUM_TESTS) {
		return 0;
	} else {
		fprintf(stderr, "%d failures %d successes\n", failure, success);
		return 1;
	}
}
