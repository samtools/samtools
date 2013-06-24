#include "../../bam_sort.c"

void dump_header(bam_header_t* hdr) {
	printf("->n_targets:(%d)\n", hdr->n_targets);
	int i;
	for (i = 0; i < hdr->n_targets; ++i) {
		printf("->target_name[%d]:(%s)\n",i,hdr->target_name[i]);
		printf("->target_len[%d]:(%d)\n",i,hdr->target_len[i]);
	}
	printf("->text:(%s)\n", hdr->text);
}

void setup_test_1(bam_header_t** translate_in, bam_header_t** out_in) {
	bam_header_t* out;
	bam_header_t* translate;

	translate = bam_header_init();
	translate->text = strdup(
							 "@HD\tVN:1.4\tSO:unknown\n"
							 "@SQ\tID:fish\tLN:133\t"
							 );
	translate->n_targets = 1;
	translate->target_name = (char**)calloc(1, sizeof(char*));
	translate->target_len = (uint32_t*)calloc(1, sizeof(uint32_t));
	translate->target_name[0] = strdup("fish");
	translate->target_len[0] = 133;
	out = bam_header_init();
	out->text = strdup(
					   "@HD\tVN:1.4\tSO:unknown\n"
					   "@SQ\tID:fish\tLN:133\tSP:frog"
					   );
	out->n_targets = 1;
	out->target_name = (char**)calloc(1, sizeof(char*));
	out->target_len = (uint32_t*)calloc(1, sizeof(uint32_t));
	out->target_name[0] = strdup("fish");
	out->target_len[0] = 133;
	
	*translate_in = translate;
	*out_in = out;
}

void setup_test_2(bam_header_t** translate_in, bam_header_t** out_in) {
	bam_header_t* out;
	bam_header_t* translate;

	translate = bam_header_init();
	translate->text = strdup(
							 "@HD\tVN:1.4\tSO:unknown\n"
							 "@SQ\tID:donkey\tLN:133\n"
							 "@SQ\tID:fish\tLN:133"
							 );
	translate->n_targets = 2;
	translate->target_name = (char**)calloc(translate->n_targets, sizeof(char*));
	translate->target_len = (uint32_t*)calloc(translate->n_targets, sizeof(uint32_t));
	translate->target_name[0] = strdup("donkey");
	translate->target_len[0] = 133;
	translate->target_name[1] = strdup("fish");
	translate->target_len[1] = 133;
	out = bam_header_init();
	out->text = strdup(
					   "@HD\tVN:1.4\tSO:unknown\n"
					   "@SQ\tID:fish\tLN:133\tSP:frog"
					   );
	out->n_targets = 1;
	out->target_name = (char**)calloc(1, sizeof(char*));
	out->target_len = (uint32_t*)calloc(1, sizeof(uint32_t));
	out->target_name[0] = strdup("fish");
	out->target_len[0] = 133;
	
	*translate_in = translate;
	*out_in = out;
}

void setup_test_3(bam_header_t** translate_in, bam_header_t** out_in) {
	bam_header_t* out;
	bam_header_t* translate;
	
	translate = bam_header_init();
	translate->text = strdup(
							 "@HD\tVN:1.4\tSO:unknown\n"
							 "@SQ\tID:donkey\tLN:133\n"
							 "@SQ\tID:fish\tLN:133\n"
							 "@RG\tID:fish\tPU:trans\n"
							 );
	translate->n_targets = 2;
	translate->target_name = (char**)calloc(translate->n_targets, sizeof(char*));
	translate->target_len = (uint32_t*)calloc(translate->n_targets, sizeof(uint32_t));
	translate->target_name[0] = strdup("donkey");
	translate->target_len[0] = 133;
	translate->target_name[1] = strdup("fish");
	translate->target_len[1] = 133;
	out = bam_header_init();
	out->text = strdup(
					   "@HD\tVN:1.4\tSO:unknown\n"
					   "@SQ\tID:fish\tLN:133\tSP:frog"
					   );
	out->n_targets = 1;
	out->target_name = (char**)calloc(1, sizeof(char*));
	out->target_len = (uint32_t*)calloc(1, sizeof(uint32_t));
	out->target_name[0] = strdup("fish");
	out->target_len[0] = 133;
	
	*translate_in = translate;
	*out_in = out;
}

void setup_test_4(bam_header_t** translate_in, bam_header_t** out_in) {
	bam_header_t* out;
	bam_header_t* translate;
	
	translate = bam_header_init();
	translate->text = strdup(
							 "@HD\tVN:1.4\tSO:unknown\n"
							 "@SQ\tID:donkey\tLN:133\n"
							 "@SQ\tID:fish\tLN:133\n"
							 "@RG\tID:fish\tPU:trans\n"
							 );
	translate->n_targets = 2;
	translate->target_name = (char**)calloc(translate->n_targets, sizeof(char*));
	translate->target_len = (uint32_t*)calloc(translate->n_targets, sizeof(uint32_t));
	translate->target_name[0] = strdup("donkey");
	translate->target_len[0] = 133;
	translate->target_name[1] = strdup("fish");
	translate->target_len[1] = 133;
	out = bam_header_init();
	out->text = strdup(
					   "@HD\tVN:1.4\tSO:unknown\n"
					   "@SQ\tID:fish\tLN:133\tSP:frog\n"
					   "@RG\tID:fish\tPU:out\n"
					   );
	out->n_targets = 1;
	out->target_name = (char**)calloc(1, sizeof(char*));
	out->target_len = (uint32_t*)calloc(1, sizeof(uint32_t));
	out->target_name[0] = strdup("fish");
	out->target_len[0] = 133;
	
	*translate_in = translate;
	*out_in = out;
}


void setup_test_5(bam_header_t** translate_in, bam_header_t** out_in) {
	bam_header_t* out;
	bam_header_t* translate;
	
	translate = bam_header_init();
	translate->text = strdup(
							 "@HD\tVN:1.4\tSO:unknown\n"
							 "@SQ\tID:donkey\tLN:133\n"
							 "@SQ\tID:fish\tLN:133\n"
							 "@RG\tID:fish\tPU:trans\n"
							 "@PG\tXX:dummy\tID:fish\tDS:trans\n"
							 "@PG\tPP:fish\tID:hook\tDS:trans\n"
							 );
	translate->n_targets = 2;
	translate->target_name = (char**)calloc(translate->n_targets, sizeof(char*));
	translate->target_len = (uint32_t*)calloc(translate->n_targets, sizeof(uint32_t));
	translate->target_name[0] = strdup("donkey");
	translate->target_len[0] = 133;
	translate->target_name[1] = strdup("fish");
	translate->target_len[1] = 133;
	out = bam_header_init();
	out->text = strdup(
					   "@HD\tVN:1.4\tSO:unknown\n"
					   "@SQ\tID:fish\tLN:133\tSP:frog\n"
					   "@RG\tID:fish\tPU:out\n"
					   "@PG\tXX:dummyx\tID:fish\tDS:out\n"
					   "@PG\tPP:fish\tID:hook\tDS:out\n"
					   );
	out->n_targets = 1;
	out->target_name = (char**)calloc(1, sizeof(char*));
	out->target_len = (uint32_t*)calloc(1, sizeof(uint32_t));
	out->target_name[0] = strdup("fish");
	out->target_len[0] = 133;
	
	*translate_in = translate;
	*out_in = out;
}

int main(int argc, char**argv)
{
	bam_header_t* out;
	bam_header_t* translate;

	printf("BEGIN test 1\n");
	// setup
	trans_tbl_t tbl;
	setup_test_1(&translate,&out);
	// test
	printf("translate\n");
	dump_header(translate);
	printf("out\n");
	dump_header(out);
	printf("RUN test 1\n");
	trans_tbl_init(out, translate, &tbl);
	printf("END RUN test 1\n");
	printf("translate\n");
	dump_header(translate);
	printf("out\n");
	dump_header(out);
	// teardown
	bam_header_destroy(translate);
	bam_header_destroy(out);
	trans_tbl_destroy(&tbl);
	printf("END test 1\n");

	// test
	printf("BEGIN test 2\n");
	// reinit
	trans_tbl_t tbl_2;
	setup_test_2(&translate,&out);
	printf("translate\n");
	dump_header(translate);
	printf("out\n");
	dump_header(out);
	printf("RUN test 2\n");
	trans_tbl_init(out, translate, &tbl_2);
	printf("END RUN test 2\n");
	printf("translate\n");
	dump_header(translate);
	printf("out\n");
	dump_header(out);
	// teardown
	bam_header_destroy(translate);
	bam_header_destroy(out);
	trans_tbl_destroy(&tbl_2);
	printf("END test 2\n");

	// test
	printf("BEGIN test 3\n");
	// reinit
	trans_tbl_t tbl_3;
	setup_test_3(&translate,&out);
	printf("translate\n");
	dump_header(translate);
	printf("out\n");
	dump_header(out);
	printf("RUN test 3\n");
	trans_tbl_init(out, translate, &tbl_3);
	printf("END RUN test 3\n");
	printf("translate\n");
	dump_header(translate);
	printf("out\n");
	dump_header(out);
	// teardown
	bam_header_destroy(translate);
	bam_header_destroy(out);
	trans_tbl_destroy(&tbl_3);
	printf("END test 3\n");

	// test
	printf("BEGIN test 4\n");
	// reinit
	trans_tbl_t tbl_4;
	setup_test_4(&translate,&out);
	printf("translate\n");
	dump_header(translate);
	printf("out\n");
	dump_header(out);
	printf("RUN test 4\n");
	trans_tbl_init(out, translate, &tbl_4);
	printf("END RUN test 4\n");
	printf("translate\n");
	dump_header(translate);
	printf("out\n");
	dump_header(out);
	// teardown
	bam_header_destroy(translate);
	bam_header_destroy(out);
	trans_tbl_destroy(&tbl_4);
	printf("END test 4\n");

	// test
	printf("BEGIN test 5\n");
	// reinit
	trans_tbl_t tbl_5;
	setup_test_5(&translate,&out);
	printf("translate\n");
	dump_header(translate);
	printf("out\n");
	dump_header(out);
	printf("RUN test 5\n");
	trans_tbl_init(out, translate, &tbl_5);
	printf("END RUN test 5\n");
	printf("translate\n");
	dump_header(translate);
	printf("out\n");
	dump_header(out);
	// teardown
	bam_header_destroy(translate);
	bam_header_destroy(out);
	trans_tbl_destroy(&tbl_5);
	printf("END test 5\n");

	return 0;
}
