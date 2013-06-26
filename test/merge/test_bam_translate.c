#include "../../bam_sort.c"

void dump_read(bam1_t* b) {
	printf("->core.tid:(%d)\n", b->core.tid);
	printf("->core.pos:(%d)\n", b->core.pos);
	printf("->core.bin:(%d)\n", b->core.bin);
	printf("->core.qual:(%d)\n", b->core.qual);
	printf("->core.l_qseq:(%d)\n", b->core.l_qseq);
	printf("->core.flag:(%d)\n", b->core.flag);
	printf("->core.n_cigar:(%d)\n", b->core.n_cigar);
	printf("->core.l_qseq:(%d)\n", b->core.l_qseq);
	printf("->core.mtid:(%d)\n", b->core.mtid);
	printf("->core.mpos:(%d)\n", b->core.mpos);
	printf("->core.isize:(%d)\n", b->core.isize);
	printf("->data:");
	if (b->data) {
		int i;
		for (i = 0; i < b->data_len; ++i) {
			printf("%x", i, b->data[i]);
		}
	}
	printf("\n");
}

void setup_test_1(bam1_t** b_in, trans_tbl_t* tbl) {
	bam1_t* b;

	b = bam_init1();
	
	tbl->tid_trans = (int*)calloc(1, sizeof(int32_t));
	tbl->rg_trans = kh_init(c2c);
	tbl->pg_trans = kh_init(c2c);
	
	*b_in = b;
}

int main(int argc, char**argv)
{
	bam1_t* b;

	printf("BEGIN test 1\n");
	// setup
	trans_tbl_t tbl;
	setup_test_1(&b,&tbl);
	// test
	printf("b\n");
	dump_read(b);
	printf("RUN test 1\n");
	bam_translate(b, &tbl);
	printf("END RUN test 1\n");
	printf("b\n");
	dump_read(b);
	// teardown
	bam_destroy1(b);
	trans_tbl_destroy(&tbl);
	printf("END test 1\n");

	return 0;
}
