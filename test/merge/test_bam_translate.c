#include "../../bam_sort.c"

void dump_read(bam1_t* b) {
	printf("->core.tid:(%d)\n", b->core.tid);
	printf("->core.pos:(%d)\n", b->core.pos);
	printf("->core.bin:(%d)\n", b->core.bin);
	printf("->core.qual:(%d)\n", b->core.qual);
	printf("->core.l_qname:(%d)\n", b->core.l_qname);
	printf("->core.flag:(%d)\n", b->core.flag);
	printf("->core.n_cigar:(%d)\n", b->core.n_cigar);
	printf("->core.l_qseq:(%d)\n", b->core.l_qseq);
	printf("->core.mtid:(%d)\n", b->core.mtid);
	printf("->core.mpos:(%d)\n", b->core.mpos);
	printf("->core.isize:(%d)\n", b->core.isize);
	if (b->data) {
		printf("->data:");
		int i;
		for (i = 0; i < b->data_len; ++i) {
			printf("%x ", b->data[i]);
		}
		printf("\n");
	}
	if (b->core.l_qname) {
		printf("qname: %s\n",bam1_qname(b));
	}
	if (b->core.l_qseq) {
		printf("qseq:");
		int i;
		for (i = 0; i < b->core.l_qseq; ++i) {
			printf("%c",bam_nt16_rev_table[bam_nt16_nt4_table[bam1_seqi(bam1_seq(b),i)]]);
		}
		printf("\n");
		printf("qual:");
		for (i = 0; i < b->core.l_qseq; ++i) {
			printf("%c",bam1_qual(b)[i]);
		}
		printf("\n");

	}

	if (b->l_aux) {
		int i;
		for (i = 0; i < b->l_aux; ++i) {
			bam1_aux(b);
		}
	}
	printf("\n");
}

void trans_tbl_test_init(trans_tbl_t* tbl, int32_t n_targets)
{
	tbl->tid_trans = (int*)calloc(n_targets, sizeof(int32_t));
	tbl->rg_trans = kh_init(c2c);
	tbl->pg_trans = kh_init(c2c);
}

void setup_test_1(bam1_t** b_in, trans_tbl_t* tbl) {
	bam1_t* b;

	b = bam_init1();
	trans_tbl_test_init(tbl, 4);
	
	tbl->tid_trans[0] = 5;
	tbl->tid_trans[1] = 6;
	tbl->tid_trans[2] = 7;
	tbl->tid_trans[3] = 8;
	
	
	b->core.tid = 0;
	b->core.pos = 1334;
	b->core.bin = 0;
	b->core.qual = 10;
	b->core.l_qname = 10;
	b->core.flag = 0;
	b->core.n_cigar = 1;
	b->core.l_qseq = 10;
	b->core.mtid = -1;
	b->core.mpos = 0;
	b->core.isize = -1;
	size_t data_len = 10 + 4 + 5 + 10 + 0;
	b->data = (uint8_t*)malloc(data_len);
	memcpy(b->data,
					 "123456789\0" // q_name
					 "\x00\x00\x00\xA0" // cigar
					 "\x00\x00\x00\x00\x00" // qseq
					 "\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF" // qual
					 "" // aux
		   , data_len
					 );
	b->l_aux = 0;
	b->m_data = b->data_len = data_len;
	
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
