#include "../../bam_sort.c"
#include <fcntl.h>
#include <poll.h>

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
		int i = 0;
		uint8_t* aux = bam1_aux(b);

		while (i < b->l_aux) {
			printf("%.2s:%c:",aux+i,*(aux+i+2));
			i += 2;
			switch (*(aux+i)) {
				case 'Z':
					while (*(aux+1+i) != '\0') { putc(*(aux+1+i), stdout); ++i; }
					break;
			}
			putc('\n',stdout);
			++i;++i;
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

void setup_test_2(bam1_t** b_in, trans_tbl_t* tbl) {
	bam1_t* b;
	
	b = bam_init1();
	trans_tbl_test_init(tbl, 4);
	
	tbl->tid_trans[0] = 5;
	tbl->tid_trans[1] = 6;
	tbl->tid_trans[2] = 7;
	tbl->tid_trans[3] = 8;
	int in_there = 0;
	khiter_t iter = kh_put(c2c, tbl->rg_trans, strdup("hello"), &in_there);
	kh_value(tbl->rg_trans, iter) = strdup("goodbye");
	
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
	size_t data_len = 10 + 4 + 5 + 10 + 9;
	b->data = (uint8_t*)malloc(data_len);
	memcpy(b->data,
		   "123456789\0" // q_name
		   "\x00\x00\x00\xA0" // cigar
		   "\x00\x00\x00\x00\x00" // qseq
		   "\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF" // qual
		   "RGZhello\0" // aux
		   , data_len
		   );
	b->l_aux = 9;
	b->m_data = b->data_len = data_len;
	
	*b_in = b;
}

void setup_test_3(bam1_t** b_in, trans_tbl_t* tbl) {
	bam1_t* b;
	
	b = bam_init1();
	trans_tbl_test_init(tbl, 4);
	
	tbl->tid_trans[0] = 5;
	tbl->tid_trans[1] = 6;
	tbl->tid_trans[2] = 7;
	tbl->tid_trans[3] = 8;
	int in_there = 0;
	khiter_t iter = kh_put(c2c, tbl->pg_trans, strdup("hello"), &in_there);
	kh_value(tbl->pg_trans,iter) = strdup("goodbye");

	
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
	size_t data_len = 10 + 4 + 5 + 10 + 9;
	b->data = (uint8_t*)malloc(data_len);
	memcpy(b->data,
		   "123456789\0" // q_name
		   "\x00\x00\x00\xA0" // cigar
		   "\x00\x00\x00\x00\x00" // qseq
		   "\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF" // qual
		   "PGZhello\0" // aux
		   , data_len
		   );
	b->l_aux = 9;
	b->m_data = b->data_len = data_len;
	
	*b_in = b;
}

void setup_test_4(bam1_t** b_in, trans_tbl_t* tbl) {
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
	size_t data_len = 10 + 4 + 5 + 10 + 9;
	b->data = (uint8_t*)malloc(data_len);
	memcpy(b->data,
		   "123456789\0" // q_name
		   "\x00\x00\x00\xA0" // cigar
		   "\x00\x00\x00\x00\x00" // qseq
		   "\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF" // qual
		   "RGZhello\0" // aux
		   , data_len
		   );
	b->l_aux = 9;
	b->m_data = b->data_len = data_len;
	
	*b_in = b;
}

void setup_test_5(bam1_t** b_in, trans_tbl_t* tbl) {
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
	size_t data_len = 10 + 4 + 5 + 10 + 9;
	b->data = (uint8_t*)malloc(data_len);
	memcpy(b->data,
		   "123456789\0" // q_name
		   "\x00\x00\x00\xA0" // cigar
		   "\x00\x00\x00\x00\x00" // qseq
		   "\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF" // qual
		   "PGZhello\0" // aux
		   , data_len
		   );
	b->l_aux = 9;
	b->m_data = b->data_len = data_len;
	
	*b_in = b;
}

void setup_test_6(bam1_t** b_in, trans_tbl_t* tbl) {
	bam1_t* b;
	
	b = bam_init1();
	trans_tbl_test_init(tbl, 4);
	
	tbl->tid_trans[0] = 5;
	tbl->tid_trans[1] = 6;
	tbl->tid_trans[2] = 7;
	tbl->tid_trans[3] = 8;
	int in_there = 0;
	khiter_t iter_rg = kh_put(c2c, tbl->rg_trans, strdup("hello"), &in_there);
	kh_value(tbl->rg_trans, iter_rg) = strdup("goodbye");
	khiter_t iter_pg = kh_put(c2c, tbl->pg_trans, strdup("quail"), &in_there);
	kh_value(tbl->pg_trans, iter_pg) = strdup("bird");

	
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
	size_t data_len = 10 + 4 + 5 + 10 + 18;
	b->data = (uint8_t*)malloc(data_len);
	memcpy(b->data,
		   "123456789\0" // q_name
		   "\x00\x00\x00\xA0" // cigar
		   "\x00\x00\x00\x00\x00" // qseq
		   "\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF" // qual
		   "RGZhello\0PGZquail\0" // aux
		   , data_len
		   );
	b->l_aux = 18;
	b->m_data = b->data_len = data_len;
	
	*b_in = b;
}


int main(int argc, char**argv)
{
	// test state
	const int NUM_TESTS = 6;
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

	bam1_t* b;

	// Setup stderr redirect
	size_t len;
	char* res;
	int restore_stderr = dup(STDERR_FILENO); // Save stderr
	char* template = strdup("test_bam_trans_XXXXXXX");
	int null_fd = mkstemp(template);
	unlink(template);
	FILE* check = fdopen(null_fd, "w+");
	
	// setup
	if (verbose) printf("BEGIN test 1\n");  // TID test
	trans_tbl_t tbl1;
	setup_test_1(&b,&tbl1);
	if (verbose) {
		printf("b\n");
		dump_read(b);
	}
	if (verbose) printf("RUN test 1\n");

	// test
	dup2(null_fd, STDERR_FILENO); // Redirect stderr to pipe
	bam_translate(b, &tbl1);
	dup2(restore_stderr, STDERR_FILENO);

	if (verbose) printf("END RUN test 1\n");
	if (verbose) {
		printf("b\n");
		dump_read(b);
	}

	// check result
	len = 0;
	rewind(check);
	res = fgetln(check,&len);
	if (feof(check) || (res && !strcmp("",res))) {
		++success;
	} else {
		++failure;
	}
	rewind(check);
	
	// teardown
	bam_destroy1(b);
	trans_tbl_destroy(&tbl1);
	if (verbose) printf("END test 1\n");
	
	// setup
	if (verbose) printf("BEGIN test 2\n");  // RG test
	trans_tbl_t tbl2;
	setup_test_2(&b,&tbl2);
	if (verbose) {
		printf("b\n");
		dump_read(b);
	}
	if (verbose) printf("RUN test 2\n");
	
	// test
	dup2(null_fd, STDERR_FILENO);
	bam_translate(b, &tbl2);
	dup2(restore_stderr, STDERR_FILENO);

	if (verbose) printf("END RUN test 2\n");
	if (verbose) {
		printf("b\n");
		dump_read(b);
	}

	// check result
	len = 0;
	rewind(check);
	res = fgetln(check,&len);
	if (feof(check) || (res && !strcmp("",res))) {
		++success;
	} else {
		++failure;
	}
	rewind(check);
	
	// teardown
	bam_destroy1(b);
	trans_tbl_destroy(&tbl2);
	if (verbose) printf("END test 2\n");
	
	if (verbose) printf("BEGIN test 3\n");  // PG test
	// setup
	trans_tbl_t tbl3;
	setup_test_3(&b,&tbl3);
	if (verbose) {
		printf("b\n");
		dump_read(b);
	}
	if (verbose) printf("RUN test 3\n");

	// test
	dup2(null_fd, 2);
	bam_translate(b, &tbl3);
	dup2(restore_stderr, 2);
	if (verbose) printf("END RUN test 3\n");
	if (verbose) {
		printf("b\n");
		dump_read(b);
	}

	// check result
	len = 0;
	rewind(check);
	res = fgetln(check,&len);
	if (feof(check) || (res && !strcmp("",res))) {
		++success;
	} else {
		++failure;
	}
	rewind(check);
	
	// teardown
	bam_destroy1(b);
	trans_tbl_destroy(&tbl3);
	if (verbose) printf("END test 3\n");
	
	if (verbose) printf("BEGIN test 4\n");  // RG test non existant
	// setup
	trans_tbl_t tbl4;
	setup_test_4(&b,&tbl4);
	if (verbose) {
		printf("b\n");
		dump_read(b);
	}
	if (verbose) printf("RUN test 4\n");
	
	// test
	dup2(null_fd, STDERR_FILENO);
	bam_translate(b, &tbl4);
	dup2(restore_stderr, STDERR_FILENO);
	
	if (verbose) printf("END RUN test 4\n");
	if (verbose) {
		printf("b\n");
		dump_read(b);
	}
	// check result
	len = 0;
	rewind(check);
	res = fgetln(check,&len);
	if (res && !strcmp("[bam_translate] RG tag \"hello\" on read \"123456789\" encountered with no corresponding entry in header, tag lost\n",res)) {
		++success;
	} else {
		++failure;
	}
	rewind(check);
	
	// teardown
	bam_destroy1(b);
	trans_tbl_destroy(&tbl4);
	if (verbose) printf("END test 4\n");
	
	if (verbose) printf("BEGIN test 5\n");  // PG test non existant
	// setup
	trans_tbl_t tbl5;
	setup_test_5(&b,&tbl5);
	if (verbose) {
		printf("b\n");
		dump_read(b);
		printf("RUN test 5\n");
	}
	// test
	dup2(null_fd, STDERR_FILENO);
	bam_translate(b, &tbl5);
	dup2(restore_stderr, STDERR_FILENO);

	if (verbose) printf("END RUN test 5\n");
	if (verbose) {
		printf("b\n");
		dump_read(b);
	}

	// check result
	len = 0;
	rewind(check);
	res = fgetln(check,&len);
	if (res && !strcmp("[bam_translate] PG tag \"hello\" on read \"123456789\" encountered with no corresponding entry in header, tag lost\n",res)) {
		++success;
	} else {
		++failure;
	}
	rewind(check);

	// teardown
	bam_destroy1(b);
	trans_tbl_destroy(&tbl5);
	if (verbose) printf("END test 5\n");
	
	if (verbose) printf("BEGIN test 6\n");  // RG and PG test
	// setup
	trans_tbl_t tbl6;
	setup_test_6(&b,&tbl6);
	if (verbose) {
		printf("b\n");
		dump_read(b);
	}
	if (verbose) printf("RUN test 6\n");

	// test
	dup2(null_fd, STDERR_FILENO);
	bam_translate(b, &tbl6);
	dup2(restore_stderr, STDERR_FILENO);

	if (verbose) printf("END RUN test 6\n");
	if (verbose) {
		printf("b\n");
		dump_read(b);
	}
	
	// check result
	len = 0;
	rewind(check);
	res = fgetln(check,&len);
	if (res && !strcmp("[bam_translate] PG tag \"hello\" on read \"123456789\" encountered with no corresponding entry in header, tag lost\n",res)) {
		res = fgetln(check,&len);
		if (feof(check) || (res && !strcmp("",res))) {
			++success;
		} else {
			++failure;
		}
	} else {
		++failure;
	}
	rewind(check);

	// teardown
	bam_destroy1(b);
	trans_tbl_destroy(&tbl6);
	if (verbose) printf("END test 6\n");

	// Cleanup
	fclose(check);
	close(restore_stderr);
	
	if (NUM_TESTS == success)
		return 0;
	else
	{
		fprintf(stderr, "%d failures %d successes\n", failure, success);
		return 1;
	}
}
