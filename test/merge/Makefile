all:
test: test_trans_tbl_init test_bam_translate
	./test_trans_tbl_init
	./test_bam_translate

test_trans_tbl_init:
	gcc -g -O0 -o test_trans_tbl_init test_trans_tbl_init.c ../../libbam.a -lz

test_bam_translate:
	gcc -g -O0 -o test_bam_translate test_bam_translate.c ../../libbam.a -lz

clean:
	rm -f test_trans_tbl_init test_bam_translate
