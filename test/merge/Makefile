all:
test: test_trans_tbl_init test_bam_translate test_rtrans_build test_pretty_header
	./test_trans_tbl_init
	./test_bam_translate
	./test_rtrans_build
	./test_pretty_header

test_trans_tbl_init: test_trans_tbl_init.c
	gcc -g -O0 -o test_trans_tbl_init test_trans_tbl_init.c ../../libbam.a ../../../htslib/htslib/libhts.a -lz -I../../../htslib -pthread

test_bam_translate: test_bam_translate.c
	gcc -g -O0 -o test_bam_translate test_bam_translate.c ../../libbam.a ../../../htslib/htslib/libhts.a -lz -I../../../htslib -pthread

test_rtrans_build: test_rtrans_build.c
	gcc -g -O0 -o test_rtrans_build test_rtrans_build.c ../../libbam.a ../../../htslib/htslib/libhts.a -lz -I../../../htslib -pthread

test_pretty_header: test_pretty_header.c
	gcc -g -O0 -o test_pretty_header test_pretty_header.c ../../libbam.a ../../../htslib/htslib/libhts.a -lz -I../../../htslib -pthread

clean:
	rm -f test_trans_tbl_init test_bam_translate test_rtrans_build test_pretty_header
