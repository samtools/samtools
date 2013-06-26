#!/bin/sh

gcc -g -O0 -o test_trans_tbl_init test_trans_tbl_init.c ../../libbam.a -lz
gcc -g -O0 -o test_bam_translate test_bam_translate.c ../../libbam.a -lz
