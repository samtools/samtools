cmake_minimum_required(VERSION 3.6)
project(samtools)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        bam.c
        bam.h
        bam2bcf.c
        bam2bcf.h
        bam2bcf_indel.c
        bam2depth.c
        bam_addrprg.c
        bam_aux.c
        bam_cat.c
        bam_color.c
        bam_endian.h
        bam_flags.c
        bam_import.c
        bam_index.c
        bam_lpileup.c
        bam_lpileup.h
        bam_mate.c
        bam_md.c
        bam_plbuf.c
        bam_plbuf.h
        bam_plcmd.c
        bam_quickcheck.c
        bam_reheader.c
        bam_rmdup.c
        bam_rmdupse.c
        bam_sort.c
        bam_split.c
        bam_stat.c
        bam_tview.c
        bam_tview.h
        bam_tview_curses.c
        bam_tview_html.c
        bamshuf.c
        bamtk.c
        bedcov.c
        bedidx.c
        cut_target.c
        dict.c
        faidx.c
        padding.c
        phase.c
        sam.c
        sam.h
        sam_header.c
        sam_header.h
        sam_opts.c
        sam_opts.h
        sam_utils.c
        sam_view.c
        sample.c
        sample.h
        samtools.h
        stats.c
        stats_isize.c
        stats_isize.h)

include_directories("../htslib")
link_directories("../htslib")
add_executable(samtools ${SOURCE_FILES})
target_link_libraries(samtools hts z)