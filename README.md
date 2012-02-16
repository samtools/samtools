## This is *NOT* the official repository of SAMtools.  
The official SAMtools repository can be found at: http://samtools.sourceforge.net/ 

## Major fixes:
    - added "samtools qa" command, to compute the mean and median coverage, as well a histogram 
    from 1 to N (defined by param) containing the number of bases covered a maximum of 1X, 2X...NX. 
    Furthermore, "other" information is also available in the output file, namely:
        - Total number of reads
        - Total number of duplicates found and ignored (duplicates are "found" based on the sam flag 
        and are ignored in the counting of coverage)
        - Percentage of unmapped reads
        - Percentage of zero quality mappings
        - Number of proper paired reads (based on sam flag of proper pairs)
        - Percentage of proper pairs.e

## Minor fixes:
    - Check the write filehandle after opening for write.
    - allow for user defined [lowercase] tags in header elements.
    - allow the maximum memory for "samtools sort" to be specified with units.
    - adjust for leading hard clip on colorspace reads.
    - catches and reports an invalid BAM header, instead of segfaulting later on.
    - fixes a small underflow/overflow bug in integer parsing.
    - checks for a lowerbound in text entry box to avoid segfault in tview.
