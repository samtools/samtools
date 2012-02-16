This is *NOT* the official repository of SAMtools.  The official SAMtools repository can be found at: http://samtools.sourceforge.net/ 

Major fixes:
- added "samtools sample" command, to sample reads from a SAM/BAM file at a given frequency.
- added "samtools qa" command, to compute the mean and median coverage, as well a histogram from 1 to N (defined by param) containing the number of bases covered a maximum of 1X, 2X...NX. Furthermore, "other" information is also available in the output file, namely:
		1. Total number of reads
		2. Total number of duplicates found and ignored (duplicates are "found" based on the sam flag and are ignored in the counting of coverage)
		3. Percentage of unmapped reads
		4. Percentage of zero quality mappings
		5. Number of proper paired reads (based on sam flag of proper_pair)
		6. Percentage of proper pairs.e

Minor fixes:
- Check the write filehandle after opening for write.
- allow for user defined [lowercase] tags in header elements.
- allow the maximum memory for "samtools sort" to be specified with units.
- adjust for leading hard clip on colorspace reads.
- catches and reports an invalid BAM header, instead of segfaulting later on.
- fixes a small underflow/overflow bug in integer parsing.
- checks for a lowerbound in text entry box to avoid segfault in tview.
