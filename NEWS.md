Release a.b
-----------

Release 1.22.1 (14th July 2025)
-------------------------------

Bug fixes:

* SECURITY: Fix use-after-free bug in `samtools mpileup -a` due to reference
  sequences being discarded too early.  This could happen when the `-a`
  option was used, and all the alignments for one of the references started
  at the same position.  This caused mpileup to try to load the next reference
  before it had finished writing data for the previous one out.
  (PR #2229, fixes #2227.  Reported by Pouya Kheradpour)

* This release bundles htslib-1.22.1 and htscodecs v1.6.4.  Fixes
  a possible buffer overflow on some invalid CRAM inputs; and a
  failure when trying to save data with very long alignment records
  with sequence '*' as CRAM version 3.1.

  See the HTSlib and htscodecs NEWS files for details.

Documentation:

* Clarify that `-@` starts additional threads in `samtools index` help,
  and fix formatting issues in the `samtools reset` manual page.
  (PR #2225)

Build changes:

* Added settings to limit memory used by `samtools sort` when running tests.
  (PR #2226)

Release 1.22 (30th May 2025)
----------------------------

Note this release changes the default output CRAM version from 3.0 to 3.1.
HTSlib and SAMtools have been able to read CRAM 3.1 since version 1.12,
however other tools may not yet be able to cope.  We know Noodles reads CRAM
3.1 and htsjdk has a draft implementation, but not yet released.

HTSlib has options for modifying the output formats, which are exposed in
SAMtools.  When specifying an output format you can explicitly change the
version via e.g. `samtools view -O cram,version=3.0 ...`.

Further documentation on this change can be found at
https://www.htslib.org/benchmarks/CRAM.html


HTSlib no longer fetches CRAM reference data from EBI's server by default. Your
organisation may wish to set up local infrastructure to supply reference
sequences, e.g., using the new ref-cache tool included in this HTSlib release.
See the REF_CACHE and REF_PATH environment variables documented in
https://www.htslib.org/doc/reference_seqs.html and the SAMtools man
page for details.


New work and changes:

* New `samtools checksum` command.  This checksums sequence, name, quality and
  barcode tags in an order and orientation agnostic manner, to facilitate
  validation of no data loss between raw fastq (or unmapped crams) through
  alignment, duplication marking, sorting, and other processing operations to
  get to the final aligned bam/cram.
  (PR#2122)

* Extend `samtools sort -M` to distinguish between mapped and unmapped files.
  (PR#2110, fixes #2105.  Reported by Armin Töpfer)

* Allow the `samtools sort` "merging from..." message to be silenced.
  Setting the verbosity to 0 or 1 will now silence this message.
  (PR#2197, resolves #2185.  Requested by Alex Predeus)

* Add `--save-counts` option to `samtools view`.  Adds an option to store counts
  of records processed, accepted and rejected by filtering to a file.
  (PR#2120, resolves #2038.  Requested by Chang Y)

* `samtools fasta` and `fastq` can now make faidx/fqidx indexes while writing
  using the `--write-index` option.
  (PR#2125, resolves #2118.  Requested by Filipe G. Vieira)

* Add a warning for `samtools fastq` on coordinate sorted data.
  (PR#2176, fixes #2169 and #2161.  Reported by wook2014)

* `samtools tview` add `-i` to hide inserts.
  (PR#2123.  Thanks to Benjamin Bræstrup Sayoc)

* Show optional headers with `samtools bedcov -H`.
  (PR#2140, fixes #2126.  Reported by biounix)

* `samtools consensus` now supports proper multi-threading.  Previously
  this was restricted to decompression only, but it should now scale better.
  (PR#2174, supersedes PR#2141)

* Add `samtools consensus -T ref.fa` functionality.  This reports the reference
  value if a consensus value cannot be calculated.
  (PR#2153, fixes an additional request in #1915)

* In `samtools consensus`, do not use consensus N for "*" (absent) calls that
  are masked due to insufficient depth.
  (PR#2204, fixes #2167.  Reported by sanschiffre)

* Improve `plot-bamstats` quality plots.
  (PR#2143 combined with PR#2116 (thanks to James Gilbert))

* Make `reheader -h` use /tmp and honour TMPDIR.
  (PR#2168, related to #2165.  Reported by Zhang Yuanfeng)

* Set sort order header tag to unsorted when ordering is lost during
  `samtools merge`.
  (PR#2173, fixes #2159.  Reported by Filipe G. Vieira)

* Protect against merging CRAM files with different headers.
  (PR#2220, fixes #2218.  Reported by Kevin Lewis)

* `samtools stats` bug-fix to checksum calculation for quality values.  This
  corrects the checksums but in turn makes the calculated value different to
  that reported by previous samtools versions.
  (PR#2193, fixes #2187)

* Clarification for `samtools stats` when used on files with different sort
  orders.
  (PR#2198, fixes #2177.  Reported by Filipe G. Vieira)

* In `samtools stats`, dovetailed (completely overlapping) read pairs are
  now always counted as inward-oriented.  Previously they could have been
  inwards or outwards depending on read ordering.
  (PR#2216, resolves #2210.  Requested by Pontus Hüer)


Documentation:

* Correct the example for 1:1 `samtools consensus` coords.
  (PR#2113, fixes #2111.  Reported by schorlton-bugseq)

* Documents the fastq format options used in SAMtools and HTSlib.
  (PR#2123, fixes #2121)

* Remove mention of threads from `samtools cat` man page.
  (PR#2162, fixes #2160.  Reported by Brandon Pickett)

* Update `samtools merge` man page to include `--template-coordinate`.
  (PR#2164.  Thanks to Nils Homer)

* Revised CRAM reference sequence documentation in the samtools man page.
  (PR#2178)

* Added fish shell completion and renamed completion for bash shell.  These
  files can be copied to appropriate directories by the user.  For full
  functionality it requires Python3.5+ and installed samtools manpages.
  (PR#2203.  Thanks to LunarEclipse363)

* Fix URL printed by the `seq_cache_populate.pl` script.
  (PR#2222.  Thanks to Charles Plessy)

Bug fixes:

* `samtools consensus` previously could give different results for BAM and
  CRAM files with the same content.  This was because MD/NM tag generation
  was disabled in CRAM, but the `decode_md=0` option did nothing with BAM.
  Note with `--no-adj-MQ` both BAM and CRAM gave identical results.
  Now use `--input-fmt-option decode_md=0` to get the old CRAM behaviour.
  Otherwise, both BAM and CRAM will be utilising MD/NM to locally modify
  mapping quality.
  (PR#2156)

* `samtools consensus` without `-a` previously still padded with leading Ns in
  some cases.  It now consistently removes both leading and trailing Ns.
  Use "-a" if you want all reference bases displayed.
  (Part of PR#2174 above)

* Change how `markdup` looks for single reads.  Due to changes to `fixmate` in
  1.21 `markdup` no longer recognised single reads that would have normally have
  been part of a pair.
  (PR#2117, fixes #2117.  Reported by Kristy Horan)

* Fix `samtools merge` crash on BAM files with malformed headers.
  (PR#2128, fixes #2127.  Reported by Frostb1te)

* Fix `faidx --write-index` invalid free.
  (PR#2147, fixes 2142.  Reported by Alex Leonard)

* Fix `samtools fastq -i` to force CRAM aux tag decoding.
  (PR#2155, fixes #2155.  Reported by Alex Leonard)

Non user-visible changes and build improvements:

* Improve htslib#PRnum support for Cirrus-CI and GitHub Actions.
  (PR#2115)

* Fix broken tests due to MSYS2 changes. Due to changes in how MSYS2 perl
  reported the identity of the OS it was built for, our tests were failing to
  adapt to the Windows style file locations.
  (PR #2196)

* Upgrade to `_XOPEN_SOURCE=700`, to match HTSlib.  Also replace `usleep()` with
  `nanosleep()`.
  (PR#2221)


Release 1.21 (12th September 2024)
----------------------------------

Notice:

* This is the last SAMtools / HTSlib release where CRAM 3.0 will be
  the default CRAM version.  From the next we will change to CRAM 3.1
  unless the version is explicitly specified, for example using
  "samtools view -O cram,version=3.0".

New work and changes:

* `samtools reset` now removes a set of predefined auxtags, as these tags are no
  longer valid after the reset operation.  This behaviour can be overridden if
  desired.
  (PR #2034, fixes #2011.  Reported by Felix Lenner)

* `samtools reset` now also removes duplicate flags.
  (PR #2047.  Reported by Kevin Lewis)

* Region and section/part filtering added to CRAM `samtools cat`.  Region
  filtering permits `samtools cat` to produce new CRAMs that only cover a
  specified region.
  (PR #2035)

* Added a report of the number of alignments for each primer to
  `samtools ampliconclip`.
  (PR #2039, PR #2101, feature request #2033.  Thanks to Brad Langhorst)

* Make `ampliconclip` primer counts output deterministic.
  (PR #2081)

* `samtools fixmate` no longer removes the PAIRED flag from reads that have no
  mate.  This is done on the understanding that the PAIRED flag is a sequencing
  technology indicator not a feature of alignment.  This is a change to
  previous `fixmate` behaviour.
  (PR #2056, fixes #2052.  Reported by John Wiedenhoeft)

* Added bgzf compressed FASTA output to `samtools faidx`.
  (PR #2067, fixes #2055. Requested by Filipe G Vieira)

* Optimise `samtools depth` histogram incrementing code.
  (PR #2078)

* In `samtools merge` zero pad unique suffix IDs.
  (PR #2087, fixes #2086.  Thanks to Chris Wright)

* `samtools idxstats` now accepts the `-X` option, making it easier
  to specify the location of the index file.
  (PR #2093, feature request #2071.  Requested by Samuel Chen)

* Improved documentation for the mpileup `--adjust-MQ` option.
  (PR #2098.  Requested by Georg Langebrake)


Bug fixes:

* Avoid `tview` buffer overflow for positions with >= 14 digits.
  (PR #2032.  Thanks to John Marshall. Reported on
  bioconda/bioconda-recipes#47137 by jmunoz94)

* Added file name and error message to 'error closing output file'
  error in `samtools sort`.
  (PR #2050, fixes #2049.  Thanks to Joshua C Randall).

* Fixed hard clip trimming issue in `ampliconclip` where right-hand side
  qualities were being removed from left-hand side trims.
  (PR #2053, fixes #2048.  Reported by Duda5)

* Fixed a bug in `samtools merge --template-coordinate` where the wrong heap
  was being tested.
  (PR #2062.  Thanks to Nils Homer.  Reported on ng-core/fastquorum#52 by
  David Mas-Ponte)

* Do not look at chr "*" for unmapped-placed reads with
  `samtools view --fetch-pairs`.  This was causing a significant slowdown when
  `--fetch-pairs` was being used.
  (PR #2070, fixes #2059.  Reported by acorvelo)

* Fixed bug which could cause `samtools view -L` to give incomplete output
  when the BED file contained nested target locations.
  (PR #2107, fixes #2104.  Reported by geertvandeweyer)

* Enable `samtools coverage` to handle alignments that do not have quality score
  data.  This was causing memory access problems.
  (PR #2083, fixes #2076.  Reported by Matthew Colpus)

* Fix undefined behaviour in `samtools fastq` with empty QUAL.
  (PR #2084)

* In `plot-bamstats` fixed read-length plot for data with limited variations in
  length. Lack of data was causing gnuplot problems.
  (PR #2085, fixes #2068.  Reported by mariyeta)

* Fixed an accidental fall-through that caused `samtools split -p` to
  also enable `--no-PG`.
  (PR #2101)

* Fixed an overflow that caused `samtools consensus -m simple` to give
  incorrect output when the input coverage exceeded several million
  reads deep.
  (PR #2099, fixes #2095.  Reported by Dylan Lawrence)

Non user-visible changes and build improvements:

* Work around address sanitizer going missing from the Cirrus CI
  ubuntu clang compiler by moving the address sanitizer build to gcc.
  Fix warnings from the new clang compiler.
  (PR #2043)

* Windows based CI has been migrated from AppVeyor to GitHub Actions.
  (PR #2072, PR #2108)

* Turn on more warning options in Cirrus-CI builds, ensure everything builds
  with `-Werror`, and add undefined behaviour checks to the address sanitizer
  test.
  (PR #2101, PR #2103, PR #2109)

* Tidy up Makefile dependencies and untracked test files.
  (PR #2106.  Thanks to John Marshall)

Release 1.20 (15th April 2024)
------------------------------

* Added a `--max-depth` option to `bedcov`, for more control over the depth
  limit used when calculating the pileup.  Previously this limit was set
  at 64000; now it is set to over 2 billion, so effectively all bases will
  be counted.
  (PR #1970, fixes #1950.  Reported by ellisjj)

* Added `mpileup --output-extra RLEN` to display the unclipped read length.
  (PR #1971, feature request #1959.  Requested by Feng Tian)

* Improved checking of symbolic flag names (e.g. UNMAP) passed to samtools.
  (PR #1981, fixes #1977.  Reported by Ilya Shlyakhter)

* The `samtools consensus --min-depth` option now works for the Bayesian
  mode as well as the simple one.
  (PR #1989, feature request #1982.  Requested by Gautier Richard)

* It's now possible to use the `samtools fastq` `-d tag:val` option multiple
  times, allowing matches on more than one tag/value.  It also gets
  a `-D` option which allows the values to be listed in a file.
  (PR #1993, feature request #1958.  Requested by Tristan Lefebure)

* Added `samtools fixmate` `-M` option to sanity check base modification
  (`ML`, `MM`, `MN`) tags, and where necessary adjust modification data
  on hard-clipped records.
  (PR #1990)

* Made `mpileup` run faster.
  (PR #1995)

* `samtools import` now adds a `@PG` header to the files it makes.
  As with other sub-commands, this can be disabled by using `--no-PG`.
  (PR #2008.  Requested by Steven Leonard)

* The `samtools split` `-d` option to split by tag value now works on
  tags with integer values.
  (PR #2005, feature request #1956.  Requested by Alex Leonard)

* Adjusted `samtools sort -n` (by name) so that primary reads are always sorted
  before secondary / supplementary.
  (PR #2012, feature request #2010.  Requested by Stijn van Dongen)

* Added `samtools bedcov` `-H` option to print column headers in the
  output.
  (PR #2025.  Thanks to Dr. K. D. Murray)

Documentation:

* Added a note that BAQ is applied before filtering and overlap removal
  during mpileup processing.
  (PR #1988, fixes #1985.  Reported by Joseph Galasso)

* Added 3.1 to the list of supported CRAM versions in the samtools manual page.
  (PR #2009.  Thanks to Andrew Thrasher)

* Made assorted improvements to ampliconclip, flagstat and markdup manual
  pages.
  (PR #2014)

Bug Fixes:

* Security fix: Fixed double free that could occur if bed file indexing failed
  due to running out of memory.  This bug first appeared in version 1.19.1.
  (PR #2026)

* Corrected error message printed when faidx fails to load the fai index.
  (PR #1987.  Thanks to Nick Moore)

* Fixed bug introduced in release 1.4 that caused incorrect reference
  bases to be printed by `samtools mpileup -a -f ref.fa` in the zero-depth
  regions at the end of each reference.
  (PR #2019, fixes #2018.  Reported by Joe Georgeson)

* Fixed a samtools view usage crash on MinGW when given invalid options.
  (PR #2030, fixes #2029.  Reported by Divon Lan)

Non user-visible changes and build improvements:

* Added tests to ensure that CRAM compression is working properly.
  (PR #1969, part of fix for #1968.  Reported by Clockris)

Release 1.19.2 (24th January 2024)
----------------------------------

Bug Fixes:

* Fixed a regression in 1.19.1 that broke BED filtering for inputs where
  the region start positions for the same reference were not sorted in
  ascending order.
  (PR #1975, fixes #1974.  Reported by Anže Starič)

Release 1.19.1 (22nd January 2024)
----------------------------------

Bug Fixes:

* Fixed a possible array bounds violation when looking up regions in
  a BED file (e.g. using `samtools view -L`).  This could lead to crashes
  or the return of incomplete results if the BED file contained a
  large number of entries all referencing low positions on a chromosome.
  (PR #1962, fixes #1961.  Reported by geertvandeweyer)

* Fixed a crash in samtools stats that occurred when trying to clean up
  after it was unable to open a CRAM reference file.
  (PR #1957, fixes crash reported in samtools/htslib#1723.  Reported by
  Alex Leonard)

Documentation:

* Fixed inverted logic in the `samtools consensus --show-del` manual page
  description.
  (PR #1955, fixes #1951.  Reported by Mikhail Schelkunov)

* Added a description of the MPC section to the `samtools stats` manual page.
  (PR #1963, fixes #1954.  Reported by litun-fkby)

Release 1.19 (12th December 2023)
---------------------------------

New work and changes:

* Samtools coverage: add a new --plot-depth option to draw depth (of
  coverage) rather than the percentage of bases covered.
  (PR #1910.  Thanks to Pierre Lindenbaum)

* Samtools merge / sort: add a lexicographical name-sort option via the -N
  option.  The "natural" alpha-numeric sort is still available via -n.
  (PR #1900, fixes #1500.  Reported by Steve Huang)

* Samtools view: add -N ^NAME_FILE and -R ^RG_FILE options.  The standard -N
  and -R options only output reads matching a specified list of names or
  read-groups.  With a caret (^) prefix these may be negated to only output
  read not matching the specified files.
  (PR #1896, fixes #1895.  Suggested by Feng Tian)

* Cope with Htslib's change to no longer close stdout on hts_close.
  Htslib companion PR is samtools/htslib#1665.
  (PR #1909.  Thanks to John Marshall)

* Plot-bamstats: add a new plot of the read lengths ("RL") from samtools
  stats output.
  (PR #1922, fixes #1824.  Thanks to @erboone, suggested by Martin Pollard)

* Samtools split: support splitting files based on the contents of
  auxiliary tags.  Also adds a -M option to limit the number of files split
  can make, to avoid accidental resource over-allocation, and fixes some
  issues with --write-index.
  (PR #1222, PR #1933, fixes #1758.  Thanks to Valeriu Ohan, suggested by
  Scott Norton)

Bug Fixes:

* Samtools stats: empty barcode tags are now treated as having no barcode.
  (PR #1929, fixes #1926.  Reported by Jukka Matilainen)

* Samtools cat: add support for non-seekable streams.  The file format
  detection blocked pipes from working before, but now files may be
  non-seekable such as stdin or a pipe.
  (PR #1930, fixes #1731.  Reported by Julian Hess)

* Samtools mpileup -aa (absolutely all positions) now produces an output
  even when given an empty input file.
  (PR #1939.  Reported by Chang Y)

* Samtools markdup: speed up optical duplicate tagging on regions with very
  deep data.
  (PR #1952)

Documentation:

* Samtools mpileup: add more usage examples to the man page.
  (PR #1913, fixes #1801)

* Samtools fastq: explicitly document the order that filters apply.
  (PR #1907)

* Samtools merge: fix example output to use an uppercase RG PL field.
  (PR #1917.  Thanks to John Marshall.  Reported by Michael Macias)

* Add hclen SAM filter documentation.
  (PR #1902.  See also samtools/htslib#1660)

* Samtools consensus: remove the -5 option from documentation.  This option
  was renamed before the consensus subcommand was merged, but accidentally
  left in the man page.
  (PR #1924)


Release 1.18 (25th July 2023)
-----------------------------

New work and changes:

* Add minimiser sort option to collate by an indexed fasta.  Expand the
  minimiser sort to arrange the minimiser values in the same order as they
  occur in the reference genome. This is acts as an extremely crude and
  simplistic read aligner that can be used to boost read compression.
  (PR #1818)

* Add a --duplicate-count option to markdup.  Adds the number of duplicates
  (including itself) to the original read in a 'dc' tag.
  (PR #1816. Thanks to wulj2)

* Make calmd handle unaligned data or empty files without throwing an error.
  This is to make pipelines work more smoothly.  A warning will still be issued.
  (PR #1841, fixes #1839.  Reported by Filipe G. Vieira)

* Consistent, more comprehensive flag filtering for fasta/fastq.  Added
  --rf/--incl[ude]-flags and long options for -F (--excl[ude]-flags and
  -f (--require-flags).
  (PR #1842.  Thanks to Devang Thakkar)

* Apply fastq --input-fmt-option settings.  Previously any options specified
  were not being applied to the input file.
  (PR #1855.  Thanks to John Marshall)

* Add fastq -d TAG[:VAL] check.  This mirrors view -d and will only output
  alignments that match TAG (and VAL if specified).
  (PR #1863, fixes #1854.  Requested by Rasmus Kirkegaard)

* Extend import --order TAG to --order TAG:length.  If length is specified, the
  tag format goes from integer to a 0-padded string format.  This is a
  workaround for BAM and CRAM that cannot encode an order tag of over 4 billion
  records.
  (PR #1850, fixes #1847.  Reported by Feng Tian)

* New -aa mode for consensus.  This works like the -aa option in depth and
  mpileup. The single 'a' reports all bases in contigs covered by alignments.
  Double 'aa' (or '-a -a') reports Ns even for the references with no alignments
  against them.
  (PR #1851, fixes #1849.  Requested by Tim Fennell)

* Add long option support to samtools index.
  (PR #1872, fixes #1869.  Reported by Jason Bacon)

* Be consistent with rounding of "average length" in samtools stats.
  (PR #1876, fixes #1867.  Reported by Jelinek-J)

* Add option to ampliconclip that marks reads as unmapped when they do not have
  enough aligned bases left after clipping.  Default is to unmap reads with zero
  aligned bases.
  (PR #1865, fixes #1856.  Requested by ces)

Bug Fixes:

* [From HTSLib] Fix a major bug when searching against a CRAM index where
  one container has start and end coordinates entirely contained within the
  previous container. This would occasionally miss data, and sometimes
  return much more than required.  The bug affected versions 1.11 to 1.17,
  although the change in 1.11 was bug-fixing multi-threaded index queries.
  This bug did not affect index building.  There is no need to reindex your
  CRAM files.
  (PR #samtools/htslib#1574, PR #samtools/htslib#1640. Fixes
  #samtools/htslib#1569, #samtools/htslib#1639, #1808, #1819.  Reported by
  xuxif, Jens Reeder and Jared Simpson)

* Fix a sort -M bug (regression) when merging sub-blocks.  Data was valid but
  in a poor order for compression.
  (PR #1812)

* Fix bug in split output format.  Now SAM and CRAM format can chosen as well
  as BAM.  Also a documentation change, see below.
  (PR #1821)

* Add error checking to view -e filter expression code.  Invalid expressions
  were not returning an error code.
  (PR #1833, fixes #1829.  Reported by Steve Huang)

* Fix reheader CRAM output version.  Sets the correct CRAM output version for
  non-3.0 CRAMs.
  (PR #1868, fixes #1866.  Reported by John Marshall)

Documentation:

* Expand the default filtering information on the mpileup man page.
  (PR #1802, fixes #1801.  Reported by gevro)

* Add an explanation of the default behaviour of split files on generating
  a file for reads with missing or unrecognised RG tags.  Also a small bug fix,
  see above.
  (PR #1821, fixes #1817.  Reported by Steve Huang)

* In the INSTALL instructions, switched back to openssl for Alpine.  This
  matches the current Alpine Linux practice.
  (PR #1837, see htslib#1591.  Reported by John Marshall)

* Fix various typos caught by lintian parsers.
  (PR #1877.  Thanks to Étienne Mollier)

* Document consensus --qual-calibration option.
  (PR #1880, fixes #1879.  Reported by John Marshall)

* Updated the page about samtools duplicate marking with more detail at
  www.htslib.org/algorithms/duplicate.html

Non user-visible changes and build improvements:

* Removed a redundant line that caused a warning in gcc-13.
  (PR #1838)


Release 1.17 (21st February 2023)
---------------------------------

New work and changes:

* New samtools reset subcommand.  Removes alignment information.  Alignment
  location, CIGAR, mate mapping and flags are updated. If the alignment was in
  reverse direction, sequence and its quality values are reversed and
  complemented and the reverse flag is reset.  Supplementary and secondary
  alignment data are discarded.
  (PR #1767, implements #1682. Requested by dkj)

* New samtools cram-size subcommand.  It writes out metrics about a CRAM file
  reporting aggregate sizes per block "Content ID" fields, the data-series
  contained within them, and the compression methods used.
  (PR #1777)

* Added a --sanitize option to fixmate and view.  This performs some sanity
  checks on the state of SAM record fields, fixing up common mistakes made by
  aligners.
  (PR #1698)

* Made better use of threading during the merge stage of samtools sort.
  It also now limits the number of temporary files present by making
  intermediate merges if the count gets too high.
  (PR #1706)

* Permit 1 thread with samtools view.  All other subcommands already allow this
  and it does provide a modest speed increase.
  (PR #1755, fixes #1743. Reported by Goran Vinterhalter)

* Add CRAM_OPT_REQUIRED_FIELDS option for view -c.  This is a big speed up for
  CRAM (maybe 5-fold), but it depends on which filtering options are being used.
  (PR #1776, fixes #1775. Reported by Chang Y)

* New filtering options in samtools depth.  The new --excl-flags option is a
  synonym for -G, with --incl-flags and --require-flags added to match view
  logic.
  (PR #1718, fixes #1702. Reported by Dario Beraldi)

* Speed up calmd's slow handling of non-position-sorted data by adding caching.
  This uses more memory but is only activated when needed.
  (PR #1723, fixes #1595. Reported by lxwgcool)

* Improve samtools consensus for platforms with instrument specific profiles,
  considerably helping for data with very different indel error models and
  providing base quality recalibration tables. On PacBio HiFi, ONT and
  Ultima Genomics consensus qualities are also redistributed within homopolymers
  and the likelihood of nearby indel errors is raised.
  (PR #1721, PR #1733)

* Consensus --mark-ins option.  This permits he consensus output to include a
  markup indicating the next base is an insertion. This is necessary as we need
  a way of outputting both consensus and also how that consensus marries up with
  the reference coordinates.
  (PR #1746)

* Consensus min depth now works on the overall depth of the data, in line with
  the documentation.  Added a --min-BQ (minimum base quality) filtering option.
  (PR #1704, fixes #1700 reported by Thiseas C. Lamnidis)

* Make faidx/fqidx output line length default to the input line length.
  (PR #1738, fixes #1734. Reported by John Marshall)

* Speed up optical duplicate checking where data has a lot of duplicates
  compared to non-duplicates.
  (PR #1779, fixes #1771. Reported by Poshi)

* Added option to add read group matching to the duplicate criteria.  Change
  stats to include read groups and add an optional JSON format.
  (PR #1699, implements #1143 requested by Jonathan Keats)

* For collate use TMPDIR environment variable, when looking for a temporary
  folder.
  (PR #1782, based on PR #1178 and fixes #1172.  Reported by Martin Pollard)

Bug Fixes:

* Fix stats breakage on long deletions when given a reference.
  (PR #1712, fixes #1707. Reported by John Didion)

* In ampliconclip, stop hard clipping from wrongly removing entire reads.
  (PR #1722, fixes #1717. Reported by Kevin Xu)

* Fix bug in ampliconstats where references mentioned in the input file headers
  but not in the bed file would cause it to complain that the SAM headers were
  inconsistent.
  (PR #1727, fixes #1650. Reported by jPontix)

* Fixed SEGV in samtools collate when no filename given.
  (PR #1724)

* Changed the default UMI barcode regex in markdup.  The old regex was too
  restrictive.  This version will at least allow the default read name UMI as
  given in the Illumina example documentation.
  (PR #1737, fixes #1730. Reported by yloemie)

* Fix samtools consensus buffer overrun with MD:Z handling.
  (PR #1745, fixes #1744. Reported by trilisser)

* Fix a buffer read-overflow in mpileup and tview on sequences with seq "*".
  (PR #1747)

* Fix view -X command line parsing that was broken in 1.15.
  (PR #1772, fixes #1720.  Reported by Francisco Rodríguez-Algarra
  and Miguel Machado)

* Stop samtools view -d from reporting meaningless system errors when
  tag validation fails.
  (PR #1796)

Documentation:

* Add a description of the samtools tview display layout to the man page.
  Documents . vs , and upper vs lowercase. Adds a -s sample example, and
  documents the -w option.
  (PR #1765, fixes #1759. Reported by Lucas Ferreira da Silva)

* Clarify intention of samtools fasta/q in man page and soft vs hard
  clipping.
  (PR #1794, fixes #1792. Reported by Ryan Lorig-Roach)

* Minor fix to wording of mpileup --rf usage and man page.
  (PR #1795, fixes #1791. Reported by Luka Pavageau)

Non user-visible changes and build improvements:

* Use POSIX grep in testing as egrep and fgrep are considered obsolete.
  (PR #1726, thanks to David Seifert)

* Switch MacOS CI tests to an ARM-based image.
  (PR #1770)


Release 1.16.1 (2nd September 2022)
-----------------------------------

Bug fixes:

 * Fixed a bug with the template-coordinate sort which caused incorrect
   ordering when using threads, or processing large files that don't
   fit completely in memory.
   (PR #1703, thanks to Nils Homer)

 * Fixed a crash that occurred when trying to use `samtools merge` in
   template-coordinate mode.
   (PR #1705, thanks to Nils Homer)

Release 1.16 (18th August 2022)
-------------------------------

New work and changes:

 * samtools reference command added.  This subcommand extracts the embedded
   reference out of a CRAM file.
   (PR #1649, addresses #723.  Requested by Torsten Seemann)

 * samtools import now adds grouped by query-name to the header.
   (PR #1633, thanks to Nils Homer)

 * Made samtools view read error messages more generic.  Former error message
   would claim that there was a "truncated file or corrupt BAM index file" with
   no real justification.  Also reset errno in stream_view which could lead to
   confusing error messages.
   (PR #1645, addresses some of the issues in #1640.  Reported by Jian-Guo Zhou)

 * Make samtools view -p also clear mqual, tlen and cigar.
   (PR #1647, fixes #1606.  Reported by eboyden)

 * Add bedcov option -c to report read count.
   (PR #1644, fixes #1629.  Reported by Natchaphon Rajudom)

 * Add UMI/barcode handling to samtools markdup.
   (PR #1630, fixes #1358 and #1514.  Reported by Gert Hulselmans and Poshi)

 * Add a new template coordinate sort order to samtools sort and
   samtools merge.  This is useful when working with unique molecular
   identifiers (UMIs).
   (PR #1605, fixes #1591.  Thanks to Nils Homer)

 * Rename mpileup --ignore-overlaps to --ignore-overlaps-removal
   or --disable-overlap-removal.  The previous name was ambiguous and was often
   read as an option to enable removal of overlapping bases, while in reality
   this is on by default and the option turns off the ability to remove
   overlapping bases.
   (PR #1666, fixes #1663.  Reported by yangdingyangding)

 * The dict command can now read BWA's .alt file and add AH:* tags
   indicating reference sequences that represent alternate loci.
   (PR #1676.  Thanks to John Marshall)

 * The "samtools index" command can now accept multiple alignment filenames
   with the new -M option, and will index each of them separately. (Specifying
   the output index filename via out.index or the new -o option is currently
   only applicable when there is only one alignment file to be indexed.)
   (PR #1674.  Reported by Abigail Ramsøe and Nicola Romanò.
    Thanks to John Marshall)

 * Allow samtools fastq -T "*".
   This allows all tags from SAM records to be written to fastq headers. This is
   a counterpart to samtools import -T "*".
   (PR #1679.  Thanks to cjw85)

Bug Fixes:

 * Re-enable --reference option for samtools depth.  The reference is not used
   but this makes the command line usage compatible with older releases.
   (PR #1646, fixes #1643.  Reported by Randy Harr)

 * Fix regex coordinate bug in samtools markdup.
   (PR #1657, fixes #1642.  Reported by Randy Harr)

 * Fix divide by zero in plot-bamstats -m, on unmapped data.
   (PR #1678, fixes #1675.  Thanks to Shane McCarthy)

 * Fix missing RG headers when using samtools merge -r.
   (PR #1683, addresses htslib#1479.  Reported by Alex Leonard)

 * Fix a possible unaligned access in samtools reference.
   (PR #1696)

Documentation:

 * Add documentation on CRAM compression profiles and some of the newer options
   that appear in CRAM 3.1 and above.
   (PR #1659, fixes #1656.  Reported by Matthias De Smet)

 * Add "sclen" filter expression keyword documentation.
   (PR #1661, see also htslib#1441)

 * Extend FILTER EXPRESSION man page section to match the changes made in
   HTSlib.
   (PR #1687, samtools/htslib#1478)

Non user-visible changes and build improvements:

 * Ensure generated test files are ignored (by git) and cleaned (by make
   testclean)
   (PR #1692, Thanks to John Marshall)

Release 1.15.1 (7th April 2022)
-------------------------------

Bug fixes:

 * A bug which prevented the samtools view --region-file (and the
   equivalent -M -L <file>) options from working in version 1.15 has
   been fixed. (#1617)

 * Fixed a crash triggered by using the samtools view -c/--count and --unmap
   options together.  The --unmap option is now ignored in counting mode.
   (#1619)

Documentation:

 * The consensus command was missing from the main samtools.1 manual page.
   It has now been added. (#1603)

 * Corrected instructions for reproducing the samtools stats "raw total
   sequences" count using samtools view -c. (#1620; reported by @krukanna)

 * Improved manual page formatting.  (#1625; thanks to John Marshall)

Non user-visible changes and build improvements:

 * Unnecessary #include lines have been removed from bam_plcmd.c. (#1607;
   thanks to John Marshall)

Release 1.15 (21st February 2022)
---------------------------------

Notice:

 * Samtools mpileup VCF and BCF output (deprecated in release 1.9) has been
   removed.  Please use bcftools mpileup instead.

New work and changes:

 * Added "--min-BQ" and "--min-MQ" options to "depth". These match the
   equivalent long options found in "samtools mpileup" and gives a consistent
   way of specifying the base and mapping quality filters.
   (#1584; fixes #1580. Reported by Chang Y)

 * Improved automatic file type detection with "view -u" or "view -1".  Setting
   either of these options would default to BAM format regardless of the usual
   automatic file type selection based on the file name.  The defaults are now
   only used when the file name does not indicate otherwise.
   (#1582)

 * For "markdup" optical duplicate marking add regex options for custom
   coordinates.  For the case of non standard read names (QNAME), add options
   to read the coordinates and, optionally, another part of the string to test
   for optical duplication.
   (#1558)

 * New "samtools consensus" subcommand for generating consensus from SAM, BAM or
   CRAM files based on the contents of the alignment records.  The consensus is
   written as FASTA, FASTQ or as a pileup oriented format.  The default
   FASTA/FASTQ output includes one base per non-gap consensus, with insertions
   with respect to the aligned reference being included and deletions removed.
   This could be used to compute a new reference from sequence assemblies to
   realign against.
   (#1557)

 * New "samtools view --fetch-pairs" option.  This options retrieves pairs even
   when the mate is outside of the requested region.  Using this option enables
   the multi-region iterator and a region to search must be specified.  The
   input file must be an indexed regular file.
   (#1542)

 * Building on #1530 below, add a tview reflist for Goto.
   (#1539, thanks to Adam Blanchet)

 * Completion of references added to tview Goto.
   (#1530; thanks to Adam Blanchet)

 * New "samtools head" subcommand for conveniently displaying the headers
   of a SAM, BAM, or CRAM file. Without options, this is equivalent to
   `samtools view --header-only --no-PG` but more succinct and memorable.
   (#1517; thanks to John Marshall)

Bug Fixes:
 * Free memory when stats fails to read the header of a file.
   (#1592; thanks to Mathias Schmitt)

 * Fixed empty field on unsupported aux tags in "mpileup --output-extra".
   Replaces the empty fields on unsupported aux tags with a '*'.
   (#1553; fixes #1544. Thanks to Adam Blanchet)

 * In mpileup, the --output-BP-5 and --output-BP are no longer mutually
   exclusive.  This fixes the problem of output columns being switched.
   (#1540; fixes 1534.  Reported by Konstantin Riege)

 * Fix for hardclip bug in ampliconclip.  Odd length sequences resulted in
   random characters appearing in sequence.
   (#1538; fixes #1527. Reported by Ivana Mihalek)

Documentation:

 * Improved mpileup documentation.
   (#1566; fixes #1564.  Reported by Chang Y)

 * Fixed "samtools depth -J" documentation, which was reversed.
   (#1552; fixes #1549.  Reported by Stephan Hutter)

 * Numerous minor man page fixes.
   (#1528, #1536, #1579, #1590.  Thanks to John Marshall for some of these)

Non user-visible changes and build improvements:

 * Replace CentOS test build with Rocky Linux.  The CentOS Docker images that
   our test build depended on has stopped working.  Switched to Rocky Linux as
   the nearest available equivalent.
   (#1589)

 * Fix missing autotools on Appveyor.  Newer versions of msys2 removed autotools
   from their base-devel package.  This is putting them back.
   (#1575)

 * Fixed bug detected by clang-13 with -Wformat-security.
   (#1553)

 * Switch to using splaysort in bam_lpileup.  Improves speed and efficiency in
   "tview".
   (#1548; thanks to Adam Blanchet)


Release 1.14 (22nd October 2021)
--------------------------------

Notice:

 * Samtools mpileup VCF and BCF output (deprecated in release 1.9)
   will be removed in the next release.  Please use bcftools mpileup
   instead.

New work and changes:

 * The legacy samtools API (libbam.a, bam_endian.h, sam.h and most of
   bam.h) has been removed.  We recommend coding against the HTSlib
   API directly. The legacy API had not been actively maintained since
   2015. (#1483)

 * New "samtools samples" command to list the samples used in a
   SAM/BAM/CRAM file. (#1432; thanks to Pierre Lindenbaum)

 * "mpileup" now supports base modifications via the SAM Mm/MM
   auxiliary tag.  Please see the "--output-mods" option. (#1311)

 * Added "mpileup --output-BP-5" option to output the BP field in 5'
   to 3' order instead of left to right. (#1484; fixes #1481)

 * Added "samtools view --rf" option as an additional FLAG filtering
   method.  This keeps records only if (FLAG & N) != 0. (#1508; fixes
   #1470)

 * New "samtools import -N" option to use the second word on a FASTQ
   header line, matching the SRA/ENA FASTQ variant. (#1485)

 * Improve "view -x" option to simplify specifying multiple tags, and
   added the reverse "--keep-tag" option to include rather than
   exclude. (#516)

 * Switched the processing order of "view" -x (tag filtering) and -e
   (expression) handling.  Expressions now happen first so we can
   filter on tags which are about to be deleted.  This is now
   consistent with the "view -d" behaviour too. (#1480; fixes
   #1476. Reported by William Rowell)

 * Added filter expression "endpos" keyword. (#1464.  Thanks to
   John Marshall)

 * "samtools view" errors now appear after any SAM output, improving
   their visibility. (#1490.  Thanks to John Marshall)

 * Improved "samtools sort" use of temporary files, both tidying up
   if it fails and recovery when facing pre-existing temporary files.
   (#1510; fixes #1035, #1503.  Reported by Vivek Rai and
   Maarten Kooyman)

 * Filtering in "samtools markdup" now sets the UNMAP BAM flag when
   given the "-p" option. (#1512; fixes #1469)

 * Make CRAM references shared during "samtools merge" so merging many
   files has a lower memory usage. (#471)

Bug fixes:

 * Prevent "samtools depth" from closing stdout when outputting to
   terminal, avoiding a bad interaction with PySam. (#1465.  Thanks
   to John Marshall)

 * In-place "samtools reheader" now works on CRAMs produced using a
   higher than default compression level. (#1479)

 * Fix setting of the dt tag in "markdup".  Optical duplicates were
   being marked too early, negating the tagging and counting elsewhere.
   (#1487; fixes #1486.  Reported by Kevin Lewis)

 * Reinstate the "samtools stats -I" option to filter by sample.
   (#1496; fixes #1489.  Reported by Matthias Bernt)

 * Fix "samtools fastq" handling of dual index tags on single-ended
   input. (#1474)

 * Improve "samtools coverage" documentation. (#1521; fixes #1504.
   Reported by Peter Menzel)

Non user-visible changes and build improvements:

 * Replace Curses mvprintw() with va_list-based
   equivalent. (#1509. Thanks to John Marshall and Andreas Tille)

 * Fixed some clang-13 warning messages. (#1506)

 * Improve quoting of options in "samtools import" tests. (#1466.
   Thanks to John Marshall)

 * Fixed a faulty test which caused test harness failures on NetBSD. (#1520)

Release 1.13 (7th July 2021)
----------------------------

 * Fixed samtools view FILE REGION, mpileup -r REGION, coverage -r REGION and
   other region queries: fixed bug introduced in 1.12, which led to region
   queries producing very few reads for some queries (especially for larger
   target regions) when unmapped reads were present.
   Thanks to @vinimfava (#1451), @JingGuo1997 (#1457) and Ramprasad Neethiraj
   (#1460) for reporting the respective issues.

 * Added options to set and clear flags to samtools view.  Along with the
   existing remove aux tags this gives the ability to remove mark duplicate
   changes (part of #1358)
   (#1441)

 * samtools view now has long option equivalents for most of its single-letter
   options. Thanks to John Marshall.
   (#1442)

 * A new tool, samtools import, has been added.  It reads one or more FASTQ
   files and converts them into unmapped SAM, BAM or CRAM.
   (#1323)

 * Fixed samtools coverage error message when the target region name is not
   present in the file header. Thanks to @Lyn16 for reporting it.
   (#1462; fixes #1461)

 * Made samtools coverage ASCII mode produce true ASCII output.  Previously it
   would produce UTF-8 characters.
   (#1423; fixes #1419)

 * samtools coverage now allows setting the maximum depth, using the -d/--depth
   option. Also, the default maximum depth has been set to 1000000.
   (#1415; fixes #1395)

 * Complete rewrite of samtools depth.  This means it is now considerably faster
   and does not need a depth limit to avoid high memory usage.  Results should
   mostly be the same as the old command with the potential exception of overlap
   removal.
   (#1428; fixes #889, helps ameliorate #1411)

 * samtools flags now accepts any number of command line arguments,
   allowing multiple SAM flag combinations to be converted at once.  Thanks to
   John Marshall.
   (#1401, fixes #749)

 * samtools ampliconclip, ampliconstats and plot-ampliconstats now support
   inputs that list more than one reference.
   (#1410 and #1417; fixes #1396 and #1418)

 * samtools ampliconclip now accepts the --tolerance option, which allows the
   user to set the number of bases within which a region is matched.  The
   default is 5.
   (#1456)

 * Updated the documentation on samtools ampliconclip to be clearer about what
   it does.  From a suggestion by Nathan S Watson-Haigh.
   (#1448)

 * Fixed negative depth values in ampliconstats output.
   (#1400)

 * samtools addreplacerg now allows for updating (replacing) an existing
   `@RG` line in the output header, if a new `@RG` line is provided in the
   command line, via the -r argument. The update still requires the user's
   approval, which can be given with the new -w option.  Thanks to Chuang Yu.
   (#1404)

 * Stopped samtools cat from outputting multiple CRAM EOF markers.
   (#1422)

 * Three new counts have been added to samtools flagstat: primary, mapped
   primary and duplicate primary.
   (#1431; fixes #1382)

 * samtools merge now accepts a `-o FILE` option specifying the output file,
   similarly to most other subcommands. The existing way of specifying it
   (as the first non-option argument, alongside the input file arguments)
   remains supported. Thanks to David McGaughey and John Marshall.
   (#1434)

 * The way samtools merge checks for existing files has been changed
   so that it does not hang when used on a named pipe.
   (#1438; fixes #1437)

 * Updated documentation on mpileup to highlight the fact that the filtering
   options on FLAGs work with ANY rules.
   (#1447; fixes #1435)

 * samtools can now be configured to use a copy of HTSlib that has been set
   up with separate build and source trees.  When this is the case, the
   `--with-htslib` configure option should be given the location of the HTSlib
   build tree.  (Note that samtools itself does not yet support out-of-tree
   builds).  Thanks to John Marshall.
   (#1427; companion change to samtools/htslib#1277)


Release 1.12 (17th March 2021)
------------------------------

 * The legacy samtools API (libbam.a, bam.h, sam.h, etc) has not been actively
   maintained since 2015. It is deprecated and will be removed entirely in a
   future SAMtools release. We recommend coding against the HTSlib API directly.

 * I/O errors and record parsing errors during the reading of SAM/BAM/CRAM
   files are now always detected. Thanks to John Marshall (#1379; fixed #101)

 * New make targets have been added: check-all, test-all, distclean-all,
   mostlyclean-all, testclean-all, which allow SAMtools installations to
   call corresponding Makefile targets from embedded HTSlib installations.

 * samtools --version now displays a summary of the compilation details and
   available features, including flags, used libraries and enabled plugins
   from HTSlib. As an alias, `samtools version` can also be used. (#1371)

 * samtools stats now displays the number of supplementary reads in the
   SN section. Also, supplementary reads are no longer considered when
   splitting read pairs by orientation (inward, outward, other). (#1363)

 * samtools stats now counts only the filtered alignments that overlap
   target regions, if any are specified. (#1363)

 * samtools view now accepts option -N, which takes a file containing
   read names of interest. This allows the output of only the reads with
   names contained in the given file. Thanks to Daniel Cameron. (#1324)

 * samtools view -d option now works without a tag associated value, which
   allows it to output all the reads with the given tag. (#1339; fixed #1317)

 * samtools view -d and -D options now accept integer and single character
   values associated with tags, not just strings. Thanks to `@dariome` and
   Keiran Raine for the suggestions. (#1357, #1392)

 * samtools view now works with the filtering expressions introduced by HTSlib.
   The filtering expression is passed to the program using the specific option
   -e or the global long option --input-fmt-option. E.g.
   `samtools view -e 'qname =~ "#49$" && mrefid != refid && refid != -1 &&
    mrefid != -1' align.bam` looks for records with query-name ending in `#49`
   that have their mate aligned in a different chromosome. More details can be
   found in the FILTER EXPRESSIONS section of the main man page. (#1346)

 * samtools markdup now benefits from an increase in performance in the
   situation when a single read has tens or hundreds of thousands of duplicates.
   Thanks to @denriquez for reporting the issue. (#1345; fixed #1325)

 * The documentation for samtools ampliconstats has been added to the samtools
   man page. (#1351)

 * A new FASTA/FASTQ sanitizer script (`fasta-sanitize.pl`) was added, which
   corrects the invalid characters in the reference names. (#1314) Thanks to
   John Marshall for the installation fix. (#1353)

 * The CI scripts have been updated to recurse the HTSlib submodules when
   cloning HTSlib, to accommodate for the CRAM codecs, which now reside in the
   htscodecs submodule. (#1359)

 * The CI integrations now include Cirrus-CI rather than Travis. (#1335; #1365)

 * Updated the Windows image used by Appveyor to 'Visual Studio 2019'.
  (#1333; fixed #1332)

 * Fixed a bug in samtools cat, which prevented the command from running in
   multi-threaded mode. Thanks to Alex Leonard for reporting the issue.
   (#1337; fixed #1336)

 * A couple of invalid CIGAR strings have been corrected in the test data.
   (#1343)

 * The documentation for `samtools depth -s` has been improved.
   Thanks to @wulj2. (#1355)

 * Fixed a `samtools merge` segmentation fault when it failed to merge
   header `@PG` records. Thanks to John Marshall.  (#1394; reported by
   Kemin Zhou in #1393)

 * Ampliconclip and ampliconstats now guard against the BED file
   containing more than one reference (chromosome) and fail when
   found.  Adding proper support for multiple references will appear
   later.  (#1398)


Release 1.11 (22nd September 2020)
----------------------------------

 * New samtools ampliconclip sub-command for removing primers from
   amplicon-based sequencing experiments, including the current
   COVID-19 projects.  The primers are listed in a BED file and can be
   either soft-clipped or hard-clipped. (#1219)

 * New samtools ampliconstats sub-command to produce a textual summary
   of primer and amplicon usage, in a similar style to "samtools
   stats". The misc/plot-ampliconstats script can generate PNG
   images based on this text report. (#1227)

 * Samtools fixmate, addreplacerg, markdup, ampliconclip and sort now
   accept a -u option to enable uncompressed output, which is useful when
   sending data over a pipe to another process.  Other subcommands which
   already support this option for the same purpose are calmd, collate,
   merge, view and depad. (#1265)

 * samtools stats has a new GCT section, where it reports ACGT content
   percentages, similar to GCC but taking into account the read orientation.
   (#1274)

 * Samtools split now supports splitting by tag content with the -d
   option (#1211)

 * samtools merge now accepts a BED file as a command line argument (-L) and
   does the merging only with reads overlapping the specified regions (#1156)

 * Samtools sort now has a minhash collation (-M) to group unmapped
   reads with similar sequence together.  This can sometimes
   significantly reduce the file size. (#1093)

 * Samtools bedcov now has -g and -G options to filter-in and
   filter-out based on the FLAG field.  Also the new -d option adds
   an extra column per file counting the number of bases with a depth
   greater than or equal to a given threshold. (#1214)

 * Fixed samtools bedcov -j option (discard deletions and ref-skips)
   with multiple input files (#1212)

 * samtools bedcov will now accept BED files with columns separated by
   spaces as well as tabs (#1246; #1188 reported by Mary Carmack)

 * samtools depth can now include deletions (D) when computing the base
   coverage depth, if the user adds the -J option to the command
   line (#1163).

 * samtools depth will count only the bases of one read, for the overlapping
   section of a read pair, if the -s option is used in the command line
   (#1241, thanks to Teng Li).

 * samtools depth will now write zeros for the entire reference length,
   when "samtools depth -aa" is run on a file with no alignments. (#1252;
   #1249 reported by Paul Donovan)

 * Stopped depth from closing stdout, which triggered test fails
   in pysam (#1208, thanks to John Marshall).

 * samtools view now accepts remote URIs for FASTA and FAI files.
   Furthermore, the reference and index file can be provided in a single
   argument, such as
   samtools view -T ftp://x.com/ref.fa##idx##ftp://y.com/index.fa.fai a.cram
   (#1176; samtools/htslib#933 reported by @uitde007)

 * samtools faidx gets new options --fai-idx and --gzi-idx to allow
   specification of the locations of the .fai and (if needed) .gzi index
   files. (#1283)

 * The samtools fasta/fastq '-T' option can now add SAM array (type 'B') tags
   to the output header lines. (#1301)

 * samtools mpileup can now display MAPQ either as ASCII characters (with
   -s/--output-MQ; column now restored to its documented order as in 1.9 and
   previous versions) or comma-separated numbers (with --output-extra MAPQ;
   in SAM column order alongside other selected --output-extra columns).

   When both -s/--output-MQ and -O/--output-BP are used, samtools 1.10 printed
   the extra columns in the opposite order.  This changes the format produced
   by 1.10's --output-extra MAPQ. (#1281, thanks to John Marshall; reported
   by Christoffer Flensburg)

 * samtools tview now accepts a -w option to set the output width in
   text mode (-d T). (#1280)

 * The dict command can now add AN tags containing alternative names with
   "chr" prefixes added to or removed from each sequence name as appropriate
   and listing both "M" and "MT" alternatives for mitochondria. (#1164,
   thanks to John Marshall)

 * The samtools import command, labelled as obsolete in May 2009 and
   removed from all help and documentation later that year, has
   finally been removed.  Use samtools view instead. (#1185)

 * Replaced the remaining usage of the Samtools 0.1 legacy API with
   htslib calls. (#1187, thanks to John Marshall)

 * Documentation / help improvements (#1154; #1168; #1191; #1199;
   #1204; #1313):
   - Fixed a few man-page table layout issues
   - Added <file>##idx##<index> filename documentation
   - Fixed usage statement for samtools addreplacerg
   - Miscellaneous spelling and grammar fixes
   - Note fixmate/markdup name collated rather than name sorted input
   - Note that fastq and fasta inputs should also be name collated
   - Reshuffled order of main man-page and added -@ to more sub-pages
   - The misc/seq_cache_populate.pl script now gives REF_CACHE guidance

 * Additional documentation improvements, thanks to John Marshall (#1181;
   #1224; #1248; #1262; #1300)
   - Emphasise that samtools index requires a position-sorted file
   - Document 2^29 chromosome length limit in BAI indexes
   - Numerous typing, spelling and formatting fixes

 * Improved the message printed when samtools view fails to read its
   input (#1296)

 * Added build support for the OpenIndiana OS (#1165, thanks to John Marshall)

 * Fixed failing tests on OpenBSD (#1151, thanks to John Marshall)

 * The samtools sort tests now use less memory so the test suite works better
   on small virtual machines. (#1159)

 * Improved markdup's calculation of insert sizes (#1161)
   Also improved tests (#1150) and made it run faster when not checking
   for optical duplicates or adding 'do' tags (#1308)

 * Fixed samtools coverage minor inconsistency vs idxstats (#1205; #1203
   reported by @calliza)

 * Fixed samtools coverage quality thresholding options which were the
   wrong way round compared to mpileup (-q is the mapping quality threshold
   and -Q is base quality). (#1279; #1278 reported by @kaspernie)

 * Fixed bug where `samtools fastq -i` would add two copies of the barcode
   in the fastq header if both reads in a pair had a "BC:Z" tag (#1309;
   #1307 reported by @mattsoup)

 * Samtools calmd no longer errors with a SEQ of "*" (#1230; #1229
   reported by Bob Harris)

 * Samtools tview now honours $COLUMNS, fixing some CI tests (#1171; #1162
   reported by @cljacobs)

 * Fixed a samtools depad overflow condition (#1200)

 * Improved curses detection in configure script (#1170, #577, #940)

 * Fixed samtools stats integer overflows and added support for long
   references (#1174; #1173)

 * Fixed a 1-byte undersized memory allocation in samtools merge. (#1302)


Release 1.10 (6th December 2019)
--------------------------------

Changes affecting the whole of samtools, or multiple sub-commands:

 * Samtools now uses the new HTSlib header API.  As this adds more checks
   for invalid headers, it is possible that some illegal files will now
   be rejected when they would have been allowed by earlier versions. (#998)

   Examples of problems that will now be rejected include `@SQ` lines with
   no SN: tag, and `@RG` or `@PG` lines with no ID: tag.

 * samtools sub-commands will now add `@PG` header lines to output sam/bam/cram
   files.  To disable this, use the '--no-PG' option. (#1087; #1097)

 * samtools now supports alignment records with reference positions greater
   than 2 gigabases.  This allows samtools to process alignments for
   species which have large chromosomes, like axolotl and lungfish.  Note that
   due to file format limitations, data with large reference positions
   must use the SAM format. (#1107; #1117)

 * Improved the efficiency of reading and writing SAM format data by 2 fold
   (single thread). This is further improved by the ability to use multiple
   threads, as previously done with BAM and CRAM.

 * samtools can now write BGZF-compressed SAM format.  To enable this,
   either save files with a '.sam.gz' suffix, or use '--output-fmt sam.gz'.

 * samtools can now index BGZF-compressed SAM files.

 * The region parsing code has been improved to handle colons in reference
   names.  Strings can be disambiguated by the use of braces, so for
   example when reference sequences called "chr1" and "chr1:100-200"
   are both present, the regions "{chr1}:100-200" and "{chr1:100-200}"
   unambiguously indicate which reference is being used. (#864)

 * samtools flags, flagstats, idxstats and stats now have aliases
   flag, flagstat, idxstat and stat. (#934)

 * A new global '--write-index' option has been added.  This allows output
   sam.gz/bam/cram files to be indexed while they are being written out.
   This should work with addreplacerg, depad, markdup, merge, sort, split,
   and view. (#1062)

 * A global '--verbosity' option has been added to enable/disable
   debugging output. (#1124, thanks to John Marshall)

 * It is now possible to have data and index files stored in different
   locations.  There are two ways to tell samtools where to find the
   index:

   1. Samtools bedcov, depth, merge, mpileup, stats, tview, and view
      accept a new option (-X).  When this is used, each input sam/bam/cram
      listed on the command line should have a corresponding index file.
      Note that all the data files should be listed first, followed by all
      the index files. (#978, thanks to Mingfei Shao)

   2. A delimiter '##idx##' can be appended to the data file name followed
      by the index file name.  This can be used both for input files and
      outputs when indexing on-the-fly.

 * HTSlib (and therefore SAMtools) now uses version 4 signatures by default
   for its s3:// plug-in.  It can also write to S3 buckets, as long as
   version 4 signatures are in use.  See HTSlib's NEWS file and
   htslib-s3-plugin manual page for more information.

 * HTSlib (and therefore SAMtools) no longer considers a zero-length empty
   file to be a valid SAM file.  This has been changed so that pipelines such
   as `somecmd | samtools ...` with `somecmd` aborting before outputting
   anything will now propagate the error to the second command.

 * The samtools manual page has been split up into one for each
   sub-command.  The main samtools.1 manual page now lists the sub-commands
   and describes the common global options. (#894)

 * The meaning of decode_md, store_md and store_nm in the fmt-option section
   of the samtools.1 man page has been clarified. (#898, thanks to Evan Benn)

 * Fixed numerous memory leaks. (#892)

 * Fixed incorrect macro definition on Windows. (#950)

 * bedcov, phase, misc/ace2sam and misc/wgsim now check for failure to open
   files.  (#1013, thanks to Julie Blommaert and John Marshall)

Changes affecting specific sub-commands:

 * A new "coverage" sub-command has been added.  This prints a tabular format
   of the average coverage and percent coverage for each reference sequence,
   as well as number of aligned reads, average mapping quality and base
   quality.  It can also (with the '-m' option) plot a histogram of
   coverage across the genome. (#992, thanks to Florian Breitwieser)

 * samtools calmd:

   - Reference bases in MD: tags are now converted to upper case. (#981, #988)

 * samtools depth:

   - Add new options to write a header to the output (-H) and to direct
     the output to a file (-o).  (#937, thanks to Pierre Lindenbaum)

   - New options '-g' and '-G' can be used to filter reads. (#953)

   - Fix memory leak when failing to set CRAM options. (#985, thanks
     to Florian Breitwieser)

   - Fix bug when using region filters where the '-a' option did not
     work for regions with no coverage. (#1113; #1112 reported by
     Paweł Sztromwasser)

 * samtools fasta and fastq:

   - '-1 FILE -2 FILE' with the same filename now works properly. (#1042)

   - '-o FILE' is added as a synonym for '-1 FILE -2 FILE'. (#1042)

   - The '-F' option now defaults to 0x900 (SECONDARY,SUPPLEMENTARY).
     Previously secondary and supplementary records were filtered internally
     in a way that could not be turned off. (#1042; #939 reported
     by @finswimmer)

   - Allow reading from a pipe without an explicit '-' on the command line.
     (#1042; #874 reported by John Marshall)

   - Turn on multi-threading for bgzf compressed output files. (#908)

   - Fixed bug where the samtools fastq -i would output incorrect information
     in the Casava tags for dual-index reads.  It also now prints the tags
     for dual indices in the same way as bcl2fastq, using a '+' sign between
     the two parts of the index. (#1059; #1047 reported by Denis Loginov)

 * samtools flagstat:

   - Samtools flagstat can now optionally write its output in JSON format or
     as a tab-separated values file. (#1106, thanks to Vivek Rai).

 * samtools markdup:

   - It can optionally tag optical duplicates (reads following Illumina
     naming conventions only).  The is enabled with the '-d' option,
     which sets the distance for duplicates to be considered as optical.
     (#1091; #1103; #1121; #1128; #1134)

   - The report stats (-s) option now outputs counts for optical and
     non-primary (supplementary / secondary) duplicates.  It also reports
     the Picard "estimate library size" statistic.  A new '-f' option
     can be used to save the statistics in a given file. (#1091)

   - The rules for calling duplicates can be changed using the new --mode
     option.  This mainly changes the position associated with each read in
     a pair.  '--mode t' (the default) is the existing behaviour where the
     position used is that of the outermost template base associated with the
     read. Alternatively '--mode s' always uses the first unclipped sequence
     base.  In practice, this only makes a difference for read pairs where the
     two reads are aligned in the same direction. (#1091)

   - A new '-c' option can be used to clear any existing duplicate tags.
     (#1091)

   - A new '--include-fails' option makes markdup include QC-failed reads.
     (#1091)

   - Fixed buffer overflow in temporary file writer when writing a mixture
     of long and short alignment records. (#911; #909)

 * samtools mpileup:

   - mpileup can now process alignments including CIGAR P (pad) operators
     correctly.  They will now also produce the correct output for alignments
     where insertions are immediately followed by deletions, or deletions by
     insertions.  Note that due to limitations in HTSlib, they are still
     unable to output sequences that have been inserted before the first
     aligned base of a read. (#847; #842 reported by Tiffany Delhomme.
     See also htslib issue #59 and pull request #699).

   - In samtools mpileup, a deletion or pad on the reverse strand is now
     marked with a different character ('#') than the one used on a forward
     strand ('*'), if the '--reverse-del' option is used. (#1070)

   - New option '--output-extra' can be used to add columns for user
     selected alignment fields or aux tags. (#1073)

   - Fixed double-counting of overlapping bases in alignment records with
     deletions or reference skips longer than twice the insert size.
     (#989; #987 reported by @dariomel)

   - Improved manual page with documentation about what each output column
     means. (#1055, thanks to John Marshall)

 * samtools quickcheck:

   - Add unmapped (-u) option, which disables the check for `@SQ` lines in
     the header. (#920, thanks to Shane McCarthy)

 * samtools reheader:

   - A new option '-c' allows the input header to be passed to a given
     command.  Samtools then takes the output of this command and uses it
     as the replacement header. (#1007)

   - Make it clear in help message that reheader --in-place only works on
     CRAM files. (#921, thanks to Julian Gehring)

   - Refuse to in-place reheader BAM files, instead of unexpectedly writing
     a BAM file to stdout. (#935)

 * samtools split:

   - In samtools split, the '-u' option no longer accepts an extra file name
     from which a replacement header was read.  The two file names were
     separated using a colon, which caused problems on Windows and prevented
     the use of URLs.  A new '-h' option has been added to allow the replacement
     header file to be specified in its own option. (#961)

   - Fixed bug where samtools split would crash if it read a SAM header that
     contained an `@RG` line with no ID tag. (#954, reported by @blue-bird1)

 * samtools stats:

   - stats will now compute base compositions for BC, CR, OX and RX tags,
     and quality histograms for QT, CY, BZ and QX tags. (#904)

   - New stats FTC and LTC showing total number of nucleotides for first and
     last fragments. (#946)

   - The rules for classifying reads as "first" or "last" fragment have been
     tightened up. (#949)

   - Fixed bug where stats could over-estimate coverage when using the
     target-regions option or when a region was specified on the command-line.
     (#1027; #1025, reported by Miguel Machado; #1029, reported by Jody Phelan).

   - Fixed error in stats GCD percentile depth calculation when the depth to be
     reported fell between two bins.  It would report the depth entirely from
     the lower bin instead of taking a weighted average of the two. (#1048)

   - Better catching and reporting of out of memory conditions. (#984;
     #982, reported by Jukka Matilainen)

   - Improved manual page. (#927)

 * samtools tview:

   - tview can now display alignments including CIGAR P operators, D followed
     by I and I followed by D correctly.  See mpileup above for more
     information. (#847; #734, reported by Ryan Lorig-Roach)

   - The "go to position" text entry box has been made wider. (#968, thanks
     to John Marshall)

   - Fixed samtools tview -s option which was not filtering reads correctly.
     It now only shows reads from the requested sample or read group. (#1089)

 * samtools view:

   - New options '-d' and '-D' to only output alignments which have a tag
     with a given type and value. (#1001, thanks to Gert Hulselmans)

 * misc/plot-bamstats script:

   - Fixed merge (-m) option. (#923, #924 both thanks to Marcus D Sherman)

   - Made the quality heatmap work with gnuplot version 5.2.7 and later.
     (#1068; #1065 reported by Martin Mokrejš)

   - Fixed --do-ref-stats bug where fasta header lines would be counted
     as part of the sequence when the --targets option was used. (#1120,
     thanks to Neil Goodgame)

 * Removed the misc/varfilter.py Python script, as it takes consensus-pileup
   as input, which was removed from samtools in release 0.1.17 in 2011. (#1125)

Release 1.9 (18th July 2018)
----------------------------

 * Samtools mpileup VCF and BCF output is now deprecated.  It is still
   functional, but will warn.  Please use bcftools mpileup instead. (#884)

 * Samtools mpileup now handles the '-d' max_depth option differently.  There
   is no longer an enforced minimum, and '-d 0' is interpreted as limitless
   (no maximum - warning this may be slow).  The default per-file depth is
   now 8000, which matches the value mpileup used to use when processing
   a single sample.  To get the previous default behaviour use the higher
   of 8000 divided by the number of samples across all input files, or 250.
   (#859)

 * Samtools stats new features:

   - The '--remove-overlaps' option discounts overlapping portions of
     templates when computing coverage and mapped base counting. (#855)

   - When a target file is in use, the number of bases inside the
     target is printed and the percentage of target bases with coverage
     above a given threshold specified by the '--cov-threshold' option. (#855)

   - Split base composition and length statistics by first and last reads.
     (#814, #816)

 * Samtools faidx new features:

   - Now takes long options. (#509, thanks to Pierre Lindenbaum)

   - Now warns about zero-length and truncated sequences due to the
     requested range being beyond the end of the sequence. (#834)

   - Gets a new option (--continue) that allows it to carry on
     when a requested sequence was not in the index. (#834)

   - It is now possible to supply the list of regions to output in a text
     file using the new '--region-file' option. (#840)

   - New '-i' option to make faidx return the reverse complement of
     the regions requested. (#878)

   - faidx now works on FASTQ (returning FASTA) and added a new
     fqidx command to index and return FASTQ. (#852)

 * Samtools collate now has a fast option '-f' that only operates on
   primary pairs, dropping secondary and supplementary.  It tries to write
   pairs to the final output file as soon as both reads have been found. (#818)

 * Samtools bedcov gets a new '-j' option to make it ignore deletions (D) and
   reference skips (N) when computing coverage. (#843)

 * Small speed up to samtools coordinate sort, by converting it to use
   radix sort. (#835, thanks to Zhuravleva Aleksandra)

 * Samtools idxstats now works on SAM and CRAM files, however this
   isn't fast due to some information lacking from indices. (#832)

 * Compression levels may now be specified with the level=N
   output-fmt-option.  E.g. with -O bam,level=3.

 * Various documentation improvements.

 * Bug-fixes:

   - Improved error reporting in several places. (#827, #834, #877, cd7197)

   - Various test improvements.

   - Fixed failures in the multi-region iterator (view -M) when regions
     provided via BED files include overlaps (#819, reported by Dave Larson).

   - Samtools stats now counts '=' and 'X' CIGAR operators when
     counting mapped bases. (#855)

   - Samtools stats has fixes for insert size filtering (-m, -i). (#845; #697
     reported by Soumitra Pal)

   - Samtools stats -F now longer negates an earlier -d option. (#830)

   - Fix samtools stats crash when using a target region. (#875, reported by
     John Marshall)

   - Samtools sort now keeps to a single thread when the -@ option is absent.
     Previously it would spawn a writer thread, which could cause the CPU
     usage to go slightly over 100%. (#833, reported by Matthias Bernt)

   - Fixed samtools phase '-A' option which was incorrectly defined to take
     a parameter. (#850; #846 reported by Dianne Velasco)

   - Fixed compilation problems when using C_INCLUDE_PATH. (#870; #817 reported
     by Robert Boissy)

   - Fixed --version when built from a Git repository. (#844, thanks to
     John Marshall)

   - Use noenhanced mode for title in plot-bamstats.  Prevents unwanted
     interpretation of characters like underscore in gnuplot version 5. (#829,
     thanks to M. Zapukhlyak)

   - blast2sam.pl now reports perfect match hits (no indels or mismatches).
     (#873, thanks to Nils Homer)

   - Fixed bug in fasta and fastq subcommands where stdout would not be flushed
     correctly if the -0 option was used.

   - Fixed invalid memory access in mpileup and depth on alignment records
     where the sequence is absent.

Release 1.8 (3rd April 2018)
----------------------------

 * samtools calmd now has a quiet mode.  This can be enabled by passing `-Q` to
   calmd. (Thanks to Colin Davenport)

 * In samtools depth `-d 0` will effectively remove the depth limit. (#764)

 * Improvements made to samtools collate's interface and documentation.  It is
   now possible to specify an output file name using `-o`, instead of deriving
   it from the prefix used for temporary files.  The prefix itself is now
   optional if `-o` or `-O` (to stdout) is used. (#780)

 * Bug-fixes:

   - Make samtools addreplacerg choose output format by file extension. (#767;
     reported by Argy Megalios)

   - Merge tests now work on ungzipped data, allowing tests to be run against
     different deflate libraries.

   - samtools markdup error messages about missing tags have been updated with
     the suggestion that samtools fixmate is run beforehand. (#765; reported by
     Yudong Cai)

   - Enables the `--reference` option for samtools fastq.  Now works like other
     programs when a reference sequence is needed for CRAM files. (#791,
     reported by Milana Kaljevic)


Release 1.7 (26th January 2018)
--------------------

* HTSlib, and so samtools, now support BAMs which include CIGARs with more
  than 65535 operations as per HTS-Specs 18th November (dab57f4 and 2f915a8).

* samtools quickcheck will now write a warning to stderr if it finds
  any problems.  These messages can be suppressed with a new `-q` option.

* samtools markdup can now mark supplementary alignments of reads where
  the primary alignment is found to be a duplicate.  Supplementary marking
  can be turned on by passing the `-S` option to markdup.  When this
  option is enabled, all the alignment data will be written to a temporary
  file so that supplementary alignments that occur before a duplicated
  primary can be correctly marked in the final output.  The location
  of this temporary file can be influenced using the new `-T` option.

* samtools view now supports HTSlib's new multi-region iterator.
  This can be enabled by passing the `-M` option to view.  When using
  this option:

  - The BED filter (`-L` option) will use the index to skip through the file
  - Reads from overlapping regions will only be output once

* samtools bedcov will now ignore BED comment and header lines (#571; thanks
  to Daniel Baker).

* samtools collate now updates the `@HD` SO: and GO: tags, and sort will
  remove a GO: tag if present.  (#757; reported by Imran Haque).

* Bug-fixes:

 - maq2sam now checks for input files that end early.  (#751; patch supplied
   by Alexandre Rebert of the Mayhem team, via Andreas Tille from Debian.)

 - Fixed incorrect check when looking up header tags that could lead
   to a crash in samtools stats. (#208; thanks to Dave Larson.)

 - Fixed bug in samtools fastq `-O` option where it would fail if
   the OQ tag in the input file had an unexpected type. (#758;
   reported by Taejeong Bae)

 - The MD5 calculations in samtools dict and md5fa did not handle
   non-alphabetic characters in the same way as the CRAM MD5 function.
   They have now been updated to match. (#704; reported by Chris Norman).

 - Fix possible infinite loop in samtools targetcut.

 - Building bam_tview_curses should no longer fail if a curses header file
   cannot be found.

Release 1.6 (28th September 2017)
--------------------

* Added new markdup sub-command and '-m' option for fixmate.  Used together,
  they allow duplicates to be marked and optionally removed.  This
  fixes a number of problems with the old 'rmdup' sub-command, for
  example samtools issue #497.  'rmdup' is kept for backwards compatibility
  but 'markdup' should be used in preference.

* Sort is now much better at keeping within the requested memory limit.  It
  should also be slightly faster and need fewer temporary files when the file
  to be sorted does not fit in memory.  (#593; thanks to Nathan Weeks.)

* Sort no longer rewrites the header when merging from files.  It can also
  now merge from memory, so fewer temporary files need to be written and
  it is better at sorting in parallel when everything fits in memory.

* Both sort and merge now resolve ties when merging based on the position
  in the input file(s).  This makes them fully stable for all ordering
  options.  (Previously position sort was stable, but name and by tag
  sorts were not).

* New --output-qname option for mpileup.

* Support for building on Windows using msys2/mingw64 or cygwin has
  been improved.

Release 1.5 [Solstice Release] (21st June 2017)
--------------------

* Samtools fastq now has a -i option to create a fastq file from an index
  tag, and a -T option (similar to -t) to add user specified aux tags to
  the fastq header line.

* Samtools fastq can now create compressed fastq files, by giving the
  output filenames an extension of .gq, .bgz, or .bgzf

* Samtools sort has a -t TAG option, that allows records to be sorted by
  the value of the specified aux tag, then by position or name.  Merge
  gets a similar option, allowing files sorted this way to be merged.
  (#675; thanks to Patrick Marks of 10xgenomics).

Release 1.4.1  (8th May 2017)
-----------------------------

* Added options to fastq to create fastq files from BC (or other)
  tags.

* Samtools view has gained a -G <flags> option to exclude on all bits
  set.  For example to discard reads where neither end has been
  mapped use "-G 12".

* Samtools cat has a -b <fofn> option to ease concatenation of many
  files.

* Added misc/samtools_tab_completion for bash auto-completion of
  samtools sub-commands. (#560)

* Samtools tview now has J and K keys for verticale movement by 20
  lines. (#257)

* Various compilation / portability improvements.

* Fixed issue with more than 65536 CIGAR operations and SAM/CRAM files.
  (#667)


Release 1.4  (13 March 2017)
----------------------------

Noteworthy changes in samtools:

* Fixed Issue #345 - out-by-one error in insert-size in samtools stats

* bam_split now add a `@PG` header to the bam file

* Added mate cigar tag support to fixmate

* Multi-threading is now supported for decoding BAM and CRAM (as well
  as the previously supported encoding).  Most commands that read BAM
  or CRAM have gained an -@ or --threads arguments, providing a
  significant speed bonus.  For commands that both read and write
  files the threads are shared between decoding and encoding tasks.

* Added -a option to samtools mpileup to show all locations, including
  sites with zero depth; repeating the option as -aa or -a -a additionally
  shows reference sequences without any reads mapped to them (#496).

* The mpileup text output no longer contains empty columns at zero coverage
  positions.  Previously it would output "...0\t\t..." in some circumstances
  (zero coverage due to being below a minimum base quality); this has been
  fixed to output as "...0\t*\t*..." with placeholder '*' characters as in
  other zero coverage circumstances (see PR #537).

* To stop it from creating too many temporary files, samtools sort
  will now not run unless its per-thread memory limit (-m) is set to
  at least 1 megabyte (#547).

* The misc/plot-bamstats script now has a -l / --log-y option to change
  various graphs to display their Y axis log-scaled.  Currently this
  affects the Insert Size graph (PR #589; thanks to Anton Kratz).

* Fixmate will now also add and update MC (mate CIGAR) tags.


Beta Release 1.3.1  (22 April 2016)
-----------------------------------

Noteworthy changes in samtools:

* The sort command creates any needed temporary files alongside the final
  output file (similarly to the pre-1.3 behaviour), and now aborts when
  it detects a collision with another sort invocation's temporary files.

  When the -T PREFIX option specified is a directory (or when sorting to
  standard output), a random component is now added to temporary filenames
  to try to avoid collisions (#432, #523, #529, #535, PR #530).

* All samtools commands now check for I/O errors more carefully, especially
  when writing output files (#111, #253, #470, PR #467).

* Build fixes for 32-bit systems; be sure to run configure on such systems
  to enable large file support and access to 2GiB+ files.

* The fasta/fastq/bam2fq command no longer ignores reads when the -s option
  is used (#532).

* The fastq -O option no longer crashes on reads that do not have an OQ tag
  field (#517).

* The merge and sort commands now handle (unusual) BAM files that have no
  textual `@SQ` headers (#548, #550).

* Sorting files containing @CO headers no longer duplicates the comment
  headers, which previously happened on large sorts for which temporary
  files were needed (#563).

* The rmdup and view -l commands no longer crash on `@RG` headers that do not
  have a LB field (#538).

* Fixed miscellaneous issues #128, #130, #131, #489, and #514.


Beta Release 1.3  (15 December 2015)
------------------------------------

Noteworthy changes in samtools:

* The obsolete "samtools sort in.bam out.prefix" usage has been removed.
  If you are still using -f, -o, or out.prefix, convert to use -T PREFIX
  and/or -o FILE instead.  (#295, #349, #356, #418, PR #441; see also
  discussions in #171, #213.)

* When writing CRAM output, samtools now defaults to writing CRAM v3.0
  rather than v2.1.

* The "bamshuf" command has been renamed to "collate" (hence the term
  bamshuf no longer appears in the documentation, though it still works
  on the command line for compatibility with existing scripts).

* The mpileup command now outputs the unseen allele in VCF/BCF as <*>
  rather than X or <X> as previously, and now has AD, ADF, ADR, INFO/AD,
  INFO/ADF, INFO/ADR --output-tags annotations that largely supersede
  the existing DV, DP4, DPR annotations.

* The mpileup command now applies BAQ calculations at all base positions,
  regardless of which -l or -r options are used (previously with -l it was
  not applied to the first few tens of bases of each chromosome, leading
  to different mpileup results with -l vs. -r; #79, #125, #286, #407).

* Samtools now has a configure script which checks your build environment
  and facilitates choosing which HTSlib to build against.  See INSTALL
  for details.

* Samtools's Makefile now fully supports the standard convention of
  allowing CC/CPPFLAGS/CFLAGS/LDFLAGS/LIBS to be overridden as needed.
  Previously it listened to $(LDLIBS) instead; if you were overriding
  that, you should now override LIBS rather than LDLIBS.

* A new addreplacerg command that adds or alters `@RG` headers and RG:Z
  record tags has been added.

* The rmdup command no longer immediately aborts (previously it always
  aborted with "bam_get_library() not yet implemented"), but remains
  not recommended for most use (#159, #252, #291, #393).

* Merging files with millions of headers now completes in a reasonable
  amount of time (#337, #373, #419, #453; thanks to Nathan Weeks,
  Chris Smowton, Martin Pollard, Rob Davies).

* Samtools index's optional index output path argument works again (#199).

* Fixed calmd, targetcut, and potential mpileup segfaults when given broken
  alignments with POS far beyond the end of their reference sequences.

* If you have source code using bam_md.c's bam_fillmd1_core(), bam_cap_mapQ(),
  or bam_prob_realn_core() functions, note that these now take an additional
  ref_len parameter.  (The versions named without "_core" are unchanged.)

* The tview command's colour scheme has been altered to be more suitable
  for users with colour blindness (#457).

* Samtools depad command now handles CIGAR N operators and accepts
  CRAM files (#201, #404).

* Samtools stats now outputs separate "N" and "other" columns in the
  ACGT content per cycle section (#376).

* Added -a option to samtools depth to show all locations, including
  zero depth sites (#374).

* New samtools dict command, which creates a sequence dictionary
  (as used by Picard) from a FASTA reference file.

* Samtools stats --target-regions option works again.

* Added legacy API sam.h functions sam_index_load() and samfetch() providing
  bam_fetch()-style iteration over either BAM or CRAM files.  (In general
  we recommend recoding against the htslib API directly, but this addition
  may help existing libbam-using programs to be CRAM-enabled easily.)

* Fixed legacy API's samopen() to write headers only with "wh" when writing
  SAM files.  Plain "w" suppresses headers for SAM file output, but this
  was broken in 1.2.

* "samtools fixmate - -" works in pipelines again; with 1.0 to 1.2,
  this failed with "[bam_mating] cannot determine output format".

* Restored previous "samtools calmd -u" behaviour of writing compression
  level 0 BAM files.  Samtools 1.0 to 1.2 incorrectly wrote raw non-BGZF
  BAM files, which cannot be read by most other tools.  (Samtools commands
  other than calmd were unaffected by this bug.)

* Restored bam_nt16_nt4_table[] to legacy API header bam.h.

* Fixed bugs #269, #305, #320, #328, #346, #353, #365, #392, #410, #445,
  #462, #475, and #495.


Beta Release 1.2  (2 February 2015)
-----------------------------------

Noteworthy changes in samtools:

* Flagstat now works on SAM, BAM, or CRAM files (rather than BAM only)
* Stats calculates mismatches per cycle for unclipped length
* Merge can now merge SAM input files
* CRAM reference files are now cached by default (see HTSlib release
  notes and samtools(1) man page)
* Tested against Intel-optimised zlib (https://github.com/jtkukunas/zlib;
  see README for details)
* Fixed bugs #302, #309, #318, and #327 and many other improvements
  and bugs fixed in HTSlib -- see the HTSlib release notes


Beta Release 1.1 (19 September, 2014)
-------------------------------------

Notable changes in samtools:

 * Sorting files with thousands of reference contigs now completes in
   a reasonable amount of time (#260)
 * Fixmate and flagstat now consider supplementary reads
 * Fixmate now only adds a template cigar tag ("ct:Z") when requested
   via a new -c option, and never adds it repeatedly (#276)
 * Mpileup DPR annotation fixed (#274)
 * Checksum added to stats command output
 * Removed view -Q option


Beta Release 1.0 (15 August, 2014)
----------------------------------

First release of HTSlib-based samtools.

Numerous changes, notably support for the CRAM sequencing file format.

The faidx command now reads either uncompressed or BGZF-compressed FASTA files
compressed with bgzip.  In previous samtools-0.1.x versions, faidx could read
either uncompressed or RAZF-compressed FASTA files, but RAZF and razip are
superseded by BGZF/bgzip and have been removed from samtools.


Beta Release 0.1.20 (15 August, 2014)
-------------------------------------

Final release of standalone samtools.


Beta Release 0.1.19 (15 March, 2013)
------------------------------------

Notable changes in samtools and bcftools:

 * The latest source code and development moved to github,
    http://github.com/samtools/samtools

 * Many important bugfixes and contributions by many people. Thanks to all!

 * Performance improvements (multi-threading)

 * Important changes in calling, see
    - samtools mpileup -p
    - bcftools view -m

 * New annotations useful for filtering (RPB, HWE, QBD, MDV)

 * New tools, bamcheck and plot-bamcheck

 * New features in samtools tview

 * And much more..

For a detailed list of commits, please see
http://github.com/samtools/samtools/commits/master

(0.1.19: 15 March 2013, commit 96b5f2294ac0054230e88913c4983d548069ea4e)


Beta Release 0.1.18 (2 September, 2011)
---------------------------------------

Notable changes in samtools:

 * Support the new =/X CIGAR operators (by Peter Cock).

 * Allow to subsample BAM while keeping the pairing intact (view -s).

 * Implemented variant distance bias as a new filter (by Petr Danecek).

 * Bugfix: huge memory usage during indexing

 * Bugfix: use of uninitialized variable in mpileup (rare)

 * Bugfix: wrong BAQ probability (rare)

Notable changes in bcftools:

 * Support indel in the contrast caller.

 * Bugfix: LRT2=nan in rare cases

(0.1.18: 2 September 2011, r982:295)



Beta Release 0.1.17 (6 July, 2011)
----------------------------------

With the maturity of `mpileup` and the lack of update in the `pileup` command,
the `pileup` command is now formally dropped. Most of the pileup functionality,
such as outputting mapping quality and read positions, have been added
`mpileup`.

Since this release, `bcftools view` is able to perform contrast SNP calling
(option -T) for discovering de novo and/or somatic mutations between a pair of
samples or in a family trio. Potential mutations are scored by a log likelihood
ratio, which is very simple in math, but should be comparable to more
sophisticated methods. Note that getting the score is only the very first step.
A lot more need to be done to reduce systematical errors due to mapping and
reference errors and structural variations.

Other notable changes in samtools:

 * Improved sorting order checking during indexing.

 * Improved region parsing. Colons in reference sequence names are parsed
   properly.

 * Fixed an issue where mpileup does not apply BAQ for the first few reads when
   a region is specified.

 * Fixed an issue where `faidx` does not work with FASTA files with long lines.

 * Bugfix: wrong SP genotype information in the BCF output.

Other notable changes in bcftools:

 * Output the ML estimate of the allele count.

 * Added the HWE plus F<0 filter to varFilter. For multiple samples, it
   effectively filters false heterozygous calls around centromeres.

 * For association mapping, perform both 1-degree and 2-degree test. The
   2-degree test is conservative but more robust to HWE violation.

(0.1.17: 6 July 2011, r973:277)



Beta Release 0.1.16 (21 April, 2011)
------------------------------------

Notable changes in samtools:

 * Support the new SAM/BAM type `B` in the latest SAM spec v1.4.

 * When the output file of `samtools merge` exists, do not overwrite it unless
   a new command-line option `-f` is applied.

 * Bugfix: BED support is not working when the input BED is not sorted.

 * Bugfix: some reads without coordinates but given on the reverse strand are
   lost in merging.

Notable changes in bcftools:

 * Code cleanup: separated max-likelihood inference and Bayesian inference.

 * Test Hardy-Weinberg equilibrium with a likelihood-ratio test.

 * Provided another association test P-value by likelihood-ratio test.

 * Use Brent's method to estimate the site allele frequency when EM converges
   slowly. The resulting ML estimate of allele frequnecy is more accurate.

 * Added the `ldpair` command, which computes r^2 between SNP pairs given in
   an input file.

Also, the `pileup` command, which has been deprecated by `mpileup` since
version 0.1.10, will be dropped in the next release. The old `pileup` command
is substandard and causing a lot of confusion.

(0.1.16: 21 April 2011, r963:234)



Beta Release 0.1.15 (10 April, 2011)
------------------------------------

Notable changes:

 * Allow to perform variant calling or to extract information in multiple
   regions specified by a BED file (`samtools mpileup -l`, `samtools view -L`
   and `bcftools view -l`).

 * Added the `depth` command to samtools to compute the per-base depth with a
   simpler interface. File `bam2depth.c`, which implements this command, is the
   recommended example on how to use the mpileup APIs.

 * Estimate genotype frequencies with ML; perform chi^2 based Hardy-Weinberg
   test using this estimate.

 * For `samtools view`, when `-R` is specified, drop read groups in the header
   that are not contained in the specified file.

 * For `samtools flagstat`, separate QC-pass and QC-fail reads.

 * Improved the command line help of `samtools mpileup` and `bcftools view`.

 * Use a global variable to control the verbose level of samtools stderr
   output. Nonetheless, it has not been full utilized.

 * Fixed an issue in association test which may report false associations,
   possibly due to floating point underflow.

(0.1.15: 10 April 2011, r949:203)



Beta release 0.1.14 (21 March, 2011)
------------------------------------

This release implements a method for testing associations for case-control
data. The method does not call genotypes but instead sums over all genotype
configurations to compute a chi^2 based test statistics. It can be potentially
applied to comparing a pair of samples (e.g. a tumor-normal pair), but this
has not been evaluated on real data.

Another new feature is to make X chromosome variant calls when female and male
samples are both present. The user needs to provide a file indicating the
ploidy of each sample (see also manual bcftools/bcftools.1).

Other notable changes:

 * Added `bcftools view -F` to parse BCF files generated by samtools r921 or
   older which encodes PL in a different way.

 * Changed the behavior of `bcftools view -s`. Now when a list of samples is
   provided, the samples in the output will be reordered to match the ordering
   in the sample list. This change is mainly designed for association test.

 * Sped up `bcftools view -v` for target sequencing given thousands of samples.
   Also added a new option `view -d` to skip loci where only a few samples are
   covered by reads.

 * Dropped HWE test. This feature has never been implemented properly. An EM
   should be much better. To be implemented in future.

 * Added the `cat` command to samtools. This command concatenate BAMs with
   identical sequence dictionaries in an efficient way. Modified from bam_cat.c
   written by Chris Saunders.

 * Added `samtools view -1` to write BAMs at a low compression level but twice
   faster to create. The `sort` command generates temporary files at a low
   compression level as well.

 * Added `samtools mpileup -6` to accept "BAM" with Illumina 1.3+ quality
   strings (strictly speaking, such a file is not BAM).

 * Added `samtools mpileup -L` to skip INDEL calling in regions with
   excessively high coverage. Such regions dramatically slow down mpileup.

 * Updated `misc/export2sam.pl`, provided by Chris Saunders from Illumina Inc.

(0.1.14: 21 March 2011, r933:170)



Beta release 0.1.13 (1 March, 2011)
-----------------------------------

The most important though largely invisible modification is the change of the
order of genotypes in the PL VCF/BCF tag. This is to conform the upcoming VCF
spec v4.1. The change means that 0.1.13 is not backward compatible with VCF/BCF
generated by samtools older than r921 inclusive.  VCF/BCF generated by the new
samtools will contain a line `##fileformat=VCFv4.1` as well as the samtools
version number.

Single Individual Haplotyping (SIH) is added as an experimental feature. It
originally aims to produce haploid consensus from fosmid pool sequencing, but
also works with short-read data. For short reads, phased blocks are usually too
short to be useful in many applications, but they can help to rule out part of
SNPs close to INDELs or between copies of CNVs.


Other notable changes in samtools:

 * Construct per-sample consensus to reduce the effect of nearby SNPs in INDEL
   calling. This reduces the power but improves specificity.

 * Improved sorting order checking in indexing. Now indexing is the preferred way
   to check if a BAM is sorted.

 * Added a switch `-E` to mpileup and calmd. This option uses an alternative way
   to apply BAQ, which increases sensistivity, especially to MNPs, at the cost of
   a little loss in specificity.

 * Added `mpileup -A` to allow to use reads in anomalous pairs in SNP calling.

 * Added `mpileup -m` to allow fine control of the collection of INDEL candidates.

 * Added `mpileup -S` to compute per-sample strand bias P-value.

 * Added `mpileup -G` to exclude read groups in variant calling.

 * Fixed segfault in indel calling related to unmapped and refskip reads.

 * Fixed an integer overflow in INDEL calling. This bug produces wrong INDEL
   genotypes for longer short INDELs, typically over 10bp.

 * Fixed a bug in tview on big-endian machines.

 * Fixed a very rare memory issue in bam_md.c

 * Fixed an out-of-boundary bug in mpileup when the read base is `N`.

 * Fixed a compiling error when the knetfile library is not used. Fixed a
   library compiling error due to the lack of bam_nt16_nt4_table[] table.
   Suppress a compiling warning related to the latest zlib.


Other notable changes in bcftools:

 * Updated the BCF spec.

 * Added the `FQ` VCF INFO field, which gives the phred-scaled probability
   of all samples being the same (identical to the reference or all homozygous
   variants). Option `view -f` has been dropped.

 * Implemented "vcfutils.pl vcf2fq" to generate a consensus sequence
   similar to "samtools.pl pileup2fq".

 * Make sure the GT FORMAT field is always the first FORMAT to conform the VCF
   spec. Drop bcf-fix.pl.

 * Output bcftools specific INFO and FORMAT in the VCF header.

 * Added `view -s` to call variants from a subset of samples.

 * Properly convert VCF to BCF with a user provided sequence dictionary. Nonetheless,
   custom fields are still unparsed and will be stored as a missing value.

 * Fixed a minor bug in Fisher's exact test; the results are rarely changed.


(0.1.13: 1 March 2011, r926:134)



Beta release 0.1.12a (2 December, 2010)
---------------------------------------

This is another bug fix release:

 * Fixed a memory violation in mpileup, which causes segfault. Release
   0.1.9 and above are affected.

 * Fixed a memory violation in the indel caller, which does not causes
   segfault, but may potentially affect deletion calls in an unexpected
   way. Release 0.1.10 and above are affected.

 * Fixed a bug in computing r-square in bcftools. Few are using this
   functionality and it only has minor effect.

 * Fixed a memory leak in bam_fetch().

 * Fixed a bug in writing meta information to the BAM index for the last
   sequence. This bug is invisible to most users, but it is a bug anyway.

 * Fixed a bug in bcftools which causes false "DP4=0,0,0,0" annotations.

(0.1.12: 2 December 2010, r862)



Beta release 0.1.11 (21 November, 2010)
---------------------------------------

This is mainly a bug fix release:

 * Fixed a bug in random retrieval (since 0.1.8). It occurs when reads
   are retrieved from a small region containing no reads.

 * Fixed a bug in pileup (since 0.1.9). The bug causes an assertion
   failure when the first CIGAR operation is a deletion.

 * Improved fault tolerance in remote access.

One minor feature has been implemented in bcftools:

 * Added a reference-free variant calling mode. In this mode, a site is
   regarded as a variat iff the sample(s) contains two or more alleles;
   the meaning of the QUAL field in the VCF output is changed
   accordingly. Effectively, the reference allele is irrelevant to the
   result in the new mode, although the reference sequence has to be
   used in realignment when SAMtools computes genotype likelihoods.

In addition, since 0.1.10, the `pileup` command has been deprecated by
`mpileup` which is more powerful and more accurate. The `pileup` command
will not be removed in the next few releases, but new features will not
be added.

(0.1.11: 21 November 2010, r851)



Beta Release 0.1.10 (16 November, 2010)
---------------------------------------

This release is featured as the first major improvement to the indel
caller. The method is similar to the old one implemented in the pileup
command, but the details are handled more carefully both in theory and
in practice. As a result, the new indel caller usually gives more
accurate indel calls, though at the cost of sensitivity. The caller is
implemented in the mpileup command and is invoked by default. It works
with multiple samples.

Other notable changes:

 * With the -r option, the calmd command writes the difference between
   the original base quality and the BAQ capped base quality at the BQ
   tag but does not modify the base quality. Please use -Ar to overwrite
   the original base quality (the 0.1.9 behavior).

 * Allow to set a maximum per-sample read depth to reduce memory. In
   0.1.9, most of memory is wasted for the ultra high read depth in some
   regions (e.g. the chr1 centromere).

 * Optionally write per-sample read depth and per-sample strand bias
   P-value.

 * Compute equal-tail (Bayesian) credible interval of site allele
   frequency at the CI95 VCF annotation.

 * Merged the vcfutils.pl varFilter and filter4vcf for better SNP/indel
   filtering.

(0.1.10: 16 November 2010, r829)



Beta Release 0.1.9 (27 October, 2010)
-------------------------------------

This release is featured as the first major improvement to the samtools'
SNP caller.  It comes with a revised MAQ error model, the support of
multi-sample SNP calling and the computation of base alignment quality
(BAQ).

The revised MAQ error model is based on the original model. It solves an
issue of miscalling SNPs in repetitive regions. Although such SNPs can
usually be filtered at a later step, they mess up unfiltered calls. This
is a theoretical flaw in the original model. The revised MAQ model
deprecates the original MAQ model and the simplified SOAPsnp model.

Multi-sample SNP calling is separated in two steps. The first is done by
samtools mpileup and the second by a new program, bcftools, which is
included in the samtools source code tree. Multi-sample SNP calling also
works for single sample and has the advantage of enabling more powerful
filtration. It is likely to deprecate pileup in future once a proper
indel calling method is implemented.

BAQ is the Phred-scaled probability of a read base being wrongly
aligned. Capping base quality by BAQ has been shown to be very effective
in suppressing false SNPs caused by misalignments around indels or in
low-complexity regions with acceptable compromise on computation
time. This strategy is highly recommended and can be used with other SNP
callers as well.

In addition to the three major improvements, other notable changes are:

 * Changes to the pileup format. A reference skip (the N CIGAR operator)
   is shown as '<' or '>' depending on the strand. Tview is also changed
   accordingly.

 * Accelerated pileup. The plain pileup is about 50% faster.

 * Regional merge. The merge command now accepts a new option to merge
   files in a specified region.

 * Fixed a bug in bgzip and razip which causes source files to be
   deleted even if option -c is applied.

 * In APIs, propagate errors to downstream callers and make samtools
   return non-zero values once errors occur.

(0.1.9: 27 October 2010, r783)



Beta Release 0.1.8 (11 July, 2010)
----------------------------------

Notable functional changes:

 * Added the `reheader` command which replaces a BAM header with a new
   header. This command is much faster than replacing header by
   BAM->SAM->BAM conversions.

 * Added the `mpileup` command which computes the pileup of multiple
   alignments.

 * The `index` command now stores the number of mapped and unmapped
   reads in the index file. This information can be retrieved quickly by
   the new `idxstats` command.

 * By default, pileup used the SOAPsnp model for SNP calling. This
   avoids the floating overflow in the MAQ model which leads to spurious
   calls in repetitive regions, although these calls will be immediately
   filtered by varFilter.

 * The `tview` command now correctly handles CIGARs like 7I10M and
   10M1P1I10M which cause assertion failure in earlier versions.

 * Tview accepts a region like `=10,000` where `=` stands for the
   current sequence name. This saves typing for long sequence names.

 * Added the `-d` option to `pileup` which avoids slow indel calling
   in ultradeep regions by subsampling reads locally.

 * Added the `-R` option to `view` which retrieves alignments in read
   groups listed in the specified file.

Performance improvements:

 * The BAM->SAM conversion is up to twice faster, depending on the
   characteristic of the input.

 * Parsing SAM headers with a lot of reference sequences is now much
   faster.

 * The number of lseek() calls per query is reduced when the query
   region contains no read alignments.

Bug fixes:

 * Fixed an issue in the indel caller that leads to miscall of indels.
   Note that this solution may not work well when the sequencing indel
   error rate is higher than the rate of SNPs.

 * Fixed another issue in the indel caller which may lead to incorrect
   genotype.

 * Fixed a bug in `sort` when option `-o` is applied.

 * Fixed a bug in `view -r`.

APIs and other changes:

 * Added iterator interfaces to random access and pileup. The callback
   interfaces directly call the iterator interfaces.

 * The BGZF blocks holding the BAM header are independent of alignment
   BGZF blocks. Alignment records shorter than 64kB is guaranteed to be
   fully contained in one BGZF block. This change is fully compatible
   with the old version of samtools/picard.

Changes in other utilities:

 * Updated export2sam.pl by Chris Saunders.

 * Improved the sam2vcf.pl script.

 * Added a Python version of varfilter.py by Aylwyn Scally.

(0.1.8: 11 July 2010, r613)



Beta Release 0.1.7 (10 November, 2009)
--------------------------------------

Notable changes:

 * Improved the indel caller in complex scenariors, in particular for
   long reads. The indel caller is now able to make reasonable indel
   calls from Craig Venter capillary reads.

 * Rewrote single-end duplicate removal with improved
   performance. Paired-end reads are not touched.

 * Duplicate removal is now library aware. Samtools remove potential
   PCR/optical duplicates inside a library rather than across libraries.

 * SAM header is now fully parsed, although this functionality is not
   used in merging and so on.

 * In samtools merge, optionally take the input file name as RG-ID and
   attach the RG tag to each alignment.

 * Added FTP support in the RAZF library. RAZF-compressed reference
   sequence can be retrieved remotely.

 * Improved network support for Win32.

 * Samtools sort and merge are now stable.

Changes in other utilities:

 * Implemented sam2vcf.pl that converts the pileup format to the VCF
   format.

 * This release of samtools is known to work with the latest
   Bio-Samtools Perl module.

(0.1.7: 10 November 2009, r510)



Beta Release 0.1.6 (2 September, 2009)
--------------------------------------

Notable changes:

 * In tview, do not show a blank screen when no reads mapped to the
   corresponding region.

 * Implemented native HTTP support in the BGZF library. Samtools is now
   able to directly open a BAM file on HTTP. HTTP proxy is also
   supported via the "http_proxy" environmental variable.

 * Samtools is now compatible with the MinGW (win32) compiler and the
   PDCurses library.

 * The calmd (or fillmd) command now calculates the NM tag and replaces
   MD tags if they are wrong.

 * The view command now recognizes and optionally prints FLAG in HEXs or
   strings to make a SAM file more friendly to human eyes. This is a
   samtools-C extension, not implemented in Picard for the time
   being. Please type `samtools view -?` for more information.

 * BAM files now have an end-of-file (EOF) marker to facilitate
   truncation detection. A warning will be given if an on-disk BAM file
   does not have this marker. The warning will be seen on BAM files
   generated by an older version of samtools. It does NO harm.

 * New key bindings in tview: 'r' to show read names and 's' to show
   reference skip (N operation) as deletions.

 * Fixed a bug in `samtools merge -n`.

 * Samtools merge now optionally copies the header of a user specified
   SAM file to the resultant BAM output.

 * Samtools pileup/tview works with a CIGAR with the first or the last
   operation is an indel.

 * Fixed a bug in bam_aux_get().


Changes in other utilities:

 * Fixed wrong FLAG in maq2sam.


(0.1.6: 2 September 2009, r453)



Beta Release 0.1.5 (7 July, 2009)
---------------------------------

Notable changes:

 * Support opening a BAM alignment on FTP. Users can now use "tview" to
   view alignments at the NCBI ftp site. Please read manual for more
   information.

 * In library, propagate errors rather than exit or complain assertion
   failure.

 * Simplified the building system and fixed compiling errors caused by
   zlib<1.2.2.1.

 * Fixed an issue about lost header information when a SAM is imported
   with "view -t".

 * Implemented "samtool.pl varFilter" which filters both SNPs and short
   indels. This command replaces "indelFilter".

 * Implemented "samtools.pl pileup2fq" to generate FASTQ consensus from
   pileup output.

 * In pileup, cap mapping quality at 60. This helps filtering when
   different aligners are in use.

 * In pileup, allow to output variant sites only.

 * Made pileup generate correct calls in repetitive region. At the same
   time, I am considering to implement a simplified model in SOAPsnp,
   although this has not happened yet.

 * In view, added '-u' option to output BAM without compression. This
   option is preferred when the output is piped to other commands.

 * In view, added '-l' and '-r' to get the alignments for one library or
   read group. The `@RG` header lines are now partially parsed.

 * Do not include command line utilities to libbam.a.

 * Fixed memory leaks in pileup and bam_view1().

 * Made faidx more tolerant to empty lines right before or after FASTA >
   lines.


Changes in other utilities:

 * Updated novo2sam.pl by Colin Hercus, the key developer of novoalign.


This release involves several modifications to the key code base which
may potentially introduce new bugs even though we have tried to minimize
this by testing on several examples. Please let us know if you catch
bugs.

(0.1.5: 7 July 2009, r373)



Beta Release 0.1.4 (21 May, 2009)
---------------------------------

Notable changes:

 * Added the 'rmdupse' command: removing duplicates for SE reads.

 * Fixed a critical bug in the indel caller: clipped alignments are not
   processed correctly.

 * Fixed a bug in the tview: gapped alignment may be incorrectly
   displayed.

 * Unified the interface to BAM and SAM I/O. This is done by
   implementing a wrapper on top of the old APIs and therefore old APIs
   are still valid. The new I/O APIs also recognize the `@SQ` header
   lines.

 * Generate the MD tag.

 * Generate "=" bases. However, the indel caller will not work when "="
   bases are present.

 * Enhanced support of color-read display (by Nils Homer).

 * Implemented the GNU building system. However, currently the building
   system does not generate libbam.a. We will improve this later. For
   the time being, `make -f Makefile.generic` is preferred.

 * Fixed a minor bug in pileup: the first read in a chromosome may be
   skipped.

 * Fixed bugs in bam_aux.c. These bugs do not affect other components as
   they were not used previously.

 * Output the 'SM' tag from maq2sam.

(0.1.4: 21 May 2009, r297)



Beta Release 0.1.3 (15 April, 2009)
-----------------------------------

Notable changes in SAMtools:

 * SAMtools is more consistent with the specification: a) '*' in the
   QUAL field is allowed; b) the field separator is TAB only and SPACE
   is treated as a character in a field; c) empty header is allowed.

 * Implemented GLFv3 support in pileup.

 * Fixed a severe bug in fixmate: strand information is wrongly
   overwritten.

 * Fixed a bug in alignment retrieval: alignments bridging n*16384bp are
   not correctly retrieved sometimes.

 * Fixed a bug in rmdup: segfault if unmapped reads are present.

 * Move indel_filter.pl to samtools.pl and improved the filtering by
   checking the actual number of alignments containing indels. The indel
   pileup line is also changed a little to make this filtration easier.

 * Fixed a minor bug in indexing: the bin number of an unmapped read is
   wrongly calculated.

 * Added `flagstat` command to show statistics on the FLAG field.

 * Improved indel caller by setting the maximum window size in local
   realignment.

Changes in other utilities:

 * Fixed a bug in maq2sam: a tag name is obsolete.

 * Improvement to wgsim: a) added support for SOLiD read simulation; b)
   show the number of substitutions/indels/errors in read name; c)
   considerable code clean up.

 * Various converters: improved functionality in general.

 * Updated the example SAM due to the previous bug in fixmate.

(0.1.3: 15 April 2009, r227)



Beta Release 0.1.2 (28 January, 2008)
-------------------------------------

Notable changes in SAMtools:

 * Implemented a Bayesian indel caller. The new caller generate scores
   and genotype and is potentially more accurate than Maq's indel
   caller. The pileup format is also changed accordingly.

 * Implemented rmdup command: remove potential PCR duplicates. Note that
   this command ONLY works for FR orientation and requires ISIZE is
   correctly set.

 * Added fixmate command: fill in mate coordinates, ISIZE and mate
   related flags from a name-sorted alignment.

 * Fixed a bug in indexing: reads bridging 16x kbp were not retrieved.

 * Allow to select reads shown in the pileup output with a mask.

 * Generate GLFv2 from pileup.

 * Added two more flags for flagging PCR/optical duplicates and for QC
   failure.

 * Fixed a bug in sort command: name sorting for large alignment did not
   work.

 * Allow to completely disable RAZF (using Makefile.lite) as some people
   have problem to compile it.

 * Fixed a bug in import command when there are reads without
   coordinates.

 * Fixed a bug in tview: clipping broke the alignment viewer.

 * Fixed a compiling error when _NO_CURSES is applied.

 * Fixed a bug in merge command.

Changes in other utilities:

 * Added wgsim, a paired-end reads simulator. Wgsim was adapted from
   maq's reads simulator. Colin Hercus further improved it to allow
   longer indels.

 * Added wgsim_eval.pl, a script that evaluates the accuracy of
   alignment on reads generated by wgsim.

 * Added soap2sam.pl, a SOAP2->SAM converter. This converter does not
   work properly when multiple hits are output.

 * Added bowtie2sam.pl, a Bowtie->SAM converter. Only the top hit will
   be retained when multiple hits are present.

 * Fixed a bug in export2sam.pl for QC reads.

 * Support RG tag at MAQ->SAM converter.

 * Added novo2sam.pl, a NovoAlign->SAM converter. Multiple hits and
   indel are not properly handled, though.

 * Added zoom2sam.pl, a ZOOM->SAM converter. It only works with the
   default Illumina output.

(0.1.2: 28 January 2008; r116)



Beta Release 0.1.1 (22 December, 2008)
--------------------------------------

The is the first public release of samtools. For more information,
please check the manual page `samtools.1` and the samtools website
http://samtools.sourceforge.net
