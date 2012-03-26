#!/usr/bin/env perl

use strict;
use warnings;

# Author: lh3
#
# This script is literally translated from the C version. It has two funtionalities:
#
# a) compute the length of the reference sequence contained in an alignment;
# b) collapse backward overlaps and generate consensus sequence and quality
#
# During the consensus construction, if two bases from two overlapping segments agree,
# the base quality is taken as the higher one of the two; if the two bases disagree,
# the base is set to the one of higher quality and the quality set to the difference
# between the two qualities.
#
# There are several special cases or errors:
#
# a) If a backward operation goes beyond the beginning of SEQ, the read is regarded to
#    be unmapped.
# b) If the CIGARs of two segments in an overlap are inconsistent (e.g. 10M3B1M1I8M)
#    the consensus CIGAR is taken as the one from the latter.
# c) If three or more segments overlap, the consensus SEQ/QUAL will be computed firstly
#    for the first two overlapping segments, and then between the two-segment consensus
#    and the 3rd segment and so on. The consensus CIGAR is always taken from the last one.

die("Usage: removeB.pl <in.sam>\n") if (@ARGV == 0 && -t STDIN);
while (<>) {
	if (/^\@/) { # do not process the header lines
		print;
		next;
	}
	my $failed = 0; # set to '1' in case of inconsistent CIGAR
	my @t = split;
	$t[5] =~ s/\d+B$//; # trim trailing 'B'
	my @cigar; # this is the old CIGAR array
	my $no_qual = ($t[10] eq '*'); # whether QUAL equals '*'

	####################################################
	# Compute the length of reference in the alignment #
	####################################################

	my $alen = 0; # length of reference in the alignment
	while ($t[5] =~ m/(\d+)([A-Z])/g) { # loop through the CIGAR string
		my ($op, $len) = ($2, $1);
		if ($op eq 'B') { # a backward operation
			my ($u, $v) = (0, 0); # $u: query length during backward lookup; $v: reference length
			my $l;
			for ($l = $#cigar; $l >= 0; --$l) { # look back
				my ($op1, $len1) = @{$cigar[$l]};
				if ($op1 =~ m/[MIS]/) { # consume the query sequence
					if ($u + $len1 >= $len) { # we have moved back enough; stop
						$v += $len - $u if ($op1 =~ m/[MDN]/);
						last;
					} else { $u += $len1; }
				}
				$v += $len1 if ($op1 =~ m/[MDN]/);
			}
			$alen = $l < 0? 0 : $alen - $v;
		} elsif ($op =~ m/[MDN]/) { # consume the reference sequence
			$alen += $len;
		}
		push(@cigar, [$op, $len]); # keep it in the @cigar array
	}
	push(@t, "XL:i:$alen"); # write $alen in a tag
	goto endloop if ($t[5] !~ /B/); # do nothing if the CIGAR does not contain 'B'

	##############################
	# Collapse backward overlaps #
	##############################

	$t[10] = '!' x length($t[9]) if $t[10] eq '*'; # when no QUAL; set all qualities to zero
	# $i: length of query that has been read; $j: length of query that has been written
	# $end_j: the right-most query position; $j may be less than $end_j after a backward operation
	my ($k, $i, $j, $end_j) = (0, 0, 0, -1);
	my @cigar2; # the new CIGAR array will be kept here; the new SEQ/QUAL will be generated in place
	for ($k = 0; $k < @cigar; ++$k) {
		my ($op, $len) = @{$cigar[$k]}; # the CIGAR operation and the operation length
		if ($op eq 'B') { # a backward operation
			my ($t, $u);  # $u: query length during backward loopup
			goto endloop if $len > $j; # an excessively long backward operation
			for ($t = $#cigar2, $u = 0; $t >= 0; --$t) { # look back along the new cigar
				my ($op1, $len1) = @{$cigar2[$t]};
				if ($op1 =~ m/[MIS]/) { # consume the query sequence
					if ($u + $len1 >= $len) { # we have looked back enough; stop
						$cigar2[$t][1] -= $len - $u;
						last;
					} else { $u += $len1; }
				}
			}
			--$t if $cigar2[$t][1] == 0; # get rid of the zero-length operation
			$#cigar2 = $t; # truncate the new CIGAR array
			$end_j = $j;
			$j -= $len;
		} else {
			push(@cigar2, $cigar[$k]); # copy the old CIGAR to the new one
			if ($op =~ m/[MIS]/) { # consume the query sequence
				my ($u, $c);
				# note that SEQ and QUAL is generated in place (i.e. overwriting the old SEQ/QUAL)
				for ($u = 0; $u < $len; ++$u) {
					$c = substr($t[9], $i + $u, 1); # the base on latter segment
					if ($j + $u < $end_j) { # we are in an backward overlap
						# for non-Perl programmers: ord() takes the ASCII of a character; chr() gets the character of an ASCII
						if ($c ne substr($t[9], $j + $u, 1)) { # a mismatch in the overlaps
							if (ord(substr($t[10], $j + $u, 1)) < ord(substr($t[10], $i + $u, 1))) { # QUAL of the 2nd segment is better
								substr($t[9], $j + $u, 1) = $c; # the consensus base is taken from the 2nd segment
								substr($t[10],$j + $u, 1) = chr(ord(substr($t[10], $i + $u, 1)) - ord(substr($t[10], $j + $u, 1)) + 33);
							} else { # then do not change the base, but reduce the quality
								substr($t[10],$j + $u, 1) = chr(ord(substr($t[10], $j + $u, 1)) - ord(substr($t[10], $i + $u, 1)) + 33);
							}
						} else { # same base; then set the quality as the higher one
							substr($t[10],$j + $u, 1) = ord(substr($t[10], $j + $u, 1)) > ord(substr($t[10], $i + $u, 1))?
														substr($t[10], $j + $u, 1) : substr($t[10], $i + $u, 1);
						}
					} else { # not in an overlap; then copy the base and quality over
						substr($t[9], $j + $u, 1) = $c;
						substr($t[10],$j + $u, 1) = substr($t[10], $i + $u, 1);
					}
				}
				$i += $len; $j += $len;
			}
		}
	}
	# merge adjacent CIGAR operations of the same type
	for ($k = 1; $k < @cigar2; ++$k) {
		if ($cigar2[$k][0] eq $cigar2[$k-1][0]) {
			$cigar2[$k][1] += $cigar2[$k-1][1];
			$cigar2[$k-1][1] = 0; # set the operation length to zero
		}
	}
	# update SEQ, QUAL and CIGAR
	$t[9] = substr($t[9], 0, $j); # SEQ
	$t[10]= substr($t[10],0, $j); # QUAL
	$t[5] = ''; # CIGAR
	for my $p (@cigar2) {
		$t[5] .= "$p->[1]$p->[0]" if ($p->[1]); # skip zero-length operations
	}

	#########
	# Print #
	#########

endloop:
	$t[1] |= 4 if $failed; # mark the read as "UNMAPPED" if something bad happens
	$t[10] = '*' if $no_qual;
	print join("\t", @t), "\n";
}
