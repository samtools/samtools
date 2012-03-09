#!/usr/bin/perl -w

# Contact: lh3
# Version: 0.1.1

use strict;
use warnings;
use Getopt::Std;

&bowtie2sam;
exit;

sub bowtie2sam {
  my %opts = ();
  die("Usage: bowtie2sam.pl <aln.bowtie>\n") if (@ARGV == 0 && -t STDIN);
  # core loop
  my (@s, $last, @staging, $k, $best_s, $subbest_s, $best_k);
  $last = '';
  while (<>) {
	my ($name, $nm) = &bowtie2sam_aux($_, \@s); # read_name, number of mismatches
	if ($name eq $last) {
	  # I do not know whether the multiple hits are ordered on the
	  # number of mismatches. I assume they are not and so I have to
	  # keep all these multiple hits in memory.
	  @{$staging[$k]} = @s;
	  if ($best_s > $nm) {
		$subbest_s = $best_s;
		$best_s = $nm;
		$best_k = $k;
	  } elsif ($subbest_s > $nm) {
		$subbest_s = $nm;
	  }
	  ++$k;
	} else {
	  if ($last) {
		if ($best_s == $subbest_s) {
		  $staging[$best_k][4] = 0;
		} elsif ($subbest_s - $best_s == 1) {
		  $staging[$best_k][4] = 15 if ($staging[$best_k][4] > 15);
		}
		print join("\t", @{$staging[$best_k]}), "\n";
	  }
	  $k = 1; $best_s = $nm; $subbest_s = 1000; $best_k = 0;
	  @{$staging[0]} = @s;
	  $last = $name;
	}
  }
  print join("\t", @{$staging[$best_k]}), "\n" if ($best_k >= 0);
}

sub bowtie2sam_aux {
  my ($line, $s) = @_;
  chomp($line);
  my @t = split("\t", $line);
  my $ret;
  @$s = ();
  # read name
  $s->[0] = $ret = $t[0];
  $s->[0] =~ s/\/[12]$//g;
  # initial flag (will be updated later)
  $s->[1] = 0;
  # read & quality
  $s->[9] = $t[4]; $s->[10] = $t[5];
  # cigar
  $s->[5] = length($s->[9]) . "M";
  # coor
  $s->[2] = $t[2]; $s->[3] = $t[3] + 1;
  $s->[1] |= 0x10 if ($t[1] eq '-');
  # mapQ
  $s->[4] = $t[6] == 0? 25 : 0;
  # mate coordinate
  $s->[6] = '*'; $s->[7] = $s->[8] = 0;
  # aux
  my $nm = @t - 7;
  push(@$s, "NM:i:" . (@t-7));
  push(@$s, "X$nm:i:" . ($t[6]+1));
  my $md = '';
  if ($t[7]) {
	$_ = $t[7];
	my $a = 0;
	while (/(\d+):[ACGTN]>([ACGTN])/gi) {
	  my ($y, $z) = ($1, $2);
	  $md .= (int($y)-$a) . $z;
	  $a += $y - $a + 1;
	}
	$md .= length($s->[9]) - $a;
  } else {
	$md = length($s->[9]);
  }
  push(@$s, "MD:Z:$md");
  return ($ret, $nm);
}
