#!/usr/bin/perl -w

# Contact: lh3
# Version: 0.1.0

use strict;
use warnings;
use Getopt::Std;

&novo2sam;
exit;

sub mating {
  my ($s1, $s2) = @_;
  my $isize = 0;
  if ($s1->[2] ne '*' && $s1->[2] eq $s2->[2]) { # then calculate $isize
	my $x1 = ($s1->[1] & 0x10)? $s1->[3] + length($s1->[9]) : $s1->[3];
	my $x2 = ($s2->[1] & 0x10)? $s2->[3] + length($s2->[9]) : $s2->[3];
	$isize = $x2 - $x1;
  }
  # update mate coordinate
  if ($s2->[2] ne '*') {
	@$s1[6..8] = (($s2->[2] eq $s1->[2])? "=" : $s2->[2], $s2->[3], $isize);
	$s1->[1] |= 0x20 if ($s2->[1] & 0x10);
  } else {
	$s1->[1] |= 0x8;
  }
  if ($s1->[2] ne '*') {
	@$s2[6..8] = (($s1->[2] eq $s2->[2])? "=" : $s1->[2], $s1->[3], -$isize);
	$s2->[1] |= 0x20 if ($s1->[1] & 0x10);
  } else {
	$s2->[1] |= 0x8;
  }
}

sub novo2sam {
  my %opts = ();
  getopts("p", \%opts);
  die("Usage: novo2sam.pl [-p] <aln.novo>\nWarning: gapped alignments are NOT converted properly!\n") if (@ARGV == 0);
  my $is_paired = defined($opts{p});
  # core loop
  my @s1 = ();
  my @s2 = ();
  my ($s_last, $s_curr) = (\@s1, \@s2);
  while (<>) {
	next if (/^#/);
	next if (/(QC|NM)\s*$/ || /(R\s+\d+)\s*$/);
	&novo2sam_aux($_, $s_curr, $is_paired);
	if (@$s_last != 0 && $s_last->[0] eq $s_curr->[0]) {
	  &mating($s_last, $s_curr);
	  print join("\t", @$s_last), "\n";
	  print join("\t", @$s_curr), "\n";
	  @$s_last = (); @$s_curr = ();
	} else {
	  print join("\t", @$s_last), "\n" if (@$s_last != 0);
	  my $s = $s_last; $s_last = $s_curr; $s_curr = $s;
	}
  }
  print join("\t", @$s_last), "\n" if (@$s_last != 0);
}

sub novo2sam_aux {
  my ($line, $s, $is_paired) = @_;
  chomp($line);
  my @t = split(/\s+/, $line);
  @$s = ();
  return if ($t[4] ne 'U');
  my $len = length($t[2]);
  # read name
  $s->[0] = substr($t[0], 1);
  $s->[0] =~ s/\/[12]$//g;
  # initial flag (will be updated later)
  $s->[1] = 0;
  $s->[1] |= 1 | 1<<($t[1] eq 'L'? 6 : 7);
  $s->[1] |= 2 if ($t[10] eq '.');
  # read & quality
  if ($t[9] eq 'R') {
	$s->[9] = reverse($t[2]);
	$s->[10] = reverse($t[3]);
	$s->[9] =~ tr/ACGTRYMKWSNacgtrymkwsn/TGCAYRKMWSNtgcayrkmwsn/;
  } else {
	$s->[9] = $t[2]; $s->[10] = $t[3];
  }
  # cigar
  $s->[5] = $len . "M"; # IMPORTANT: this cigar is not correct for gapped alignment
  # coor
  $s->[2] = substr($t[7], 1); $s->[3] = $t[8];
  $s->[1] |= 0x10 if ($t[9] eq 'R');
  # mapQ
  $s->[4] = $t[5] > $t[6]? $t[5] : $t[6];
  # mate coordinate
  $s->[6] = '*'; $s->[7] = $s->[8] = 0;
  # aux
  push(@$s, "NM:i:".(@t-13));
  my $md = '';
  if ($t[13]) {
	my @x;
	for (13 .. $#t) {
	  push(@x, sprintf("%.4d,$2", $1-1)) if ($t[$_] =~ /^(\d+)([ACGT])>([ACGT])/i);
	}
	#@x = sort(@x);
	my $a = 0;
	for (@x) {
	  my ($y, $z) = split(",");
	  $md .= (int($y)-$a? (int($y)-$a) : '') . $z;
	  $a += $y - $a + 1;
	}
	$md .= $len - $a if ($len - $a);
  } else {
	$md = $len;
  }
  push(@$s, "MD:Z:$md");
}
