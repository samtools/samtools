#!/usr/bin/perl -w

# Contact: lh3
# Version: 0.1.1

use strict;
use warnings;
use Getopt::Std;

&soap2sam;
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

sub soap2sam {
  my %opts = ();
  getopts("p", \%opts);
  die("Usage: soap2sam.pl [-p] <aln.soap>\n") if (@ARGV == 0 && -t STDIN);
  my $is_paired = defined($opts{p});
  # core loop
  my @s1 = ();
  my @s2 = ();
  my ($s_last, $s_curr) = (\@s1, \@s2);
  while (<>) {
	s/[\177-\377]|[\000-\010]|[\012-\040]//g;
	next if (&soap2sam_aux($_, $s_curr, $is_paired) < 0);
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

sub soap2sam_aux {
  my ($line, $s, $is_paired) = @_;
  chomp($line);
  my @t = split(/\s+/, $line);
  return -1 if (@t < 9 || $line =~ /^\s/ || !$t[0]);
  @$s = ();
  # fix SOAP-2.1.x bugs
  @t = @t[0..2,4..$#t] unless ($t[3] =~ /^\d+$/);
  # read name
  $s->[0] = $t[0];
  $s->[0] =~ s/\/[12]$//g;
  # initial flag (will be updated later)
  $s->[1] = 0;
  $s->[1] |= 1 | 1<<($t[4] eq 'a'? 6 : 7);
  $s->[1] |= 2 if ($is_paired);
  # read & quality
  $s->[9] = $t[1];
  $s->[10] = (length($t[2]) > length($t[1]))? substr($t[2], 0, length($t[1])) : $t[2];
  # cigar
  $s->[5] = length($s->[9]) . "M";
  # coor
  $s->[2] = $t[7]; $s->[3] = $t[8];
  $s->[1] |= 0x10 if ($t[6] eq '-');
  # mapQ
  $s->[4] = $t[3] == 1? 30 : 0;
  # mate coordinate
  $s->[6] = '*'; $s->[7] = $s->[8] = 0;
  # aux
  push(@$s, "NM:i:$t[9]");
  my $md = '';
  if ($t[9]) {
	my @x;
	for (10 .. $#t) {
	  push(@x, sprintf("%.3d,$1", $2)) if ($t[$_] =~ /^([ACGT])->(\d+)/i);
	}
	@x = sort(@x);
	my $a = 0;
	for (@x) {
	  my ($y, $z) = split(",");
	  $md .= (int($y)-$a) . $z;
	  $a += $y - $a + 1;
	}
	$md .= length($t[1]) - $a;
  } else {
	$md = length($t[1]);
  }
  push(@$s, "MD:Z:$md");
  return 0;
}
