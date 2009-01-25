#!/usr/bin/perl -w

# Contact: lh3
# Version: 0.1.0

use strict;
use warnings;
use Getopt::Std;

&zoom2sam;
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

sub zoom2sam {
  my %opts = ();
  getopts("p", \%opts);
  die("Usage: zoom2sam.pl [-p] <readLen> <aln.zoom>
Warnings: This script only supports the default Illumina outputs.\n") if (@ARGV < 2);
  my $is_paired = defined($opts{p});
  my $len = shift(@ARGV);
  # core loop
  my @s1 = ();
  my @s2 = ();
  my ($s_last, $s_curr) = (\@s1, \@s2);
  while (<>) {
	&zoom2sam_aux($_, $s_curr, $is_paired, $len);
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

sub zoom2sam_aux {
  my ($line, $s, $is_paired, $len) = @_;
  chomp($line);
  my @t = split("\t", $line);
  @$s = ();
  # read name
  $s->[0] = $t[0];
  # initial flag (will be updated later)
  $s->[1] = 0;
  $s->[1] |= 1 | 1<<6 if ($s->[0] =~ /_F$/);
  $s->[1] |= 1 | 1<<7 if ($s->[0] =~ /_R$/);
  $s->[1] |= 2 if ($is_paired);
  # read & quality
  $s->[9] = "*"; $s->[10] = "*";
  # cigar
  $s->[5] = $len . "M";
  # coor
  my @s = split(/\s+/, $t[1]);
  $s->[2] = $s[0];
  $t[1] =~ /:(\d+)$/;
  $s->[3] = $1 + 1;
  if ($s->[0] =~ /_[FR]$/) {
	my $u = ($s->[0] =~ /_F$/)? 1 : 0;
	my $w = ($t[2] eq '+')? 1 : 0;
	$s->[1] |= 0x10 if ($u ^ $w);
	$s->[0] =~ s/_[FR]$//;
  } else {
	$s->[1] |= 0x10 if ($t[2] eq '-');
  }
  # mapQ
  $s->[4] = 30;
  # mate coordinate
  $s->[6] = '*'; $s->[7] = $s->[8] = 0;
  # aux
  push(@$s, "NM:i:$t[3]");
}
