#!/usr/bin/perl -w

# Author: lh3

use strict;
use warnings;
use Getopt::Std;

&main;
exit;

sub main {
  my $version = '0.1.0';
  &usage if (@ARGV < 1);
  my $command = shift(@ARGV);
  my %func = (subsam=>\&subsam);
  die("Unknown command \"$command\".\n") if (!defined($func{$command}));
  &{$func{$command}};
}

sub subsam {
  die(qq/Usage: vcfutils.pl subsam <in.vcf> [samples]\n/) if (@ARGV == 0);
  my ($fh, %h);
  my $fn = shift(@ARGV);
  my @col;
  open($fh, ($fn =~ /\.gz$/)? "gzip -dc $fn |" : $fn) || die;
  $h{$_} = 1 for (@ARGV);
  while (<$fh>) {
	if (/^##/) {
	  print;
	} elsif (/^#/) {
	  my @t = split;
	  my @s = @t[0..8]; # all fixed fields + FORMAT
	  for (9 .. $#t) {
		if ($h{$t[$_]}) {
		  push(@s, $t[$_]);
		  push(@col, $_);
		}
	  }
	  pop(@s) if (@s == 9); # no sample selected; remove the FORMAT field
	  print join("\t", @s), "\n";
	} else {
	  my @t = split;
	  if (@col == 0) {
		print join("\t", @t[0..7]), "\n";
	  } else {
		print join("\t", @t[0..8], map {$t[$_]} @col), "\n";
	  }
	}
  }
  close($fh);
}

sub usage {
  die(qq/
Usage:   vcfutils.pl <command> [<arguments>]\n
Command: subsam       get a subset of samples
\n/);
}
