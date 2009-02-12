#!/usr/bin/perl -w

# Contact: lh3

use strict;
use warnings;
use Getopt::Std;

my %opts = (D=>100, m=>10, r=>undef);
getopts('D:m:r', \%opts);

die(qq/
Usage:   indel_filter.pl [options] <in.indel>\n
Options: -D INT    maximum read depth [$opts{D}]
         -m INT    minimum distance between two adjacent indels [$opts{m}]
\n/) if (@ARGV == 0 && -t STDIN);

my (@arr1, @arr2);
my ($curr, $last) = (\@arr1, \@arr2);
my $is_ref = defined($opts{r})? 1 : 0;
while (<>) {
  my @t = split;
  next if ($t[2] ne '*');
  if (!$is_ref) {
	next if ($t[3] eq '*/*');
	next if ($t[5] == 0);
  }
  next if ($t[7] > $opts{D});
  @$curr = ($t[0], $t[1], $t[5], $_);
  my $do_swap = 1;
  if (defined $last->[0]) {
	if ($curr->[0] eq $last->[0] && $last->[1] + $opts{m} > $curr->[1]) {
	  $do_swap = 0 if ($last->[2] > $curr->[2]);
	} else { # then print
	  print $last->[3];
	}
  }
  if ($do_swap) {
	my $tmp = $curr; $curr = $last; $last = $tmp;
  }
}
print $last->[3] if (defined $last->[0]);
