#!/usr/bin/env perl
#
#    Copyright (C) 2009 Genome Research Ltd.
#
#    Author: Heng Li <lh3@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# This script calculates a score using the BLAST scoring
# system. However, I am not sure how to count gap opens and gap
# extensions. It seems to me that column 5-8 are not what I am
# after. This script counts gaps from the last three columns. It does
# not generate reference skip (N) in the CIGAR as it is not easy to
# directly tell which gaps correspond to introns.

use strict;
use warnings;
use Getopt::Std;

my %opts = (a=>1, b=>3, q=>5, r=>2);
getopts('a:b:q:r:', \%opts);
die("Usage: psl2sam.pl [-a $opts{a}] [-b $opts{b}] [-q $opts{q}] [-r $opts{r}] <in.psl>\n") if (@ARGV == 0 && -t STDIN);

my @stack;
my $last = '';
my ($a, $b, $q, $r) = ($opts{a}, $opts{b}, $opts{q}, $opts{r});
while (<>) {
  next unless (/^\d/);
  my @t = split;
  my @s;
  my $cigar = '';
  if ($t[8] eq '-') {
    my $tmp = $t[11];
    $t[11] = $t[10] - $t[12];
    $t[12] = $t[10] - $tmp;
  }
  @s[0..4] = ($t[9], (($t[8] eq '+')? 0 : 16), $t[13], $t[15]+1, 0);
  @s[6..10] = ('*', 0, 0, '*', '*');
  $cigar .= $t[11].'H' if ($t[11]); # 5'-end clipping
  my @x = split(',', $t[18]);
  my @y = split(',', $t[19]);
  my @z = split(',', $t[20]);
  my ($y0, $z0) = ($y[0], $z[0]);
  my ($gap_open, $gap_ext) = (0, 0, 0);
  for (1 .. $t[17]-1) {
    my $ly = $y[$_] - $y[$_-1] - $x[$_-1];
    my $lz = $z[$_] - $z[$_-1] - $x[$_-1];
    if ($ly < $lz) { # del: the reference gap is longer
      ++$gap_open;
      $gap_ext += $lz - $ly;
      $cigar .= ($y[$_] - $y0) . 'M';
      $cigar .= ($lz - $ly) . 'D';
      ($y0, $z0) = ($y[$_], $z[$_]);
    } elsif ($lz < $ly) { # ins: the query gap is longer
      ++$gap_open;
      $gap_ext += $ly - $lz;
      $cigar .= ($z[$_] - $z0) . 'M';
      $cigar .= ($ly - $lz) . 'I';
      ($y0, $z0) = ($y[$_], $z[$_]);
    }
  }
  $cigar .= ($t[12] - $y0) . 'M';
  $cigar .= ($t[10] - $t[12]).'H' if ($t[10] != $t[12]); # 3'-end clipping
  $s[5] = $cigar;
  my $score = $a * $t[0] - $b * $t[1] - $q * $gap_open - $r * $gap_ext;
  $score = 0 if ($score < 0);
  $s[11] = "AS:i:$score";
  print join("\t", @s), "\n";
}
