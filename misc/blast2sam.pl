#!/usr/bin/env perl
#
#    Copyright (C) 2009 Genome Research Ltd.
#    Portions copyright (C) 2014 Ontario Institute for Cancer Research.
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

use strict;
use warnings;
use Getopt::Std;

my $dummy_score;

&blast2sam;

sub blast2sam {
  my %opts = ();
  getopts('sd', \%opts);
  die("Usage: blast2sam.pl <in.blastn>\n") if (-t STDIN && @ARGV == 0);
  my ($qlen, $slen, $q, $s, $qbeg, $qend, @sam, @cigar, @cmaux, $show_seq);
  $show_seq = defined($opts{s});
  $dummy_score = defined($opts{d});
  @sam = (); @sam[0,4,6..8,10] = ('', 255, '*', 0, 0, '*');
  while (<>) {
    if (@cigar && (/^Query=/ || /Score =.*bits.*Expect/ || /^>\S+/)) { # print
      &blast_print_sam(\@sam, \@cigar, \@cmaux, $qlen - $qend);
      @cigar = ();
    }
    if (/^Query=\s(\S+)/) {
      $sam[2] = undef;
      $sam[0] = $1;
      my $next_line = <>;
      if ($next_line=~/^(\S+)$/) {
        $sam[0] .= $1;
      }
    } elsif (/(\S+)\s+total letters/) {
      $qlen = $1; $qlen =~ s/,//g;
    } elsif (/^>(\S+)/) {
      $sam[2] = $1;
    } elsif (/Length\s*=\s*(\d+)/) {
      $slen = $1;
    } elsif (/Score\s+=\s+(\S+) bits.+Expect(\(\d+\))?\s+=\s+(\S+)/) { # the start of an alignment block
      my ($as, $ev) = (int($1 + .499), $3);
      $ev = "1$ev" if ($ev =~ /^e/);
      @sam[1,3,9,11,12] = (0, 0, '', "AS:i:$as", "EV:Z:$ev");
      @cigar = (); $qbeg = 0;
      @cmaux = (0, 0, 0, '');
    } elsif (/Strand=(\S+)\/(\S+)/) {
      $sam[1] |= 0x10 if ($2 eq 'Minus');
    } elsif (/Query\s+(\d+)\s*(\S+)\s+(\d+)/) {
      $q = $2;
      unless ($qbeg) {
        $qbeg = $1;
        push(@cigar, ($1-1) . "H") if ($1 > 1);
      }
      $qend = $3;
      if ($show_seq) {
        my $x = $q;
        $x =~ s/-//g; $sam[9] .= $x;
      }
    } elsif (/Sbjct\:*\s+(\d+)\s*(\S+)\s+(\d+)/) {
     $s = $2;
      if ($sam[1] & 0x10) {
        $sam[3] = $3;
      } else {
        $sam[3] = $1 unless ($sam[3]);
      }
      &aln2cm(\@cigar, \$q, \$s, \@cmaux);
    }
  }
  if ($sam[2]) {
   &blast_print_sam(\@sam, \@cigar, \@cmaux, $qlen - $qend); # the last argument may be a problem
  }
}

sub blast_print_sam {
  my ($sam, $cigar, $cmaux, $qrest) = @_;
  push(@$cigar, $cmaux->[1] . substr("MDI", $cmaux->[0], 1));
  #push(@$cigar, $qrest . 'H') if ($qrest);
  if ($sam->[1] & 0x10) {
    @$cigar = reverse(@$cigar);
    $sam->[9] = reverse($sam->[9]);
    $sam->[9] =~ tr/atgcrymkswATGCRYMKSW/tacgyrkmswTACGYRKMSW/;
  }
  if ($sam->[9]) {
    if ($dummy_score) {
      $sam->[10] = "";
      map {$sam->[10].="I"} (1..length($sam->[9]));
    }
  } else {
    $sam->[9] = '*';
  }
  $sam->[5] = join('', @$cigar);
  print join("\t", @$sam), "\n";
}

sub aln2cm {
  my ($cigar, $q, $s, $cmaux) = @_;
  my $l = length($$q);
  for (my $i = 0; $i < $l; ++$i) {
    my $op;
    # set $op
    if (substr($$q, $i, 1) eq '-') { $op = 1; }
    elsif (substr($$s, $i, 1) eq '-') { $op = 2; }
    else { $op = 0; }
    # for CIGAR
    if ($cmaux->[0] == $op) {
      ++$cmaux->[1];
    } else {
      push(@$cigar, $cmaux->[1] . substr("MDI", $cmaux->[0], 1));
      $cmaux->[0] = $op; $cmaux->[1] = 1;
    }
  }
}

=head2 SYNOPSIS

blast2sam.pl is a script for parsing output of NCBI's blastn output (default format) into sam format

=over

blast2sam.pl out.blast > out.blast.sam

=back

=head2 OPTIONS

The script has some (hopefully) useful options for tweaking the output sam

B<-s> Print out sequence of the query.

Note that the current implementation prints out
the sequence of aligned query which may be trimmed or otherwise
different from the sequence of raw read in the input fastq. The CIGAR string
is also calculated for this query sequence, not the original read

B<-d> Dummy base quality score will be printed as field #11 in sam file.

Blast output does not have base quality information for a read, so this option
allows to have some fake value instead, may help when using sam file with some
programs. Hardcoded to be a string of 'I' that corresponds to Phred score 40
according to Sanger format.

Using both options:

=over

blast2sam.pl -sd out.blast > out.blast.sam

=back

Note that there is no header generated, so you will need to run

=over

samtools -hT your_ref.fasta your_file.sam > your_file_with_header.sam

=back

=cut
