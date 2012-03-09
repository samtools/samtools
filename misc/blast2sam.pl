#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

&blast2sam;

sub blast2sam {
  my %opts = ();
  getopts('s', \%opts);
  die("Usage: blast2sam.pl <in.blastn>\n") if (-t STDIN && @ARGV == 0);
  my ($qlen, $slen, $q, $s, $qbeg, $qend, @sam, @cigar, @cmaux, $show_seq);
  $show_seq = defined($opts{s});
  @sam = (); @sam[0,4,6..8,10] = ('', 255, '*', 0, 0, '*');
  while (<>) {
	if (@cigar && (/^Query=/ || /Score =.*bits.*Expect/)) { # print
	  &blast_print_sam(\@sam, \@cigar, \@cmaux, $qlen - $qend);
	  @cigar = ();
	}
	if (/^Query= (\S+)/) {
	  $sam[0] = $1;
	} elsif (/\((\S+)\s+letters\)/) {
	  $qlen = $1; $qlen =~ s/,//g;
	} elsif (/^>(\S+)/) {
	  $sam[2] = $1;
	} elsif (/Length = (\d+)/) {
	  $slen = $1;
	} elsif (/Score =\s+(\S+) bits.+Expect(\(\d+\))? = (\S+)/) { # the start of an alignment block
	  my ($as, $ev) = (int($1 + .499), $3);
	  $ev = "1$ev" if ($ev =~ /^e/);
	  @sam[1,3,9,11,12] = (0, 0, '', "AS:i:$as", "EV:Z:$ev");
	  @cigar = (); $qbeg = 0;
	  @cmaux = (0, 0, 0, '');
	} elsif (/Strand = (\S+) \/ (\S+)/) {
	  $sam[1] |= 0x10 if ($2 eq 'Minus');
	} elsif (/Query\:\s(\d+)\s*(\S+)\s(\d+)/) {
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
	} elsif (/Sbjct\:\s(\d+)\s*(\S+)\s(\d+)/) {
	  $s = $2;
	  if ($sam[1] & 0x10) {
		$sam[3] = $3;
	  } else {
		$sam[3] = $1 unless ($sam[3]);
	  }
	  &aln2cm(\@cigar, \$q, \$s, \@cmaux);
	}
  }
  &blast_print_sam(\@sam, \@cigar, \@cmaux, $qlen - $qend);
}

sub blast_print_sam {
  my ($sam, $cigar, $cmaux, $qrest) = @_;
  push(@$cigar, $cmaux->[1] . substr("MDI", $cmaux->[0], 1));
  push(@$cigar, $qrest . 'H') if ($qrest);
  if ($sam->[1] & 0x10) {
	@$cigar = reverse(@$cigar);
	$sam->[9] = reverse($sam->[9]);
	$sam->[9] =~ tr/atgcrymkswATGCRYMKSW/tacgyrkmswTACGYRKMSW/;
  }
  $sam->[9] = '*' if (!$sam->[9]);
  $sam->[5] = join('', @$cigar);
  print join("\t", @$sam), "\n";
}

sub aln2cm {
  my ($cigar, $q, $s, $cmaux) = @_;
  my $l = length($$q);
  for (my $i = 0; $i < $l; ++$i) {
	my $op;
	# set $op
	if (substr($$q, $i, 1) eq '-') { $op = 2; }
	elsif (substr($$s, $i, 1) eq '-') { $op = 1; }
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
