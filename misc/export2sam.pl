#!/usr/bin/perl -w

# Contact: lh3
# Version: 0.1.1 (03JAN2009)

use strict;
use warnings;
use Getopt::Std;

&export2sam;
exit;

sub export2sam {
  my ($fh1, $fh2, $is_paired);
  $is_paired = (@ARGV >= 2);
  die("export2sam.pl <read1.export> [<read2.export>]\n") if (@ARGV == 0);
  open($fh1, $ARGV[0]) || die;
  if ($is_paired) {
	open($fh2, $ARGV[1]) || die;
  }
  # conversion table
  my @conv_table;
  for (-64..64) {
	$conv_table[$_+64] = chr(int(33 + 10*log(1+10**($_/10.0))/log(10)+.499));
  }
  # core loop
  while (<$fh1>) {
	my (@s1, @s2);
	&export2sam_aux($_, \@s1, \@conv_table, $is_paired);
	if ($is_paired) {
	  $_ = <$fh2>;
	  &export2sam_aux($_, \@s2, \@conv_table, $is_paired);
	  if (@s1 && @s2) { # then set mate coordinate
		my $isize = 0;
		if ($s1[2] ne '*' && $s1[2] eq $s2[2]) { # then calculate $isize
		  my $x1 = ($s1[1] & 0x10)? $s1[3] + length($s1[9]) : $s1[3];
		  my $x2 = ($s2[1] & 0x10)? $s2[3] + length($s2[9]) : $s2[3];
		  $isize = $x2 - $x1;
		}
		# update mate coordinate
		if ($s2[2] ne '*') {
		  @s1[6..8] = (($s2[2] eq $s1[2])? "=" : $s2[2], $s2[3], $isize);
		  $s1[1] |= 0x20 if ($s2[1] & 0x10);
		} else {
		  $s1[1] |= 0x8;
		}
		if ($s1[2] ne '*') {
		  @s2[6..8] = (($s1[2] eq $s2[2])? "=" : $s1[2], $s1[3], -$isize);
		  $s2[1] |= 0x20 if ($s1[1] & 0x10);
		} else {
		  $s2[1] |= 0x8;
		}
	  }
	}
	print join("\t", @s1), "\n" if (@s1);
	print join("\t", @s2), "\n" if (@s2 && $is_paired);
  }
  close($fh1);
  close($fh2) if ($is_paired);
}

sub export2sam_aux {
  my ($line, $s, $ct, $is_paired) = @_;
  chomp($line);
  my @t = split("\t", $line);
  @$s = ();
  return if ($t[21] ne 'Y');
  # read name
  $s->[0] = $t[1]? "$t[0]_$t[1]:$t[2]:$t[3]:$t[4]:$t[5]" : "$t[0]:$t[2]:$t[3]:$t[4]:$t[5]";
  # initial flag (will be updated later)
  $s->[1] = 0;
  $s->[1] |= 1 | 1<<(5 + $t[7]) if ($is_paired);
  # read & quality
  $s->[9] = $t[8]; $s->[10] = $t[9];
  if ($t[13] eq 'R') { # then reverse the sequence and quality
	$s->[9] = reverse($t[8]);
	$s->[9] =~ tr/ACGTacgt/TGCAtgca/;
	$s->[10] = reverse($t[9]);
  }
  $s->[10] =~ s/(.)/$ct->[ord($1)]/eg; # change coding
  # cigar
  $s->[5] = length($s->[9]) . "M";
  # coor
  my $has_coor = 0;
  $s->[2] = "*";
  if ($t[10] eq 'NM' || $t[10] eq 'QC') {
	$s->[1] |= 0x8; # unmapped
  } elsif ($t[10] =~ /(\d+):(\d+):(\d+)/) {
	$s->[1] |= 0x8; # TODO: should I set BAM_FUNMAP in this case?
	push(@$s, "H0:i:$1", "H1:i:$2", "H2:i:$3")
  } else {
	$s->[2] = $t[10];
	$has_coor = 1;
  }
  $s->[3] = $has_coor? $t[12] : 0;
  $s->[1] |= 0x10 if ($has_coor && $t[13] eq 'R');
  # mapQ (TODO: should I choose the larger between $t[15] and $t[16]?)
  $s->[4] = 0;
  $s->[4] = $t[15] if ($t[15] ne '');
  $s->[4] = $t[16] if ($t[16] ne '' && $s->[4] < $t[16]);
  # mate coordinate
  $s->[6] = '*'; $s->[7] = $s->[8] = 0;
  # aux
  push(@$s, "BC:Z:$t[6]") if ($t[6]);
  push(@$s, "MD:Z:$t[14]") if ($has_coor);
  push(@$s, "SM:i:$t[15]") if ($is_paired && $has_coor);
}
