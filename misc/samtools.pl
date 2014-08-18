#!/usr/bin/env perl
#
#    Copyright (C) 2009, 2010 Genome Research Ltd.
#    Portions copyright (C) 2010 Broad Institute.
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


# Author: lh3

use strict;
use warnings;
use Getopt::Std;

my $version = '0.3.3';
&usage if (@ARGV < 1);

my $command = shift(@ARGV);
my %func = (showALEN=>\&showALEN, pileup2fq=>\&pileup2fq, varFilter=>\&varFilter, plp2vcf=>\&plp2vcf,
            unique=>\&unique, uniqcmp=>\&uniqcmp, sra2hdr=>\&sra2hdr, sam2fq=>\&sam2fq);

die("Unknown command \"$command\".\n") if (!defined($func{$command}));
&{$func{$command}};
exit(0);

#
# showALEN
#

sub showALEN {
  die(qq/Usage: samtools.pl showALEN <in.sam>\n/) if (@ARGV == 0 && -t STDIN);
  while (<>) {
    my @t = split;
    next if (/^\@/ || @t < 11);
    my $l = 0;
    $_ = $t[5];
    s/(\d+)[MI]/$l+=$1/eg;
    print join("\t", @t[0..5]), "\t$l\t", join("\t", @t[6..$#t]), "\n";
  }
}

#
# varFilter
#

#
# Filtration code:
#
# d low depth
# D high depth
# W too many SNPs in a window (SNP only)
# G close to a high-quality indel (SNP only)
# Q low RMS mapping quality (SNP only)
# g close to another indel with higher quality (indel only)
# s low SNP quality (SNP only)
# i low indel quality (indel only)

sub varFilter {
  my %opts = (d=>3, D=>100, l=>30, Q=>25, q=>10, G=>25, s=>100, w=>10, W=>10, N=>2, p=>undef, S=>'', i=>'');
  getopts('pq:d:D:l:Q:w:W:N:G:S:i:', \%opts);
  die(qq/
Usage:   samtools.pl varFilter [options] <in.cns-pileup>

Options: -Q INT    minimum RMS mapping quality for SNPs [$opts{Q}]
         -q INT    minimum RMS mapping quality for gaps [$opts{q}]
         -d INT    minimum read depth [$opts{d}]
         -D INT    maximum read depth [$opts{D}]
         -S INT    minimum SNP quality [$opts{S}]
         -i INT    minimum indel quality [$opts{i}]

         -G INT    min indel score for nearby SNP filtering [$opts{G}]
         -w INT    SNP within INT bp around a gap to be filtered [$opts{w}]

         -W INT    window size for filtering dense SNPs [$opts{W}]
         -N INT    max number of SNPs in a window [$opts{N}]

         -l INT    window size for filtering adjacent gaps [$opts{l}]

         -p        print filtered variants
\n/) if (@ARGV == 0 && -t STDIN);

  # calculate the window size
  my ($ol, $ow, $oW) = ($opts{l}, $opts{w}, $opts{W});
  my $max_dist = $ol > $ow? $ol : $ow;
  $max_dist = $oW if ($max_dist < $oW);
  # the core loop
  my @staging; # (indel_filtering_score, flt_tag)
  while (<>) {
    my @t = split;
    next if (uc($t[2]) eq uc($t[3]) || $t[3] eq '*/*'); # skip non-var sites
    # clear the out-of-range elements
    while (@staging) {
      # Still on the same chromosome and the first element's window still affects this position?
      last if ($staging[0][3] eq $t[0] && $staging[0][4] + $staging[0][2] + $max_dist >= $t[1]);
      varFilter_aux(shift(@staging), $opts{p}); # calling a function is a bit slower, not much
    }
    my ($flt, $score) = (0, -1);
    # first a simple filter
    if ($t[7] < $opts{d}) {
      $flt = 2;
    } elsif ($t[7] > $opts{D}) {
      $flt = 3;
    }
    if ($t[2] eq '*') { # an indel
        if ($opts{i} && $opts{i}>$t[5]) { $flt = 8; }
    }
    elsif ($opts{S} && $opts{S}>$t[5]) { $flt = 7; }    # SNP

    # site dependent filters
    my $len=0;
    if ($flt == 0) {
      if ($t[2] eq '*') { # an indel
        # If deletion, remember the length of the deletion
        my ($a,$b) = split(m{/},$t[3]);
        my $alen = length($a) - 1;
        my $blen = length($b) - 1;
        if ( $alen>$blen )
        {
            if ( substr($a,0,1) eq '-' ) { $len=$alen; }
        }
        elsif ( substr($b,0,1) eq '-' ) { $len=$blen; }

        $flt = 1 if ($t[6] < $opts{q});
        # filtering SNPs
        if ($t[5] >= $opts{G}) {
          for my $x (@staging) {
            # Is it a SNP and is it outside the SNP filter window?
            next if ($x->[0] >= 0 || $x->[4] + $x->[2] + $ow < $t[1]);
            $x->[1] = 5 if ($x->[1] == 0);
          }
        }
        # calculate the filtering score (different from indel quality)
        $score = $t[5];
        $score += $opts{s} * $t[10] if ($t[8] ne '*');
        $score += $opts{s} * $t[11] if ($t[9] ne '*');
        # check the staging list for indel filtering
        for my $x (@staging) {
          # Is it a SNP and is it outside the gap filter window
          next if ($x->[0] < 0 || $x->[4] + $x->[2] + $ol < $t[1]);
          if ($x->[0] < $score) {
            $x->[1] = 6;
          } else {
            $flt = 6; last;
          }
        }
      } else { # a SNP
        $flt = 1 if ($t[6] < $opts{Q});
        # check adjacent SNPs
        my $k = 1;
        for my $x (@staging) {
          ++$k if ($x->[0] < 0 && $x->[4] + $x->[2] + $oW >= $t[1] && ($x->[1] == 0 || $x->[1] == 4 || $x->[1] == 5));
        }
        # filtering is necessary
        if ($k > $opts{N}) {
          $flt = 4;
          for my $x (@staging) {
             $x->[1] = 4 if ($x->[0] < 0 && $x->[4] + $x->[2] + $oW >= $t[1] && $x->[1] == 0);
          }
        } else { # then check gap filter
          for my $x (@staging) {
            next if ($x->[0] < 0 || $x->[4] + $x->[2] + $ow < $t[1]);
            if ($x->[0] >= $opts{G}) {
              $flt = 5; last;
            }
          }
        }
      }
    }
    push(@staging, [$score, $flt, $len, @t]);
  }
  # output the last few elements in the staging list
  while (@staging) {
    varFilter_aux(shift @staging, $opts{p});
  }
}

sub varFilter_aux {
  my ($first, $is_print) = @_;
  if ($first->[1] == 0) {
    print join("\t", @$first[3 .. @$first-1]), "\n";
  } elsif ($is_print) {
    print STDERR join("\t", substr("UQdDWGgsiX", $first->[1], 1), @$first[3 .. @$first-1]), "\n";
  }
}

#
# pileup2fq
#

sub pileup2fq {
  my %opts = (d=>3, D=>255, Q=>25, G=>25, l=>10);
  getopts('d:D:Q:G:l:', \%opts);
  die(qq/
Usage:   samtools.pl pileup2fq [options] <in.cns-pileup>

Options: -d INT    minimum depth        [$opts{d}]
         -D INT    maximum depth        [$opts{D}]
         -Q INT    min RMS mapQ         [$opts{Q}]
         -G INT    minimum indel score  [$opts{G}]
         -l INT    indel filter winsize [$opts{l}]\n
/) if (@ARGV == 0 && -t STDIN);

  my ($last_chr, $seq, $qual, @gaps, $last_pos);
  my $_Q = $opts{Q};
  my $_d = $opts{d};
  my $_D = $opts{D};

  $last_chr = '';
  while (<>) {
    my @t = split;
    if ($last_chr ne $t[0]) {
      &p2q_post_process($last_chr, \$seq, \$qual, \@gaps, $opts{l}) if ($last_chr);
      $last_chr = $t[0];
      $last_pos = 0;
      $seq = ''; $qual = '';
      @gaps = ();
    }
    if ($t[1] - $last_pos != 1) {
      $seq .= 'n' x ($t[1] - $last_pos - 1);
      $qual .= '!' x ($t[1] - $last_pos - 1);
    }
    if ($t[2] eq '*') {
      push(@gaps, $t[1]) if ($t[5] >= $opts{G});
    } else {
      $seq .= ($t[6] >= $_Q && $t[7] >= $_d && $t[7] <= $_D)? uc($t[3]) : lc($t[3]);
      my $q = $t[4] + 33;
      $q = 126 if ($q > 126);
      $qual .= chr($q);
    }
    $last_pos = $t[1];
  }
  &p2q_post_process($last_chr, \$seq, \$qual, \@gaps, $opts{l});
}

sub p2q_post_process {
  my ($chr, $seq, $qual, $gaps, $l) = @_;
  &p2q_filter_gaps($seq, $gaps, $l);
  print "\@$chr\n"; &p2q_print_str($seq);
  print "+\n"; &p2q_print_str($qual);
}

sub p2q_filter_gaps {
  my ($seq, $gaps, $l) = @_;
  for my $g (@$gaps) {
    my $x = $g > $l? $g - $l : 0;
    substr($$seq, $x, $l + $l) = lc(substr($$seq, $x, $l + $l));
  }
}

sub p2q_print_str {
  my ($s) = @_;
  my $l = length($$s);
  for (my $i = 0; $i < $l; $i += 60) {
    print substr($$s, $i, 60), "\n";
  }
}

#
# sam2fq
#

sub sam2fq {
  my %opts = (n=>20, p=>'');
  getopts('n:p:', \%opts);
  die("Usage: samtools.pl sam2fq [-n 20] [-p <prefix>] <inp.sam>\n") if (@ARGV == 0 && -t STDIN);
  if ($opts{p} && $opts{n} > 1) {
    my $pre = $opts{p};
    my @fh;
    for (0 .. $opts{n}-1) {
      open($fh[$_], sprintf("| gzip > $pre.%.3d.fq.gz", $_)) || die;
    }
    my $i = 0;
    while (<>) {
      next if (/^@/);
      chomp;
      my @t = split("\t");
      next if ($t[9] eq '*');
      my ($name, $seq, $qual);
      if ($t[1] & 16) { # reverse strand
        $seq = reverse($t[9]);
        $qual = reverse($t[10]);
        $seq =~ tr/ACGTacgt/TGCAtgca/;
      } else {
        ($seq, $qual) = @t[9,10];
      }
      $name = $t[0];
      $name .= "/1" if ($t[1] & 0x40);
      $name .= "/2" if ($t[1] & 0x80);
      print {$fh[$i]} "\@$name\n$seq\n";
      if ($qual ne '*') {
        print {$fh[$i]} "+\n$qual\n";
      }
      $i = 0 if (++$i == $opts{n});
    }
    close($fh[$_]) for (0 .. $opts{n}-1);
  } else {
    die("To be implemented.\n");
  }
}

#
# sra2hdr
#

# This subroutine does not use an XML parser. It requires that the SRA
# XML files are properly formated.
sub sra2hdr {
  my %opts = ();
  getopts('', \%opts);
  die("Usage: samtools.pl sra2hdr <SRA.prefix>\n") if (@ARGV == 0);
  my $pre = $ARGV[0];
  my $fh;
  # read sample
  my $sample = 'UNKNOWN';
  open($fh, "$pre.sample.xml") || die;
  while (<$fh>) {
    $sample = $1 if (/<SAMPLE.*alias="([^"]+)"/i);
  }
  close($fh);
  # read experiment
  my (%exp2lib, $exp);
  open($fh, "$pre.experiment.xml") || die;
  while (<$fh>) {
    if (/<EXPERIMENT.*accession="([^\s"]+)"/i) {
      $exp = $1;
    } elsif (/<LIBRARY_NAME>\s*(\S+)\s*<\/LIBRARY_NAME>/i) {
      $exp2lib{$exp} = $1;
    }
  }
  close($fh);
  # read run
  my ($run, @fn);
  open($fh, "$pre.run.xml") || die;
  while (<$fh>) {
    if (/<RUN.*accession="([^\s"]+)"/i) {
      $run = $1; @fn = ();
    } elsif (/<EXPERIMENT_REF.*accession="([^\s"]+)"/i) {
      print "\@RG\tID:$run\tSM:$sample\tLB:$exp2lib{$1}\n";
    } elsif (/<FILE.*filename="([^\s"]+)"/i) {
      push(@fn, $1);
    } elsif (/<\/RUN>/i) {
      if (@fn == 1) {
        print STDERR "$fn[0]\t$run\n";
      } else {
        for (0 .. $#fn) {
          print STDERR "$fn[$_]\t$run", "_", $_+1, "\n";
        }
      }
    }
  }
  close($fh);
}

#
# unique
#

sub unique {
  my %opts = (f=>250.0, q=>5, r=>2, a=>1, b=>3);
  getopts('Qf:q:r:a:b:m', \%opts);
  die("Usage: samtools.pl unique [-f $opts{f}] <in.sam>\n") if (@ARGV == 0 && -t STDIN);
  my $last = '';
  my $recal_Q = !defined($opts{Q});
  my $multi_only = defined($opts{m});
  my @a;
  while (<>) {
    my $score = -1;
    print $_ if (/^\@/);
    $score = $1 if (/AS:i:(\d+)/);
    my @t = split("\t");
    next if (@t < 11);
    if ($score < 0) { # AS tag is unavailable
      my $cigar = $t[5];
      my ($mm, $go, $ge) = (0, 0, 0);
      $cigar =~ s/(\d+)[ID]/++$go,$ge+=$1/eg;
      $cigar = $t[5];
      $cigar =~ s/(\d+)M/$mm+=$1/eg;
      $score = $mm * $opts{a} - $go * $opts{q} - $ge * $opts{r}; # no mismatches...
    }
    $score = 1 if ($score < 1);
    if ($t[0] ne $last) {
      &unique_aux(\@a, $opts{f}, $recal_Q, $multi_only) if (@a);
      $last = $t[0];
    }
    push(@a, [$score, \@t]);
  }
  &unique_aux(\@a, $opts{f}, $recal_Q, $multi_only) if (@a);
}

sub unique_aux {
  my ($a, $fac, $is_recal, $multi_only) = @_;
  my ($max, $max2, $max_i) = (0, 0, -1);
  for (my $i = 0; $i < @$a; ++$i) {
    if ($a->[$i][0] > $max) {
      $max2 = $max; $max = $a->[$i][0]; $max_i = $i;
    } elsif ($a->[$i][0] > $max2) {
      $max2 = $a->[$i][0];
    }
  }
  if ($is_recal) {
    if (!$multi_only || @$a > 1) {
      my $q = int($fac * ($max - $max2) / $max + .499);
      $q = 250 if ($q > 250);
      $a->[$max_i][1][4] = $q < 250? $q : 250;
    }
  }
  print join("\t", @{$a->[$max_i][1]});
  @$a = ();
}

#
# uniqcmp: compare two SAM files
#

sub uniqcmp {
  my %opts = (q=>10, s=>100);
  getopts('pq:s:', \%opts);
  die("Usage: samtools.pl uniqcmp <in1.sam> <in2.sam>\n") if (@ARGV < 2);
  my ($fh, %a);
  warn("[uniqcmp] read the first file...\n");
  &uniqcmp_aux($ARGV[0], \%a, 0);
  warn("[uniqcmp] read the second file...\n");
  &uniqcmp_aux($ARGV[1], \%a, 1);
  warn("[uniqcmp] stats...\n");
  my @cnt;
  $cnt[$_] = 0 for (0..9);
  for my $x (keys %a) {
    my $p = $a{$x};
    my $z;
    if (defined($p->[0]) && defined($p->[1])) {
      $z = ($p->[0][0] == $p->[1][0] && $p->[0][1] eq $p->[1][1] && abs($p->[0][2] - $p->[1][2]) < $opts{s})? 0 : 1;
      if ($p->[0][3] >= $opts{q} && $p->[1][3] >= $opts{q}) {
        ++$cnt[$z*3+0];
      } elsif ($p->[0][3] >= $opts{q}) {
        ++$cnt[$z*3+1];
      } elsif ($p->[1][3] >= $opts{q}) {
        ++$cnt[$z*3+2];
      }
      print STDERR "$x\t$p->[0][1]:$p->[0][2]\t$p->[0][3]\t$p->[0][4]\t$p->[1][1]:$p->[1][2]\t$p->[1][3]\t$p->[1][4]\t",
        $p->[0][5]-$p->[1][5], "\n" if ($z && defined($opts{p}) && ($p->[0][3] >= $opts{q} || $p->[1][3] >= $opts{q}));
    } elsif (defined($p->[0])) {
      ++$cnt[$p->[0][3]>=$opts{q}? 6 : 7];
      print STDERR "$x\t$p->[0][1]:$p->[0][2]\t$p->[0][3]\t$p->[0][4]\t*\t0\t*\t",
        $p->[0][5], "\n" if (defined($opts{p}) && $p->[0][3] >= $opts{q});
    } else {
      print STDERR "$x\t*\t0\t*\t$p->[1][1]:$p->[1][2]\t$p->[1][3]\t$p->[1][4]\t",
        -$p->[1][5], "\n" if (defined($opts{p}) && $p->[1][3] >= $opts{q});
      ++$cnt[$p->[1][3]>=$opts{q}? 8 : 9];
    }
  }
  print "Consistent (high, high):   $cnt[0]\n";
  print "Consistent (high, low ):   $cnt[1]\n";
  print "Consistent (low , high):   $cnt[2]\n";
  print "Inconsistent (high, high): $cnt[3]\n";
  print "Inconsistent (high, low ): $cnt[4]\n";
  print "Inconsistent (low , high): $cnt[5]\n";
  print "Second missing (high):     $cnt[6]\n";
  print "Second missing (low ):     $cnt[7]\n";
  print "First  missing (high):     $cnt[8]\n";
  print "First  missing (low ):     $cnt[9]\n";
}

sub uniqcmp_aux {
  my ($fn, $a, $which) = @_;
  my $fh;
  $fn = "samtools view $fn |" if ($fn =~ /\.bam/);
  open($fh, $fn) || die;
  while (<$fh>) {
    my @t = split;
    next if (@t < 11);
#   my $l = ($t[5] =~ /^(\d+)S/)? $1 : 0;
    my $l = 0;
    my ($x, $nm) = (0, 0);
    $nm = $1 if (/NM:i:(\d+)/);
    $_ = $t[5];
    s/(\d+)[MI]/$x+=$1/eg;
    @{$a->{$t[0]}[$which]} = (($t[1]&0x10)? 1 : 0, $t[2], $t[3]-$l, $t[4], "$x:$nm", $x - 4 * $nm);
  }
  close($fh);
}

sub plp2vcf {
  while (<>) {
    my @t = split;
    next if ($t[3] eq '*/*');
    if ($t[2] eq '*') { # indel
      my @s = split("/", $t[3]);
      my (@a, @b);
      my ($ref, $alt);
      for (@s) {
        next if ($_ eq '*');
        if (/^-/) {
          push(@a, 'N'.substr($_, 1));
          push(@b, 'N');
        } elsif (/^\+/) {
          push(@a, 'N');
          push(@b, 'N'.substr($_, 1));
        }
      }
      if ($a[0] && $a[1]) {
        if (length($a[0]) < length($a[1])) {
          $ref = $a[1];
          $alt = ($b[0] . ('N' x (length($a[1]) - length($a[0])))) . ",$b[1]";
        } elsif (length($a[0]) > length($a[1])) {
          $ref = $a[0];
          $alt = ($b[1] . ('N' x (length($a[0]) - length($a[1])))) . ",$b[0]";
        } else {
          $ref = $a[0];
          $alt = ($b[0] eq $b[1])? $b[0] : "$b[0],$b[1]";
        }
      } else {
        $ref = $a[0]; $alt = $b[0];
      }
      print join("\t", @t[0,1], '.', $ref, $alt, $t[5], '.', '.'), "\n";
    } else { # SNP
    }
  }
}

#
# Usage
#

sub usage {
  die(qq/
Program: samtools.pl (helper script for SAMtools)
Version: $version
Contact: Heng Li <lh3\@sanger.ac.uk>\n
Usage:   samtools.pl <command> [<arguments>]\n
Command: varFilter     filtering SNPs and short indels
         pileup2fq     generate fastq from `pileup -c'
         showALEN      print alignment length (ALEN) following CIGAR
\n/);
}
