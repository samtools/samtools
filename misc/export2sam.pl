#!/usr/bin/env perl
#
#
# export2sam.pl converts GERALD export files to SAM format.
#
#
#
########## License:
#
# The MIT License
#
# Copyright (c) 2008-2009 Genome Research Ltd.
# Modifications Copyright (c) 2010 Illumina, Inc.
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
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
#
#
#
########## ChangeLog:
#
# Version: 2.3.1 (18MAR2011)
#
#   - Restore file '-' as stdin input.
#
# Version: 2.3.0 (24JAN2011)
#
#   - Add support for export reserved chromosome name "CONTROL",
#       which is translated to optional field "XC:Z:CONTROL".
#   - Check for ".gz" file extension on export files and open
#       these as gzip pipes when the extension is found.
#
# Version: 2.2.0 (16NOV2010)
#
#   - Remove any leading zeros in export fields: RUNNO,LANE,TILE,X,Y
#   - For export records with reserved chromosome name identifiers
#       "QC" and "RM", add the optional field "XC:Z:QC" or "XC:Z:RM"
#       to the SAM record, so that these cases can be distinguished
#       from other unmatched reads.
#
# Version: 2.1.0 (21SEP2010)
#
#   - Additional export record error checking.
#   - Convert export records with chromomsome value of "RM" to unmapped
#       SAM records.
#
# Version: 2.0.0 (15FEB2010)
#
#   Script updated by Illumina in conjunction with CASAVA 1.7.0
#   release.
#
#   Major changes are as follows:
#   - The CIGAR string has been updated to include all gaps from
#       ELANDv2 alignments.
#   - The ELAND single read alignment score is always stored in the
#       optional "SM" field and the ELAND paired read alignment score
#       is stored in the optional "AS" field when it exists.
#   - The MAPQ value is set to the higher of the two alignment scores,
#       but no greater than 254, i.e. min(254,max(SM,AS))
#   - The SAM "proper pair" bit (0x0002) is now set for read pairs
#       meeting ELAND's expected orientation and insert size criteria.
#   - The default quality score translation is set for export files
#       which contain Phread+64 quality values. An option,
#       "--qlogodds", has been added to translate quality values from
#       the Solexa+64 format used in export files prior to Pipeline
#       1.3
#   - The export match descriptor is now reverse-complemented when
#       necessary such that it always corresponds to the forward
#       strand of the reference, to be consistent with other
#       information in the SAM record. It is now written to the
#       optional 'XD' field (rather than 'MD') to acknowledge its
#       minor differences from the samtools match descriptor (see
#       additional detail below).
#   - An option, "--nofilter", has been added to include reads which
#       have failed primary analysis quality filtration. Such reads
#       will have the corresponding SAM flag bit (0x0200) set.
#   - Labels in the export 'contig' field are preserved by setting
#       RNAME to "$export_chromosome/$export_contig" when the contig
#       label exists.
#
#
# Contact: lh3
# Version: 0.1.2 (03JAN2009)
#
#
#
########## Known Conversion Limitations:
#
# - Export records for reads that map to a position < 1 (allowed
#     in export format), are converted to unmapped reads in the SAM
#     record.
# - Export records contain the reserved chromosome names: "NM",
#     "QC","RM" and "CONTROL". "NM" indicates that the aligner could
#     not map the read to the reference sequence set. "QC" means that
#     the aligner did not attempt to map the read due to some
#     technical limitation. "RM" means that the read mapped to a set
#     of 'contaminant' sequences specified in GERALD's RNA-seq
#     workflow. "CONTROL" means that the read is a control. All of
#     these alignment types are collapsed to the single unmapped
#     alignment state in the SAM record, but the optional SAM "XC"
#     field is used to record the original reserved chromosome name of
#     the read for all but the "NM" case.
# - The export match descriptor is slightly different than the
#     samtools match descriptor. For this reason it is stored in the
#     optional SAM field 'XD' (and not 'MD'). Note that the export
#     match descriptor differs from the samtools version in two
#     respects: (1) indels are explicitly closed with the '$'
#     character and (2) insertions must be enumerated in the match
#     descriptor. For example a 35-base read with a two-base insertion
#     is described as: 20^2$14
#
#
#

my $version = "2.3.1";

use strict;
use warnings;

use Getopt::Long;
use File::Spec;
use List::Util qw(min max);


use constant {
  EXPORT_MACHINE => 0,
  EXPORT_RUNNO => 1,
  EXPORT_LANE => 2,
  EXPORT_TILE => 3,
  EXPORT_X => 4,
  EXPORT_Y => 5,
  EXPORT_INDEX => 6,
  EXPORT_READNO => 7,
  EXPORT_READ => 8,
  EXPORT_QUAL => 9,
  EXPORT_CHROM => 10,
  EXPORT_CONTIG => 11,
  EXPORT_POS => 12,
  EXPORT_STRAND => 13,
  EXPORT_MD => 14,
  EXPORT_SEMAP => 15,
  EXPORT_PEMAP => 16,
  EXPORT_PASSFILT => 21,
  EXPORT_SIZE => 22,
};


use constant {
  SAM_QNAME => 0,
  SAM_FLAG => 1,
  SAM_RNAME => 2,
  SAM_POS => 3,
  SAM_MAPQ => 4,
  SAM_CIGAR => 5,
  SAM_MRNM => 6,
  SAM_MPOS => 7,
  SAM_ISIZE => 8,
  SAM_SEQ => 9,
  SAM_QUAL => 10,
};


# function prototypes for Richard's code
sub match_desc_to_cigar($);
sub match_desc_frag_length($);
sub reverse_compl_match_descriptor($);
sub write_header($;$;$);


&export2sam;
exit;




sub export2sam {

  my $cmdline = $0 . " " . join(" ",@ARGV);
  my $arg_count = scalar @ARGV;
  my $progname = (File::Spec->splitpath($0))[2];

  my $is_logodds_qvals = 0; # if true, assume files contain logodds (i.e. "solexa") quality values
  my $is_nofilter = 0;
  my $read1file;
  my $read2file;
  my $print_version = 0;
  my $help = 0;

  my $result = GetOptions( "qlogodds" => \$is_logodds_qvals,
                           "nofilter" => \$is_nofilter,
                           "read1=s"  => \$read1file,
                           "read2=s"  => \$read2file,
                           "version"  => \$print_version,
                           "help"     => \$help );

  my $usage = <<END;

$progname converts GERALD export files to SAM format.

Usage: $progname --read1=FILENAME [ options ] | --version | --help

  --read1=FILENAME  read1 export file or '-' for stdin (mandatory)
                      (file may be gzipped with ".gz" extension)
  --read2=FILENAME  read2 export file or '-' for stdin
                      (file may be gzipped with ".gz" extension)
  --nofilter        include reads that failed the basecaller
                      purity filter
  --qlogodds        assume export file(s) use logodds quality values
                      as reported by OLB (Pipeline) prior to v1.3
                      (default: phred quality values)

END

  my $version_msg = <<END;

$progname version: $version

END

  if((not $result) or $help or ($arg_count==0)) {
    die($usage);
  }

  if(@ARGV) {
    print STDERR "\nERROR: Unrecognized arguments: " . join(" ",@ARGV) . "\n\n";
    die($usage);
  }

  if($print_version) {
    die($version_msg);
  }

  if(not defined($read1file)) {
    print STDERR "\nERROR: read1 export file must be specified\n\n";
    die($usage);
  }

  unless((-f $read1file) or ($read1file eq '-')) {
    die("\nERROR: Can't find read1 export file: '$read1file'\n\n");
  }

  if (defined $read2file) {
    unless((-f $read2file) or ($read2file eq '-')) {
      die("\nERROR: Can't find read2 export file: '$read2file'\n\n");
    }
    if($read1file eq $read2file) {
      die("\nERROR: read1 and read2 export filenames are the same: '$read1file'\n\n");
    }
  }

  my ($fh1, $fh2, $is_paired);

  my $read1cmd="$read1file";
  $read1cmd = "gzip -dc $read1file |" if($read1file =~ /\.gz$/);
  open($fh1, $read1cmd)
      or die("\nERROR: Can't open read1 process: '$read1cmd'\n\n");
  $is_paired = defined $read2file;
  if ($is_paired) {
    my $read2cmd="$read2file";
    $read2cmd = "gzip -dc $read2file |" if($read2file =~ /\.gz$/);
    open($fh2, $read2cmd)
        or die("\nERROR: Can't open read2 process: '$read2cmd'\n\n");
  }
  # quality value conversion table
  my @conv_table;
  if($is_logodds_qvals){ # convert from solexa+64 quality values (pipeline pre-v1.3):
    for (-64..64) {
      $conv_table[$_+64] = int(33 + 10*log(1+10**($_/10.0))/log(10)+.499);
    }
  } else {               # convert from phred+64 quality values (pipeline v1.3+):
    for (-64..-1) {
      $conv_table[$_+64] = undef;
    }
    for (0..64) {
      $conv_table[$_+64] = int(33 + $_);
    }
  }
  # write the header
  print write_header( $progname, $version, $cmdline );
  # core loop
  my $export_line_count = 0;
  while (<$fh1>) {
    $export_line_count++;
    my (@s1, @s2);
    &export2sam_aux($_, $export_line_count, \@s1, \@conv_table, $is_paired, 1, $is_nofilter);
    if ($is_paired) {
      my $read2line = <$fh2>;
      if(not $read2line){
        die("\nERROR: read1 and read2 export files do not contain the same number of reads.\n  Extra reads observed in read1 file at line no: $export_line_count.\n\n");
      }
      &export2sam_aux($read2line, $export_line_count, \@s2, \@conv_table, $is_paired, 2, $is_nofilter);

      if (@s1 && @s2) { # then set mate coordinate
        if($s1[SAM_QNAME] ne $s2[SAM_QNAME]){
          die("\nERROR: Non-paired reads in export files on line: $export_line_count.\n  Read1: $_  Read2: $read2line\n");
        }

        my $isize = 0;
        if ($s1[SAM_RNAME] ne '*' && $s1[SAM_RNAME] eq $s2[SAM_RNAME]) { # then calculate $isize
          my $x1 = ($s1[SAM_FLAG] & 0x10)? $s1[SAM_POS] + length($s1[SAM_SEQ]) : $s1[SAM_POS];
          my $x2 = ($s2[SAM_FLAG] & 0x10)? $s2[SAM_POS] + length($s2[SAM_SEQ]) : $s2[SAM_POS];
          $isize = $x2 - $x1;
        }

        foreach ([\@s1,\@s2,$isize],[\@s2,\@s1,-$isize]){
          my ($sa,$sb,$is) = @{$_};
          if ($sb->[SAM_RNAME] ne '*') {
            $sa->[SAM_MRNM] = ($sb->[SAM_RNAME] eq $sa->[SAM_RNAME]) ? "=" : $sb->[SAM_RNAME];
            $sa->[SAM_MPOS] = $sb->[SAM_POS];
            $sa->[SAM_ISIZE] = $is;
            $sa->[SAM_FLAG] |= 0x20 if ($sb->[SAM_FLAG] & 0x10);
          } else {
            $sa->[SAM_FLAG] |= 0x8;
          }
        }
      }
    }
    print join("\t", @s1), "\n" if (@s1);
    print join("\t", @s2), "\n" if (@s2 && $is_paired);
  }
  close($fh1);
  if($is_paired) {
    while(my $read2line = <$fh2>){
      $export_line_count++;
      die("\nERROR: read1 and read2 export files do not contain the same number of reads.\n  Extra reads observed in read2 file at line no: $export_line_count.\n\n");
    }
    close($fh2);
  }
}

sub export2sam_aux {
  my ($line, $line_no, $s, $ct, $is_paired, $read_no, $is_nofilter) = @_;
  chomp($line);
  my @t = split("\t", $line);
  if(scalar(@t) < EXPORT_SIZE) {
    my $msg="\nERROR: Unexpected number of fields in export record on line $line_no of read$read_no export file. Found " . scalar(@t) . " fields but expected " . EXPORT_SIZE . ".\n";
    $msg.="\t...erroneous export record:\n" . $line . "\n\n";
    die($msg);
  }
  @$s = ();
  my $isPassFilt = ($t[EXPORT_PASSFILT] eq 'Y');
  return if(not ($isPassFilt or $is_nofilter));
  # read name
  my $samQnamePrefix = $t[EXPORT_MACHINE] . (($t[EXPORT_RUNNO] ne "") ? "_" .  int($t[EXPORT_RUNNO]) : "");
  $s->[SAM_QNAME] = join(':', $samQnamePrefix, int($t[EXPORT_LANE]), int($t[EXPORT_TILE]),
                         int($t[EXPORT_X]), int($t[EXPORT_Y]));
  # initial flag (will be updated later)
  $s->[SAM_FLAG] = 0;
  if($is_paired) {
    if($t[EXPORT_READNO] != $read_no){
      die("\nERROR: read$read_no export file contains record with read number: " .$t[EXPORT_READNO] . " on line: $line_no\n\n");
    }
    $s->[SAM_FLAG] |= 1 | 1<<(5 + $read_no);
  }
  $s->[SAM_FLAG] |= 0x200 if (not $isPassFilt);

  # read & quality
  my $is_export_rev = ($t[EXPORT_STRAND] eq 'R');
  if ($is_export_rev) { # then reverse the sequence and quality
    $s->[SAM_SEQ] = reverse($t[EXPORT_READ]);
    $s->[SAM_SEQ] =~ tr/ACGTacgt/TGCAtgca/;
    $s->[SAM_QUAL] = reverse($t[EXPORT_QUAL]);
  } else {
    $s->[SAM_SEQ] = $t[EXPORT_READ];
    $s->[SAM_QUAL] = $t[EXPORT_QUAL];
  }
  my @convqual = ();
  foreach (unpack('C*', $s->[SAM_QUAL])){
    my $val=$ct->[$_];
    if(not defined $val){
      my $msg="\nERROR: can't interpret export quality value: " . $_ . " in read$read_no export file, line: $line_no\n";
      if( $_ < 64 ) { $msg .= "  Use --qlogodds flag to translate logodds (solexa) quality values.\n"; }
      die($msg . "\n");
    }
    push @convqual,$val;
  }

  $s->[SAM_QUAL] = pack('C*',@convqual); # change coding


  # coor
  my $has_coor = 0;
  $s->[SAM_RNAME] = "*";
  if (($t[EXPORT_CHROM] eq 'NM') or
      ($t[EXPORT_CHROM] eq 'QC') or
      ($t[EXPORT_CHROM] eq 'RM') or
      ($t[EXPORT_CHROM] eq 'CONTROL')) {
    $s->[SAM_FLAG] |= 0x4; # unmapped
    push(@$s,"XC:Z:".$t[EXPORT_CHROM]) if($t[EXPORT_CHROM] ne 'NM');
  } elsif ($t[EXPORT_CHROM] =~ /(\d+):(\d+):(\d+)/) {
    $s->[SAM_FLAG] |= 0x4; # TODO: should I set BAM_FUNMAP in this case?
    push(@$s, "H0:i:$1", "H1:i:$2", "H2:i:$3")
  } elsif ($t[EXPORT_POS] < 1) {
    $s->[SAM_FLAG] |= 0x4; # unmapped
  } else {
    $s->[SAM_RNAME] = $t[EXPORT_CHROM];
    $s->[SAM_RNAME] .= "/" . $t[EXPORT_CONTIG] if($t[EXPORT_CONTIG] ne '');
    $has_coor = 1;
  }
  $s->[SAM_POS] = $has_coor? $t[EXPORT_POS] : 0;

#  print STDERR "t[14] = " . $t[14] . "\n";
  my $matchDesc = '';
  $s->[SAM_CIGAR] = "*";
  if($has_coor){
    $matchDesc = ($is_export_rev) ? reverse_compl_match_descriptor($t[EXPORT_MD]) : $t[EXPORT_MD];

    if($matchDesc =~ /\^/){
      # construct CIGAR string using Richard's function
      $s->[SAM_CIGAR] = match_desc_to_cigar($matchDesc); # indel processing
    } else {
      $s->[SAM_CIGAR] = length($s->[SAM_SEQ]) . "M";
    }
  }

#  print STDERR "cigar_string = $cigar_string\n";

  $s->[SAM_FLAG] |= 0x10 if ($has_coor && $is_export_rev);
  if($has_coor){
    my $semap = ($t[EXPORT_SEMAP] ne '') ? $t[EXPORT_SEMAP] : 0;
    my $pemap = 0;
    if($is_paired) {
      $pemap = ($t[EXPORT_PEMAP] ne '') ? $t[EXPORT_PEMAP] : 0;

      # set `proper pair' bit if non-blank, non-zero PE alignment score:
      $s->[SAM_FLAG] |= 0x02 if ($pemap > 0);
    }
    $s->[SAM_MAPQ] = min(254,max($semap,$pemap));
  } else {
    $s->[SAM_MAPQ] = 0;
  }
  # mate coordinate
  $s->[SAM_MRNM] = '*';
  $s->[SAM_MPOS] = 0;
  $s->[SAM_ISIZE] = 0;
  # aux
  push(@$s, "BC:Z:$t[EXPORT_INDEX]") if ($t[EXPORT_INDEX]);
  if($has_coor){
    # The export match descriptor differs slightly from the samtools match descriptor.
    # In order for the converted SAM files to be as compliant as possible,
    # we put the export match descriptor in optional field 'XD' rather than 'MD':
    push(@$s, "XD:Z:$matchDesc");
    push(@$s, "SM:i:$t[EXPORT_SEMAP]") if ($t[EXPORT_SEMAP] ne '');
    push(@$s, "AS:i:$t[EXPORT_PEMAP]") if ($is_paired and ($t[EXPORT_PEMAP] ne ''));
  }
}



#
# the following code is taken from Richard Shaw's sorted2sam.pl file
#
sub reverse_compl_match_descriptor($)
{
#    print "\nREVERSING THE MATCH DESCRIPTOR!\n";
    my ($match_desc) = @_;
    my $rev_compl_match_desc = reverse($match_desc);
    $rev_compl_match_desc =~ tr/ACGT\^\$/TGCA\$\^/;

    # Unreverse the digits of numbers.
    $rev_compl_match_desc = join('',
                                 map {($_ =~ /\d+/)
                                      ? join('', reverse(split('', $_)))
                                      : $_} split(/(\d+)/,
                                                  $rev_compl_match_desc));

    return $rev_compl_match_desc;
}



sub match_desc_to_cigar($)
{
    my ($match_desc) = @_;

    my @match_desc_parts = split(/(\^.*?\$)/, $match_desc);
    my $cigar_str = '';
    my $cigar_del_ch = 'D';
    my $cigar_ins_ch = 'I';
    my $cigar_match_ch = 'M';

    foreach my $match_desc_part (@match_desc_parts) {
        next if (!$match_desc_part);

        if ($match_desc_part =~ /^\^([ACGTN]+)\$$/) {
            # Deletion
            $cigar_str .= (length($1) . $cigar_del_ch);
        } elsif ($match_desc_part =~ /^\^(\d+)\$$/) {
            # Insertion
            $cigar_str .= ($1 . $cigar_ins_ch);
        } else {
            $cigar_str .= (match_desc_frag_length($match_desc_part)
                           . $cigar_match_ch);
        }
    }

    return $cigar_str;
}


#------------------------------------------------------------------------------

sub match_desc_frag_length($)
                           {
    my ($match_desc_str) = @_;
    my $len = 0;

    my @match_desc_fields = split(/([ACGTN]+)/, $match_desc_str);

    foreach my $match_desc_field (@match_desc_fields) {
        next if ($match_desc_field eq '');

        $len += (($match_desc_field =~ /(\d+)/)
                 ? $1 : length($match_desc_field));
    }

    return $len;
}


# argument holds the command line
sub write_header($;$;$)
{
    my ($progname,$version,$cl) = @_;
    my $complete_header = "";
    $complete_header .= "\@PG\tID:$progname\tVN:$version\tCL:$cl\n";

    return $complete_header;
}
