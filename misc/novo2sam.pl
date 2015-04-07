#!/usr/bin/env perl

# Contact: lh3
# Version: 0.1.3

#Modified by Zayed Albertyn(zayed.albertyn@gmail.com) & Colin Hercus(colin@novocraft.com)

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;

&novo2sam;
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

sub novo2sam {
  my %opts = ();
  getopts("p", \%opts);
  die("Usage: novo2sam.pl [-p] <aln.novo>\n") if (@ARGV == 0);
  my $is_paired = defined($opts{p});
  # core loop
  my @s1 = ();
  my @s2 = ();
  my ($s_last, $s_curr) = (\@s1, \@s2);
  while (<>) {
    next if (/^#/);
    next if (/(QC|NM)\s*$/ || /(R\s+\d+)\s*$/);
    &novo2sam_aux($_, $s_curr, $is_paired);
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

sub novo2sam_aux {
  my ($line, $s, $is_paired) = @_;

  chomp($line);
  my @t = split(/\s+/, $line);
  my @variations =  @t[13 .. $#t];
  @$s = ();
  return if ($t[4] ne 'U');
  my $len = length($t[2]);
  # read name
  $s->[0] = substr($t[0], 1);
  $s->[0] =~ s/\/[12]$//g;
  # initial flag (will be updated later)
  $s->[1] = 0;
  $s->[1] |= 1 | 1<<($t[1] eq 'L'? 6 : 7);
  $s->[1] |= 2 if ($t[10] eq '.');
  # read & quality
  if ($t[9] eq 'R') {
    $s->[9] = reverse($t[2]);
    $s->[10] = reverse($t[3]);
    $s->[9] =~ tr/ACGTRYMKWSNacgtrymkwsn/TGCAYRKMWSNtgcayrkmwsn/;
  } else {
    $s->[9] = $t[2]; $s->[10] = $t[3];
  }
  # cigar
   my $cigarstring ="";
  if (scalar @variations ==0 ) {
     $s->[5] = $len . "M"; # IMPORTANT: this cigar is not correct for gapped alignment
  } else {
    #convert to correct CIGAR
    my $tmpstr =  join" ",@variations ;
    if ( $tmpstr=~ /\+|\-/ ) {
        $cigarstring  = cigar_method($line,\@variations,$len);
        $s->[5]=$cigarstring;
    } else {
        $s->[5]=$len. "M";
    }
}

# coor
  $s->[2] = substr($t[7], 1); $s->[3] = $t[8];
  $s->[1] |= 0x10 if ($t[9] eq 'R');
  # mapQ
  $s->[4] = $t[5] > $t[6]? $t[5] : $t[6];
  # mate coordinate
  $s->[6] = '*'; $s->[7] = $s->[8] = 0;
  # aux
  push(@$s, "NM:i:".(@t-13));
  my $md = '';
  $md = mdtag($md,$line,\@variations,$len);
  push(@$s, "MD:Z:$md");

}

sub mdtag {
    my $oldmd = shift;
    my $line = shift;
    my $ref =shift;
    my $rdlen  = shift;
    my @variations = @$ref;
    my $string="";
    my $mdtag="";
    my $t=1;
    my $q=1;
    my $deleteflag=0;
    my $len =0;
    foreach $string (@variations) {
        my ($indeltype,$insert) = indeltype($string);
        if ($indeltype eq "+") {
            $len = length ($insert);
            $q+=$len;
                  next;
        }
        my $pos = $1 if $string =~ /^(\d+)/;
        $len = $pos - $t;
        if ($len !=0 || ($deleteflag eq 1 && $indeltype eq ">")) {
            $mdtag.=$len;
        }
        $t+=$len;
        $q+=$len;
        if ($indeltype eq ">") {
            $mdtag.=$insert;
            $deleteflag=0;
              $t+=1;
              $q+=1;
        }
        if ($indeltype eq "-") {
            my $deletedbase = $2 if $string =~ /(\d+)\-([A-Za-z]+)/;
            if ($deleteflag == 0 ) {
                $mdtag.="^";
            }
            $mdtag.=$deletedbase;
            $deleteflag=1;
            $t+=1;
        }
    }
    $len = $rdlen - $q + 1;
    if ($len > 0) {
        $mdtag.="$len";
    }
#   print "In:$line\n";
#   print "MD: OLD => NEW\nMD: $oldmd => $mdtag\n\n";

    return $mdtag;
}

sub indeltype {
    my $string =  shift;
    my $insert="";
    my $indeltype;
               if ($string =~ /([A-Za-z]+)\>/) {
                        $indeltype=">";
                        $insert=$1;
                } elsif ($string =~ /\-/) {
                        $indeltype="-";
                } elsif ($string =~ /\+([A-Za-z]+)/) {
                        $indeltype="+";
                        $insert=$1;
                }
     return ($indeltype,$insert);

}


sub cigar_method {
    my $line = shift;
    my $ref =shift;
    my $rdlen  = shift;
    my @variations = @$ref;
    my $string="";
    my $type="";
    my $t =1;
    my $q=1;
    my $indeltype="";
    my $cigar=  "";
    my $insert = "";
    my $len=0;
    my @cig=();
    foreach $string (@variations) {
        next if $string =~  />/;
        my $pos = $1 if $string =~ /^(\d+)/;

        if ($string =~ /\+([A-Za-z]+)/) {
            $indeltype="+";
            $insert = $1;
        }elsif ($string =~ /\-([A-Za-z]+)/) {
            $indeltype="-";
            $insert = $1;
        }
#print "$pos $indeltype $insert $t $q\n";
        $len = $pos - $t;
        if ( $len > 0) {
            $cigar.=$len."M";
            push(@cig,$len."M");
        }
        $t+=$len;
        $q+=$len;

        if ($indeltype eq "-") {
            $cigar.="D";
            push(@cig,"D");
            $t++;
        }
        if ($indeltype eq "+") {
            $len = length ($insert);
            if ($len == 1) {
                $cigar.="I";
                push(@cig,"I");
            }
            if ($len > 1) {
                $cigar.=$len."I";
                push(@cig,$len."I")
            }
            $q+=$len;
        }
        $insert="";
    }
    $len= $rdlen - $q + 1;
    if ($len > 0) {
        $cigar.=$len."M";
        push(@cig,$len."M");
    }

      $cigar = newcigar($cigar,'D');
      $cigar = newcigar($cigar,'I');

    #print "$line\n";
    #print "c CIGAR:\t$cigar\n\n";
    return $cigar;

}



sub newcigar {
    my $cigar = shift;
    my $char = shift;
    my $new = "";
    my $copy = $cigar;
#print "$cigar\n";
    $copy =~ s/^($char+)/$1;/g;
#print "$copy\n";
    $copy =~ s/([^0-9$char])($char+)/$1;$2;/g;
#print "$copy\n";
    my @parts = split(/;/,$copy);
    my $el="";
    foreach $el (@parts) {
#print "$el\n";
        if ($el =~ /^$char+$/) {
            $new.=length($el).$char;
        }else {
            $new.=$el;
        }

    }
     return  $new;
}
