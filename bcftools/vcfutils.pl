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
  my %func = (subsam=>\&subsam, listsam=>\&listsam, fillac=>\&fillac, qstats=>\&qstats, varFilter=>\&varFilter);
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

sub listsam {
  die(qq/Usage: vcfutils.pl listsam <in.vcf>\n/) if (@ARGV == 0 && -t STDIN);
  while (<>) {
	if (/^#/ && !/^##/) {
	  my @t = split;
	  print join("\n", @t[9..$#t]), "\n";
	  exit;
	}
  }
}

sub fillac {
  die(qq/Usage: vcfutils.pl fillac <in.vcf>\n\nNote: The GT field MUST BE present and always appear as the first field.\n/) if (@ARGV == 0 && -t STDIN);
  while (<>) {
	if (/^#/) {
	  print;
	} else {
	  my @t = split;
	  my @c;
	  my $n = 0;
	  $c[1] = 0;
	  for (9 .. $#t) {
		if ($t[$_] =~ /^(\d+).(\d+)/) {
		  ++$c[$1]; ++$c[$2];
		  $n += 2;
		}
	  }
	  my $AC = "AC=" . join("\t", @c[1..$#c]) . ";AN=$n";
	  my $info = $t[7];
	  $info =~ s/(;?)AC=(\d+)//;
	  $info =~ s/(;?)AN=(\d+)//;
	  if ($info eq '.') {
		$info = $AC;
	  } else {
		$info .= ";$AC";
	  }
	  $t[7] = $info;
	  print join("\t", @t), "\n";
	}
  }
}

sub qstats {
  my %opts = (r=>'', s=>0.01);
  getopts('r:s:', \%opts);
  die("Usage: vcfutils.pl qstats [-r ref.vcf] <in.vcf>\n
Note: This command discards indels. Output: QUAL #non-indel #SNPs #transitions #joint ts/tv #joint/#ref #joint/#non-indel \n") if (@ARGV == 0 && -t STDIN);
  my %ts = (AG=>1, GA=>1, CT=>1, TC=>1);
  my %h = ();
  if ($opts{r}) { # read the reference positions
	my $fh;
	open($fh, $opts{r}) || die;
	while (<$fh>) {
	  next if (/^#/);
	  $h{$1,$2} = 1 if (/^(\S+)\s+(\d+)/);
	}
	close($fh);
  }
  my $hsize = scalar(keys %h);
  my @a;
  while (<>) {
	next if (/^#/);
	my @t = split;
	next if (length($t[3]) != 1 || uc($t[3]) eq 'N');
	$t[3] = uc($t[3]); $t[4] = uc($t[4]);
	my @s = split(',', $t[4]);
	$t[5] = 3 if ($t[5] < 0);
	next if (length($s[0]) != 1);
	push(@a, [$t[5], ($t[4] eq '.' || $t[4] eq $t[3])? 0 : 1, $ts{$t[3].$s[0]}? 1 : 0, $h{$t[0],$t[1]}? 1 : 0]);
  }
  push(@a, [-1, 0, 0, 0]); # end marker
  die("[qstats] No SNP data!\n") if (@a == 0);
  @a = sort {$b->[0]<=>$a->[0]} @a;
  my $next = $opts{s};
  my $last = $a[0];
  my @c = (0, 0, 0, 0);
  for my $p (@a) {
	if ($p->[0] == -1 || ($p->[0] != $last && $c[0]/@a > $next)) {
	  my @x;
	  $x[0] = sprintf("%.3f", $c[1]-$c[2]? $c[2] / ($c[1] - $c[2]) : 100);
	  $x[1] = sprintf("%.3f", $hsize? $c[3] / $hsize : 0);
	  $x[2] = sprintf("%.3f", $c[3] / $c[1]);
	  print join("\t", $last, @c, @x), "\n";
	  $next = $c[0]/@a + $opts{s};
	}
	++$c[0]; $c[1] += $p->[1]; $c[2] += $p->[2]; $c[3] += $p->[3];
	$last = $p->[0];
  }
}

sub varFilter {
  my %opts = (d=>1, D=>10000, l=>30, Q=>25, q=>10, G=>25, s=>100, w=>10, W=>10, N=>2, p=>undef, F=>.001);
  getopts('pq:d:D:l:Q:w:W:N:G:F:', \%opts);
  die(qq/
Usage:   vcfutils.pl varFilter [options] <in.vcf>

Options: -Q INT    minimum RMS mapping quality for SNPs [$opts{Q}]
         -q INT    minimum RMS mapping quality for gaps [$opts{q}]
         -d INT    minimum read depth [$opts{d}]
         -D INT    maximum read depth [$opts{D}]

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
	next if (/^#/);
	next if ($t[4] eq '.'); # skip non-var sites
	my $is_snp = 1;
	if (length($t[3]) > 1) {
	  $is_snp = 0;
	} else {
	  my @s = split(',', $t[4]);
	  for (@s) {
		$is_snp = 0 if (length > 1);
	  }
	}
	# clear the out-of-range elements
	while (@staging) {
      # Still on the same chromosome and the first element's window still affects this position?
	  last if ($staging[0][3] eq $t[0] && $staging[0][4] + $staging[0][2] + $max_dist >= $t[1]);
	  varFilter_aux(shift(@staging), $opts{p}); # calling a function is a bit slower, not much
	}
	my ($flt, $score) = (0, -1);

	# collect key annotations
	my ($dp, $mq, $af) = (-1, -1, 1);
	if ($t[7] =~ /DP=(\d+)/i) {
	  $dp = $1;
	} elsif ($t[7] =~ /DP4=(\d+),(\d+),(\d+),(\d+)/i) {
	  $dp = $1 + $2 + $3 + $4;
	}
	if ($t[7] =~ /MQ=(\d+)/i) {
	  $mq = $1;
	}
	if ($t[7] =~ /AF=([^\s;=]+)/i) {
	  $af = $1;
	} elsif ($t[7] =~ /AF1=([^\s;=]+)/i) {
	  $af = $1;
	}
	# the depth filter
	if ($dp >= 0) {
	  if ($dp < $opts{d}) {
		$flt = 2;
	  } elsif ($dp > $opts{D}) {
		$flt = 3;
	  }
	}

	# site dependent filters
	my $dlen = 0;
	if ($flt == 0) {
	  if (!$is_snp) { # an indel
        # If deletion, remember the length of the deletion
		$dlen = length($t[3]) - 1;
		$flt = 1 if ($mq < $opts{q});
		# filtering SNPs
		if ($t[5] >= $opts{G}) {
		  for my $x (@staging) {
            # Is it a SNP and is it outside the SNP filter window?
			next if ($x->[0] >= 0 || $x->[4] + $x->[2] + $ow < $t[1]);
			$x->[1] = 5 if ($x->[1] == 0);
		  }
		}
		# the indel filtering score
		$score = $t[5];
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
		$flt = 1 if ($mq < $opts{Q});
		# check adjacent SNPs
		my $k = 1;
		for my $x (@staging) {
		  ++$k if ($x->[0] < 0 && -($x->[0] + 1) > $opts{F} && $x->[4] + $x->[2] + $oW >= $t[1] && ($x->[1] == 0 || $x->[1] == 4 || $x->[1] == 5));
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
	push(@staging, [$score < 0? -$af-1 : $score, $flt, $dlen, @t]);
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

sub usage {
  die(qq/
Usage:   vcfutils.pl <command> [<arguments>]\n
Command: subsam       get a subset of samples
         listsam      list the samples
         fillac       fill the allele count field
         qstats       SNP stats stratified by QUAL
         varFilter    filtering short variants
\n/);
}
