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
  my %func = (subsam=>\&subsam, listsam=>\&listsam, fillac=>\&fillac, qstats=>\&qstats);
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
	my @s = split(',', $t[4]);
	$t[3] = uc($t[3]); $t[4] = uc($t[4]);
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
	  $x[2] = sprintf("%.3f", $c[3] / $c[0]);
	  print join("\t", $last, @c, @x), "\n";
	  $next = $c[0]/@a + $opts{s};
	}
	++$c[0]; $c[1] += $p->[1]; $c[2] += $p->[2]; $c[3] += $p->[3];
	$last = $p->[0];
  }
}

sub usage {
  die(qq/
Usage:   vcfutils.pl <command> [<arguments>]\n
Command: subsam       get a subset of samples
         listsam      list the samples
         fillac       fill the allele count field
         qstats       SNP stats stratified by QUAL
\n/);
}
