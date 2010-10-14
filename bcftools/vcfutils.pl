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
  my %func = (subsam=>\&subsam, listsam=>\&listsam, fillac=>\&fillac, qstats=>\&qstats, varFilter=>\&varFilter,
			  hapmap2vcf=>\&hapmap2vcf, ucscsnp2vcf=>\&ucscsnp2vcf, filter4vcf=>\&filter4vcf, ldstats=>\&ldstats);
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
	  my @c = (0);
	  my $n = 0;
	  my $s = -1;
	  @_ = split(":", $t[8]);
	  for (0 .. $#_) {
		if ($_[$_] eq 'GT') { $s = $_; last; }
	  }
	  if ($s < 0) {
		print join("\t", @t), "\n";
		next;
	  }
	  for (9 .. $#t) {
		if ($t[$_] =~ /^0,0,0/) {
		} elsif ($t[$_] =~ /^([^\s:]+:){$s}(\d+).(\d+)/) {
		  ++$c[$2]; ++$c[$3];
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

sub ldstats {
  my %opts = (t=>0.9);
  getopts('t:', \%opts);
  die("Usage: vcfutils.pl ldstats [-t $opts{t}] <in.vcf>\n") if (@ARGV == 0 && -t STDIN);
  my $cutoff = $opts{t};
  my ($last, $lastchr) = (0x7fffffff, '');
  my ($x, $y, $n) = (0, 0, 0);
  while (<>) {
	if (/^([^#\s]+)\s(\d+)/) {
	  my ($chr, $pos) = ($1, $2);
	  if (/NEIR=([\d\.]+)/) {
		++$n;
		++$y, $x += $pos - $last if ($lastchr eq $chr && $pos > $last && $1 > $cutoff);
	  }
	  $last = $pos; $lastchr = $chr;
	}
  }
  print "Number of SNP intervals in strong LD (r > $opts{t}): $y\n";
  print "Fraction: ", $y/$n, "\n";
  print "Length: $x\n";
}

sub qstats {
  my %opts = (r=>'', s=>0.02, v=>undef);
  getopts('r:s:v', \%opts);
  die("Usage: vcfutils.pl qstats [-r ref.vcf] <in.vcf>\n
Note: This command discards indels. Output: QUAL #non-indel #SNPs #transitions #joint ts/tv #joint/#ref #joint/#non-indel \n") if (@ARGV == 0 && -t STDIN);
  my %ts = (AG=>1, GA=>1, CT=>1, TC=>1);
  my %h = ();
  my $is_vcf = defined($opts{v})? 1 : 0;
  if ($opts{r}) { # read the reference positions
	my $fh;
	open($fh, $opts{r}) || die;
	while (<$fh>) {
	  next if (/^#/);
	  if ($is_vcf) {
		my @t = split;
		$h{$t[0],$t[1]} = $t[4];
	  } else {
		$h{$1,$2} = 1 if (/^(\S+)\s+(\d+)/);
	  }
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
	my $hit;
	if ($is_vcf) {
	  $hit = 0;
	  my $aa = $h{$t[0],$t[1]};
	  if (defined($aa)) {
		my @aaa = split(",", $aa);
		for (@aaa) {
		  $hit = 1 if ($_ eq $s[0]);
		}
	  }
	} else {
	  $hit = defined($h{$t[0],$t[1]})? 1 : 0;
	}
	push(@a, [$t[5], ($t[4] eq '.' || $t[4] eq $t[3])? 0 : 1, $ts{$t[3].$s[0]}? 1 : 0, $hit]);
  }
  push(@a, [-1, 0, 0, 0]); # end marker
  die("[qstats] No SNP data!\n") if (@a == 0);
  @a = sort {$b->[0]<=>$a->[0]} @a;
  my $next = $opts{s};
  my $last = $a[0];
  my @c = (0, 0, 0, 0);
  my @lc;
  $lc[1] = $lc[2] = 0;
  for my $p (@a) {
	if ($p->[0] == -1 || ($p->[0] != $last && $c[0]/@a > $next)) {
	  my @x;
	  $x[0] = sprintf("%.4f", $c[1]-$c[2]? $c[2] / ($c[1] - $c[2]) : 100);
	  $x[1] = sprintf("%.4f", $hsize? $c[3] / $hsize : 0);
	  $x[2] = sprintf("%.4f", $c[3] / $c[1]);
	  my $a = $c[1] - $lc[1];
	  my $b = $c[2] - $lc[2];
	  $x[3] = sprintf("%.4f", $a-$b? $b / ($a-$b) : 100);
	  print join("\t", $last, @c, @x), "\n";
	  $next = $c[0]/@a + $opts{s};
	  $lc[1] = $c[1]; $lc[2] = $c[2];
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

sub filter4vcf {
  my %opts = (d=>3, D=>2000, 1=>1e-4, 2=>1e-100, 3=>0, 4=>1e-4, Q=>10, q=>3);
  getopts('d:D:1:2:3:4:Q:q:', \%opts);
  die(qq/
Usage:   vcfutils.pl filter4vcf [options] <in.vcf>

Options: -d INT     min total depth (given DP or DP4) [$opts{d}]
         -D INT     max total depth [$opts{D}]
         -q INT     min SNP quality [$opts{q}]
         -Q INT     min RMS mapQ (given MQ) [$opts{Q}]
         -1 FLOAT   min P-value for strand bias (given PV4) [$opts{1}]
         -2 FLOAT   min P-value for baseQ bias [$opts{2}]
         -3 FLOAT   min P-value for mapQ bias [$opts{3}]
         -4 FLOAT   min P-value for end distance bias [$opts{4}]\n
/) if (@ARGV == 0 && -t STDIN);

  my %ts = (AG=>1, GA=>1, CT=>1, TC=>1);

  my @n = (0, 0);
  while (<>) {
	if (/^#/) {
	  print;
	  next;
	}
	next if (/PV4=([^,]+),([^,]+),([^,]+),([^,;\t]+)/ && ($1<$opts{1} || $2<$opts{2} || $3<$opts{3} || $4<$opts{4}));
	my $depth = -1;
	$depth = $1 if (/DP=(\d+)/);
	$depth = $1+$2+$3+$4 if (/DP4=(\d+),(\d+),(\d+),(\d+)/);
	next if ($depth > 0 && ($depth < $opts{d} || $depth > $opts{D}));
	next if (/MQ=(\d+)/ && $1 < $opts{Q});
	my @t = split;
	next if ($t[5] >= 0 && $t[5] < $opts{q});
	++$n[0];
	my @s = split(',', $t[4]);
	++$n[1] if ($ts{$t[3].$s[0]});
	print;
  }
}

sub ucscsnp2vcf {
  die("Usage: vcfutils.pl <in.ucsc.snp>\n") if (@ARGV == 0 && -t STDIN);
  print "##fileformat=VCFv4.0\n";
  print join("\t", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), "\n";
  while (<>) {
	my @t = split("\t");
	my $indel = ($t[9] =~ /^[ACGT](\/[ACGT])+$/)? 0 : 1;
	my $pos = $t[2] + 1;
	my @alt;
	push(@alt, $t[7]);
	if ($t[6] eq '-') {
	  $t[9] = reverse($t[9]);
	  $t[9] =~ tr/ACGTRYMKWSNacgtrymkwsn/TGCAYRKMWSNtgcayrkmwsn/;
	}
	my @a = split("/", $t[9]);
	for (@a) {
	  push(@alt, $_) if ($_ ne $alt[0]);
	}
	if ($indel) {
	  --$pos;
	  for (0 .. $#alt) {
		$alt[$_] =~ tr/-//d;
		$alt[$_] = "N$alt[$_]";
	  }
	}
	my $ref = shift(@alt);
	my $af = $t[13] > 0? ";AF=$t[13]" : '';
	my $valid = ($t[12] eq 'unknown')? '' : ";valid=$t[12]";
	my $info = "molType=$t[10];class=$t[11]$valid$af";
	print join("\t", $t[1], $pos, $t[4], $ref, join(",", @alt), 0, '.', $info), "\n";
  }
}

sub hapmap2vcf {
  die("Usage: vcfutils.pl <in.ucsc.snp> <in.hapmap>\n") if (@ARGV == 0);
  my $fn = shift(@ARGV);
  # parse UCSC SNP
  warn("Parsing UCSC SNPs...\n");
  my ($fh, %map);
  open($fh, ($fn =~ /\.gz$/)? "gzip -dc $fn |" : $fn) || die;
  while (<$fh>) {
	my @t = split;
	next if ($t[3] - $t[2] != 1); # not SNP
	@{$map{$t[4]}} = @t[1,3,7];
  }
  close($fh);
  # write VCF
  warn("Writing VCF...\n");
  print "##fileformat=VCFv4.0\n";
  while (<>) {
	my @t = split;
	if ($t[0] eq 'rs#') { # the first line
	  print join("\t", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", @t[11..$#t]), "\n";
	} else {
	  next unless ($map{$t[0]});
	  next if (length($t[1]) != 3); # skip non-SNPs
	  my $a = \@{$map{$t[0]}};
	  my $ref = $a->[2];
	  my @u = split('/', $t[1]);
	  if ($u[1] eq $ref) {
		$u[1] = $u[0]; $u[0] = $ref;
	  } elsif ($u[0] ne $ref) { next; }
	  my $alt = $u[1];
	  my %w;
	  $w{$u[0]} = 0; $w{$u[1]} = 1;
	  my @s = (@$a[0,1], $t[0], $ref, $alt, 0, '.', '.', 'GT');
	  my $is_tri = 0;
	  for (@t[11..$#t]) {
		if ($_ eq 'NN') {
		  push(@s, './.');
		} else {
		  my @a = ($w{substr($_,0,1)}, $w{substr($_,1,1)});
		  if (!defined($a[0]) || !defined($a[1])) {
			$is_tri = 1;
			last;
		  }
		  push(@s, "$a[0]/$a[1]");
		}
	  }
	  next if ($is_tri);
	  print join("\t", @s), "\n";
	}
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
         filter4vcf   filtering VCFs produced by samtools+bcftools
         hapmap2vcf   convert the hapmap format to VCF
         ucscsnp2vcf  convert UCSC SNP SQL dump to VCF
\n/);
}
