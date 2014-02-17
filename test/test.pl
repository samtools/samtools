#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;

my $opts = parse_params();

test_bgzip($opts);
test_faidx($opts);
test_mpileup($opts);
test_usage($opts, cmd=>'samtools');
test_view($opts);

print "\nNumber of tests:\n";
printf "    total   .. %d\n", $$opts{nok}+$$opts{nfailed};
printf "    passed  .. %d\n", $$opts{nok};
printf "    failed  .. %d\n", $$opts{nfailed};
print "\n";

exit ($$opts{nfailed} > 0);

#--------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    print 
        "About: samtools/htslib consistency test script\n",
        "Usage: test.pl [OPTIONS]\n",
        "Options:\n",
        "   -r, --redo-outputs              Recreate expected output files.\n",
        "   -t, --temp-dir <path>           When given, temporary files will not be removed.\n",
        "   -h, -?, --help                  This help message.\n",
        "\n";
    exit 1;
}
sub parse_params
{
    my $opts = { keep_files=>0, nok=>0, nfailed=>0 };
    my $help;
    Getopt::Long::Configure('bundling');
    my $ret = GetOptions (
            't|temp-dir:s' => \$$opts{keep_files}, 
            'r|redo-outputs' => \$$opts{redo_outputs}, 
            'h|?|help' => \$help
            );
    if ( !$ret or $help ) { error(); }
    $$opts{tmp}  = $$opts{keep_files} ? $$opts{keep_files} : tempdir(CLEANUP=>1);
    $$opts{path} = $FindBin::RealBin;
    $$opts{bin}  = $FindBin::RealBin;
    $$opts{bin}  =~ s{/test/?$}{};
    if ( $$opts{keep_files} ) { cmd("mkdir -p $$opts{keep_files}"); }
    else 
    { 
        $SIG{TERM} = $SIG{INT} = sub { clean_files($opts); }; 
    }
    return $opts;
}
sub clean_files
{
    my ($opts) = @_;
    print STDERR "Signal caught, cleaning and exiting...\n";
    `rm -rf $$opts{tmp}`;
}
sub _cmd
{
    my ($cmd) = @_;
    my $kid_io;
    my @out;
    my $pid = open($kid_io, "-|");
    if ( !defined $pid ) { error("Cannot fork: $!"); }
    if ($pid) 
    {
        # parent
        @out = <$kid_io>;
        close($kid_io);
    } 
    else 
    {      
        # child
        exec('/bin/bash', '-o','pipefail','-c', $cmd) or error("Cannot execute the command [/bin/sh -o pipefail -c $cmd]: $!");
    }
    return ($? >> 8, join('',@out));
}
sub cmd
{
    my ($cmd) = @_;
    my ($ret,$out) = _cmd($cmd);
    if ( $ret ) { error("The command failed [$ret]: $cmd\n", $out); }
    return $out;
}
sub test_cmd
{
    my ($opts,%args) = @_;
    if ( !exists($args{out}) )
    {
        if ( !exists($args{in}) ) { error("FIXME: expected out or in key\n"); }
        $args{out} = "$args{in}.out";
    }
    my ($package, $filename, $line, $test)=caller(1);
    $test =~ s/^.+:://;

    print "$test:\n";
    print "\t$args{cmd}\n";

    my ($ret,$out) = _cmd("$args{cmd} 2>&1");
    if ( $ret ) { failed($opts,$test); return; }
    if ( $$opts{redo_outputs} && -e "$$opts{path}/$args{out}" )
    {
        rename("$$opts{path}/$args{out}","$$opts{path}/$args{out}.old");
        open(my $fh,'>',"$$opts{path}/$args{out}") or error("$$opts{path}/$args{out}: $!");
        print $fh $out;
        close($fh);
        my ($ret,$out) = _cmd("diff -q $$opts{path}/$args{out} $$opts{path}/$args{out}.old");
        if ( !$ret && $out eq '' ) { unlink("$$opts{path}/$args{out}.old"); }
        else
        {
            print "\tthe expected output changed, saving:\n";
            print "\t  old .. $$opts{path}/$args{out}.old\n";
            print "\t  new .. $$opts{path}/$args{out}\n";
        }
    }
    my $exp = '';
    if ( open(my $fh,'<',"$$opts{path}/$args{out}") )
    {
        my @exp = <$fh>;
        $exp = join('',@exp);
        close($fh);
    }
    elsif ( !$$opts{redo_outputs} ) { failed($opts,$test,"$$opts{path}/$args{out}: $!"); return; }

    if ( $exp ne $out ) 
    { 
        open(my $fh,'>',"$$opts{path}/$args{out}.new") or error("$$opts{path}/$args{out}.new");
        print $fh $out;
        close($fh);
        if ( !-e "$$opts{path}/$args{out}" )
        {
            rename("$$opts{path}/$args{out}.new","$$opts{path}/$args{out}") or error("rename $$opts{path}/$args{out}.new $$opts{path}/$args{out}: $!");
            print "\tthe file with expected output does not exist, creating new one:\n";
            print "\t\t$$opts{path}/$args{out}\n";
        }
        else
        {
            failed($opts,$test,"The outputs differ:\n\t\t$$opts{path}/$args{out}\n\t\t$$opts{path}/$args{out}.new"); 
        }
        return; 
    }
    passed($opts,$test);
}
sub failed
{
    my ($opts,$test,$reason) = @_;
    $$opts{nfailed}++;
    if ( defined $reason ) { print "\n\t$reason"; }
    print "\n.. failed ...\n\n";
}
sub passed
{
    my ($opts,$test) = @_;
    $$opts{nok}++;
    print ".. ok\n\n";
}
sub is_file_newer
{
    my ($afile,$bfile) = @_;
    my (@astat) = stat($afile) or return 0;
    my (@bstat) = stat($bfile) or return 0;
    if ( $astat[9]>$bstat[9] ) { return 1 }
    return 0;
}


# The tests --------------------------

sub test_bgzip
{
    my ($opts,%args) = @_;

    # Create test data: The beggining of each line gives the 0-based offset (including '\n's)
    #
    open(my $fh,'>',"$$opts{tmp}/bgzip.dat") or error("$$opts{tmp}/bgzip.dat: $!");
    my $ntot = 1_000_000;
    my $nwr  = 0;
    while ($nwr < $ntot)
    {
        my $out = sprintf("%d\n", $nwr);
        $nwr += length($out);
        print $fh $out;
    }
    close($fh);
    cmd("cat $$opts{tmp}/bgzip.dat | $$opts{bin}/bgzip -ci -I$$opts{tmp}/bgzip.dat.gz.gzi > $$opts{tmp}/bgzip.dat.gz");

    # Run tests
    my ($test,$out);

    $test = "$$opts{bin}/bgzip -c -b 65272 -s 5 $$opts{tmp}/bgzip.dat.gz";
    print "$test\n";
    $out = cmd($test);
    if ( $out ne '65272' ) { failed($opts,$test,"Expected \"65272\" got \"$out\"\n"); } 
    else { passed($opts,$test); }

    $test = "$$opts{bin}/bgzip -c -b 979200 -s 6 $$opts{tmp}/bgzip.dat.gz";
    print "$test\n";
    $out = cmd($test);
    if ( $out ne '979200' ) { failed($opts,$test,"Expected \"979200\" got \"$out\"\n"); } 
    else { passed($opts,$test); }

    $test = "$$opts{bin}/bgzip -c -b 652804 -s 6 $$opts{tmp}/bgzip.dat.gz";
    print "$test\n";
    $out = cmd($test);
    if ( $out ne '652804' ) { failed($opts,$test,"Expected \"652804\" got \"$out\"\n"); } 
    else { passed($opts,$test); }
}


# For faidx testing: code numbers as base 4 ACGT strings
sub faidx_num_to_seq
{
    my ($dec) = @_;
    my $out = '';
    my @base = qw(A C G T);
    while ( $dec>=0 )
    {
        my $r = $dec % 4;
        $out  = $base[$r] . $out;
        $dec  = int( ($dec - $r) / 4 );
        if ( !$dec ) { last; }
    }
    return $out;
}
sub faidx_seq_to_num
{
    my ($seq) = @_;
    my $out = 0;
    my $len = length($seq);
    my %base = ( A=>0, C=>1, G=>2, T=>3 );
    for (my $i=0; $i<$len; $i++)
    {
        my $b = substr($seq,$i,1);
        $out += $base{$b} * 4**($len-$i-1);
    }
    return $out;
}
sub faidx_wrap
{
    my ($seq,$width) = @_;
    my $out;
    while ( length($seq) )
    {
        $out .= substr($seq,0,60,'');
        $out .= "\n";
    }
    return $out;
}
sub test_faidx
{
    my ($opts,%args) = @_;

    # Create test data: The fake sequence consists of sequence offsets coded
    # into A,C,G,T and separated with Ns. The offsets are 1-based.
    #
    open(my $fh,'>',"$$opts{tmp}/faidx.fa") or error("$$opts{tmp}/faidx.fa: $!");
    my $ntot = 100_000;
    my @dat  = qw(A C G T);
    for (my $seq=1; $seq<=3; $seq++)
    {
        my $nwr = 1;
        my $out = '';
        while ($nwr < $ntot)
        {
            my $tmp = faidx_num_to_seq($nwr) . 'N';
            $out .= $tmp;
            $nwr += length($tmp);
        }
        print $fh ">$seq\n";
        print $fh faidx_wrap($out);
    }
    close($fh);

    # Run tests: index and retrieval from plain text and compressed files
    cmd("$$opts{bin}/samtools faidx $$opts{tmp}/faidx.fa");
    cmd("cat $$opts{tmp}/faidx.fa | $$opts{bin}/bgzip -ci -I $$opts{tmp}/faidx.fa.gz.gzi > $$opts{tmp}/faidx.fa.gz");
    cmd("$$opts{bin}/samtools faidx $$opts{tmp}/faidx.fa.gz");

    for my $reg ('3:11-13','2:998-1003','1:100-104','1:99998-100007')
    {
        for my $file ("$$opts{tmp}/faidx.fa","$$opts{tmp}/faidx.fa.gz")
        {
            my $test = "$$opts{bin}/samtools faidx $file $reg | grep -v '^>'";
            print "$test\n";
            my $seq = cmd($test);
            chomp($seq);
            $seq =~ s/N.*$//;
            my $num = faidx_seq_to_num($seq);
            my $xreg = $reg;
            $xreg =~ s/^[^:]*://;
            $xreg =~ s/-.*$//;
            if ( $num ne $xreg ) { failed($opts,$test,"Expected \"". faidx_num_to_seq($xreg) ."\" got \"$seq\"\n"); } 
            else { passed($opts,$test); }
        }
    }
}

sub test_mpileup
{
    my ($opts,%args) = @_;

    my @bams = ('mpileup.1.bam','mpileup.2.bam','mpileup.3.bam');
    my $ref  = 'mpileup.ref.fa';

    # make a local copy, index the bams and the reference
    open(my $fh,'>',"$$opts{tmp}/mpileup.list") or error("$$opts{tmp}/mpileup.list: $!");
    for my $bam (@bams)
    {
        cmd("cp $$opts{path}/dat/$bam $$opts{tmp}/$bam");
        cmd("$$opts{bin}/samtools index $$opts{tmp}/$bam");
        print $fh "$$opts{tmp}/$bam\n";
    }
    close($fh);
    cmd("cp $$opts{path}/dat/$ref $$opts{tmp}/$ref");
    cmd("$$opts{bin}/bgzip -fi $$opts{tmp}/$ref");
    cmd("$$opts{bin}/samtools faidx $$opts{tmp}/$ref.gz");

    # print "$$opts{bin}samtools mpileup -gb $$opts{tmp}/mpileup.list -f $$opts{tmp}/$args{ref}.gz > $$opts{tmp}/mpileup.bcf\n";
    test_cmd($opts,out=>'dat/mpileup.out.1',cmd=>"$$opts{bin}/samtools mpileup -b $$opts{tmp}/mpileup.list -f $$opts{tmp}/mpileup.ref.fa.gz -r17:100-150 2>/dev/null");
    test_cmd($opts,out=>'dat/mpileup.out.2',cmd=>"$$opts{bin}/samtools mpileup -uvDV -b $$opts{tmp}/mpileup.list -f $$opts{tmp}/mpileup.ref.fa.gz -r17:100-600 2>/dev/null| grep -v ^##samtools | grep -v ^##ref");
    # test that filter mask replaces (not just adds to) default mask
    test_cmd($opts,out=>'dat/mpileup.out.3',cmd=>"$$opts{bin}/samtools mpileup -B --ff 0x14 -f $$opts{tmp}/mpileup.ref.fa.gz -r17:1050-1060 $$opts{tmp}/mpileup.1.bam 2>/dev/null | grep -v mpileup");
}

sub test_usage
{
    my ($opts,%args) = @_;

    my $test = "test_usage";
    print "$test:\n";
    print "\t$args{cmd}\n";
    
    my $command = $args{cmd};
    my $commandpath = $$opts{bin}."/".$command;
    my ($ret,$out) = _cmd("$commandpath 2>&1");
    if ( $out =~ m/\/bin\/bash.*no.*such/i ) { failed($opts,$test,"could not run $commandpath: $out"); return; }

    my @sections = ($out =~ m/(^[A-Za-z]+.*?)(?:(?=^[A-Za-z]+:)|\z)/msg);
    
    my $have_usage = 0;
    my $have_version = 0;
    my $have_subcommands = 0;
    my $usage = "";
    my @subcommands = ();
    foreach my $section (@sections) {
	if ( $section =~ m/^usage/i ) {
	    $have_usage = 1;
	    $section =~ s/^[[:word:]]+[[:punct:][:space:]]*//;
	    $usage = $section;
	} elsif ( $section =~ m/^version/i ) {
	    $have_version = 1;
	} elsif ( $section =~ m/^command/i ) {
	    $have_subcommands = 1;
	    $section =~ s/^[[:word:]]+[[:punct:]]?[[:space:]]*//;
	    $section =~ s/^[[:space:]]+//mg;
	    $section =~ s/^[[:punct:]]+.*?\n//msg;
	    @subcommands = ($section =~ m/^([[:word:]]+)[[:space:]].*/mg);
	}
    }
    
    if ( !$have_usage ) { failed($opts,$test,"did not have Usage:"); return; }
    if ( !$have_version ) { failed($opts,$test,"did not have Version:"); return; }
    if ( !$have_subcommands ) { failed($opts,$test,"did not have Commands:"); return; }

    if ( !($usage =~ m/$command/) ) { failed($opts,$test,"usage did not mention $command"); return; } 
    
    if ( scalar(@subcommands) < 1 ) { failed($opts,$test,"could not parse subcommands"); return; }
    print "\t$command has subcommands: ".join(", ", @subcommands)."\n";

    passed($opts,$test);
    
    # now test subcommand usage as well
    foreach my $subcommand (@subcommands) {
	test_usage_subcommand($opts,%args,subcmd=>$subcommand);
    }
}

sub test_usage_subcommand
{
    my ($opts,%args) = @_;

    my $test = "test_usage_subcommand";
    print "$test:\n";
    print "\t$args{cmd} $args{subcmd}\n";

    my $command = $args{cmd};
    my $subcommand = $args{subcmd};
    my $commandpath = $$opts{bin}."/".$command;
    my ($ret,$out) = _cmd("$commandpath $subcommand 2>&1");
    if ( $out =~ m/\/bin\/bash.*no.*such/i ) { failed($opts,$test,"could not run $commandpath $subcommand: $out"); return; }

    if ( $out =~ m/not.*implemented/is ) { passed($opts,$test,"subcommand indicates it is not implemented"); return; }

    my @sections = ($out =~ m/(^[A-Za-z]+.*?)(?:(?=^[A-Za-z]+:)|\z)/msg);
    
    my $have_usage = 0;
    my $usage = "";
    foreach my $section (@sections) {
	if ( $section =~ m/^usage/i ) {
	    $have_usage = 1;
	    $section =~ s/^[[:word:]]+[[:punct:]]?[[:space:]]*//;
	    $usage = $section;
	}
    }
    
    if ( !$have_usage ) { failed($opts,$test,"did not have Usage:"); return; }

    if ( !($usage =~ m/$command[[:space:]]+$subcommand/) ) { failed($opts,$test,"usage did not mention $command $subcommand"); return; } 
    
    passed($opts,$test);
}

sub add_ur_tags
{
    my ($in, $out, $fasta_location) = @_;

    # Add suitable UR: tag to the source SAM file @SQ lines (for CRAM to use)
    open(my $sam_in, '<', $in) || die "Couldn't open $in : $!\n";
    open(my $sam_out, '>', $out) || die "Couldn't open $out for writing : $!\n";
    while (<$sam_in>) {
	if (/^\@SQ/) {
	    chomp;
	    print $sam_out "$_\tUR:$fasta_location\n"
		|| die "Error writing to $out : $!\n";
	} else {
	    print $sam_out $_ || die "Error writing to $out : $!\n";
	}
    }
    close($sam_in) || die "Error reading $in : $!\n";
    close($sam_out) || die "Error writing to $out : $!\n";    
}

sub convert_flags
{
    my ($in, $out, $style) = @_;

    my @letters = qw(p P u U r R 1 2 s f d);

    open(my $sam_in, '<', $in) || die "Couldn't open $in : $!\n";
    open(my $sam_out, '>', $out) || die "Couldn't open $out for writing : $!\n";
    while (<$sam_in>) {
	if (/^@/) {
	    print $sam_out $_ || die "Error writing to $out : $!\n";
	} else {
	    chomp;
	    my @sam = split(/\t/, $_);
	    if ($style eq 'h') {
		$sam[1] = sprintf("0x%x", $sam[1]);
	    } else {
		my $str;
		for (my $i = 0; $i < @letters; $i++) {
		    if ($sam[1] & (1 << $i)) { $str .= $letters[$i]; }
		}
		$sam[1] = $str;
	    }
	    print $sam_out join("\t", @sam), "\n"
		|| die "Error writing to $out : $!\n";
	}
    }
    close($sam_in) || die "Error reading $in : $!\n";
    close($sam_out) || die "Error writing to $out : $!\n";
}

sub reflen
{
    my ($cigar) = @_;
    
    my $len = 0;
    my %m = ( M => 1, D => 1, N => 1, '=' => 1, X => 1 );
    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
	if (exists($m{$2})) { $len += $1; }
    }
    return $len;
}

sub filter_sam
{
    my ($in, $out, $args) = @_;

    my $no_body   = exists($args->{no_body});
    my $no_header = exists($args->{no_header});
    my $no_m5     = exists($args->{no_m5});
    my $no_sq     = exists($args->{no_sq});
    my $min_map_qual   = $args->{min_map_qual} || 0;
    my $flags_required = $args->{flags_required} || 0;
    my $flags_rejected = $args->{flags_rejected} || 0;
    my $read_groups    = $args->{read_groups};
    my $libraries      = $args->{libraries};
    my $region         = $args->{region};
    my $body_filter = ($flags_required || $flags_rejected || $read_groups
		       || $min_map_qual || $region);

    open(my $sam_in, '<', $in) || die "Couldn't open $in : $!\n";
    open(my $sam_out, '>', $out) || die "Couldn't open $out for writing : $!\n";
    while (<$sam_in>) {
	if (/^@/) {
	    next if ($no_header);
	    if ($libraries && /^\@RG/) {
		my ($id) = /\tID:([^\t]+)/;
		my ($lib) = /\tLB:([^\t]+)/;
		if (exists($libraries->{$lib})) { $read_groups->{$id} = 1; }
	    }
	    next if ($no_sq && /^\@SQ/);
	    if ($no_m5 && /^\@SQ/) {
		s/\tM5:[^\t\n]+//;
	    }
	    print $sam_out $_ || die "Error writing to $out : $!\n";
	} else {
	    next if ($no_body);
	    if ($body_filter) {
		my @sam = split(/\t/, $_);
		next if ($flags_required
			 && ($sam[1] & $flags_required) != $flags_required);
		next if ($flags_rejected && ($sam[1] & $flags_rejected) != 0);
		next if ($min_map_qual && $sam[4] < $min_map_qual);
		if ($read_groups) {
		    my $group = '';
		    for my $i (11 .. $#sam) {
			last if (($group) = $sam[$i] =~ /^RG:Z:(.*)/);
		    }
		    next unless (exists($read_groups->{$group}));
		}
		if ($region) {
		    my $in_range = 0;
		    foreach my $r (@$region) {
			next if ($r->[0] ne $sam[2]);
			next if (@$r > 1
				 && $sam[3] + reflen($sam[5]) - 1 < $r->[1]);
			next if (@$r > 2 && $sam[3] > $r->[2]);
			$in_range = 1;
			last;
		    }
		    next if (!$in_range);
		}
	    }
	    print $sam_out $_ || die "Error writing to $out : $!\n";
	}
    }
    close($sam_in) || die "Error reading $in : $!\n";
    close($sam_out) || die "Error writing to $out : $!\n";    
}

sub run_view_test
{
    my ($opts, %args) = @_;

    printf "\t%-60s", $args{msg};
    my @cmd = ("$$opts{bin}/samtools", 'view');
    if ($args{out} && !$args{redirect}) { push(@cmd, '-o', $args{out}); }
    if ($args{args}) { push(@cmd, @{$args{args}}); }
    
    my $pid = fork();
    unless (defined($pid)) { die "Couldn't fork : $!\n"; }
    unless ($pid) {
	if ($args{stdin}) {
	    open(STDIN, '-|', 'cat', $args{stdin})
		|| die "Couldn't pipe cat $args{stdin} to STDIN : $!\n";
	}
	if ($args{redirect} && $args{out}) {
	    open(STDOUT, '>', $args{out})
		|| die "Couldn't redirect STDOUT to $args{out} : $!\n";
	}
	exec(@cmd) || die "Couldn't exec @cmd : $!\n";
    }
    my $reaped = waitpid($pid, 0);
    my $res = $reaped == $pid && $? == 0 ? 0 : 1;
    
    if (!$res && $args{compare_sam} && $args{out}) {
	# Convert output back to sam and compare
	my $sam_out = "$args{out}.sam";
	my @cmd2 = ("$$opts{bin}/samtools",
		    'view', '-h', '-o', $sam_out, $args{out});
	push(@cmd, '&&', @cmd2); # For the 'Failed command' message below
	$res = system(@cmd2) == 0 ? 0 : 1;
	# Hack $args so the comparison gets done
	$args{compare} = $args{compare_sam};
	$args{out}     = $sam_out;
    }

    if ($res) {
	failed($opts, $args{msg});
    } else {
	if ($args{compare} && $args{out}) {
	    $res = sam_compare($opts, $args{msg}, $args{out}, $args{compare});
	} elsif ($args{compare_bam} && $args{out}) {
	    $res = bam_compare($opts, $args{msg},
			       $args{out}, $args{compare_bam});
	} elsif ($args{compare_count}) {
	    $res = count_compare($opts, $args{msg},
				 $args{out}, $args{compare_count});
	} else {
	    passed($opts, $args{msg});
	}
    }
    print "\tFailed command:\n\t@cmd\n\n" if ($res);
}

sub sam_compare
{
    my ($opts, $msg, $sam1, $sam2) = @_;

    unless (-e $sam1 && -e $sam2) {
	failed($opts, $msg);
	return 1;
    }

    my %hdr1;
    my %hdr2;

    my ($lno1, $lno2) = (0, 0);
    my ($l1, $l2, $ht1, $ht2);

    open(my $f1, '<', $sam1) || die "Couldn't open $sam1: $!\n";
    while ($l1 = <$f1>) {
	$lno1++;
	if (($ht1) = $l1 =~ /^(@\S+)/) {
	    push(@{$hdr1{$ht1}}, $l1);
	} else {
	    last;
	}
    }

    open(my $f2, '<', $sam2) || die "Couldn't open $sam2: $!\n";
    while ($l2 = <$f2>) {
	$lno2++;
	if (($ht2) = $l2 =~ /^(@\S+)/) {
	    push(@{$hdr2{$ht2}}, $l2);
	} else {
	    last;
	}
    }

    while (my ($ht, $h1) = each %hdr1) {
	my $h2 = $hdr2{$ht};
	my $same = $h2 && @$h1 == @$h2;
	if ($same) {
	    for (my $i = 0; $i < @$h1; $i++) {
		if ($h1->[$i] ne $h2->[$i]) {
		    $same = 0;
		    last;
		}
	    }
	}
	if (!$same) {
	    print "\n\tHeader type $ht differs.\n";
	    failed($opts, $msg);
	    close($f1);
	    close($f2);
	    return 1;
	}
    }

    while ($l1 && $l2) {
	last if ($l1 ne $l2);
	$l1 = <$f1>;
	$l2 = <$f2>;
	$lno1++;
	$lno2++;
    }

    close($f1) || die "Error reading $sam1: $!\n";
    close($f2) || die "Error reading $sam2: $!\n";

    if ($l1 || $l2) {
	print "\n\tSAM files differ at $sam1:$lno1 / $sam2:$lno2\n";
	failed($opts, $msg);
	return 1;
    } else {
	passed($opts, $msg);
	return 0;
    }
}

sub bam_compare
{
    my ($opts, $msg, $bam1, $bam2) = @_;

    my $buffer1;
    my $buffer2;
    my $bytes1;
    my $bytes2;
    my $fail = 0;

    open(my $b1, '-|', 'gunzip', '-c', $bam1)
	|| die "Couldn't open pipe to gunzip -c $bam1 : $!\n";
    open(my $b2, '-|', 'gunzip', '-c', $bam2)
	|| die "Couldn't open pipe to gunzip -c $bam2 : $!\n";
    do {
	$bytes1 = read($b1, $buffer1, 65536);
	$bytes2 = read($b2, $buffer2, 65536);
	if (!defined($bytes1)) { die "Error reading $bam1 : $!\n"; }
	if (!defined($bytes2)) { die "Error reading $bam2 : $!\n"; }
	if ($bytes1 != $bytes2 || $buffer1 ne $buffer2) {
	    $fail = 1;
	    last;
	}
    } while ($bytes1 && $bytes2);
    close($b1) || die "Error running gunzip -c $bam1\n";
    close($b2) || die "Error running gunzip -c $bam2\n";
    if ($fail) {
	print "\n\tBAM files $bam1 and $bam2 differ.\n";
	failed($opts, $msg);
	return 1;
    } else {
	passed($opts, $msg);
	return 0;
    }
}

sub count_compare
{
    my ($opts, $msg, $count1, $count2) = @_;
    
    open(my $c1, '<', $count1) || die "Couldn't open $count1 : $!\n";
    my $number1 = <$c1>;
    chomp($number1);
    close($c1) || die "Error reading $count1 : $!\n";

    unless ($number1 =~ /^\d+$/) {
	print "\n\tExpected a number in $count1 but got '$number1'\n";
	failed($opts, $msg);
	return 1;
    }

    my $number2 = 0;
    if ($count2 =~ /^\d+$/) {
	$number2 = $count2;
    } else {
	open(my $c2, '<', $count2) || die "Couldn't open $count2 : $!\n";
	while (<$c2>) {
	    if (!/^@/) { $number2++; }
	}
	close($c2) || die "Error reading $count2 : $!\n";
    }
    if ($number1 != $number2) {
	print "\n\tIncorrect count: Got $number1; expected $number2.\n";
	failed($opts, $msg);
	return 1;
    } else {
	passed($opts, $msg);
	return 0;
    }
}

sub test_view
{
    my ($opts) = @_;

    my $test_name = "test_view";
    print "$test_name:\n";
    
    my $sam_no_ur   = "$$opts{path}/dat/view.001.sam";
    my $sam_with_ur = "$$opts{tmp}/view.001.sam";
    add_ur_tags($sam_no_ur, $sam_with_ur, "$$opts{path}/dat/view.001.fa");

    my $test = 1;

    my $out = "$$opts{tmp}/view_out.001";

    # SAM -> BAM -> SAM
    my $bam_with_ur_out = sprintf("%s.test%02d.bam", $out, $test);
    run_view_test($opts,
		  msg => "$test: SAM -> BAM -> SAM",
		  args => ['-b', $sam_with_ur],
		  out => $bam_with_ur_out,
		  compare_sam => $sam_with_ur);
    $test++;

    # SAM -> uncompressed BAM -> SAM
    my $uncomp_bam = sprintf("%s.test%02d.bam", $out, $test);
    run_view_test($opts,
		  msg => "$test: SAM -> uncompressed BAM",
		  args => ['-u', $sam_with_ur],
		  out => $uncomp_bam,
		  compare_bam => $bam_with_ur_out);
    $test++;
    run_view_test($opts,
		  msg => "$test: uncompressed BAM -> SAM and compare",
		  args => ['-h', $uncomp_bam],
		  out => sprintf("%s.test%02d.sam", $out, $test),
		  compare => $sam_with_ur);
    $test++;

    # SAM -> fast compressed BAM -> SAM
    my $fastcomp_bam = sprintf("%s.test%02d.bam", $out, $test);
    run_view_test($opts,
		  msg => "$test: SAM -> fast compressed BAM",
		  args => ['-1', $sam_with_ur],
		  out => $fastcomp_bam,
		  compare_bam => $bam_with_ur_out);
    $test++;
    run_view_test($opts,
		  msg => "$test: fast compressed BAM -> SAM and compare",
		  args => ['-h', $fastcomp_bam],
		  out => sprintf("%s.test%02d.sam", $out, $test),
		  compare => $sam_with_ur);
    $test++;

    # SAM -> CRAM -> SAM with UR tags
    my $cram_with_ur_out = sprintf("%s.test%02d.cram", $out, $test);
    run_view_test($opts,
		  msg => "$test: SAM -> CRAM -> SAM (UR header tags)",
		  args => ['-C', $sam_with_ur],
		  out => $cram_with_ur_out,
		  compare_sam => $sam_with_ur);
    $test++;

    # SAM -> BAM -> CRAM -> SAM with UR tags
    my $cram_from_bam = sprintf("%s.test%02d.cram", $out, $test);
    run_view_test($opts,
		  msg => "$test: BAM -> CRAM with UR -> SAM",
		  args => ['-C', $bam_with_ur_out],
		  out => $cram_from_bam,
		  compare_sam => $sam_with_ur);
    $test++;

    # SAM -> BAM -> CRAM -> BAM -> SAM with UR tags
    my $bam_from_cram = sprintf("%s.test%02d.bam", $out, $test);
    run_view_test($opts,
		  msg => "$test: CRAM -> BAM with UR",
		  args => ['-b', $cram_from_bam],
		  out => $bam_from_cram,
		  compare_sam => $sam_with_ur);
    $test++;

    # Write to stdout
    run_view_test($opts,
		  msg => "$test: SAM -> SAM via stdout",
		  args => ['-h', $sam_with_ur],
		  out => sprintf("%s.test%02d.sam", $out, $test),
		  redirect => 1,
		  compare => $sam_with_ur);
    $test++;
    run_view_test($opts,
		  msg => "$test: SAM -> BAM via stdout",
		  args => ['-b', $sam_with_ur],
		  out => sprintf("%s.test%02d.bam", $out, $test),
		  redirect => 1,
		  compare_bam => $bam_with_ur_out);
    $test++;
    my $cram_via_stdout = sprintf("%s.test%02d.cram", $out, $test);
    run_view_test($opts,
		  msg => "$test: SAM -> CRAM via stdout",
		  args => ['-C', $sam_with_ur],
		  out => $cram_via_stdout,
		  redirect => 1,
		  compare_sam => $sam_with_ur);
    $test++;

    # Read from stdin
    run_view_test($opts,
		  msg => "$test: SAM from stdin -> SAM",
		  args => ['-h', '-'],
		  out => sprintf("%s.test%02d.sam", $out, $test),
		  stdin => $sam_with_ur,
		  compare => $sam_with_ur);
    $test++;
    run_view_test($opts,
		  msg => "$test: BAM from stdin -> SAM",
		  args => ['-h', '-'],
		  out => sprintf("%s.test%02d.sam", $out, $test),
		  stdin => $bam_with_ur_out,
		  compare => $sam_with_ur);
    $test++;
    run_view_test($opts,
		  msg => "$test: CRAM from stdin -> SAM",
		  args => ['-h', '-'],
		  out => sprintf("%s.test%02d.sam", $out, $test),
		  stdin => $cram_with_ur_out,
		  compare => $sam_with_ur);
    $test++;
    

    # Header only options
    my $sam_header = "$$opts{tmp}/view.001.header.sam";
    filter_sam($sam_with_ur, $sam_header, {no_body => 1});

    run_view_test($opts,
		  msg => "$test: samtools view -H (SAM input)",
		  args => ['-H', $sam_with_ur],
		  out => sprintf("%s.test%02d.header", $out, $test),
		  compare => $sam_header);
    $test++;
    run_view_test($opts,
		  msg => "$test: samtools view -H (BAM input)",
		  args => ['-H', $bam_with_ur_out],
		  out => sprintf("%s.test%02d.header", $out, $test),
		  compare => $sam_header);
    $test++;
    run_view_test($opts,
		  msg => "$test: samtools view -H (CRAM input)",
		  args => ['-H', $cram_with_ur_out],
		  out => sprintf("%s.test%02d.header", $out, $test),
		  compare => $sam_header);
    $test++;
    
    # Body only
    my $sam_body = "$$opts{tmp}/view.001.body.sam";
    filter_sam($sam_with_ur, $sam_body, {no_header => 1});

    run_view_test($opts,
		  msg => "$test: No headers (SAM input)",
		  args => [$sam_with_ur],
		  out => sprintf("%s.test%02d.body", $out, $test),
		  compare => $sam_body);
    $test++;
    run_view_test($opts,
		  msg => "$test: No headers (BAM input)",
		  args => [$bam_with_ur_out],
		  out => sprintf("%s.test%02d.body", $out, $test),
		  compare => $sam_body);
    $test++;
    run_view_test($opts,
		  msg => "$test: No headers (CRAM input)",
		  args => [$cram_with_ur_out],
		  out => sprintf("%s.test%02d.body", $out, $test),
		  compare => $sam_body);
    $test++;

    # -x / -X options
    my $sam_hexflag = "$$opts{tmp}/view.001.hexflag.sam";
    my $sam_strflag = "$$opts{tmp}/view.001.strflag.sam";
    convert_flags($sam_with_ur, $sam_hexflag, 'h');
    convert_flags($sam_with_ur, $sam_strflag, 's');

    run_view_test($opts,
		  msg => "$test: -x hex flags option (SAM input)",
		  args => ['-h', '-x', $sam_with_ur],
		  out => sprintf("%s.test%02d.sam", $out, $test),
		  compare => $sam_hexflag);
    $test++;
    run_view_test($opts,
		  msg => "$test: -X string flags option (SAM input)",
		  args => ['-h', '-X', $sam_with_ur],
		  out => sprintf("%s.test%02d.sam", $out, $test),
		  compare => $sam_strflag);
    $test++;

    # -x / -X without headers
    my $sam_hexflag_body = "$$opts{tmp}/view.001.hexflag.body.sam";
    my $sam_strflag_body = "$$opts{tmp}/view.001.strflag.body.sam";
    filter_sam($sam_hexflag, $sam_hexflag_body, {no_header => 1});
    filter_sam($sam_strflag, $sam_strflag_body, {no_header => 1});
    run_view_test($opts,
		  msg => "$test: -x hex flags option, no headers (SAM input)",
		  args => ['-x', $sam_with_ur],
		  out => sprintf("%s.test%02d.body", $out, $test),
		  compare => $sam_hexflag_body);
    $test++;
    run_view_test($opts,
		  msg => "$test: -X string flags option, no headers (SAM input)",
		  args => ['-X', $sam_with_ur],
		  out => sprintf("%s.test%02d.body", $out, $test),
		  compare => $sam_strflag_body);
    $test++;
    
    # Filter and counting tests.
    
    # Group names file for -R test
    my $fogn = "$$opts{tmp}/view.001.fogn";
    open(my $f, '>', $fogn) || die "Couldn't open $fogn : $!\n";
    print $f "grp1\ngrp3\n" || die "Error writing to $fogn : $!\n";
    close($f) || die "Error writing to $fogn : $!\n";


    my @filter_tests = (
	['req128', {flags_required => 128}, ['-f', 128]],
	['rej128', {flags_rejected => 128}, ['-F', 128]],
	['rej128req2', { flags_rejected => 128, flags_required => 2 },
	 ['-F', 128, '-f', 2]],
	['rg_grp2', { read_groups => { grp2 => 1 }}, ['-r', 'grp2']],
	['rg_fogn', { read_groups => { grp1 => 1, grp3 => 1 }}, ['-R', $fogn]],
	['lib2', { libraries => { 'Library 2' => 1 }}, ['-l', 'Library 2']],
	['mq50',  { min_map_qual => 50 },  ['-q', 50]],
	['mq99',  { min_map_qual => 99 },  ['-q', 99]],
	['mq100', { min_map_qual => 100 }, ['-q', 100]],
	);

    my @filter_inputs = ([SAM  => $sam_with_ur],
			 [BAM  => $bam_with_ur_out],
			 [CRAM => $cram_with_ur_out]);

    foreach my $filter (@filter_tests) {
	my $sam_file = "$$opts{tmp}/view.001.$$filter[0].sam";
	filter_sam($sam_with_ur, $sam_file, $$filter[1]);

	foreach my $ip (@filter_inputs) {

	    # Filter test
	    run_view_test($opts,
			  msg => "$test: Filter @{$$filter[2]} ($$ip[0] input)",
			  args => ['-h', @{$$filter[2]}, $$ip[1]],
			  out => sprintf("%s.test%02d.sam", $out, $test),
			  compare => $sam_file);
	    $test++;
	    
	    # Count test
	    run_view_test($opts,
			  msg => "$test: Count @{$$filter[2]} ($$ip[0] input)",
			  args => ['-c', @{$$filter[2]}, $$ip[1]],
			  out => sprintf("%s.test%02d.sam", $out, $test),
			  redirect => 1,
			  compare_count => $sam_file);
	    $test++;
	}
    }

    # Region query tests
    my $sam_no_ur2 = "$$opts{path}/dat/view.002.sam";
    my $sam_with_ur2 = "$$opts{tmp}/view.002.sam";
    add_ur_tags($sam_no_ur2, $sam_with_ur2, "$$opts{path}/dat/view.002.fa");

    my $bam_with_ur_out2 = sprintf("%s.test%02d.bam", $out, $test);
    run_view_test($opts,
		  msg => "$test: 1bp reads file SAM -> BAM -> SAM",
		  args => ['-b', $sam_with_ur2],
		  out => $bam_with_ur_out2,
		  compare_sam => $sam_with_ur2);
    $test++;
    my $cram_with_ur_out2 = sprintf("%s.test%02d.cram", $out, $test);
    run_view_test($opts,
		  msg => "$test: 1bp reads file SAM -> CRAM -> SAM",
		  args => ['-C', $sam_with_ur2],
		  out => $cram_with_ur_out2,
		  compare_sam => $sam_with_ur2);

    my @region_sams   = ($sam_with_ur2, $sam_with_ur);
    my @region_inputs = ([BAM  => [$bam_with_ur_out2, $bam_with_ur_out]],
			 [CRAM => [$cram_with_ur_out2, $cram_with_ur_out]]);
    # Add indicies
    cmd("'$$opts{bin}/samtools' index '$bam_with_ur_out'");
    cmd("'$$opts{bin}/samtools' index '$cram_with_ur_out'");
    cmd("'$$opts{bin}/samtools' index '$bam_with_ur_out2'");
    cmd("'$$opts{bin}/samtools' index '$cram_with_ur_out2'");

    my $bed1 = "$$opts{path}/dat/view.001.01.bed";
    my $bed2 = "$$opts{path}/dat/view.001.02.bed";
    my $bed3 = "$$opts{path}/dat/view.002.01.bed";
    my $bed4 = "$$opts{path}/dat/view.002.02.bed";
    my $bed1reg = [['ref1', 11, 24], ['ref1', 45, 45], ['ref2', 17, 17]];

    my @region_tests = (
	['1bp_reg1', 0, { region => [['Z']]}, [], ['Z']],
	['1bp_reg2', 0, { region => [['Z', 30]]}, [], ['Z:30']],
	['1bp_reg3', 0, { region => [['Z', 20, 40]]}, [], ['Z:20-40']],
	['reg1', 1, { region => [['ref1']] }, [], ['ref1']],
	['reg2', 1, { region => [['ref1', 15]] }, [], ['ref1:15']],
	['reg3', 1, { region => [['ref1', 15, 45]] }, [], ['ref1:15-45']],
	['reg4', 1, { region => [['ref1', 15, 45], ['ref2']] },
	 [], ['ref1:15-45', 'ref2']],
	['reg5', 1, { region => [['ref1', 15, 45], ['ref2', 16]] },
	 [], ['ref1:15-45', 'ref2:16']],
	['reg6', 1, { region => [['ref1', 15, 45], ['ref2', 16, 31]] },
	 [], ['ref1:15-45', 'ref2:16-31']],
	['reg7', 1, { region => [['ref1', 15, 15], ['ref1', 45, 45]] },
	 [], ['ref1:15-15', 'ref1:45-45']],
	# Regions combined with other filters.
	['reg8', 1, { region => [['ref1']], flags_required => 128 },
	 ['-f', 128], ['ref1']],
	['reg9', 1, { region => [['ref1', 15, 45]], min_map_qual => 50 },
	 ['-q', 50], ['ref1:15-45']],
	['reg10', 1,
	 { region => [['ref1', 15, 45]], read_groups => { grp2 => 1 }},
	 ['-r', 'grp2'], ['ref1:15-45']],

	# Regions from BED files.  Regions here need to be kept in synch.
	# with the .bed files in test/dat.  Note that BED counts from
	# base 0 and ends ranges one past the last base.  Note also that
	# samtools also understands a two-column format that describes a
	# single point, and in that case the position is 1-based.  And
	# even better, you can mix the two in the same file.  Confusing? Yes.
	#
	# $bed1reg is defined above as it's used a lot of times.

	['1bp_bed3', 0, { region => [['Z', 20, 40]]}, ['-L', $bed3], []],
	['1bp_bed4', 0,
	 { region => [['Z', 10, 15], ['Z', 20, 20], ['Z', 25, 30],
		      ['Z', 35, 35], ['Z', 40, 50]]}, ['-L', $bed4], []],
	['bed1', 1, { region => $bed1reg }, ['-L', $bed1], []],
	['bed2', 1, { region => [['ref1', 6, 20]] }, ['-L', $bed2], []],

	# BED file plus region specification.

	['1bp_bed3r1', 0, { region => [['Z', 20, 40]]},
	 ['-L', $bed3], ['Z']],
	['1bp_bed3r2', 0, { region => [['Z', 30, 40]]},
	 ['-L', $bed3], ['Z:30']],
	['1bp_bed3r3', 0, { region => [['Z', 20, 30]]},
	 ['-L', $bed3], ['Z:10-30']],
	['1bp_bed3r4', 0, { region => [['Z', 25, 35]]},
	 ['-L', $bed3], ['Z:25-35']],
	['bed1r1', 1, { region => [['ref1', 11, 16]] },
	 ['-L', $bed1], ['ref1:11-16']],
	['bed1r2', 1,
	 { region => [['ref2', 17, 17]] }, ['-L', $bed1], ['ref2']],

	# BED file plus other filters

	['bed1f1', 1, { region => $bed1reg, read_groups => { grp1 => 1} },
	 ['-L', $bed1, '-r', 'grp1'], []],
	['bed1f2', 1, { region => $bed1reg, flags_required => 128 },
	 ['-L', $bed1, '-f', 128], []],
	['bed1f3', 1, { region => $bed1reg, min_map_qual => 5 },
	 ['-L', $bed1, '-q', 5], []],

	# BED file, region and filters
	['bed1f1', 1,
	 { region => [['ref1', 11, 16]], read_groups => { grp1 => 1}},
	 ['-L', $bed1, '-r', 'grp1'], ['ref1:11-16']],
	);
    foreach my $rt (@region_tests) {
	my $input_sam = $region_sams[$$rt[1]];
	my $sam_file = "$$opts{tmp}/view.$$rt[0].sam";
	filter_sam($input_sam, $sam_file, $$rt[2]);

	foreach my $ip (@region_inputs) {
	    my $input_file = $$ip[1]->[$$rt[1]];
	    # Region test
	    run_view_test($opts,
			  msg => "$test: Region @{$$rt[3]} @{$$rt[4]} ($$ip[0] input)",
			  args => ['-h', @{$$rt[3]}, $input_file, @{$$rt[4]}],
			  out => sprintf("%s.test%02d.sam", $out, $test),
			  compare => $sam_file);
	    $test++;
	    # Count test
	    run_view_test($opts,
			  msg => "$test: Count @{$$rt[3]} @{$$rt[4]} ($$ip[0] input)",
			  args => ['-c', @{$$rt[3]}, $input_file, @{$$rt[4]}],
			  out => sprintf("%s.test%02d.sam", $out, $test),
			  redirect => 1,
			  compare_count => $sam_file);
	    $test++;
	}
    }

    # -T / -t options
    my $sam_no_sq = "$$opts{tmp}/view.001.no_sq.sam";
    filter_sam($sam_no_ur, $sam_no_sq, {no_sq => 1});
    my $sam_no_m5 = "$$opts{tmp}/view.001.no_m5.sam";
    filter_sam($sam_no_ur, $sam_no_m5, {no_m5 => 1});
    
    my $ref_file = "$$opts{path}/dat/view.001.fa";
    my $ref_idx  = "$$opts{path}/dat/view.001.fa.fai";

    # Test SAM output
    foreach my $topt (['-t', $ref_idx], ['-T', $ref_file]) {
	run_view_test($opts,
		      msg => "$test: Add \@SQ with $topt->[0] (SAM output)",
		      args => ['-h', @$topt, $sam_no_sq],
		      out => sprintf("%s.test%02d.sam", $out, $test),
		      compare => $sam_no_m5);
	$test++;
    }

    # Test BAM output.
    foreach my $topt (['-t', $ref_idx], ['-T', $ref_file]) {
	my $bam = sprintf("%s.test%02d.bam", $out, $test);
	run_view_test($opts,
		      msg => "$test: Add \@SQ with $topt->[0] (BAM output)",
		      args => ['-b', @$topt, $sam_no_sq],
		      out => $bam,
		      compare_sam => $sam_no_m5);
	$test++;
    }

    # Don't bother testing CRAM for the moment

    # CIGAR B-operator removal tests.
    my $b_op_sam      = "$$opts{path}/dat/view.003.sam";
    my $b_op_expected = "$$opts{path}/dat/view.003.expected.sam";
    run_view_test($opts,
		  msg => "$test: CIGAR B-operator removal",
		  args => ['-h', '-B', $b_op_sam],
		  out => sprintf("%s.test%02d.sam", $out, $test),
		  compare => $b_op_expected);
    $test++;
    
}
