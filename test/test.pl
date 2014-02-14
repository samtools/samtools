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

    my @files = ('mpileup.1','mpileup.2','mpileup.3');
    my $ref   = 'mpileup.ref.fa';

    # temporary hack: annotate SAMs with absolute reference path, otherwise SAM to CRAM conversion is not working
    for my $file (@files)
    {
        open(my $fh,'>',"$$opts{tmp}/$file.sam") or error("$$opts{tmp}/$file.sam");
        print $fh join("\t", qw(@HD VN:1.0 SO:coordinate)), "\n";
        print $fh join("\t", qw(@SQ SN:17 LN:81195210 M5:a9a06ca09c111789d92723fbf39820f6 AS:NCBI37 SP:Human), "UR:$$opts{path}/dat/$ref"), "\n";
        close($fh) or error("close failed: $$opts{tmp}/$file.sam");
        cmd("cat $$opts{path}/dat/$file.sam >> $$opts{tmp}/$file.sam");
    }

    # make a local copy, create bams, index the bams and the reference
    open(my $fh1,'>',"$$opts{tmp}/mpileup.bam.list") or error("$$opts{tmp}/mpileup.bam.list: $!");
    open(my $fh2,'>',"$$opts{tmp}/mpileup.cram.list") or error("$$opts{tmp}/mpileup.cram.list: $!");
    for my $file (@files)
    {
        cmd("$$opts{bin}/samtools view -b $$opts{tmp}/$file.sam > $$opts{tmp}/$file.bam");
        cmd("$$opts{bin}/samtools view -C $$opts{tmp}/$file.sam > $$opts{tmp}/$file.cram");
        cmd("$$opts{bin}/samtools index $$opts{tmp}/$file.bam");
        cmd("$$opts{bin}/samtools index $$opts{tmp}/$file.cram");
        print $fh1 "$$opts{tmp}/$file.bam\n";
        print $fh2 "$$opts{tmp}/$file.cram\n";
    }
    close($fh1);
    close($fh2);
    cmd("cp $$opts{path}/dat/$ref $$opts{tmp}/$ref");
    cmd("$$opts{bin}/bgzip -fi $$opts{tmp}/$ref");
    cmd("$$opts{bin}/samtools faidx $$opts{tmp}/$ref.gz");

    # print "$$opts{bin}samtools mpileup -gb $$opts{tmp}/mpileup.list -f $$opts{tmp}/$args{ref}.gz > $$opts{tmp}/mpileup.bcf\n";
    test_cmd($opts,out=>'dat/mpileup.out.1',cmd=>"$$opts{bin}/samtools mpileup -b $$opts{tmp}/mpileup.bam.list -f $$opts{tmp}/mpileup.ref.fa.gz -r17:100-150 2>/dev/null");
    test_cmd($opts,out=>'dat/mpileup.out.1',cmd=>"$$opts{bin}/samtools mpileup -b $$opts{tmp}/mpileup.cram.list -f $$opts{tmp}/mpileup.ref.fa.gz -r17:100-150 2>/dev/null");
    test_cmd($opts,out=>'dat/mpileup.out.2',cmd=>"$$opts{bin}/samtools mpileup -uvDV -b $$opts{tmp}/mpileup.bam.list -f $$opts{tmp}/mpileup.ref.fa.gz -r17:100-600 2>/dev/null| grep -v ^##samtools | grep -v ^##ref");
    test_cmd($opts,out=>'dat/mpileup.out.2',cmd=>"$$opts{bin}/samtools mpileup -uvDV -b $$opts{tmp}/mpileup.cram.list -f $$opts{tmp}/mpileup.ref.fa.gz -r17:100-600 2>/dev/null| grep -v ^##samtools | grep -v ^##ref");
    # test that filter mask replaces (not just adds to) default mask
    test_cmd($opts,out=>'dat/mpileup.out.3',cmd=>"$$opts{bin}/samtools mpileup -B --ff 0x14 -f $$opts{tmp}/mpileup.ref.fa.gz -r17:1050-1060 $$opts{tmp}/mpileup.1.bam 2>/dev/null | grep -v mpileup");
    test_cmd($opts,out=>'dat/mpileup.out.3',cmd=>"$$opts{bin}/samtools mpileup -B --ff 0x14 -f $$opts{tmp}/mpileup.ref.fa.gz -r17:1050-1060 $$opts{tmp}/mpileup.1.cram 2>/dev/null | grep -v mpileup");
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
