#!/usr/bin/env perl
#
#    Copyright (C) 2013-2017 Genome Research Ltd.
#
#    Author: Petr Danecek <pd3@sanger.ac.uk>
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

use strict;
use warnings;
use Carp;
use Cwd qw/ abs_path /;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use File::Temp;
use IO::Handle;

my $opts = parse_params();

test_bgzip($opts);
test_faidx($opts);
test_dict($opts);
test_index($opts);
test_index($opts, threads=>2);
test_mpileup($opts);
test_usage($opts, cmd=>'samtools');
test_view($opts);
test_cat($opts);
test_bam2fq($opts);
test_bam2fq($opts, threads=>2);
test_depad($opts);
test_stats($opts);
test_merge($opts);
test_merge($opts, threads=>2);
test_sort($opts);
test_sort($opts, threads=>2);
test_fixmate($opts);
test_fixmate($opts, threads=>2);
test_calmd($opts);
test_calmd($opts, threads=>2);
test_idxstat($opts);
test_quickcheck($opts);
test_reheader($opts);
test_addrprg($opts);
test_addrprg($opts, threads=>2);
test_markdup($opts);
test_markdup($opts, threads=>2);


print "\nNumber of tests:\n";
printf "    total            .. %d\n", $$opts{nok}+$$opts{nfailed}+$$opts{nxfail}+$$opts{nxpass};
printf "    passed           .. %d\n", $$opts{nok};
printf "    failed           .. %d\n", $$opts{nfailed};
printf "    expected failure .. %d\n", $$opts{nxfail};
printf "    unexpected pass  .. %d\n", $$opts{nxpass};
print "\n";

# Exit non-zero if there is a failure or an unexpected pass.  In the case
# of an unexpected pass, the test script itself is at fault and should
# be updated to expect a pass instead of failure.
exit ($$opts{nfailed} > 0 || $$opts{nxpass} > 0);

#--------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    print
        "About: samtools/htslib consistency test script\n",
        "Usage: test.pl [OPTIONS]\n",
        "Options:\n",
        "   -e, --exec <tool>=[<path>]      Path to use for specified tool executable.\n",
        "   -r, --redo-outputs              Recreate expected output files.\n",
        "   -t, --temp-dir <path>           When given, temporary files will not be removed.\n",
        "   -h, -?, --help                  This help message.\n",
        "\n";
    exit 1;
}

sub tempfile {
    my ($fh, $name) = File::Temp::tempfile(@_);
    if (wantarray) {
	if ($^O =~ /^(?:msys|MSWin32)/) {
	    $name = abs_path($name);
	}
	return ($fh, $name);
    }
    return $fh;
}

sub cygpath {
    my ($path) = @_;
    $path = `cygpath -m $path`;
    $path =~ s/\r?\n//;
    return $path
}

sub tempdir
{
    my $dir = File::Temp::tempdir(@_);
    if ($^O =~ /^msys/) {
        $dir = cygpath($dir);
    } elsif ($^O eq 'MSWin32') {
	$dir =~ s/\\/\//g;
    }
    return $dir;
}

sub parse_params
{
    my $opts = { bgzip=>"bgzip", keep_files=>0, nok=>0, nfailed=>0, nxfail => 0, nxpass => 0 };
    my $help;
    Getopt::Long::Configure('bundling');
    my $ret = GetOptions (
            'e|exec=s' => sub {
		my ($tool, $path) = split /=/, $_[1];
		if ($^O eq 'MSWin32' && $path !~ /\.exe$/) {
		    $path .= '.exe';
		}
		$$opts{$tool} = abs_path($path) if $path;
	},
            't|temp-dir:s' => \$$opts{keep_files},
            'r|redo-outputs' => \$$opts{redo_outputs},
            'h|?|help' => \$help
            );
    if ( !$ret or $help ) { error(); }
    $$opts{tmp}  = $$opts{keep_files} ? $$opts{keep_files} : tempdir(CLEANUP => 1);
    $$opts{path} = $FindBin::RealBin;
    $$opts{bin}  = $FindBin::RealBin;
    if ($^O =~ /^msys/) {
	$$opts{path} = cygpath($$opts{path});
	$$opts{bin}  = cygpath($$opts{bin});
    }
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
    my ($cmd, $args) = @_;
    my $kid_io;
    my $out;
    my $err;
    my ($err_fh, $err_filename) = tempfile(UNLINK => 1);
    my $pid = open($kid_io, "-|", 'bash', '-o','pipefail','-c', "($cmd) 2> $err_filename");
    if ( !defined $pid ) { error("Cannot fork: $!"); }
    if ($pid)
    {
        # parent
	binmode($kid_io);
	binmode($err_fh);
	local $/;
        $out = <$kid_io>;
        close($kid_io);
        my $child_retval = $?;
        $err = <$err_fh>;
        close ($err_fh);
	$out =~ s/\x0d\x0a/\x0a/g if (!$args->{binary});
	$err =~ s/\x0d\x0a/\x0a/g if (!$args->{binary_err});
	$out =~ s/e([-+])0(\d\d)/e$1$2/g if ($args->{exp_fix});

        return ($child_retval, $out, $err);
    }
}
sub cmd
{
    my ($cmd, $args) = @_;
    my ($ret,$out,$err) = _cmd($cmd, $args);
    if ( $ret ) { error("The command failed [$ret]: $cmd\n", "out:$out\n", "err:$err\n"); }
    return $out;
}
# test harness for a command
# %args out=> expected output
#       err=> expected stderr output (optional)
#       cmd=> command to test
#       expect_fail=> as per passed()/failed()
#       want_fail=> consider passed() if cmd() returns non-zero
#       out_map => map output filenames to their expected result file (can be used alongside out)
#       hskip => number of header lines to ignore during diff
sub test_cmd
{
    my ($opts,%args) = @_;
    if ( !exists($args{out}) && !exists($args{out_map}) )
    {
        if ( !exists($args{in}) ) { error("FIXME: expected out or in key\n"); }
        $args{out} = "$args{in}.out";
    }
    my ($package, $filename, $line, $test)=caller(1);
    $test =~ s/^.+:://;

    print "$test:\n";
    print "\t$args{cmd}\n";

    my ($ret,$out,$err) = _cmd("$args{cmd}", \%args);
    if ( $args{want_fail}? ($ret == 0) : ($ret != 0) ) { failed($opts,%args,msg=>$test,reason=>"ERR: $err"); return; }
    if ( $$opts{redo_outputs} && -e "$$opts{path}/$args{out}" )
    {
        rename("$$opts{path}/$args{out}","$$opts{path}/$args{out}.old");
        open(my $fh,'>',"$$opts{path}/$args{out}") or error("$$opts{path}/$args{out}: $!");
	binmode($fh);
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
    if ( $$opts{redo_outputs} && -e "$$opts{path}/$args{err}" )
    {
        rename("$$opts{path}/$args{err}","$$opts{path}/$args{err}.old");
        open(my $fh,'>',"$$opts{path}/$args{err}") or error("$$opts{path}/$args{err}: $!");
	binmode($fh);
        print $fh $err;
        close($fh);
        my ($ret,$out) = _cmd("diff -q $$opts{path}/$args{err} $$opts{path}/$args{err}.old");
        if ( !$ret && $err eq '' ) { unlink("$$opts{path}/$args{err}.old"); }
        else
        {
            print "\tthe expected stderr output changed, saving:\n";
            print "\t  old .. $$opts{path}/$args{err}.old\n";
            print "\t  new .. $$opts{path}/$args{err}\n";
        }
    }

    # check stdout
    my $exp = '';
    if ( open(my $fh,'<',"$$opts{path}/$args{out}") )
    {
	binmode($fh);
	local $/;
        $exp = <$fh>;
        close($fh);
    }
    elsif ( !$$opts{redo_outputs} ) { failed($opts,%args,msg=>$test,reason=>"$$opts{path}/$args{out}: $!"); return; }

    if ( $exp ne $out )
    {
        open(my $fh,'>',"$$opts{path}/$args{out}.new") or error("$$opts{path}/$args{out}.new");
	binmode($fh);
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
            failed($opts,%args,msg=>$test,reason=>"The outputs stdout differ:\n\t\t$$opts{path}/$args{out}\n\t\t$$opts{path}/$args{out}.new");
        }
        return;
    }
    # check stderr
    if ( exists($args{err}) ) {
        my $exp_err = '';
        if ( open(my $fh,'<',"$$opts{path}/$args{err}") )
        {
	    binmode($fh);
            my @exp = <$fh>;
            $exp_err = join('',@exp);
            close($fh);
        }
        elsif ( !$$opts{redo_outputs} ) { failed($opts,%args,msg=>$test,reason=>"$$opts{path}/$args{err}: $!"); return; }

        if ( $exp_err ne $err )
        {
            open(my $fh,'>',"$$opts{path}/$args{err}.new") or error("$$opts{path}/$args{err}.new");
	    binmode($fh);
            print $fh $err;
            close($fh);
            if ( !-e "$$opts{path}/$args{err}" )
            {
                rename("$$opts{path}/$args{err}.new","$$opts{path}/$args{err}") or error("rename $$opts{path}/$args{err}.new $$opts{path}/$args{err}: $!");
                print "\tthe file with expected output does not exist, creating new one:\n";
                print "\t\t$$opts{path}/$args{err}\n";
            }
            else
            {
                failed($opts,%args,msg=>$test,reason=>"The outputs stderr differ:\n\t\t$$opts{path}/$args{err}\n\t\t$$opts{path}/$args{err}.new");
            }
            return;
        }
    }

    # check other output files
    if( exists($args{out_map}) )
    {
        while ( my($out_actual, $out_expected) = each %{$args{out_map}} )
        {
            my $exp = '';
            if ( open(my $fh,'<',"$$opts{path}/$out_expected") )
            {
		binmode($fh);
		local $/;
                $exp = <$fh>;
                close($fh);
            }
            elsif ( !$$opts{redo_outputs} ) { failed($opts,%args,msg=>$test,reason=>"$$opts{path}/$out_expected: $!"); return; }

            my $out = '';
            if ( open(my $fh,'<',"$$opts{path}/$out_actual") )
            {
		binmode($fh);
                if( exists($args{hskip}) ){
                    # Strip hskip lines to allow match to the expected output
                    for (my $i = 0; $i < $args{hskip}; $i++) {
			my $ignore = <$fh>;
		    }
                }
		local $/;
                $out = <$fh>;
		$out =~ s/\x0d\x0a/\x0a/g if (!$args{binary});
		$out =~ s/e\+0(\d\d)/e+$1/g if ($args{exp_fix});
                close($fh);
            }
            elsif ( !$$opts{redo_outputs} ) { failed($opts,%args,msg=>$test,reason=>"$$opts{path}/$out_actual: $!"); return; }

            if ( $exp ne $out )
            {
                open(my $fh,'>',"$$opts{path}/$out_expected.new") or error("$$opts{path}/$out_expected.new");
		binmode($fh);
                print $fh $out;
                close($fh);
                if ( !-e "$$opts{path}/$out_expected" )
                {
                    rename("$$opts{path}/$out_expected.new","$$opts{path}/$out_expected") or error("rename $$opts{path}/$out_expected.new $$opts{path}/$out_expected: $!");
                    print "\tthe file with expected output does not exist, creating new one:\n";
                    print "\t\t$$opts{path}/$out_expected\n";
                }
                else
                {
                    failed($opts,%args,msg=>$test,reason=>"The output files differ:\n\t\t$$opts{path}/$out_expected\n\t\t$$opts{path}/$out_expected.new");
                }
                return;
            }
	    unlink("$$opts{path}/$out_actual");
        }
    }
    passed($opts,%args,msg=>$test);
}

# Record the success or failure of a test.  $opts is the global settings hash;
# %args is a list of key => value pairs that may include:
#  msg         => Message describing the test (currently unused)
#  reason      => Description of why the test failed (only for failed())
#  expect_fail => If set, records failed as xfail and passed as xpass
sub failed
{
    my ($opts, %args) = @_;
    my $reason = exists $args{reason}? "\t$args{reason}\n" : "";
    if (!$args{expect_fail}) {
        $$opts{nfailed}++;
        print "\n"; STDOUT->flush();
        print STDERR "$reason.. failed ...\n"; STDERR->flush();
        print "\n";
    }
    else { $$opts{nxfail}++; print "\n$reason.. expected failure\n\n"; }
}
sub passed
{
    my ($opts, %args) = @_;
    if (!$args{expect_fail}) { $$opts{nok}++; print ".. ok\n\n"; }
    else {
        $$opts{nxpass}++;
        STDOUT->flush();
        print STDERR ".. unexpected pass\n"; STDERR->flush();
        print "\n";
    }
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

    # Create test data: The beginning of each line gives the 0-based offset (including '\n's)
    #
    open(my $fh,'>',"$$opts{tmp}/bgzip.dat") or error("$$opts{tmp}/bgzip.dat: $!");
    binmode($fh);
    my $ntot = 1_000_000;
    my $nwr  = 0;
    while ($nwr < $ntot)
    {
        my $out = sprintf("%d\n", $nwr);
        $nwr += length($out);
        print $fh $out;
    }
    close($fh);
    cmd("cat $$opts{tmp}/bgzip.dat | $$opts{bgzip} -ci -I$$opts{tmp}/bgzip.dat.gz.gzi > $$opts{tmp}/bgzip.dat.gz");

    # Run tests
    my ($test,$out);

    $test = "$$opts{bgzip} -c -b 65272 -s 5 $$opts{tmp}/bgzip.dat.gz";
    print "$test\n";
    $out = cmd($test);
    if ( $out ne '65272' ) { failed($opts,msg=>$test,reason=>"Expected \"65272\" got \"$out\"\n"); }
    else { passed($opts,msg=>$test); }

    $test = "$$opts{bgzip} -c -b 979200 -s 6 $$opts{tmp}/bgzip.dat.gz";
    print "$test\n";
    $out = cmd($test);
    if ( $out ne '979200' ) { failed($opts,msg=>$test,reason=>"Expected \"979200\" got \"$out\"\n"); }
    else { passed($opts,msg=>$test); }

    $test = "$$opts{bgzip} -c -b 652804 -s 6 $$opts{tmp}/bgzip.dat.gz";
    print "$test\n";
    $out = cmd($test);
    if ( $out ne '652804' ) { failed($opts,msg=>$test,reason=>"Expected \"652804\" got \"$out\"\n"); }
    else { passed($opts,msg=>$test); }
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
    cmd("cat $$opts{tmp}/faidx.fa | $$opts{bgzip} -ci -I $$opts{tmp}/faidx.fa.gz.gzi > $$opts{tmp}/faidx.fa.gz");
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
            if ( $num ne $xreg ) { failed($opts,msg=>$test,reason=>"Expected \"". faidx_num_to_seq($xreg) ."\" got \"$seq\"\n"); }
            else { passed($opts,msg=>$test); }
        }
    }
}

sub test_dict
{
    my ($opts,%args) = @_;
    cmd("cat $$opts{path}/dat/dict.fa | $$opts{bgzip} -c > $$opts{tmp}/dict.fa.gz");
    test_cmd($opts,out=>'dat/dict.out',cmd=>"$$opts{bin}/samtools dict -a hf37d5 -s 'Homo floresiensis' -u ftp://orthanc.org/hf37d5.fa.gz $$opts{path}/dat/dict.fa");
    test_cmd($opts,out=>'dat/dict.out',cmd=>"$$opts{bin}/samtools dict -a hf37d5 -s 'Homo floresiensis' -u ftp://orthanc.org/hf37d5.fa.gz $$opts{tmp}/dict.fa.gz");
    test_cmd($opts,out=>'dat/dict.out',cmd=>"cat $$opts{path}/dat/dict.fa | $$opts{bin}/samtools dict -a hf37d5 -s 'Homo floresiensis' -u ftp://orthanc.org/hf37d5.fa.gz");
}

sub test_index
{
    my ($opts,%args) = @_;
    my $threads = exists($args{threads}) ? " -@ $args{threads}" : "";

    cmd("$$opts{bin}/samtools view${threads} -b $$opts{path}/dat/large_chrom.sam > $$opts{tmp}/large_chrom.bam");
    # test_cmd($opts,out=>'dat/empty.expected',err=>'dat/large_chrom_bai_index.err',cmd=>"$$opts{bin}/samtools index${threads} $$opts{tmp}/large_chrom.bam",want_fail=>1,expect_fail=>1); # command should fail and give an error message, but isn't at the moment
    cmd("$$opts{bin}/samtools index${threads} -c $$opts{tmp}/large_chrom.bam");
    test_cmd($opts,out=>'dat/large_chrom.out',cmd=>"$$opts{bin}/samtools view${threads} $$opts{tmp}/large_chrom.bam ref2");
    test_cmd($opts,out=>'dat/large_chrom.out',cmd=>"$$opts{bin}/samtools view${threads} $$opts{tmp}/large_chrom.bam ref2:1-541556283");
    test_cmd($opts,out=>'dat/test_input_1_a.bam.bai.expected',cmd=>"$$opts{bin}/samtools index${threads} $$opts{path}/dat/test_input_1_a.bam && cat $$opts{path}/dat/test_input_1_a.bam.bai",binary=>1);
}

sub test_mpileup
{
    my ($opts,%args) = @_;

    my @files = ('mpileup.1','mpileup.2','mpileup.3');
    my $ref   = 'mpileup.ref.fa';

    # make a local copy, create bams, index the bams and the reference
    open(my $fh1,'>',"$$opts{tmp}/mpileup.bam.list") or error("$$opts{tmp}/mpileup.bam.list: $!");
    open(my $fh2,'>',"$$opts{tmp}/mpileup.cram.list") or error("$$opts{tmp}/mpileup.cram.list: $!");
    open(my $fh3,'>',"$$opts{tmp}/mpileup.bam.urllist") or error("$$opts{tmp}/mpileup.bam.urllist: $!");
    open(my $fh4,'>',"$$opts{tmp}/mpileup.cram.urllist") or error("$$opts{tmp}/mpileup.cram.urllist: $!");
    for my $file (@files)
    {
        cmd("$$opts{bin}/samtools view -b $$opts{path}/dat/$file.sam > $$opts{tmp}/$file.bam");
        cmd("$$opts{bin}/samtools view -C -T $$opts{path}/dat/mpileup.ref.fa $$opts{path}/dat/$file.sam > $$opts{tmp}/$file.cram");
        cmd("$$opts{bin}/samtools index $$opts{tmp}/$file.bam");
        cmd("$$opts{bin}/samtools index $$opts{tmp}/$file.cram");
        print $fh1 "$$opts{tmp}/$file.bam\n";
        print $fh2 "$$opts{tmp}/$file.cram\n";
	my $atmp = $^O =~ /^msys/ ? cygpath($$opts{tmp}) : abs_path($$opts{tmp});
	unless ($atmp =~ /^\//) { $atmp = "/$atmp"; }
	print $fh3 "file://$atmp/$file.bam\n";
	print $fh4 "file://$atmp/$file.cram\n";
   }
    close($fh1);
    close($fh2);
    close($fh3);
    close($fh4);
    cmd("cp $$opts{path}/dat/$ref $$opts{tmp}/$ref");
    cmd("$$opts{bgzip} -fi $$opts{tmp}/$ref");
    cmd("$$opts{bin}/samtools faidx $$opts{tmp}/$ref.gz");

    # print "$$opts{bin}samtools mpileup -gb $$opts{tmp}/mpileup.list -f $$opts{tmp}/$args{ref}.gz > $$opts{tmp}/mpileup.bcf\n";
    test_cmd($opts,out=>'dat/mpileup.out.1',err=>'dat/mpileup.err.1',cmd=>"$$opts{bin}/samtools mpileup -b $$opts{tmp}/mpileup.bam.list -f $$opts{tmp}/mpileup.ref.fa.gz -r17:100-150");
    test_cmd($opts,out=>'dat/mpileup.out.1',err=>'dat/mpileup.err.1',cmd=>"$$opts{bin}/samtools mpileup -b $$opts{tmp}/mpileup.cram.list -f $$opts{tmp}/mpileup.ref.fa.gz -r17:100-150");
    test_cmd($opts,out=>'dat/mpileup.out.1',err=>'dat/mpileup.err.1',cmd=>"$$opts{bin}/samtools mpileup -b $$opts{tmp}/mpileup.bam.urllist -f $$opts{tmp}/mpileup.ref.fa.gz -r17:100-150");
    test_cmd($opts,out=>'dat/mpileup.out.1',err=>'dat/mpileup.err.1',cmd=>"$$opts{bin}/samtools mpileup -b $$opts{tmp}/mpileup.cram.urllist -f $$opts{tmp}/mpileup.ref.fa.gz -r17:100-150");
    test_cmd($opts,out=>'dat/mpileup.out.2',cmd=>"$$opts{bin}/samtools mpileup -uvDV -b $$opts{tmp}/mpileup.bam.list -f $$opts{tmp}/mpileup.ref.fa.gz -r17:100-600| grep -v ^##samtools | grep -v ^##ref");
    test_cmd($opts,out=>'dat/mpileup.out.2',cmd=>"$$opts{bin}/samtools mpileup -uvDV -b $$opts{tmp}/mpileup.cram.list -f $$opts{tmp}/mpileup.ref.fa.gz -r17:100-600| grep -v ^##samtools | grep -v ^##ref");
    test_cmd($opts,out=>'dat/mpileup.out.4',cmd=>"$$opts{bin}/samtools mpileup -uv -t DP,DPR,DV,DP4,INFO/DPR,SP -b $$opts{tmp}/mpileup.cram.list -f $$opts{tmp}/mpileup.ref.fa.gz -r17:100-600| grep -v ^##samtools | grep -v ^##ref");
    test_cmd($opts,out=>'dat/mpileup.out.4',cmd=>"$$opts{bin}/samtools mpileup -uv -t DP,DPR,DV,DP4,INFO/DPR,SP -b $$opts{tmp}/mpileup.cram.list -f $$opts{tmp}/mpileup.ref.fa.gz -r17:100-600| grep -v ^##samtools | grep -v ^##ref");
    # test that filter mask replaces (not just adds to) default mask
    test_cmd($opts,out=>'dat/mpileup.out.3',cmd=>"$$opts{bin}/samtools mpileup -B --ff 0x14 -f $$opts{tmp}/mpileup.ref.fa.gz -r17:1050-1060 $$opts{tmp}/mpileup.1.bam | grep -v mpileup");
    test_cmd($opts,out=>'dat/mpileup.out.3',cmd=>"$$opts{bin}/samtools mpileup -B --ff 0x14 -f $$opts{tmp}/mpileup.ref.fa.gz -r17:1050-1060 $$opts{tmp}/mpileup.1.cram | grep -v mpileup");
}

sub test_usage
{
    my ($opts,%args) = @_;

    my $test = "test_usage";
    print "$test:\n";
    print "\t$args{cmd}\n";

    my $tty_input;
    if (-t) {
        $args{redirection} = "";  # no redirection necessary
    }
    elsif (eval { require IO::Pty; $tty_input = new IO::Pty; }) {
        # ensure stdin is a terminal, so that subcommands display their usage
        $args{redirection} = "<'" . $tty_input->ttyname . "'";
    }
    else {
        warn "$0: no IO::Pty module or can't open pty; skipping usage tests\n";
        return;
    }

    my $command = $args{cmd};
    my $commandpath = $$opts{bin}."/".$command;
    my ($ret,$out,$err) = _cmd("$commandpath $args{redirection}");
    if ( $err =~ m/\/bin\/bash.*no.*such/i ) { failed($opts,msg=>$test,reason=>"could not run $commandpath: $out"); return; }

    my @sections = ($err =~ m/(^[A-Za-z]+.*?)(?:(?=^[A-Za-z]+:)|\z)/msg);

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

    if ( !$have_usage ) { failed($opts,msg=>$test,reason=>"did not have Usage:"); return; }
    if ( !$have_version ) { failed($opts,msg=>$test,reason=>"did not have Version:"); return; }
    if ( !$have_subcommands ) { failed($opts,msg=>$test,reason=>"did not have Commands:"); return; }

    if ( !($usage =~ m/$command/) ) { failed($opts,msg=>$test,reason=>"usage did not mention $command"); return; }

    if ( scalar(@subcommands) < 1 ) { failed($opts,msg=>$test,reason=>"could not parse subcommands"); return; }
    print "\t$command has subcommands: ".join(", ", @subcommands)."\n";

    passed($opts,msg=>$test);

    # now test subcommand usage as well
    foreach my $subcommand (@subcommands) {
	# Under msys the isatty function fails to recognise the terminal.
	# Skip these tests for now.
	next if ($^O =~ /^msys/ && $subcommand =~ /^(dict|sort|stats|view)$/);
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
    my ($ret,$out,$err) = _cmd("$commandpath $subcommand $args{redirection}");

    if ( $err =~ m/\/bin\/bash.*no.*such/i ) { failed($opts,msg=>$test,reason=>"could not run $commandpath $subcommand: $out"); return; }

    if ( $err =~ m/not.*implemented/is ) { failed($opts,msg=>$test,reason=>"subcommand indicates it is not implemented",expect_fail=>1); return; }

    my $have_usage = 0;
    my $usage = "";
    my @sections = ($out =~ m/(^[A-Za-z]+.*?)(?:(?=^[A-Za-z]+:)|\z)/msg);
    foreach my $section (@sections) {
        if ( $section =~ m/^usage/i ) {
            $have_usage = 1;
            $section =~ s/^[[:word:]]+[[:punct:]]?[[:space:]]*//;
            $usage = $section;
        }
    }
    @sections = ($err =~ m/(^[A-Za-z]+.*?)(?:(?=^[A-Za-z]+:)|\z)/msg);
    foreach my $section (@sections) {
        if ( $section =~ m/^usage/i ) {
            $have_usage = 2;
            $section =~ s/^[[:word:]]+[[:punct:]]?[[:space:]]*//;
            $usage = $section;
        }
    }

    my $fail = 0;
    if ( !$have_usage ) { failed($opts,msg=>$test,reason=>"did not have Usage:"); $fail = 1; }
    elsif ( $have_usage == 2 ) { failed($opts,msg=>$test,reason=>"Usage on stderr rather than stdout",expect_fail=>1); $fail = 1; }

    if ( !($usage =~ m/$command[[:space:]]+$subcommand/) ) { failed($opts,msg=>$test,reason=>"usage did not mention $command $subcommand"); $fail = 1; }

    passed($opts,msg=>$test) if $fail == 0;
}

# Add UR tags to @SQ lines in a SAM file.
# UR points to a local file.  As this could be anywhere, it is not possible
# to put UR tags into the test SAM files that are checked in to version
# control.  Instead this is used to add them when running the tests.
#
# $in  is the SAM file to read.
# $out is the SAM file to write.
# $fasta_location is the filename to put in the UR tag.

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

# Calculates the number of reference bases consumed by a CIGAR string.
# $cigar is the CIGAR string.
# returns the number of reference bases.

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

# Calculates the number of query bases consumed by a CIGAR string.
# $cigar is the CIGAR string.
# returns the number of query bases.

sub querylen
{
    my ($cigar) = @_;

    my $len = 0;
    my %m = ( M => 1, I => 1, S => 1, '=' => 1, X => 1 );
    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
        if (exists($m{$2})) { $len += $1; }
    }
    return $len;
}

# Filter SAM files.
# This is used to convert a test SAM file into one that has been modified
# in a number of ways.  This is then used to compare with the output of
# samtools view, to see if they are the same.
#
# $in  is the SAM file to read.
# $out is the SAM file to write.
# $args is a hash ref indicating which filters to apply:
#   $args->{no_body}        outputs only the header (-H option)
#   $args->{no_header}      outputs only the alignments (no -h)
#   $args->{no_m5}          removes M5: tags from header @SQ lines
#   $args->{no_sq}          removes all @SQ lines from the header
#   $args->{min_map_qual}   minimum mapping quality it output (-q option)
#   $args->{flags_required} bits which must be set in flags (-f option)
#   $args->{flags_rejected} bits which must not be set in flags (-F option)
#   $args->{read_groups}    hash of read groups to output (-r or -R)
#   $args->{libraries}      hash of libraries to output (-l)
#   $args->{region}         region list to output (-L and region list)
#   $args->{strip_tags}     hash of tags to strip from alignments (-x)
#   $args->{min_qlen}       minimum query length to output (-m)
#
# The region list is a reference to an array of region definitions.  Each
# region definition is itself a reference to an array with between one and
# three elements - reference, start position and end position.  The positions
# are inclusive and count from base 1.
#
# *  If one element is present, then all the reads aligned to the given
#    reference are output.
# *  If two elements are present, then all the reads on the given reference
#    that end on or after the start position are output.
# *  If three elements are present, then all the reads on the given reference
#    that end on or after the start position and start on or before the
#    end position are output.
#
# N.B. The region list is just a simple filter, so does not work in
# exactly the same way as samtools view which treats the regions as index
# queries.  This means samtools view will output the same alignment more than
# once if given overlapping regions, while filter_sam will not.  It will
# work in a similar manner to samtools for BED file input though as samtools
# implements BED file regions as a filter rather than a range request.
#
# The BED file support in samtools view is a bit complicated as it actually
# allows two different formats, which can both be mixed in the same file.
#
# The standard BED format has three or more whitespace-separated elements
# on each line, the first three corresponding to reference name, start
# and end positions.  Unlike samtools, BED counts bases from zero, and the
# end position is exclusive, i.e. the base after the desired range.
# To convert BED format to a range definition:
# ref S E        becomes:   ['ref', S + 1, E]
#
# The other format has only two elements (not legal in a normal BED file).
# They define a reference and a single base position.  Unfortunately this
# format counts from base one.  To convert this to a filter_sam range:
# ref P          becomes:   ['ref', P, P]

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
    my $strip_tags     = $args->{strip_tags};
    my $min_qlen       = $args->{min_qlen} || 0;
    my $body_filter = ($flags_required || $flags_rejected || $read_groups
                       || $min_map_qual || $libraries || $region
                       || $strip_tags || $min_qlen);
    my $lib_read_groups = $libraries ? {} : undef;

    open(my $sam_in, '<', $in) || die "Couldn't open $in : $!\n";
    open(my $sam_out, '>', $out) || die "Couldn't open $out for writing : $!\n";
    while (<$sam_in>) {
        chomp;
        if (/^@/) {
            next if ($no_header);
            if ($libraries && /^\@RG/) {
                my ($id) = /\tID:([^\t]+)/;
                my ($lib) = /\tLB:([^\t]+)/;
                if (exists($libraries->{$lib||""})) {
                    $lib_read_groups->{$id} = 1;
                }
            }
            if ($read_groups && /^\@RG/) {
                my ($id) = /\tID:([^\t]+)/;
                next if (!exists($read_groups->{$id||""}));
            }
            next if ($no_sq && /^\@SQ/);
            if ($no_m5 && /^\@SQ/) {
                s/\tM5:[^\t\n]+//;
            }
            print $sam_out "$_\n" || die "Error writing to $out : $!\n";
        } else {
            next if ($no_body);
            if ($body_filter) {
                my @sam = split(/\t/, $_);
                next if ($flags_required
                         && ($sam[1] & $flags_required) != $flags_required);
                next if ($flags_rejected && ($sam[1] & $flags_rejected) != 0);
                next if ($min_map_qual && $sam[4] < $min_map_qual);
                if ($read_groups || $lib_read_groups) {
                    my $group = '';
                    for my $i (11 .. $#sam) {
                        last if (($group) = $sam[$i] =~ /^RG:Z:(.*)/);
                    }
                    next if ($read_groups && !exists($read_groups->{$group||""}));
                    next if ($lib_read_groups
                             && !exists($lib_read_groups->{$group||""}));
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
                if ($min_qlen > 0) {
                    next if (querylen($sam[5]) < $min_qlen);
                }
                if ($strip_tags) {
                    my $stripped = 0;
                    for (my $i = $#sam; $i >= 11; --$i) {
                        if (exists($strip_tags->{substr($sam[$i], 0, 2)})) {
                            $stripped = 1;
                            splice(@sam, $i, 1);
                        }
                    }
                    if ($stripped) { $_ = join("\t", @sam); }
                }
            }
            print $sam_out "$_\n" || die "Error writing to $out : $!\n";
        }
    }
    close($sam_in) || die "Error reading $in : $!\n";
    close($sam_out) || die "Error writing to $out : $!\n";
}

# Run a samtools subcommand test.  As well as running samtools, various
# options are available to compare the output with a known version.  The
# default subcommand run is 'view', unless overridden by the 'cmd' argument.
# This will call passed() or failed() depending on the outcome of the test.
#
# $opts is the global settings hash
# %args is a list of key => value pairs giving options for run_view_test.
# The options for running samtools view are:
#  msg         => Message describing the test.
#  cmd         => samtools subcommand to run (default is view)
#  args        => array ref of arguments to add to the samtools view command
#                 Does not include the output file, which is in $args{out}
#  stdin       => Redirect STDIN to a pipe running "cat $args{stdin}"
#                 Otherwise the input file is included in $args{args}
#  out         => Name of the output file to make.
#  redirect    => If set, redirect STDOUT, otherwise use "-o $args{out}"
#  ref_path    => Setting for the REF_PATH environment variable
#  expect_fail => Expected failure, convert failed to xfail and passed to xpass
#
# One of the compare* options can be used to compare to an existing file:
#  compare     => Compare $args{out} with $args{compare} using sam_compare()
#  compare_sam => As compare, but run "samtools view -h $args{out}" first to
#                 convert an alternate format to SAM.  The resulting SAM file
#                 is saved to "$args{out}.sam" unless $args{pipe} is set.
#  pipe        => When used with compare_sam, the "samtools view -h" is read
#                 from a pipe instead of being stored in a file.  This can
#                 be used to avoid making very big SAM files for some tests.
#  compare_bam => Compare $args{out} with $args{compare_bam} using bam_compare()
#                 This is a binary comparison of the gunzipped BAM files.
#  compare_count => Compare count stored in $args{out} and $args{count_compare}
#                   The latter can either be a number, or a file containing a
#                   number.
#  compare_text => Compare $args{out} with $args{compare_text} using
#                  text_compare().  This is an exact line-by-line comparison.

sub run_view_test
{
    my ($opts, %args) = @_;

    printf "\t%-60s", $args{msg};

    local $ENV{REF_PATH} =  $args{ref_path} if ($args{ref_path});

    my @cmd = ("$$opts{bin}/samtools", $args{cmd} ? $args{cmd} : 'view');
    if ($args{out} && !$args{redirect}) { push(@cmd, '-o', $args{out}); }
    if ($args{args}) { push(@cmd, @{$args{args}}); }

    my $pid = fork();
    unless (defined($pid)) { die "Couldn't fork : $!\n"; }
    my $save_stdout;
    if ($^O eq 'MSWin32') {
	# Ensure we can restore STDOUT properly on Windows.  Not doing this
	# causes output to be lost if STDOUT is redirected below.
	open($save_stdout, '>&STDOUT') || die "Couldn't dup STDOUT: $!\n";
    }
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
    if ($^O eq 'MSWin32') {
	open(STDOUT, '>&', $save_stdout)
	    || die "Couldn't restore STDOUT : $!\n";
    }
    
    if (!$res && $args{compare_sam} && $args{out}) {
        # Convert output back to sam and compare
        my $sam_name = "$args{out}.sam";
        my @cmd2 = ("$$opts{bin}/samtools", 'view', '-h');
        if (!$args{pipe}) {
            push(@cmd2, '-o', $sam_name);
        }
        push(@cmd2, $args{out});

        push(@cmd, '&&', @cmd2); # For the 'Failed command' message below

        my $sam_out;
        if (!$args{pipe}) {
            $sam_out = $sam_name;
            $res = system(@cmd2) == 0 ? 0 : 1;
        } else {
            open($sam_out, '-|', @cmd2)
                || die "Couldn't open pipe from @cmd2: $!\n";
	    binmode($sam_out);
        }
        # Hack $args so the comparison gets done
        $args{compare} = $args{compare_sam};
        $args{out}     = $sam_out;
    }

    if (!$res) {
        if ($args{compare} && $args{out}) {
            $res = sam_compare($opts, $args{out}, $args{compare});
        } elsif ($args{compare_bam} && $args{out}) {
            $res = bam_compare($opts, $args{out}, $args{compare_bam});
        } elsif ($args{compare_count}) {
            $res = count_compare($args{out}, $args{compare_count});
        } elsif ($args{compare_text}) {
            $res = text_compare($args{out}, $args{compare_text});
        }
    }
    if ($res) {
        print "\tFailed command:\n\t@cmd\n\n";
        failed($opts, %args);
    } else {
        passed($opts, %args);
    }
}

# Runs a test of the samtools view -s subsampling option.
# The subsampling is pseudo-random, so multiple tests are run with different
# seeds, and the resulting read counts are tested against an acceptable
# range.  As the input file is paired, this also checks that both reads from
# each pair get returned.
#
# $opts is the global settings hash
# %args is a list of key => value pairs giving options
# The options are:
#  msg    => Message describing the test.
#  trials => The number of trials to run.  The trial number is used as the seed.
#  input  => The input file to read.
#  frac   => The fraction of reads to return (must be < 1).
#  min    => The minimum acceptable read count
#  max    => The maximum acceptable read count

sub run_view_subsample_test
{
    my ($opts, %args) = @_;

    printf "\t%s ", $args{msg};

    my @counts;
    my $res = 0;
    for (my $try = 0; $res == 0 && $try < $args{trials}; $try++) {
        my %reads;
        my @cmd = ("$$opts{bin}/samtools", 'view',
                   '-s', $try + $args{frac}, $args{input});
        open(my $samp, '-|', @cmd) || die "Couldn't open pipe from @cmd: $!\n";
        while (<$samp>) {
            my ($name) = split(/\t/);
            $reads{$name}++;
        }
        close($samp) || die "Error running @cmd\n";
        my $count = 0;
        while (my ($name, $num) = each %reads) {
            if ($num != 2) {
                print "\n\tGot one of read $name, expected two.\n";
                $res = 1;
                last;
            }
            $count += $num;
        }
        print ".";
        push(@counts, $count);
    }
    if (0 == @counts) {
        print "samtools view -s returned no counts\n";
        $res = 1;
    }

    @counts = sort { $a <=> $b } @counts;
    if ($counts[0] < $args{min} || $counts[-1] > $args{max}) {
        printf("\n\tOutput out of range: target (%d..%d) got (%d..%d)\n",
               $args{min}, $args{max}, $counts[0], $counts[-1]);
        $res = 1;
    }

    if (!$res) {
        passed($opts, %args);
        return;
    }
    failed($opts, %args);
}

# Open a pipe from bgzip to decompress a gzipped file.  bgzip is annoying
# as it checks the suffix of the file it is decompressing even when
# writing to STDOUT.  Hence we get it to read from a redirected STDIN
# instead so it doesn't care about the filename any more.
#
# $opts is the global settings hash
# $in is the compressed file to read
#
# Returns a file handle reference to the pipe from bgzip.

sub open_bgunzip
{
    my ($opts, $in) = @_;

    my $bgzip;
    open($bgzip, "$$opts{bgzip} -c -d < $in |")
	|| die "Couldn't open pipe to bgzip -c -d $!\n";
    binmode($bgzip);
    return $bgzip;
}

# Compare two SAM files.  Headers are collated by type to allow for reordering
# although all headers of a particular type should be in the same order.
# Optional tags on alignment lines can also be reordered.
#
# $opts is the global settings hash
# $sam1 is the first SAM file
# $sam2 is the second SAM file.
#
# Returns 0 if the files were the same, 1 if different.
#
# If $sam1 is a reference, it is assumed to be a pipe from which the data
# should be read.  The pipe will be closed after being read.  Otherwise it
# is treated as a filename.
#
# If $sam2 ends in '.gz' it is uncompressed with bgzip.

sub sam_compare
{
    my ($opts, $sam1, $sam2) = @_;

    unless (-e $sam2) {
        print "\n\tMissing SAM file $sam2.\n";
        return 1;
    }

    unless ((ref($sam1) || -e $sam1)) {
        print "\n\tMissing SAM file $sam1.\n";
        return 1;
    }

    my %hdr1;
    my %hdr2;

    my ($lno1, $lno2) = (0, 0);
    my ($l1, $l2, $ht1, $ht2);

    my $f1;
    if (ref($sam1)) {
        $f1 = $sam1;
    } else {
        open($f1, '<', $sam1) || die "Couldn't open $sam1: $!\n";
    }
    while ($l1 = <$f1>) {
        $lno1++;
        if (($ht1) = $l1 =~ /^(@\S+)/) {
	    $l1 =~ s/\x0d\x0a$/\x0a/;
            push(@{$hdr1{$ht1}}, $l1);
        } else {
            last;
        }
    }

    my $f2;
    if ($sam2 =~ /\.gz$/) {
        $f2 = open_bgunzip($opts, $sam2);
    } else {
        open($f2, '<', $sam2) || die "Couldn't open $sam2: $!\n";
    }
    while ($l2 = <$f2>) {
        $lno2++;
        if (($ht2) = $l2 =~ /^(@\S+)/) {
	    $l2 =~ s/\x0d\x0a$/\x0a/;
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
                next if ($h1->[$i] eq $h2->[$i]);

                # Hack to deal with CRAM adding M5 tags
                if ($ht eq '@SQ' && $h1->[$i] =~ /\tM5/ && $h2->[$i] !~ /\tM5/){
                    $h1->[$i] =~ s/\tM5:[0-9a-f]+//;
                    next if ($h1->[$i] eq $h2->[$i]);
                }

                $same = 0;
                last;
            }
        }
        if (!$same) {
            print "\n\tHeader type $ht differs.\n";
            print "\t$sam1 has:\n";
            foreach my $t (@{$hdr1{$ht}}) {
                print "\t$t\n";
            }
            print "\t$sam2 has:\n";
            foreach my $t (@{$hdr2{$ht}}) {
                print "\t$t\n";
            }
            close($f1);
            close($f2);
            return 1;
        }
    }

    while ($l1 && $l2) {
        chomp($l1); $l1 =~ s/\r$//;
        chomp($l2); $l2 =~ s/\r$//;
        my @sam1 = split(/\t/, $l1);
        my @sam2 = split(/\t/, $l2);
        my @tags1 = sort(splice(@sam1, 11));
        my @tags2 = sort(splice(@sam2, 11));
	# Windows uses e.g. 1.9e+009 vs 1.9e+09 on Linux
	@tags1 = map {s/(e[-+])0(\d\d)/$1$2/g;$_} @tags1;
	@tags2 = map {s/(e[-+])0(\d\d)/$1$2/g;$_} @tags2;
        if (join("\t", @sam1, @tags1) ne join("\t", @sam2, @tags2)) {
            last;
        }
        $l1 = <$f1>;
        $l2 = <$f2>;
        $lno1++;
        $lno2++;
    }

    close($f1) || die "Error reading $sam1: $!\n";
    close($f2) || die "Error reading $sam2: $!\n";

    if ($l1 || $l2) {
        print "\n"; STDOUT->flush();
        print STDERR "\tSAM files differ at $sam1 line $lno1 / $sam2 line $lno2\n";
        print STDERR "$l1\n$l2\n";
        return 1;
    }

    return 0;
}

# Compare two BAM files.  This is done by uncompressing them and performing
# a binary comparison.  Used to check that (for example) the uncompressed and
# compressed versions are the same.
#
# $opts is the global settings hash
# $bam1 is the first BAM file
# $bam2 is the second BAM file.
#
# Returns 0 if the files are the same, 1 if they are different.

sub bam_compare
{
    my ($opts, $bam1, $bam2) = @_;

    my $buffer1;
    my $buffer2;
    my $bytes1;
    my $bytes2;
    my $fail = 0;

    my $b1 = open_bgunzip($opts, $bam1);
    my $b2 = open_bgunzip($opts, $bam2);
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
    close($b1) || die "Error running bgzip -c -d $bam1\n";
    close($b2) || die "Error running bgzip -c -d $bam2\n";
    if ($fail) {
        print "\n\tBAM files $bam1 and $bam2 differ.\n";
        return 1;
    }

    return 0;
}

# Compare two text files.
#
# $txt1 is the first file
# $txt2 is the second file.
#
# Returns 0 if the files are the same, 1 if they are different.

sub text_compare
{
    my ($txt1, $txt2) = @_;

    open(my $t1, '<', $txt1) || die "Couldn't open $txt1 : $!\n";
    open(my $t2, '<', $txt2) || die "Couldn't open $txt2 : $!\n";
    my $line = 0;
    my $diff = 0;
    while (!$diff) {
        $line++;
        my $l1 = <$t1>;
        my $l2 = <$t2>;
        last if (!defined($l1) && !defined($l2));
	$l1 =~ s/\x0d\x0a$/\x0a/;
	$l2 =~ s/\x0d\x0a$/\x0a/;
        if (($l1 || '') ne ($l2 || '')) {
            $diff = 1;
            if (defined($l1)) { chomp($l1); }
            if (defined($l2)) { chomp($l2); }
            print "\n\tFiles differ at line $line:\n";
            print "\t$txt1 ", defined($l1) ? " : $l1\n" : "End of file\n";
            print "\t$txt2 ", defined($l2) ? " : $l2\n" : "End of file\n";
        }
    }
    close($t1) || die "Error reading $t1 : $!\n";
    close($t2) || die "Error reading $t2 : $!\n";

    return $diff;
}

# Compare two counts to see if they are the same.  The first is in a file
# created by, for example, 'samtools view -c'.  The second can either be
# a number, or the name of a file to read.
#
# $count1 is the first file containing a count
# $count2 is either a number or the second file containing a count
#
# Returns 0 if the counts were the same, 1 if they differ.

sub count_compare
{
    my ($count1, $count2) = @_;

    open(my $c1, '<', $count1) || die "Couldn't open $count1 : $!\n";
    my $number1 = <$c1>;
    chomp($number1); $number1 =~ s/\r$//;
    close($c1) || die "Error reading $count1 : $!\n";

    unless ($number1 =~ /^\d+$/) {
        print "\n\tExpected a number in $count1 but got '$number1'\n";
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
        return 1;
    }

    return 0;
}

# Generate a pair of reads, for use in gen_file.
#
# $reads is a hash that the reads are put into.
# $seq is the sequence (passed by reference for speed)
# $qual is the quality (passed by reference for speed)
# $pos1 is the position of the leftmost read.
# $size is the length of $$seq
# $rnum is the read number, used to generate the read name.

sub gen_pair
{
    my ($reads, $seq, $qual, $pos1, $size, $rnum) = @_;

    my $l1 = int(rand(50) + 75);
    my $l2 = int(rand(50) + 75);
    my $pos2 = $pos1 + int(rand(50) + 275);
    if ($pos2 + $l2 > $size) { return; }
    my $orient = int(rand(2));
    my $name = "ERR123456.${rnum}";

    my $rd1 = sprintf("%s\t%d\tref1\t%d\t40\t%dM\t=\t%d\t%d\t%s\t%s\tNM:i:0\tMD:Z:$l1\tRG:Z:g1",
                      $name, $orient ? 99 : 163, $pos1 + 1, $l1, $pos2 + 1,
                      $pos2 + $l2 - $pos1, substr($$seq, $pos1, $l1),
                      substr($$qual, $pos1, $l1));
    my $rd2 = sprintf("%s\t%d\tref1\t%d\t40\t%dM\t=\t%d\t%d\t%s\t%s\tNM:i:0\tMD:Z:$l2\tRG:Z:g1",
                      $name, $orient ? 147 : 83, $pos2 + 1, $l2, $pos1 + 1,
                      -($pos2 + $l2 - $pos1), substr($$seq, $pos2, $l2),
                      substr($$qual, $pos2, $l2));
    push(@{$reads->{$pos1}}, $rd1);
    push(@{$reads->{$pos2}}, $rd2);
}

# Randomly generate a test SAM file.  Used to create the big SAM file for
# testing the -s and -@ options of samtools view.  As well as the SAM file
# it writes a fasta file with the reference sequence in it and a .fai file
# which is an index of the fasta file.  Both of these are used for CRAM.
# The SAM file is compressed with bgzf to keep the size down.
#
# To generate random sequence aligned against identical referenes,
# pass in the 4th argument (ref_seed) with a constant value.  Otherwise a
# random reference is generated too.
#
# $opts is the global setting hash.
# $prefix is the prefix for the generated file names.
# $size is the length of the random reference sequence.
#
# Returns a list containing the name of the SAM file and the read count.


sub gen_file
{
    my ($opts, $prefix, $size, $ref_seed) = @_;

    local $| = 1;

    $ref_seed = rand(1<<31) unless defined($ref_seed);
    my $pair_seed = rand(1<<31);

    print "\tGenerating test data file ";
    my $dot_interval = $size / 10;

    my $seq  = "!" x $size;
    my $qual = $seq;

    srand($ref_seed);
    my $next_dot = $dot_interval;
    for (my $b = 0; $b < $size; $b++) {
        if ($b == $next_dot) {
            print ".";
            $next_dot += $dot_interval;
        }
        substr($seq, $b, 1) = (qw(A C G T))[int(rand(4))];
        substr($qual, $b, 1) = chr(33 + int(rand(40)));
    }
    my $fasta = "$prefix.fa";
    open(my $fa, '>', $fasta) || die "Couldn't open $fasta for writing : $!\n";
    binmode($fa);
    print $fa ">ref1\n";
    for (my $i = 0; $i < $size; $i += 60) {
        print $fa substr($seq, $i, 60), "\n";
    }
    close($fa) || die "Error writing to $fasta : $!\n";

    my $fai = "$prefix.fa.fai";
    open($fa, '>', $fai) || die "Couldn't open $fai for writing : $!\n";
    binmode($fa);
    print $fa "ref1\t$size\t6\t60\t61\n";
    close($fa) || die "Error writing to $fai : $!\n";

    my $sam = "$prefix.sam.gz";
    open(my $s, "| $$opts{bgzip} -c > $sam")
	|| die "Couldn't open pipe to bgzip $!";
    binmode($s);
    print $s "\@HD\tVN:1.4\tSO:coordinate\n";
    print $s "\@RG\tID:g1\tDS:Group 1\tLB:Lib1\tSM:Sample1\n";
    print $s "\@SQ\tSN:ref1\tLN:$size\tUR:$fasta\n";
    my %read_store;
    my $rnum = 0;
    $next_dot = $dot_interval;

    srand($pair_seed);
    for (my $i = 0; $i < $size; $i++) {
        if ($i == $next_dot) {
            print ".";
            $next_dot += $dot_interval;
        }
        if ($i % 20 == 0) {
            gen_pair(\%read_store, \$seq, \$qual, $i, $size, ++$rnum);
        }
        if (exists($read_store{$i})) {
            foreach my $read (@{$read_store{$i}}) {
                print $s $read, "\n";
            }
            delete($read_store{$i});
        }
    }
    close($s) || die "Error running bgzip -c > $sam : $!\n";
    print " done\n";
    return ($sam, $rnum * 2);
}

# Run samtools view tests.

sub test_view
{
    my ($opts) = @_;

    my $test_name = "test_view";
    print "$test_name:\n";

    my $ref_path = "$$opts{path}/dat/cram_md5/%s";

    # Add @SQ UR: tags to the stored test file.

    my $sam_no_ur   = "$$opts{path}/dat/view.001.sam";
    my $sam_with_ur = "$$opts{tmp}/view.001.sam";
    add_ur_tags($sam_no_ur, $sam_with_ur, "$$opts{path}/dat/view.001.fa");

    # Generate a big SAM file for use later.

    my ($big_sam, $big_sam_count)
        = gen_file($opts, "$$opts{tmp}/view.big", 1000000);

    my $test = 1;

    my $out = "$$opts{tmp}/view";

    # SAM -> BAM -> SAM
    my $bam_with_ur_out = sprintf("%s.test%03d.bam", $out, $test);
    run_view_test($opts,
                  msg => "$test: SAM -> BAM -> SAM",
                  args => ['-b', $sam_with_ur],
                  out => $bam_with_ur_out,
                  compare_sam => $sam_with_ur);
    $test++;

    # SAM -> uncompressed BAM -> SAM
    my $uncomp_bam = sprintf("%s.test%03d.bam", $out, $test);
    run_view_test($opts,
                  msg => "$test: SAM -> uncompressed BAM",
                  args => ['-u', $sam_with_ur],
                  out => $uncomp_bam,
                  compare_bam => $bam_with_ur_out);
    $test++;
    run_view_test($opts,
                  msg => "$test: uncompressed BAM -> SAM and compare",
                  args => ['-h', $uncomp_bam],
                  out => sprintf("%s.test%03d.sam", $out, $test),
                  compare => $sam_with_ur);
    $test++;

    # SAM -> fast compressed BAM -> SAM
    my $fastcomp_bam = sprintf("%s.test%03d.bam", $out, $test);
    run_view_test($opts,
                  msg => "$test: SAM -> fast compressed BAM",
                  args => ['-1', $sam_with_ur],
                  out => $fastcomp_bam,
                  compare_bam => $bam_with_ur_out);
    $test++;
    run_view_test($opts,
                  msg => "$test: fast compressed BAM -> SAM and compare",
                  args => ['-h', $fastcomp_bam],
                  out => sprintf("%s.test%03d.sam", $out, $test),
                  compare => $sam_with_ur);
    $test++;

    # SAM -> CRAM -> SAM with UR tags
    my $cram_with_ur_out = sprintf("%s.test%03d.cram", $out, $test);
    run_view_test($opts,
                  msg => "$test: SAM -> CRAM -> SAM (UR header tags)",
                  args => ['-C', $sam_with_ur],
                  out => $cram_with_ur_out,
                  compare_sam => $sam_with_ur);
    $test++;

    # SAM -> CRAM -> SAM with M5 tags
    run_view_test($opts,
                  msg => "$test: SAM -> CRAM -> SAM (M5 header tags)",
                  args => ['-C', $sam_no_ur],
                  out => sprintf("%s.test%03d.cram", $out, $test),
                  compare_sam => $sam_no_ur,
                  ref_path => $ref_path);
    $test++;

    # SAM -> BAM -> CRAM -> SAM with UR tags
    my $cram_from_bam = sprintf("%s.test%03d.cram", $out, $test);
    run_view_test($opts,
                  msg => "$test: BAM -> CRAM with UR -> SAM",
                  args => ['-C', $bam_with_ur_out],
                  out => $cram_from_bam,
                  compare_sam => $sam_with_ur);
    $test++;

    # SAM -> BAM -> CRAM -> SAM with M5 tags
    my $bam_no_ur = sprintf("%s.test%03d.bam", $out, $test);
    run_view_test($opts,
                  msg => "$test: SAM -> BAM (M5 header tags)",
                  args => ['-b', $sam_no_ur],
                  out => $bam_no_ur,
                  compare_sam => $sam_no_ur);
    $test++;
    my $cram_no_ur = sprintf("%s.test%03d.cram", $out, $test);
    run_view_test($opts,
                  msg => "$test: SAM -> BAM -> CRAM -> SAM (M5 header tags)",
                  args => ['-C', $bam_no_ur],
                  out => $cram_no_ur,
                  compare_sam => $sam_no_ur,
                  ref_path => $ref_path);
    $test++;

    # SAM -> BAM -> CRAM -> BAM -> SAM with UR tags
    my $bam_from_cram = sprintf("%s.test%03d.bam", $out, $test);
    run_view_test($opts,
                  msg => "$test: CRAM -> BAM with UR",
                  args => ['-b', $cram_from_bam],
                  out => $bam_from_cram,
                  compare_sam => $sam_with_ur);
    $test++;

    # SAM -> BAM -> CRAM -> BAM -> SAM with M5 tags
    run_view_test($opts,
                  msg => "$test: CRAM -> BAM with M5",
                  args => ['-b', $cram_no_ur],
                  out => sprintf("%s.test%03d.bam", $out, $test),
                  compare_sam => $sam_no_ur,
                  ref_path => $ref_path);
    $test++;

    # Write to stdout
    run_view_test($opts,
                  msg => "$test: SAM -> SAM via stdout",
                  args => ['-h', $sam_with_ur],
                  out => sprintf("%s.test%03d.sam", $out, $test),
                  redirect => 1,
                  compare => $sam_with_ur);
    $test++;
    run_view_test($opts,
                  msg => "$test: SAM -> BAM via stdout",
                  args => ['-b', $sam_with_ur],
                  out => sprintf("%s.test%03d.bam", $out, $test),
                  redirect => 1,
                  compare_bam => $bam_with_ur_out);
    $test++;
    my $cram_via_stdout = sprintf("%s.test%03d.cram", $out, $test);
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
                  out => sprintf("%s.test%03d.sam", $out, $test),
                  stdin => $sam_with_ur,
                  compare => $sam_with_ur);
    $test++;
    run_view_test($opts,
                  msg => "$test: BAM from stdin -> SAM",
                  args => ['-h', '-'],
                  out => sprintf("%s.test%03d.sam", $out, $test),
                  stdin => $bam_with_ur_out,
                  compare => $sam_with_ur);
    $test++;
    run_view_test($opts,
                  msg => "$test: CRAM from stdin -> SAM",
                  args => ['-h', '-'],
                  out => sprintf("%s.test%03d.sam", $out, $test),
                  stdin => $cram_with_ur_out,
                  compare => $sam_with_ur);
    $test++;


    # Header only options
    my $sam_header = "$$opts{tmp}/view.001.header.sam";
    filter_sam($sam_with_ur, $sam_header, {no_body => 1});

    run_view_test($opts,
                  msg => "$test: samtools view -H (SAM input)",
                  args => ['-H', $sam_with_ur],
                  out => sprintf("%s.test%03d.header", $out, $test),
                  compare => $sam_header);
    $test++;
    run_view_test($opts,
                  msg => "$test: samtools view -H (BAM input)",
                  args => ['-H', $bam_with_ur_out],
                  out => sprintf("%s.test%03d.header", $out, $test),
                  compare => $sam_header);
    $test++;
    run_view_test($opts,
                  msg => "$test: samtools view -H (CRAM input)",
                  args => ['-H', $cram_with_ur_out],
                  out => sprintf("%s.test%03d.header", $out, $test),
                  compare => $sam_header);
    $test++;

    # Body only
    my $sam_body = "$$opts{tmp}/view.001.body.sam";
    filter_sam($sam_with_ur, $sam_body, {no_header => 1});

    run_view_test($opts,
                  msg => "$test: No headers (SAM input)",
                  args => [$sam_with_ur],
                  out => sprintf("%s.test%03d.body", $out, $test),
                  compare => $sam_body);
    $test++;
    run_view_test($opts,
                  msg => "$test: No headers (BAM input)",
                  args => [$bam_with_ur_out],
                  out => sprintf("%s.test%03d.body", $out, $test),
                  compare => $sam_body);
    $test++;
    run_view_test($opts,
                  msg => "$test: No headers (CRAM input)",
                  args => [$cram_with_ur_out],
                  out => sprintf("%s.test%03d.body", $out, $test),
                  compare => $sam_body);
    $test++;

    # Filter and counting tests.

    # Group names file for -R test
    my $fogn = "$$opts{tmp}/view.001.fogn";
    open(my $f, '>', $fogn) || die "Couldn't open $fogn : $!\n";
    print $f "grp1\ngrp3\n" || die "Error writing to $fogn : $!\n";
    close($f) || die "Error writing to $fogn : $!\n";


    my @filter_tests = (
        # [test_name, {filter_sam options}, [samtools options], expect_fail]
        # Flags
        ['req128', {flags_required => 128}, ['-f', 128], 0],
        ['rej128', {flags_rejected => 128}, ['-F', 128], 0],
        ['rej128req2', { flags_rejected => 128, flags_required => 2 },
         ['-F', 128, '-f', 2], 0],
        # Read groups
        ['rg_grp2', { read_groups => { grp2 => 1 }}, ['-r', 'grp2'], 0],
        ['rg_fogn', { read_groups => { grp1 => 1, grp3 => 1 }},
         ['-R', $fogn], 0],
        ['rg_both', { read_groups => { grp1 => 1, grp2 => 1, grp3 => 1 }},
         ['-R', $fogn, '-r', 'grp2'], 0],
        ['rg_both2', { read_groups => { grp1 => 1, grp2 => 1, grp3 => 1 }},
         ['-r', 'grp2', '-R', $fogn], 0],
        # Libraries
        ['lib2', { libraries => { 'Library 2' => 1 }}, ['-l', 'Library 2'], 0],
        ['lib3', { libraries => { 'Library 3' => 1 }}, ['-l', 'Library 3'], 0],
        # Mapping qualities
        ['mq50',  { min_map_qual => 50 },  ['-q', 50], 0],
        ['mq99',  { min_map_qual => 99 },  ['-q', 99], 0],
        ['mq100', { min_map_qual => 100 }, ['-q', 100], 0],
        # Tag stripping
        ['tags1', { strip_tags => { fa => 1 } }, ['-x', 'fa'], 0],
        ['tags2', { strip_tags => { fa => 1, ha => 1 } },
         ['-x', 'fa', '-x', 'ha'], 0],
        # Tag strip plus read group
        ['tags_rg1', { strip_tags => { fa => 1 }, read_groups => { grp2 => 1 }},
         ['-x', 'fa', '-r', 'grp2'], 0],
        ['tags_rg2', { strip_tags => { RG => 1 }, read_groups => { grp2 => 1 }},
         ['-x', 'RG', '-r', 'grp2'], 0],
        # Minimum query length
        ['qlen10', { min_qlen => 10 }, ['-m', 10], 0],
        ['qlen11', { min_qlen => 11 }, ['-m', 11], 0],
        ['qlen15', { min_qlen => 15 }, ['-m', 15], 0],
        ['qlen16', { min_qlen => 16 }, ['-m', 16], 0],
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
                          out => sprintf("%s.test%03d.sam", $out, $test),
                          compare => $sam_file,
                          expect_fail => $$filter[3]);
            $test++;

            # Count test
            run_view_test($opts,
                          msg => "$test: Count @{$$filter[2]} ($$ip[0] input)",
                          args => ['-c', @{$$filter[2]}, $$ip[1]],
                          out => sprintf("%s.test%03d.count", $out, $test),
                          redirect => 1,
                          compare_count => $sam_file,
                          expect_fail => $$filter[3]);
            $test++;
        }
    }

    # Region query tests
    my $sam_no_ur2 = "$$opts{path}/dat/view.002.sam";
    my $sam_with_ur2 = "$$opts{tmp}/view.002.sam";
    add_ur_tags($sam_no_ur2, $sam_with_ur2, "$$opts{path}/dat/view.002.fa");

    my $bam_with_ur_out2 = sprintf("%s.test%03d.bam", $out, $test);
    run_view_test($opts,
                  msg => "$test: 1bp reads file SAM -> BAM -> SAM",
                  args => ['-b', $sam_with_ur2],
                  out => $bam_with_ur_out2,
                  compare_sam => $sam_with_ur2);
    $test++;
    my $cram_with_ur_out2 = sprintf("%s.test%03d.cram", $out, $test);
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
        # Unmapped reads
        ['reg_unmapped1', 1, { region => [['*']] }, [], ['*']],

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
                          out => sprintf("%s.test%03d.sam", $out, $test),
                          compare => $sam_file);
            $test++;
            # Count test
            run_view_test($opts,
                          msg => "$test: Count @{$$rt[3]} @{$$rt[4]} ($$ip[0] input)",
                          args => ['-c', @{$$rt[3]}, $input_file, @{$$rt[4]}],
                          out => sprintf("%s.test%03d.count", $out, $test),
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

    # We can't make a BAM without @SQ lines.  Instead make one
    # with no M5/UR in @SQ and see if we can still use it to make CRAM.
    my $bam_no_m5 = "$$opts{tmp}/view.001.no_sq.bam";
    run_view_test($opts,
                  msg => "$test: Make BAM with no M5/UR tags",
                  args => ['-b', $sam_no_m5],
                  out => $bam_no_m5,
                  compare_sam => $sam_no_m5);
    $test++;

    my $ref_file = "$$opts{path}/dat/view.001.fa";
    my $ref_idx  = "$$opts{path}/dat/view.001.fa.fai";

    # Test SAM output
    foreach my $in ([SAM => $sam_no_sq], [BAM => $bam_no_m5]) {
        foreach my $topt (['-t', $ref_idx], ['-T', $ref_file]) {
            run_view_test($opts,
                          msg => "$test: Add \@SQ with $topt->[0] ($in->[0] -> SAM)",
                          args => ['-h', @$topt, $in->[1]],
                          out => sprintf("%s.test%03d.sam", $out, $test),
                          compare => $sam_no_m5);
            $test++;
        }
    }

    # Test BAM output.
    foreach my $in ([SAM => $sam_no_sq], [BAM => $bam_no_m5]) {
        foreach my $topt (['-t', $ref_idx], ['-T', $ref_file]) {
            my $bam = sprintf("%s.test%03d.bam", $out, $test);
            run_view_test($opts,
                          msg => "$test: Add \@SQ with $topt->[0] ($in->[0] -> BAM)",
                          args => ['-b', @$topt, $in->[1]],
                          out => $bam,
                          compare_sam => $sam_no_m5);
            $test++;
        }
    }

    # Test CRAM output
    foreach my $in ([SAM => $sam_no_sq], [BAM => $bam_no_m5]) {
        foreach my $topt (['-t', $ref_idx], ['-T', $ref_file]) {
            my $cram = sprintf("%s.test%03d.cram", $out, $test);
            run_view_test($opts,
                          msg => "$test: Add \@SQ with $topt->[0] ($in->[0] -> CRAM)",
                          args => ['-C', @$topt, $in->[1]],
                          out => $cram,
                          compare_sam => $sam_with_ur);
            $test++;
        }
    }

    # Test CRAM with explicit -T
    my $cram_no_ur_t = sprintf("%s.test%03d.cram", $out, $test);
    run_view_test($opts,
                  msg => "$test: Make CRAM with no UR field",
                  args => ['-C', $sam_no_ur],
                  ref_path => "$$opts{path}/dat/cram_md5",
                  out => $cram_no_ur_t);
    $test++;

    run_view_test($opts,
                  msg => "$test: Decoding CRAM with no UR field via -T",
                  args => ['-T', $ref_file, $cram_no_ur_t],
                  out => sprintf("%s.test%03d.sam", $out, $test),
                  compare => $sam_no_ur);
    $test++;


    # CIGAR B-operator removal tests.
    my $b_op_sam      = "$$opts{path}/dat/view.003.sam";
    my $b_op_expected = "$$opts{path}/dat/view.003.expected.sam";
    run_view_test($opts,
                  msg => "$test: CIGAR B-operator removal",
                  args => ['-h', '-B', $b_op_sam],
                  out => sprintf("%s.test%03d.sam", $out, $test),
                  compare => $b_op_expected);
    $test++;

    # Threads
    # big SAM -> BAM
    my $big_bam = sprintf("%s.test%03d.bam", $out, $test);
    run_view_test($opts,
                  msg => "$test: Big SAM -> BAM",
                  args => ['-b', $big_sam],
                  out => $big_bam);
    $test++;

    foreach my $threads (2, 4) {
        run_view_test($opts,
                      msg => "$test: Big SAM -> BAM ($threads threads)",
                      args => ['-b', '-@', $threads, $big_sam],
                      out => sprintf("%s.test%03d.bam", $out, $test),
                      compare_bam => $big_bam);
        $test++;
    }

    # big SAM -> CRAM
    my $big_cram = sprintf("%s.test%03d.cram", $out, $test);
    run_view_test($opts,
                  msg => "$test: Big SAM -> CRAM",
                  args => ['-C', $big_sam],
                  out => $big_cram,
                  compare_sam => $big_sam,
                  pipe => 1);
    $test++;

    foreach my $threads (2, 4) {
        run_view_test($opts,
                      msg => "$test: Big SAM -> CRAM ($threads threads)",
                      args => ['-C', '-@', $threads, $big_sam],
                      out => sprintf("%s.test%03d.cram", $out, $test),
                      compare_sam => $big_sam,
                      pipe => 1);
        $test++;
    }


    # Subsampling (-s) option.  This is random so accept if within +/- 5%
    # of the expected number of reads.

    foreach my $ip ([SAM => $big_sam], [BAM => $big_bam], [CRAM => $big_cram]) {
        foreach my $frac (0.2, 0.5, 0.8) {

            run_view_subsample_test($opts,
                                    msg => "$test: Subsample $frac ($ip->[0] input)",
                                    frac => $frac,
                                    input => $ip->[1],
                                    trials => 10,
                                    min => $big_sam_count * $frac * 0.95,
                                    max => $big_sam_count * $frac * 1.05);
            $test++;
        }
    }
}

# cat SAM files in the same way as samtools cat does with BAMs
#
# $sam_out is the name of the output file
# $sam_in is the name of the input file.  More than one of these can be
#   passed in.  The header is taken from the first file.

sub cat_sams
{
    my $sam_out = shift(@_);
    my $first = 1;

    open(my $out, '>', $sam_out)
        || die "Couldn't open $sam_out for writing: $!\n";

    while (my $sam_in = shift(@_)) {
        my $in;
        if ($sam_in =~ /\.gz$/) {
            $in = open_bgunzip($opts, $sam_in);
        } else {
            open($in, '<', $sam_in) || die "Couldn't open $sam_in : $!\n";
        }
        while (<$in>) {
            next if (/^@/ && !$first);
            print $out $_ || die "Error writing to $sam_out : $!\n";
        }
        close($in) || die "Error reading $sam_in : $!\n";
        $first = 0;
    }
    close($out) || die "Error writing to $sam_out : $!\n";
}

# Test samtools cat.

sub test_cat
{
    my ($opts) = @_;

    my $test_name = "test_cat";
    print "$test_name:\n";

    my @sams;
    my @bams;
    my @crams;
    my @bgbams;
    my $nfiles = 4;
    my $test = 1;
    my $out = "$$opts{tmp}/cat";

    # Generate some files big enough to include a few bgzf blocks
    for (my $i = 0; $i < $nfiles; $i++) {
        ($sams[$i]) = gen_file($opts, sprintf("%s.%d", $out, $i + 1), 10000, 15551);

        # Convert to BAM
        $bams[$i] = sprintf("%s.%d.bam", $out, $i + 1);
        run_view_test($opts,
                      msg =>  sprintf("Generate BAM file #%d", $i + 1),
                      args => ['-b', $sams[$i]],
                      out => $bams[$i],
                      compare_sam => $sams[$i],
                      pipe => 1);

        # Recompress with bgzip to alter the location of the bgzf boundaries.
        $bgbams[$i] = sprintf("%s.%d.bgzip.bam", $out, $i + 1);
        cmd("'$$opts{bgzip}' -c -d < '$bams[$i]' | '$$opts{bgzip}' -c > '$bgbams[$i]'");

        # Create CRAMs
        $crams[$i] = sprintf("%s.%d.cram", $out, $i + 1);
        run_view_test($opts,
                      msg =>  sprintf("Generate CRAM file #%d", $i + 1),
                      args => ['-C', $sams[$i]],
                      out => $crams[$i],
                      compare_sam => $sams[$i],
                      pipe => 1);
    }

    # Make a concatenated SAM file to compare
    my $catsam1 = "$out.all1.sam";
    cat_sams($catsam1, @sams);

    foreach my $redirect (0, 1) {
        my $to_stdout = $redirect ? ' to stdout' : '';

        # Basic test
        run_view_test($opts,
                      msg =>  "$test: cat BAM files$to_stdout",
                      cmd => 'cat',
                      args => [@bams],
                      out => sprintf("%s.test%03d.bam", $out, $test),
                      redirect => $redirect,
                      compare_sam => $catsam1);
        $test++;

        # Test BAM files recompressed with bgzip
        run_view_test($opts,
                      msg =>  "$test: cat recompressed BAM files$to_stdout",
                      cmd => 'cat',
                      args => [@bgbams],
                      out => sprintf("%s.test%03d.bam", $out, $test),
                      redirect => $redirect,
                      compare_sam => $catsam1);
        $test++;

        # Test CRAM files
        run_view_test($opts,
                      msg =>  "$test: cat CRAM files$to_stdout",
                      cmd => 'cat',
                      args => [@crams],
                      out => sprintf("%s.test%03d.cram", $out, $test),
                      redirect => $redirect,
                      compare_sam => $catsam1);
        $test++;
    }

    # Test reheader option
    my $hdr_no_ur   = "$$opts{path}/dat/cat.hdr";
    my $header      = "$$opts{tmp}/cat.hdr";
    add_ur_tags($hdr_no_ur, $header, "$$opts{tmp}/cat.1.fa");
    my $catsam2 = "$out.all2.sam";
    cat_sams($catsam2, $header, @sams);

    run_view_test($opts,
                  msg =>  "$test: cat BAM files with new header",
                  cmd => 'cat',
                  args => ['-h', $header, @bams],
                  out => sprintf("%s.test%03d.bam", $out, $test),
                  compare_sam => $catsam2);
    $test++;

    run_view_test($opts,
                  msg =>  "$test: cat recompressed BAM files with new header",
                  cmd => 'cat',
                  args => ['-h', $header, @bgbams],
                  out => sprintf("%s.test%03d.bam", $out, $test),
                  compare_sam => $catsam2);
    $test++;

    run_view_test($opts,
                  msg =>  "$test: cat CRAM files with new header",
                  cmd => 'cat',
                  args => ['-h', $header, @crams],
                  out => sprintf("%s.test%03d.cram", $out, $test),
                  compare_sam => $catsam2);
    $test++;
}

sub sam2fq
{
    my ($sam_in, $fq_out, $suffixes) = @_;

    open(my $in, '<', $sam_in) || die "Couldn't open $sam_in : $!\n";
    open(my $out, '>', $fq_out)
        || die "Couldn't open $fq_out for writing : $!\n";
    while (<$in>) {
        next if (/^@/);
        my @s = split(/\t/, $_);
        next if ($s[1] & (256|2048));
        my $dirn = ($s[1] & 0xc0) >> 6;
        my $suff = $suffixes ? ('', '/1', '/2', '')[$dirn] : '';
        if (($s[1] & 0x10) != 0) { # reverse complement
            $s[9] =~ tr/ACGTMRWSYKVHDBN/TGCAKYWSRMBDHVN/;
            $s[9] = reverse($s[9]);
            $s[10] = reverse($s[10]);
        }
        print $out "\@$s[0]$suff\n$s[9]\n+\n$s[10]\n";
    }
    close($out) || die "Error writing $fq_out : $!\n";
    close($in) || die "Error reading $sam_in : $!\n";
}

sub test_bam2fq
{
    my ($opts, %args) = @_;

    my $threads = exists($args{threads}) ? ['-@', $args{threads}] : [];
    my $nthreads = exists($args{threads}) ? " ($args{threads} threads)" : "";

    my $test_name = "test_bam2fq";
    print "$test_name:\n";

    my $out = "$$opts{tmp}/bam2fq";

    my $sam = "$$opts{path}/dat/bam2fq.001.sam";
    my $bam = "$out.001.bam";
    my $cram = "$out.001.cram";
    my $fqsuffix   = "$out.001.fq";
    my $fqnosuffix = "$out.001.nosuff.fq";

    sam2fq($sam, $fqsuffix,   1);
    sam2fq($sam, $fqnosuffix, 0);

    # Make a BAM file with the test data
    run_view_test($opts,
                  msg =>  "Generate BAM file$nthreads",
                  args => ['-b', $sam, @$threads],
                  out => $bam,
                  compare_sam => $sam,
                  pipe => 1);

    # Make a CRAM file with the test data
    run_view_test($opts,
                  msg =>  "Generate CRAM file$nthreads",
                  args => ['-C', $sam, @$threads],
                  out => $cram,
                  ref_path => "$$opts{path}/dat/cram_md5",
                  compare_sam => $sam,
                  pipe => 1);

    my $test = 1;
    my @inputs = ([SAM => $sam], [BAM => $bam], [CRAM => $cram]);
    foreach my $input (@inputs) {
        foreach my $nosuffix (0, 1) {
            my @n = $nosuffix ? ('-n') : ();

            run_view_test($opts,
                          msg => "$test: fastq @n ($input->[0] input)$nthreads",
                          cmd => 'fastq',
                          args => [@n, @$threads, $input->[1]],
                          out => sprintf("%s.test%03d.fq", $out, $test),
                          ref_path => "$$opts{path}/dat/cram_md5",
                          redirect => 1,
                          compare_text => $nosuffix ? $fqnosuffix : $fqsuffix);
            $test++;
        }
    }
    # basic 2 output test without singleton tracking
    test_cmd($opts, out=>'bam2fq/1.stdout.expected', out_map=>{'1.fq' => 'bam2fq/1.1.fq.expected', '2.fq' => 'bam2fq/1.2.fq.expected'},cmd=>"$$opts{bin}/samtools fastq @$threads -1 $$opts{path}/1.fq -2 $$opts{path}/2.fq $$opts{path}/dat/bam2fq.001.sam");
    # basic 2 output test with singleton tracking but no singleton
    test_cmd($opts, out=>'bam2fq/2.stdout.expected', out_map=>{'1.fq' => 'bam2fq/2.1.fq.expected', '2.fq' => 'bam2fq/2.2.fq.expected', 's.fq' => 'bam2fq/2.s.fq.expected'}, cmd=>"$$opts{bin}/samtools fastq @$threads -s $$opts{path}/s.fq -1 $$opts{path}/1.fq -2 $$opts{path}/2.fq $$opts{path}/dat/bam2fq.001.sam");
    # basic 2 output test with singleton tracking with a singleton in the middle
    test_cmd($opts, out=>'bam2fq/2.stdout.expected', out_map=>{'1.fq' => 'bam2fq/3.1.fq.expected', '2.fq' => 'bam2fq/3.2.fq.expected', 's.fq' => 'bam2fq/3.s.fq.expected'}, cmd=>"$$opts{bin}/samtools fastq @$threads -s $$opts{path}/s.fq -1 $$opts{path}/1.fq -2 $$opts{path}/2.fq $$opts{path}/dat/bam2fq.002.sam");
    # basic 2 output test with singleton tracking with a singleton as last read
    test_cmd($opts, out=>'bam2fq/2.stdout.expected', out_map=>{'1.fq' => 'bam2fq/4.1.fq.expected', '2.fq' => 'bam2fq/4.2.fq.expected', 's.fq' => 'bam2fq/4.s.fq.expected'}, cmd=>"$$opts{bin}/samtools fastq @$threads -s $$opts{path}/s.fq -1 $$opts{path}/1.fq -2 $$opts{path}/2.fq $$opts{path}/dat/bam2fq.003.sam");
    # tag output test with singleton tracking with a singleton as last read
    test_cmd($opts, out=>'bam2fq/2.stdout.expected', out_map=>{'1.fq' => 'bam2fq/4.1.fq.expected', '2.fq' => 'bam2fq/4.2.fq.expected', 's.fq' => 'bam2fq/4.s.fq.expected', 'bc.fq' => 'bam2fq/bc.fq.expected'}, cmd=>"$$opts{bin}/samtools fastq @$threads --barcode-tag BC --index-format 'n2i2' --i1 $$opts{path}/bc.fq -s $$opts{path}/s.fq -1 $$opts{path}/1.fq -2 $$opts{path}/2.fq $$opts{path}/dat/bam2fq.004.sam");
    # tag output test with separators and -N flag
    test_cmd($opts, out=>'bam2fq/2.stdout.expected', out_map=>{'1.fq' => 'bam2fq/5.1.fq.expected', '2.fq' => 'bam2fq/5.2.fq.expected', 's.fq' => 'bam2fq/5.s.fq.expected', 'bc_split.fq' => 'bam2fq/bc_split.fq.expected'}, cmd=>"$$opts{bin}/samtools fastq @$threads --barcode-tag BC -N --index-format 'n*i*' --i1 $$opts{path}/bc_split.fq -s $$opts{path}/s.fq -1 $$opts{path}/1.fq -2 $$opts{path}/2.fq $$opts{path}/dat/bam2fq.005.sam");
    # -t flag
    test_cmd($opts, out=>'bam2fq/2.stdout.expected', out_map=>{'1.fq' => 'bam2fq/6.1.fq.expected', '2.fq' => 'bam2fq/6.2.fq.expected', 's.fq' => 'bam2fq/6.s.fq.expected'}, cmd=>"$$opts{bin}/samtools fastq @$threads -N -t -s $$opts{path}/s.fq -1 $$opts{path}/1.fq -2 $$opts{path}/2.fq $$opts{path}/dat/bam2fq.005.sam");
    # -T flag
    test_cmd($opts, out=>'bam2fq/2.stdout.expected', out_map=>{'1.fq' => 'bam2fq/7.1.fq.expected', '2.fq' => 'bam2fq/7.2.fq.expected', 's.fq' => 'bam2fq/7.s.fq.expected'}, cmd=>"$$opts{bin}/samtools fastq @$threads -N -t -T MD,ia -s $$opts{path}/s.fq -1 $$opts{path}/1.fq -2 $$opts{path}/2.fq $$opts{path}/dat/bam2fq.005.sam");
    # -i flag with no index
    test_cmd($opts, out=>'bam2fq/2.stdout.expected', out_map=>{'1.fq' => 'bam2fq/7.1.fq.expected', '2.fq' => 'bam2fq/7.2.fq.expected', 's.fq' => 'bam2fq/7.s.fq.expected'}, cmd=>"$$opts{bin}/samtools fastq @$threads -N -t -i -T MD,ia -s $$opts{path}/s.fq -1 $$opts{path}/1.fq -2 $$opts{path}/2.fq $$opts{path}/dat/bam2fq.005.sam");
    # -i flag with index
    test_cmd($opts, out=>'bam2fq/2.stdout.expected', out_map=>{'1.fq' => 'bam2fq/8.1.fq.expected', '2.fq' => 'bam2fq/8.2.fq.expected', 's.fq' => 'bam2fq/8.s.fq.expected', 'i.fq' => 'bam2fq/8.i.fq.expected'}, cmd=>"$$opts{bin}/samtools fastq @$threads --barcode-tag BC -i --index-format 'n2i2' --i1 $$opts{path}/i.fq -s $$opts{path}/s.fq -1 $$opts{path}/1.fq -2 $$opts{path}/2.fq $$opts{path}/dat/bam2fq.004.sam");
}

sub test_depad
{
    my ($opts) = @_;

    my $test_name = "test_depad";
    print "$test_name:\n";

    my $pad_sam   = "$$opts{path}/dat/depad.001p.sam";
    my $unpad_sam = "$$opts{path}/dat/depad.001u.sam";
    my $ref       = "$$opts{path}/dat/depad.001.fa";

    my $out = "$$opts{tmp}/depad";
    my $pad_bam = "$out.001p.bam";
    my $pad_cram = "$out.001p.cram";

    # Make a BAM file with the test data
    run_view_test($opts,
                  msg =>  "Generate BAM file",
                  args => ['-b', $pad_sam],
                  out => $pad_bam,
                  compare_sam => $pad_sam);

    # Don't try CRAM for now as it loses the CIGAR string from the embedded
    # reference.

    # # Make a CRAM file with the test data
    # run_view_test($opts,
    #               msg =>  "Generate CRAM file",
    #               args => ['-C', $pad_sam],
    #               out => $pad_cram,
    #               ref_path => "$$opts{path}/dat/cram_md5",
    #               compare_sam => $pad_sam);

    my $test = 1;
#    my @inputs = ([SAM => $pad_sam], [BAM => $pad_bam], [CRAM => $pad_cram]);
    my @inputs = ([SAM => $pad_sam], [BAM => $pad_bam]);
    my @formats = (['', 'bam'], ['-s', 'sam'], ['-u', 'uncomp.bam'],
                   ['-1', 'fast.bam']);
    foreach my $input (@inputs) {
        foreach my $format (@formats) {
            # Only test -T option for now due to problems with reformatting
            # the @SQ lines in the headers with embedded references.
            foreach my $use_t (1) {
                my @args = $use_t ? ('-T', $ref) : ();
                if ($format->[0]) { push(@args, $format->[0]); }
                push(@args, $input->[1]);
                my $compare = $format->[1] eq 'sam' ? 'compare' : 'compare_sam';
                run_view_test($opts,
                              msg => "depad $format->[0] ($input->[0] input)",
                              cmd => 'depad',
                              args => \@args,
                              out => sprintf("%s.test%03d.%s",
                                             $out, $test, $format->[1]),
                              ref_path => "$$opts{path}/dat/cram_md5",
                              redirect => 1,
                              $compare => $unpad_sam);
            }
        }
    }
}

sub test_stats
{
    my ($opts,%args) = @_;

    my $efix = ($^O =~ /^(?:msys|MSWin32)$/) ? 1 : 0;

    test_cmd($opts,out=>'stat/1.stats.expected',cmd=>"$$opts{bin}/samtools stats -r $$opts{path}/stat/test.fa $$opts{path}/stat/1_map_cigar.sam | tail -n+4", exp_fix=>$efix);
    test_cmd($opts,out=>'stat/2.stats.expected',cmd=>"$$opts{bin}/samtools stats -r $$opts{path}/stat/test.fa $$opts{path}/stat/2_equal_cigar_full_seq.sam | tail -n+4", exp_fix=>$efix);
    test_cmd($opts,out=>'stat/3.stats.expected',cmd=>"$$opts{bin}/samtools stats -r $$opts{path}/stat/test.fa $$opts{path}/stat/3_map_cigar_equal_seq.sam | tail -n+4", exp_fix=>$efix);
    test_cmd($opts,out=>'stat/4.stats.expected',cmd=>"$$opts{bin}/samtools stats -r $$opts{path}/stat/test.fa $$opts{path}/stat/4_X_cigar_full_seq.sam | tail -n+4", exp_fix=>$efix);
    test_cmd($opts,out=>'stat/5.stats.expected',cmd=>"$$opts{bin}/samtools stats -r $$opts{path}/stat/test.fa $$opts{path}/stat/5_insert_cigar.sam | tail -n+4", exp_fix=>$efix);
    test_cmd($opts,out=>'stat/6.stats.expected',cmd=>"$$opts{bin}/samtools stats -r $$opts{path}/stat/test.fa -i 0 $$opts{path}/stat/5_insert_cigar.sam | tail -n+4", exp_fix=>$efix);
    test_cmd($opts,out=>'stat/7.stats.expected',cmd=>"$$opts{bin}/samtools stats -r $$opts{path}/stat/test.fa $$opts{path}/stat/7_supp.sam | tail -n+4", exp_fix=>$efix);
    test_cmd($opts,out=>'stat/8.stats.expected',cmd=>"$$opts{bin}/samtools stats -r $$opts{path}/stat/test.fa $$opts{path}/stat/8_secondary.sam | tail -n+4", exp_fix=>$efix);

    test_cmd($opts,out=>'stat/9.stats.expected',cmd=>"$$opts{bin}/samtools stats -S RG -r $$opts{path}/stat/test.fa $$opts{path}/stat/1_map_cigar.sam | tail -n+4", exp_fix=>$efix,out_map=>{"stat/1_map_cigar.sam_s1_a_1.bamstat"=>"stat/1_map_cigar.sam_s1_a_1.expected.bamstat"},hskip=>3);
    test_cmd($opts,out=>'stat/10.stats.expected',cmd=>"$$opts{bin}/samtools stats -S RG -r $$opts{path}/stat/test.fa $$opts{path}/stat/10_map_cigar.sam | tail -n+4", exp_fix=>$efix,out_map=>{"stat/10_map_cigar.sam_s1_a_1.bamstat"=>"stat/10_map_cigar.sam_s1_a_1.expected.bamstat", "stat/10_map_cigar.sam_s1_b_1.bamstat"=>"stat/10_map_cigar.sam_s1_b_1.expected.bamstat"},hskip=>3);
}

sub test_merge
{
    my ($opts,%args) = @_;

    my $threads = exists($args{threads}) ? " -@ $args{threads}" : "";

    # Note the use of -s 1 to fix the random seed in place

    # Merge 1 - Standard 3 file SAM merge all presented on the command line
    test_cmd($opts,out=>'merge/2.merge.expected.bam',cmd=>"$$opts{bin}/samtools merge${threads} -s 1 - $$opts{path}/dat/test_input_1_a.sam $$opts{path}/dat/test_input_1_b.sam $$opts{path}/dat/test_input_1_c.sam",binary=>1);
    # Merge 2 - Standard 3 file BAM merge all files presented on the command line
    test_cmd($opts,out=>'merge/2.merge.expected.bam',cmd=>"$$opts{bin}/samtools merge${threads} -s 1 - $$opts{path}/dat/test_input_1_a.bam $$opts{path}/dat/test_input_1_b.bam $$opts{path}/dat/test_input_1_c.bam",binary=>1);
    # Merge 3 - Standard 3 file BAM merge 2 files in fofn 1 on command line
    open(my $fofn, "$$opts{path}/merge/test_3.fofn");
    my ($tmpfile_fh, $tmpfile_filename) = tempfile(UNLINK => 1);

    while (<$fofn>) {
        print $tmpfile_fh "$$opts{path}/$_";
    }
    close($tmpfile_fh);
    test_cmd($opts,out=>'merge/3.merge.expected.bam', err=>'merge/3.merge.expected.err',cmd=>"$$opts{bin}/samtools merge${threads} -s 1 -b $tmpfile_filename - $$opts{path}/dat/test_input_1_a.bam",binary=>1);
    # Merge 4 - 1 file BAM merge with file presented on the command line
    test_cmd($opts,out=>'merge/4.merge.expected.bam',cmd=>"$$opts{bin}/samtools merge${threads} -s 1 - $$opts{path}/dat/test_input_1_b.bam",binary=>1);
    # Merge 5 - 3 file SAM merge all presented on the command line override IDs to file names
    test_cmd($opts,out=>'merge/5.merge.expected.bam',cmd=>"$$opts{bin}/samtools merge${threads} -r -s 1 - $$opts{path}/dat/test_input_1_a.sam $$opts{path}/dat/test_input_1_b.sam $$opts{path}/dat/test_input_1_c.sam",binary=>1);
    # Merge 6 - merge all presented on the command line, combine PG and RG rather than dedup
    test_cmd($opts,out=>'merge/6.merge.expected.bam',cmd=>"$$opts{bin}/samtools merge${threads} -cp -s 1 - $$opts{path}/dat/test_input_1_a.sam $$opts{path}/dat/test_input_1_b.sam",binary=>1);
    # Merge 7 - ID and SN with regex in them
    test_cmd($opts,out=>'merge/7.merge.expected.bam',cmd=>"$$opts{bin}/samtools merge${threads} -s 1 - $$opts{path}/dat/test_input_1_a_regex.sam $$opts{path}/dat/test_input_1_b_regex.sam",binary=>1);

    # Sort inputs by PG, then merge
    system("$$opts{bin}/samtools sort -o $$opts{tmp}/merge.tag.1.bam -t PG  $$opts{path}/dat/test_input_1_b.sam") == 0 or die "failed to create sort BAM: $?";
    system("$$opts{bin}/samtools sort -o $$opts{tmp}/merge.tag.2.bam -t PG  $$opts{path}/dat/test_input_1_d.sam") == 0 or die "failed to create sort BAM: $?";
    test_cmd($opts,out=>'merge/tag.pg.merge.expected.sam',cmd=>"$$opts{bin}/samtools merge${threads} -s 1 -p -c -t PG -O SAM - $$opts{tmp}/merge.tag.1.bam $$opts{tmp}/merge.tag.2.bam");

    # Sort inputs by PG, then merge (name sorted)
    system("$$opts{bin}/samtools sort -o $$opts{tmp}/merge.tag.3.bam -n -t PG  $$opts{path}/dat/test_input_1_c.sam") == 0 or die "failed to create sort BAM: $?";
    system("$$opts{bin}/samtools sort -o $$opts{tmp}/merge.tag.4.bam -n -t PG  $$opts{path}/dat/test_input_1_d.sam") == 0 or die "failed to create sort BAM: $?";
    test_cmd($opts,out=>'merge/tag.pg.n.merge.expected.sam',cmd=>"$$opts{bin}/samtools merge${threads} -s 1 -p -c -n -t PG -O SAM - $$opts{tmp}/merge.tag.3.bam $$opts{tmp}/merge.tag.4.bam");
}

sub test_sort
{
    my ($opts, %args) = @_;

    my $threads = exists($args{threads}) ? " -@ $args{threads}" : "";

    # TODO Sort test cases

    # Check obsolete invocation is detected
    test_cmd($opts, out=>"dat/empty.expected", cmd=>"$$opts{bin}/samtools sort${threads} $$opts{path}/dat/test_input_1_a.bam $$opts{tmp}/sortout", want_fail=>1);
    test_cmd($opts, out=>"dat/empty.expected", cmd=>"$$opts{bin}/samtools sort${threads} -f $$opts{path}/dat/test_input_1_a.bam $$opts{tmp}/sortout.bam", want_fail=>1);
    test_cmd($opts, out=>"dat/empty.expected", cmd=>"$$opts{bin}/samtools sort${threads} -o $$opts{path}/dat/test_input_1_a.bam $$opts{tmp}/sorttmp", want_fail=>1);


    # Pos sort
    test_cmd($opts, out=>"sort/pos.sort.expected.sam", cmd=>"$$opts{bin}/samtools sort${threads}  $$opts{path}/dat/test_input_1_a.bam -O SAM -o -");

    # Name sort
    test_cmd($opts, out=>"sort/name.sort.expected.sam", cmd=>"$$opts{bin}/samtools sort${threads} -n  $$opts{path}/dat/test_input_1_a.bam -O SAM -o -");

    # Tag sort (RG)
    test_cmd($opts, out=>"sort/tag.rg.sort.expected.sam", cmd=>"$$opts{bin}/samtools sort${threads} -t RG  $$opts{path}/dat/test_input_1_a.bam -O SAM -o -");

    # Tag sort (RG); secondary by name
    test_cmd($opts, out=>"sort/tag.rg.n.sort.expected.sam", cmd=>"$$opts{bin}/samtools sort${threads} -n -t RG  $$opts{path}/dat/test_input_1_a.bam -O SAM -o -");

    # Tag sort (AS)
    test_cmd($opts, out=>"sort/tag.as.sort.expected.sam", cmd=>"$$opts{bin}/samtools sort${threads} -t AS $$opts{path}/dat/test_input_1_d.sam -O SAM -o -");

    # Tag sort (FI)
    test_cmd($opts, out=>"sort/tag.fi.sort.expected.sam", cmd=>"$$opts{bin}/samtools sort${threads} -t FI $$opts{path}/dat/test_input_1_d.sam -O SAM -o -");
}

sub test_fixmate
{
    my ($opts,%args) = @_;

    my $threads = exists($args{threads}) ? " -@ $args{threads}" : "";
    test_cmd($opts,out=>'fixmate/1_coord_sort.sam.expected', err=>'fixmate/1_coord_sort.sam.expected.err', cmd=>"$$opts{bin}/samtools fixmate${threads} -O sam $$opts{path}/fixmate/1_coord_sort.sam -", expect_fail=>1);
    test_cmd($opts,out=>'fixmate/2_isize_overflow.sam.expected', cmd=>"$$opts{bin}/samtools fixmate${threads} -O sam $$opts{path}/fixmate/2_isize_overflow.sam -");
    test_cmd($opts,out=>'fixmate/3_reverse_read_pp_lt.sam.expected', cmd=>"$$opts{bin}/samtools fixmate${threads} -O sam $$opts{path}/fixmate/3_reverse_read_pp_lt.sam -");
    test_cmd($opts,out=>'fixmate/4_reverse_read_pp_equal.sam.expected', cmd=>"$$opts{bin}/samtools fixmate${threads} -O sam $$opts{path}/fixmate/4_reverse_read_pp_equal.sam -");
    test_cmd($opts,out=>'fixmate/5_ct.sam.expected', cmd=>"$$opts{bin}/samtools fixmate${threads} -cO sam $$opts{path}/fixmate/5_ct.sam -");
    test_cmd($opts,out=>'fixmate/6_ct_replace.sam.expected', cmd=>"$$opts{bin}/samtools fixmate${threads} -cO sam $$opts{path}/fixmate/6_ct_replace.sam -");
    test_cmd($opts,out=>'fixmate/7_two_read_mapped.sam.expected', cmd=>"$$opts{bin}/samtools fixmate${threads} -O sam $$opts{path}/fixmate/7_two_read_mapped.sam -");
}

sub test_calmd
{
    my ($opts, %args) = @_;

    my $threads = exists($args{threads}) ? " -@ $args{threads}" : "";
    my $test = "$$opts{bin}/samtools calmd${threads} -uAr $$opts{path}/dat/mpileup.1.sam $$opts{path}/dat/mpileup.ref.fa";
    print "$test\n";
    my $out = cmd($test);
    if (substr($out, 0, 2) eq "\x1f\x8b") { passed($opts,msg=>$test); }
    else { failed($opts,msg=>$test,reason=>"Expected BGZF-compressed output"); }
}

sub test_idxstat
{
    my ($opts,%args) = @_;

    test_cmd($opts,out=>'idxstats/test_input_1_a.bam.expected', err=>'idxstats/test_input_1_a.bam.expected.err', cmd=>"$$opts{bin}/samtools idxstats $$opts{path}/dat/test_input_1_a.bam", expect_fail=>0);
}

sub test_quickcheck
{
    my ($opts,%args) = @_;

    my @testfiles = (
        'quickcheck/1.quickcheck.badeof.bam',
        'quickcheck/2.quickcheck.badheader.bam',
        'quickcheck/3.quickcheck.ok.bam',
        'quickcheck/4.quickcheck.ok.bam',
        'quickcheck/5.quickcheck.scramble30.truncated.cram',
        'quickcheck/6.quickcheck.cram21.ok.cram',
        'quickcheck/7.quickcheck.cram30.ok.cram',
        'quickcheck/8.quickcheck.cram21.truncated.cram',
        'quickcheck/9.quickcheck.cram30.truncated.cram',
        );

    my $all_testfiles;

    foreach my $fn (@testfiles) {
        $all_testfiles .= " $$opts{path}/$fn";
        test_cmd($opts, out => 'dat/empty.expected',
            want_fail => ($fn !~ /[.]ok[.]/),
            cmd => "$$opts{bin}/samtools quickcheck $$opts{path}/$fn");
    }

    test_cmd($opts, out => 'quickcheck/all.expected', want_fail => 1,
        cmd => "$$opts{bin}/samtools quickcheck -v $all_testfiles | sed 's,.*/quickcheck/,,'");
}

sub test_reheader
{
    my ($opts,%args) = @_;

    local $ENV{REF_PATH} = "$$opts{path}/dat/cram_md5/%s";

    my $fn = "$$opts{path}/dat/view.001";

    # Create local BAM and CRAM inputs
    system("$$opts{bin}/samtools view -b $fn.sam > $fn.tmp.bam")  == 0 or die "failed to create bam: $?";
    system("$$opts{bin}/samtools view -C --output-fmt-option version=2.1 $fn.sam > $fn.tmp.v21.cram") == 0 or die "failed to create cram: $?";
    system("$$opts{bin}/samtools view -C --output-fmt-option version=3.0 $fn.sam > $fn.tmp.v30.cram") == 0 or die "failed to create cram: $?";

    # Fudge @PG lines.  The version number will differ each commit.
    # Also the pathname will differ for each install. We'll take it on faith
    # that these bits work.
    test_cmd($opts,
             out=>'reheader/1_view1.sam.expected',
             err=>'reheader/1_view1.sam.expected.err',
             cmd=>"$$opts{bin}/samtools reheader $$opts{path}/reheader/hdr.sam $fn.tmp.bam | $$opts{bin}/samtools view -h | perl -pe 's/\tVN:.*//'",
	     exp_fix=>1);

    test_cmd($opts,
             out=>'reheader/2_view1.sam.expected',
             err=>'reheader/2_view1.sam.expected.err',
             cmd=>"$$opts{bin}/samtools reheader $$opts{path}/reheader/hdr.sam $fn.tmp.v21.cram | $$opts{bin}/samtools view -h | perl -pe 's/\tVN:.*//'",
	     exp_fix=>1);

    test_cmd($opts,
             out=>'reheader/2_view1.sam.expected',
             err=>'reheader/2_view1.sam.expected.err',
             cmd=>"$$opts{bin}/samtools reheader $$opts{path}/reheader/hdr.sam $fn.tmp.v30.cram | $$opts{bin}/samtools view -h | perl -pe 's/\tVN:.*//'",
	     exp_fix=>1);

    # In-place testing
    test_cmd($opts,
             out=>'reheader/3_view1.sam.expected',
             err=>'reheader/3_view1.sam.expected.err',
             cmd=>"$$opts{bin}/samtools reheader --in-place $$opts{path}/reheader/hdr.sam $fn.tmp.v21.cram && $$opts{bin}/samtools view -h $fn.tmp.v21.cram | perl -pe 's/\tVN:.*//'",
	     exp_fix=>1);

    test_cmd($opts,
             out=>'reheader/3_view1.sam.expected',
             err=>'reheader/3_view1.sam.expected.err',
             cmd=>"$$opts{bin}/samtools reheader --in-place $$opts{path}/reheader/hdr.sam $fn.tmp.v30.cram && $$opts{bin}/samtools view -h $fn.tmp.v30.cram | perl -pe 's/\tVN:.*//'",
	     exp_fix=>1);
}

sub test_addrprg
{
    my ($opts,%args) = @_;

    my $threads = exists($args{threads}) ? " -@ $args{threads}" : "";
    test_cmd($opts,out=>'addrprg/1_fixup.sam.expected', err=>'addrprg/1_fixup.sam.expected.err', cmd=>"$$opts{bin}/samtools addreplacerg${threads} -O sam -m overwrite_all $$opts{path}/addrprg/1_fixup.sam");
    test_cmd($opts,out=>'addrprg/2_fixup_orphan.sam.expected', err=>'addrprg/2_fixup_orphan.sam.expected.err', cmd=>"$$opts{bin}/samtools addreplacerg${threads} -O sam -m orphan_only $$opts{path}/addrprg/2_fixup_orphan.sam");
    test_cmd($opts,out=>'addrprg/3_fixup.sam.expected', err=>'addrprg/3_fixup.sam.expected.err', cmd=>"$$opts{bin}/samtools addreplacerg${threads} -O sam -R '1#7' $$opts{path}/addrprg/1_fixup.sam", want_fail=>1);
    test_cmd($opts,out=>'addrprg/4_fixup_norg.sam.expected', err=>'addrprg/4_fixup_norg.sam.expected.err', cmd=>"$$opts{bin}/samtools addreplacerg${threads} -O sam -r '\@RG\\tID:1#8\\tCN:SC' $$opts{path}/addrprg/4_fixup_norg.sam");
    test_cmd($opts,out=>'addrprg/1_fixup.sam.expected', err=>'addrprg/1_fixup.sam.expected.err', cmd=>"$$opts{bin}/samtools addreplacerg${threads} -O sam -m overwrite_all -R '1#8' $$opts{path}/addrprg/1_fixup.sam");
    test_cmd($opts,out=>'addrprg/4_fixup_norg.sam.expected', err=>'addrprg/4_fixup_norg.sam.expected.err', cmd=>"$$opts{bin}/samtools addreplacerg${threads} -O sam -r 'ID:1#8' -r 'CN:SC' $$opts{path}/addrprg/4_fixup_norg.sam");
}

sub test_markdup
{
    my ($opts,%args) = @_;

    my $threads = exists($args{threads}) ? " -@ $args{threads}" : "";
    test_cmd($opts, out=>'markdup/1_name_sort.expected.sam', err=>'1_name_sort.expected.sam.err', cmd=>"$$opts{bin}/samtools markdup${threads} -O sam $$opts{path}/markdup/1_name_sort.sam -", expect_fail=>1);
    test_cmd($opts, out=>'markdup/2_bad_order.expected.sam', err=>'2_bad_order.expected.sam.err', cmd=>"$$opts{bin}/samtools markdup${threads} -O sam $$opts{path}/markdup/2_bad_order.sam -", expect_fail=>1);
    test_cmd($opts, out=>'markdup/3_missing_mc.expected.sam', err=>'3_missing_mc.expected.sam.err', cmd=>"$$opts{bin}/samtools markdup${threads} -O sam $$opts{path}/markdup/3_missing_mc.sam -", expect_fail=>1);
    test_cmd($opts, out=>'markdup/4_missing_ms.expected.sam', err=>'4_missing_ms.expected.sam.err', cmd=>"$$opts{bin}/samtools markdup${threads} -O sam $$opts{path}/markdup/4_missing_ms.sam -", expect_fail=>1);
    test_cmd($opts, out=>'markdup/5_markdup.expected.sam', cmd=>"$$opts{bin}/samtools markdup${threads} -O sam $$opts{path}/markdup/5_markdup.sam -");
    test_cmd($opts, out=>'markdup/6_remove_dups.expected.sam', cmd=>"$$opts{bin}/samtools markdup${threads} -O sam -r $$opts{path}/markdup/6_remove_dups.sam -");
}
