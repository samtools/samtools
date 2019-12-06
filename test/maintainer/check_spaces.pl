#!/usr/bin/perl
# check_spaces.pl : Check source files for tabs and trailing spaces
#
#     Author : Rob Davies <rmd@sanger.ac.uk>
#
#     Copyright (C) 2018-2019 Genome Research Ltd.
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
use File::Find;
use Getopt::Long;

my $verbose = 0;
GetOptions('v' => \$verbose);

my ($root) = @ARGV;
if (!$root) {
    die "Usage: $0 [-v] <directory>\n";
}
my $errors = 0;
find({ wanted => \&check, no_chdir=>1}, $root);
exit($errors ? 1 : 0);

sub check {
    # Ignore embedded copied of htslib
    if (m#/[^/]*htslib[^/]*$# && -d $_) {
        $File::Find::prune = 1;
        return;
    }

    # Only check C, perl and shell files
    return unless (/(?:\.[ch]|\.pl|\.sh)$/);

    my %allow_tabs = map { ("$root/$_", 1) } (
    );

    my $check_tabs = !exists($allow_tabs{$_});

    my $in;
    if (!open($in, '<', $_)) {
        print STDERR "Couldn't open $_ : $!\n";
        $errors++;
        return;
    }
    my $tab = 0;
    my $trailing = 0;
    while (my $line = <$in>) {
        chomp($line);
        if ($check_tabs && $line =~ /\t/)  { $tab = 1; }
        if ($line =~ /\s$/) { $trailing = 1; }
    }
    if (!close($in)) {
        print STDERR "Error on closing $_ : $!\n";
        $errors++;
        return;
    }
    my $failed = ($tab || $trailing);
    if ($verbose || $failed) {
        my $msg = ($failed ? join(" ",
                                  $tab ? ("includes_tabs") : (),
                                  $trailing ? "trailing_spaces" : ())
                   : "ok");
        print "$_ : $msg\n";
    }
    if ($failed) {
        $errors++;
    }
}
