#!/usr/bin/perl
# check_copyright.pl : Basic source file checks for copyright boilerplate
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

    # Exclusions:
    my %exclude = map { ("$root/$_", 1) } (
        'config.h',         # Auto-generated
        'version.h',        # Auto-generated
    );
    return if exists($exclude{$_});

    my $remove_left = /\.[ch]$/ ? qr/\s*\*?\s*/ : qr/\s*#\s*/;

    return unless (-f $_);      # Only check plain files
    my $in;
    if (!open($in, '<', $_)) {
        print STDERR "Couldn't open $_ : $!\n";
        $errors++;
        return;
    }
    my $count = 0;
    my $copyright_found = 0;
    my $license_found = "";
    my $line;
    while ($count < 100 && ($line = <$in>)) {
        $count++;
        $line =~ s/^$remove_left//;
        $line =~ s/\s+/ /g;
        if ($line =~ /^Copyright\s+\([cC]\)\s+(?:19|20)\d\d[-, ]/) {
            $copyright_found = 1;
        } elsif ($line =~ /^Redistribution and use in source and binary forms/) {
            $license_found = "BSD";
        } elsif ($line =~ /^Permission is hereby granted, free of charge/) {
            $license_found = "MIT";
        }
        last if ($copyright_found && $license_found);
    }
    if (!close($in)) {
        print STDERR "Error on closing $_ : $!\n";
        $errors++;
        return;
    }
    my $failed = (!$copyright_found || !$license_found);
    if ($verbose || $failed) {
        printf("$_ : %s%s\n",
               $license_found ? $license_found : "no_license",
               $copyright_found ? "" : " no_copyright_line");
    }
    if ($failed) {
        $errors++;
    }
}
