#!/usr/bin/env perl
#
#    Copyright (C) 2020 Genome Research Ltd.
#
#    Author: James Bonfield <jkb@sanger.ac.uk>
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


# Usage: fasta-sanitize.pl file.in > file.out
# Or via stdin, e.g. zcat file.in.gz | fasta-sanitize.pl > file.out
#
# Also supports and autodetects fastq.

# This tool sanitizes the reference names as per the SAM specification.
# See SAM pull request https://github.com/samtools/hts-specs/pull/333

# It is important that this is run prior to aligning.  This ensures the
# SAM file contains @SQ lines which adhere to the specification, and that
# VCF files produced from those SAM files also match the VCF specification.
#
# Furthermore, doing this early rather than via later with "samtools reheader"
# means that the fasta file matches the SAM file.
# Several samtools and bcftools sub-commands also require a fasta
# reference, which must match the primary SQ lines, plus it helps CRAM
# if using "view -T" (although it'll still work if using the
# recommended practice of an MD5sum based ref cache).

# The regexp permitted is:
# [0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*

use strict;

my $re = qr/^[0-9A-Za-z!#$%&+.\/:;?@^_|~-][0-9A-Za-z!#$%&*+.\/:;=?@^_|~-]*$/;

my $fastq = 0;
my $in_qual = 0;
my $seq_len = 0;

my $name_re = qr/^([>@])\s*(\S*)(.*)/;

while (<>) {
    # Name
    if (/$name_re/ && !$in_qual) {
        my ($prefix, $name, $other) = ($1,$2,$3);
        $fastq = ($prefix eq "@") ? 1 : 0;

        if ($name !~ /$re/) {
            my ($l,$r)=($name=~/(.)(.*)/);
            $l =~ tr/[0-9A-Za-z!#$%&+.\/:;?@^_|~\-]/_/c;
            $r =~ tr/[0-9A-Za-z!#$%&*+.\/:;=?@^_|~\-]/_/c;
            my $new_name = $l.$r;

            print STDERR "Renaming reference $name to $new_name\n";
            $name = $new_name;
            $seq_len = 0;
        }

        print "$prefix$name$other\n";
        next;
    }

    if (!$in_qual) {
        # FASTQ separator between seq and qual
        if ($fastq && /^\+/) {
            print;
            $in_qual = 1;
            next;
        }

        # Seq
        print;
        chomp($_);
        $seq_len += length($_);
    } else {
        # Qual
        print;
        chomp($_);
        $in_qual = 0 if (($seq_len -= length($_)) <= 0);
    }
}
