#!/usr/bin/env perl

# The MIT License

# Copyright (c) 2014 Genome Research Ltd.
# Author: Rob Davies <rmd+sam@sanger.ac.uk>

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# Import references into a cram reference cache from fasta files.
# See below __END__ for POD documentation.

use strict;
use warnings;
use Digest::MD5;
use Getopt::Long;
use File::Find;
use File::Temp qw(tempfile);
use File::Spec::Functions;
use File::Path 'make_path';
use IO::Handle;

$| = 1;

# Directory where the cache will be built
my $root_dir;

# Number of subdirectories to make below $root_dir
# Each subdir will eat up two hex digits of the file MD5
my $subdirs = 2;

# Directory tree to search when using the -find option
my $find = '';

# How much data to read before spilling to a file
my $max_acc = 256 * 1024 * 1024;

my $usage = "Usage: $0 -root <dir> [-subdirs <n>] input1.fasta ...\n       $0 -root <dir> [-subdirs <n>] -find <dir>\n";

# Deal with options
GetOptions("root=s" => \$root_dir, "subdirs=s" => \$subdirs,
       "find=s" => \$find) || die $usage;

unless ($root_dir && $subdirs =~ /^\d+$/) { die $usage; }
if ($subdirs >= 16) {
    die "$0: Error: -subdirs should be less than 15.\n";
}

# Regexp to convert a hex MD5 to a list of $subdirs subdirectory names, the
# remainder making the filename in the leaf directory
my $dest_regexp = "(..)" x $subdirs . "(" . ".." x (16 - $subdirs) . ")";
my $dest_re = qr/$dest_regexp/;

# Ensure $root_dir exists
unless (-e $root_dir) {
    make_path($root_dir);
}

if ($find) {
    # Find mode - search a directory tree for anything that looks like a
    # fasta file.  Any that are found will be put into the new cache, if
    # they are not already there.
    find({
    wanted => sub {
        find_files($File::Find::name, $root_dir, $dest_re, $max_acc);
    },
    no_chdir => 1,
     },
     $find);
} elsif (@ARGV) {
    # If a list of files was given on the command line, go through them
    # and try to add each one.
    foreach my $name (@ARGV) {
    open(my $fh, '<', $name) || die "Couldn't open $name: $!\n";
    process_file($name, $fh, $root_dir, $dest_re, $max_acc);
    close($fh) || die "Error closing $name: $!\n";
    }
} else {
    # Otherwise read from STDIN
    process_file('STDIN', \*STDIN, $root_dir, $dest_re, $max_acc);
}
exit;

sub find_files {
    my ($name, $root_dir, $dest_re, $max_acc) = @_;

    # See if $name is a candidate file

    my $fh;
    return if ($name =~ /~$/);        # Ignore backup files
    return unless (-f $name && -r _); # Ignore non-regular and unreadable files

    # Inspect the first two lines of the candidate
    my $buffer;
    open($fh, '<', $name) || die "Couldn't open $name: $!\n";
    read($fh, $buffer, 8192); # Should be enough to find the header & sequence
    close($fh) || die "Error closing $name: $!\n";
    my ($l1, $l2) = split(/\n/, $buffer);

    # Check for fasta-like content
    return unless ($l1 && $l1 =~ /^>\S+/);
    return unless ($l2 && $l2 =~ /^[ACGTMRWSYKVHDBNacgtmrwsykvhdbn]+$/);

    # Looks like a fasta file, so process it
    open($fh, '<', $name) || die "Couldn't open $name: $!\n";
    process_file($name, $fh, $root_dir, $dest_re, $max_acc);
    close($fh) || die "Error closing $name: $!\n";
}

sub process_file {
    my ($name, $in_fh, $root_dir, $dest_re, $max_acc) = @_;

    # Process the fasta file $in_fh.  Each entry in the file is read, and
    # the MD5 calculated as described in the SAM specification (i.e.
    # all uppercased with whitespace stripped out).  The MD5 is then
    # split using $dest_re to convert it into the path name for the entry
    # in the cache.  If the path is not already present, the entry (in
    # uppercased and stripped form) is saved into the cache.

    # Entries shorter that $max_acc will be kept in memory.  For fasta files
    # with lots of short entries this can save a lot of unnecessary writing
    # if the data is already in the cache.  Anything longer
    # gets written out to a file to keep memory consumption under control.
    # The temporary files have to be made in $root_dir, as the final
    # destination is not known until the entry has been completely read.

    my $id;        # Name of current fasta entry
    my $ctx;       # MD5 context
    my $acc = '';  # The accumulated sequence
    my $tmpfile;   # Temporary file name
    my $tmpfh;     # Temporary file handle
    my $extra = 1024; # Extra space to pre-allocate to account for reading
                      # 1 line past $max_acc
    vec($acc, $max_acc + $extra, 8) = 1; # Pre-allocate some space
    $acc = '';

    # Use an eval block so any stray temporary file can be cleaned up before
    # exiting.
    eval {
    print "Reading $name ...\n";
    for (;;) { # Error catching form of while (<>) {...}
        undef($!);
        last if (eof($in_fh)); # Needed if last line isn't terminated
        unless (defined($_ = readline($in_fh))) {
        die "Error reading $name: $!" if $!;
        last; # EOF
        }

        if (/^>(\S+)/) {
        # Found a fasta header
        if ($ctx) { # Finish previous entry, if there is one
            finish_entry($id, $ctx, \$acc, $tmpfh, $tmpfile,
                 $root_dir, $dest_re);
            undef($tmpfile);
            $acc = '';
        }
        $id = $1;
        $ctx = Digest::MD5->new();
        } else {
        unless ($id) { die "Found sequence with no header\n"; }
        # Read some sequence
        chomp;
        s/\s+//g;
        if ($_) {
            $_ = uc($_);
            $acc .= $_;
            $ctx->add($_);

            if (length($acc) > $max_acc) {
            # Spill long sequences out to a temporary file in
            # $root_dir.
            unless ($tmpfile) {
                ($tmpfh, $tmpfile) = tempfile(DIR => $root_dir,
                              SUFFIX => '.tmp');
            }
            print $tmpfh $acc
                || die "Error writing to $tmpfile: $!\n";
            $acc = '';
            }
        }
        }
    }
    if ($ctx) {
        # Finish off the last entry
        finish_entry($id, $ctx, \$acc, $tmpfh, $tmpfile,
             $root_dir, $dest_re);
        undef($tmpfile);
    }
    };
    my $err = $@;
    if ($tmpfile) { unlink($tmpfile); }
    if ($err) { die $err; }
}

sub finish_entry {
    my ($id, $ctx, $acc_ref, $tmpfh, $tmpfile, $root_dir, $dest_re) = @_;

    # Finish writing an entry

    my $digest = $ctx->hexdigest;

    # Get the destination directory and filename
    my @segs = $digest =~ /$dest_re/;
    my $dest_dir = (@segs > 1
            ? catdir($root_dir, @segs[0..($#segs - 1)])
            : $root_dir);
    my $dest = catfile($dest_dir, $segs[-1]);

    # Make the destination dir if necessary
    unless (-e $dest_dir) {
    make_path($dest_dir);
    }

    if (-e $dest) {
    # If the file is already present, there's nothing to do apart from
    # remove the temporary file if it was made.
    print "Already exists: $digest $id\n";
    if ($tmpfile) {
        close($tmpfh) || die "Error closing $tmpfile: $!\n";
        unlink($tmpfile) || die "Couldn't remove $tmpfile: $!\n";
    }
    } else {
    # Need to add the data to the cache.
    unless ($tmpfile) {
        # If the data hasn't been written already, it needs to be done
        # now.  Write to a temp file in $dest_dir so if it goes wrong
        # we won't leave a file with the right name but half-written
        # content.
        ($tmpfh, $tmpfile) = tempfile(DIR => $dest_dir,
                      SUFFIX => '.tmp');
    }

    # Assert that the $tmpfile is now set
    unless ($tmpfile) { die "Error: Didn't make a temp file"; }

    eval {
        # Flush out any remaining data
        if ($$acc_ref) {
        print $tmpfh $$acc_ref || die "Error writing to $tmpfile: $!\n";
        }
        # Paranoid file close
        $tmpfh->flush() || die "Error flushing to $tmpfile: $!\n";
        $tmpfh->sync() || die "Error syncing $tmpfile: $!\n";
        close($tmpfh) || die "Error writing to $tmpfile: $!\n";
    };
    if ($@) {
        # Attempt to clean up if writing failed
        my $save = $@;
        unlink($tmpfile) || warn "Couldn't remove $tmpfile: $!";
        die $save;
    }

    # Finished writing, and everything is flushed as far as possible
    # so rename the temp file
    print "$dest $id\n";
    rename($tmpfile, $dest)
        || die "Error moving $tmpfile to $dest: $!\n";
    }
}

__END__

=head1 NAME

seq_cache_populate.pl

=head1 SYNOPSIS

 seq_cache_populate.pl -root <dir> [-subdirs <n>] input1.fasta ...

 seq_cache_populate.pl -root <dir> [-subdirs <n>]  -find <dir>

=head1 DESCRIPTION

Import references into a cram reference cache from fasta files.

When run with a list of fasta files, this program reads the files and stores
the sequences within it in the reference cache under directory <dir>.  The
sequences in the cache are stored with names based on the MD5 checksum
of the sequence.

By default, sequences are stored in a hierarchy two directories deep, to
keep the number of items in a single directory to a reasonable number.  This
depth can be chaged using the -subdirs option.

If the -find option is used, the program will scan the given directory tree.
Any files that appear to be fasta (by looking for a first line starting with
'>' followed by somthing that looks like DNA sequence) will be read and
added to the reference cache.  The traversal will ignore symbolic links.

Samtools/htslib can be made to use the cache by appropriate setting of the
REF_PATH environment varaiable.  For example, if seq_cache_populate was run
using options '-root /tmp/ref_cache -subdirs 2', setting REF_PATH to
'/tmp/ref_cache/%2s/%2s/%s' should allow samtools to find the references that
it stored.

Note that if no REF_PATH is specified, htslib will default to downloading from
the EBI reference server and caching locally (see the samtools(1) man page for
details), defaulting to $HOME/.cache/hts-ref/%2s/%2s/%s.  This is functionally
equivalent to running this tool with '-root $HOME/.cache/hts-ref -subdirs 2'.

=head1 AUTHOR

Rob Davies.

=cut
