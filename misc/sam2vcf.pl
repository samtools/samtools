#!/usr/bin/env perl
#
#    Copyright (C) 2009, 2010 Genome Research Ltd.
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

# VCF specs: http://www.1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcf3.3

use strict;
use warnings;
use Carp;

my $opts = parse_params();
do_pileup_to_vcf($opts);

exit;

#---------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { croak(@msg); }
    die
        "Usage: sam2vcf.pl [OPTIONS] < in.pileup > out.vcf\n",
        "Options:\n",
        "   -h, -?, --help                  This help message.\n",
        "   -i, --indels-only               Ignore SNPs.\n",
        "   -r, --refseq <file.fa>          The reference sequence, required when indels are present.\n",
        "   -R, --keep-ref                  Print reference alleles as well.\n",
        "   -s, --snps-only                 Ignore indels.\n",
        "   -t, --column-title <string>     The column title.\n",
        "\n";
}


sub parse_params
{
    my %opts = ();

    $opts{fh_in}  = *STDIN;
    $opts{fh_out} = *STDOUT;

    while (my $arg=shift(@ARGV))
    {
        if ( $arg eq '-R' || $arg eq '--keep-ref' ) { $opts{keep_ref}=1; next; }
        if ( $arg eq '-r' || $arg eq '--refseq' ) { $opts{refseq}=shift(@ARGV); next; }
        if ( $arg eq '-t' || $arg eq '--column-title' ) { $opts{title}=shift(@ARGV); next; }
        if ( $arg eq '-s' || $arg eq '--snps-only' ) { $opts{snps_only}=1; next; }
        if ( $arg eq '-i' || $arg eq '--indels-only' ) { $opts{indels_only}=1; next; }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }

        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    return \%opts;
}

sub iupac_to_gtype
{
    my ($ref,$base) = @_;
    my %iupac = (
            'K' => ['G','T'],
            'M' => ['A','C'],
            'S' => ['C','G'],
            'R' => ['A','G'],
            'W' => ['A','T'],
            'Y' => ['C','T'],
            );
    if ( !exists($iupac{$base}) )
    {
        if ( $base ne 'A' && $base ne 'C' && $base ne 'G' && $base ne 'T' ) { error("FIXME: what is this [$base]?\n"); }
        if ( $ref eq $base ) { return ('.','0/0'); }
        return ($base,'1/1');
    }
    my $gt = $iupac{$base};
    if ( $$gt[0] eq $ref  ) { return ($$gt[1],'0/1'); }
    elsif ( $$gt[1] eq $ref ) { return ($$gt[0],'0/1'); }
    return ("$$gt[0],$$gt[1]",'1/2');
}


sub parse_indel
{
    my ($cons) = @_;
    if ( $cons=~/^-/ )
    {
        my $len = length($');
        return "D$len";
    }
    elsif ( $cons=~/^\+/ ) { return "I$'"; }
    elsif ( $cons eq '*' ) { return undef; }
    error("FIXME: could not parse [$cons]\n");
}


# An example of the pileup format:
#   1       3000011 C       C       32      0       98      1       ^~,     A
#   1       3002155 *       +T/+T   53      119     52      5       +T      *       4       1       0
#   1       3003094 *       -TT/-TT 31      164     60      11      -TT     *       5       6       0
#   1       3073986 *       */-AAAAAAAAAAAAAA       3       3       45      9       *       -AAAAAAAAAAAAAA 7       2       0
#
sub do_pileup_to_vcf
{
    my ($opts) = @_;

    my $fh_in  = $$opts{fh_in};
    my $fh_out = $$opts{fh_out};
    my ($prev_chr,$prev_pos,$prev_ref);
    my $refseq;
    my $ignore_indels = $$opts{snps_only} ? 1 : 0;
    my $ignore_snps   = $$opts{indels_only} ? 1 : 0;
    my $keep_ref      = $$opts{keep_ref} ? 1 : 0;
    my $title = exists($$opts{title}) ? $$opts{title} : 'data';

    print $fh_out
        qq[##fileformat=VCFv3.3\n],
        qq[##INFO=DP,1,Integer,"Total Depth"\n],
        qq[##FORMAT=GT,1,String,"Genotype"\n],
        qq[##FORMAT=GQ,1,Integer,"Genotype Quality"\n],
        qq[##FORMAT=DP,1,Integer,"Read Depth"\n],
        qq[#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$title\n]
        ;

    while (my $line=<$fh_in>)
    {
        chomp($line);
        my (@items) = split(/\t/,$line);
        if ( scalar @items<8 )
        {
            error("\nToo few columns, does not look like output of 'samtools pileup -c': $line\n");
        }
        my ($chr,$pos,$ref,$cons,$cons_qual,$snp_qual,$rms_qual,$depth,$a1,$a2) = @items;
        $ref  = uc($ref);
        $cons = uc($cons);

        my ($alt,$gt);
        if ( $ref eq '*' )
        {
            # An indel is involved.
            if ( $ignore_indels )
            {
                $prev_ref = $ref;
                $prev_pos = $pos;
                $prev_chr = $chr;
                next;
            }

            if (!defined $prev_chr || $chr ne $prev_chr || $pos ne $prev_pos)
            {
                if ( !$$opts{refseq} ) { error("Cannot do indels without the reference.\n"); }
                if ( !$refseq ) { $refseq = Fasta->new(file=>$$opts{refseq}); }
                $ref = $refseq->get_base($chr,$pos);
                $ref = uc($ref);
            }
            else { $ref = $prev_ref; }

            # One of the alleles can be a reference and it can come in arbitrary order. In some
            #   cases */* can be encountered. In such a case, look in the additional columns.
            my ($al1,$al2) = split(m{/},$cons);
            if ( $al1 eq $al2 && $al1 eq '*' ) { $al1=$a1; $al2=$a2; }
            my $alt1 = parse_indel($al1);
            my $alt2 = parse_indel($al2);
            if ( !$alt1 && !$alt2 ) { error("FIXME: could not parse indel:\n", $line); }
            if ( !$alt1 )
            {
                $alt=$alt2;
                $gt='0/1';
            }
            elsif ( !$alt2 )
            {
                $alt=$alt1;
                $gt='0/1';
            }
            elsif ( $alt1 eq $alt2 )
            {
                $alt="$alt1";
                $gt='1/1';
            }
            else
            {
                $alt="$alt1,$alt2";
                $gt='1/2';
            }
        }
        else
        {
            if ( $ignore_snps || (!$keep_ref && $ref eq $cons) )
            {
                $prev_ref = $ref;
                $prev_pos = $pos;
                $prev_chr = $chr;
                next;
            }

            # SNP
            ($alt,$gt) = iupac_to_gtype($ref,$cons);
        }

        print $fh_out "$chr\t$pos\t.\t$ref\t$alt\t$snp_qual\t0\tDP=$depth\tGT:GQ:DP\t$gt:$cons_qual:$depth\n";

        $prev_ref = $ref;
        $prev_pos = $pos;
        $prev_chr = $chr;
    }
}


#------------- Fasta --------------------
#
# Uses samtools to get a requested base from a fasta file. For efficiency, preloads
#   a chunk to memory. The size of the cached sequence can be controlled by the 'size'
#   parameter.
#
package Fasta;

use strict;
use warnings;
use Carp;

sub Fasta::new
{
    my ($class,@args) = @_;
    my $self = {@args};
    bless $self, ref($class) || $class;
    if ( !$$self{file} ) { $self->throw(qq[Missing the parameter "file"\n]); }
    $$self{chr}  = undef;
    $$self{from} = undef;
    $$self{to}   = undef;
    if ( !$$self{size} ) { $$self{size}=10_000_000; }
    bless $self, ref($class) || $class;
    return $self;
}

sub read_chunk
{
    my ($self,$chr,$pos) = @_;
    my $to = $pos + $$self{size};
    my $cmd = "samtools faidx $$self{file} $chr:$pos-$to";
    my @out = `$cmd`;
    if ( $? ) { $self->throw("$cmd: $!"); }
    my $line = shift(@out);
    if ( !($line=~/^>$chr:(\d+)-(\d+)/) ) { $self->throw("Could not parse: $line"); }
    $$self{chr}  = $chr;
    $$self{from} = $1;
    $$self{to}   = $2;
    my $chunk = '';
    while ($line=shift(@out))
    {
        chomp($line);
        $chunk .= $line;
    }
    $$self{chunk} = $chunk;
    return;
}

sub get_base
{
    my ($self,$chr,$pos) = @_;
    if ( !$$self{chr} || $chr ne $$self{chr} || $pos<$$self{from} || $pos>$$self{to} )
    {
        $self->read_chunk($chr,$pos);
    }
    my $idx = $pos - $$self{from};
    return substr($$self{chunk},$idx,1);
}

sub throw
{
    my ($self,@msg) = @_;
    croak(@msg);
}
