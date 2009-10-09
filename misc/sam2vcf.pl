#!/usr/bin/perl -w
# 
# VCF specs: http://www.1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcfv3.2

# Contact: pd3@sanger
# Version: 2009-10-08

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
        "   -r, -refseq <file.fa>            The reference sequence, required when indels are present.\n",
        "   -h, -?, --help                   This help message.\n",
        "\n";
}


sub parse_params
{
    my %opts = ();

    $opts{fh_in}  = *STDIN;
    $opts{fh_out} = *STDOUT;

    while (my $arg=shift(@ARGV))
    {
        if ( $arg eq '-r' || $arg eq '--refseq' ) { $opts{refseq}=shift(@ARGV); next; }
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
        if ( $ref eq $base ) { return ('.','0|0'); }
        return ($base,'1|1');
    }
    my $gt = $iupac{$base};
    if ( $$gt[0] eq $ref  ) { return ($$gt[1],'0|1'); }
    elsif ( $$gt[1] eq $ref ) { return ($$gt[0],'0|1'); }
    return ("$$gt[0],$$gt[1]",'1|2');
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

    while (my $line=<$fh_in>)
    {
        chomp($line);
        my ($chr,$pos,$ref,$cons,$cons_qual,$snp_qual,$rms_qual,$depth,@items) = split(/\t/,$line);

        my ($alt,$gt);
        if ( $ref eq '*' )
        {
            # An indel is involved.
            if ($chr ne $prev_chr || $pos ne $prev_pos) 
            {
                if ( !$$opts{refseq} ) { error("Cannot do indels without the reference.\n"); }
                if ( !$refseq ) { $refseq = Fasta->new(file=>$$opts{refseq}); }
                $ref = $refseq->get_base($chr,$pos);
            }
            else { $ref = $prev_ref; }

            # One of the alleles can be a reference and it can come in arbitrary order
            my ($al1,$al2) = split(m{/},$cons);
            my $alt1 = parse_indel($al1);
            my $alt2 = parse_indel($al2);
            if ( !$alt1 && !$alt2 ) { error("FIXME: could not parse indel:\n", $line); }
            if ( $alt1 && $alt2 && $alt1 eq $alt2 ) { $alt2=''; }
            if ( !$alt1 ) 
            { 
                $alt=$alt2; 
                $gt='0|1'; 
            }
            elsif ( !$alt2 ) 
            { 
                $alt=$alt1; 
                $gt='0|1'; 
            }
            else 
            { 
                $alt="$alt1,$alt2"; 
                $gt='1|2'; 
            }
        }
        else
        {
            # SNP
            ($alt,$gt) = iupac_to_gtype($ref,$cons);
        }

        print $fh_out "$chr\t$pos\t.\t$ref\t$alt\t$snp_qual\t0\t\tGT:GQ:DP\t$gt:$cons_qual:$depth\n";

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
    my $self = @args ? {@args} : {};
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
