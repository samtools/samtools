#!/usr/bin/perl -w

use strict;
use warnings;
use Carp;

my $opts = parse_params();
bcf_fix();

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    die
        "Usage: bcftools view test.bcf | bcf-fix.pl > test.vcf\n",
        "Options:\n",
        "   -h, -?, --help                  This help message.\n",
        "\n";
}


sub parse_params
{
    my $opts = {};
    while (my $arg=shift(@ARGV))
    {
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    return $opts;
}

sub bcf_fix
{
    while (my $line=<STDIN>)
    {
        if ( $line=~/^#CHROM/ )
        {
            print 
qq[##INFO=<ID=DP4,Number=4,Type=Integer,Description="Read depth for 1) forward reference bases; 2) reverse ref; 3) forward non-ref; 4) reverse non-ref">
##INFO=<ID=PV4,Number=4,Type=Float,Description="P-values for 1) strand bias (exact test); 2) baseQ bias (t-test); 3) mapQ bias (t); 4) tail distance bias (t)">
##INFO=<ID=AF1,Number=1,Type=Float,Description="EM estimate of site allele frequency without prior">
##INFO=<ID=AFE,Number=1,Type=Float,Description="Posterior expectation of site allele frequency (with prior)">
##INFO=<ID=HWE,Number=1,Type=Float,Description="P-value for Hardy-Weinberg equillibrium (chi-square test)">
##INFO=<ID=MQ,Number=1,Type=Integer,Descriptin="RMS mapping quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GL,Number=3,Type=Float,Description="Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)">
];
            print $line;
        }
        elsif ( $line=~/^#/ )
        {
            print $line;
        }
        else
        {
            my @items = split(/\t/,$line);
            my @tags = split(/:/,$items[8]);    # FORMAT tags

            my $nidx=2;
            my @idxs;   # Mapping which defines new ordering: $idxs[$inew]=$iold; GT comes first, PL second
            for (my $i=0; $i<@tags; $i++)
            {
                if ( $tags[$i] eq 'GT' ) { $idxs[0]=$i; }
                elsif ( $tags[$i] eq 'PL' ) { $idxs[1]=$i; }
                else { $idxs[$nidx++]=$i; }
            }
            if ( !exists($tags[0]) or !exists($tags[1]) ) { error("FIXME: expected GT and PL in the format field.\n"); }

            # First fix the FORMAT column
            $items[8] = 'GT:GL';
            for (my $i=2; $i<@tags; $i++)
            {
                $items[8] .= ':'.$tags[$idxs[$i]];
            }

            # Now all the genotype columns
            for (my $iitem=9; $iitem<@items; $iitem++)
            {
                @tags = split(/:/,$items[$iitem]);
                $items[$iitem] = $tags[$idxs[0]] .':';

                # GL=-PL/10
                my ($a,$b,$c) = split(/,/,$tags[$idxs[1]]);
                $items[$iitem] .= sprintf "%.2f,%.2f,%.2f",-$a/10.,-$b/10.,-$c/10.;

                for (my $itag=2; $itag<@tags; $itag++)
                {
                    $items[$iitem] .= ':'.$tags[$idxs[$itag]];
                }
            }
            print join("\t",@items);
        }
    }
}

