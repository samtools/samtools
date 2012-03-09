#!/usr/bin/perl
use strict;

###Builds interpolated pileup from SAM file
##@description counts bases between paired ends and piles up single end reads.
##@output, uses a #header for the RNAME and then the number of reads per base
##@author sm8@sanger.ac.uk, Stephen B. Montgomery

##@caveats
##Requires RNAME to have format as per example
##      chromosome:NCBI36:18:1:76117153:1
##      supercontig::NT_113883:1:137703:1
##      clone::AC138827.3:1:149397:1
##Expects simple CIGAR characters, M, I and D
##Expects SAM file to be sorted.
##Expects 0x0010 to mark second read in PE file (as has been the observed case from MAQ output) (important for line 77)

##Verify and read in SAM file
my $sam_file = $ARGV[0];
if(!defined($sam_file)) { die("No sam file defined on arg 1"); }
unless(-f $sam_file) { die("Sam file does not exist: $sam_file"); }
open(SAM, $sam_file) || die("Cannot open sam file"); 

##Globals
my $current_location = ""; ##Current RNAME being processed
my $current_size = 0; ##Size of sequence region being processed
my $current_position = 1; ##Current base being processed
my $open = 0; ##Number of open reads (PE reads that have not been closed)
my %close = (); ##Hash of closing positions, when the current_position gets to this position it subtracts the
    ##contained value from those open and deletes the indexed position from the hash

while (my $line = <SAM>) {
    my @tokens = split /\t/, $line;
    
    if ($current_location ne $tokens[2]) { ##Start a new sequence region 
        for (my $i = $current_position; $i <= $current_size; $i++) { ##Close the previous sequence region
            if (defined($close{$i})) {
                $open = $open - $close{$i};
                delete $close{$i};
            }
            print $open . "\n";
        }
        if ($current_location ne "") {
            print "\n";
        }
        
        ##Initiate a new sequence region
        my @location_tokens = split /:/, $tokens[2]; 
        $current_position = 1;
        $current_location = $tokens[2];
        $current_size = $location_tokens[4];
        $open = 0;
        %close = (); 
        print "#" . $tokens[2] . "\n";
        
        ##Print pileup to just before the first read (will be 0)
        for (my $current_position = 1; $current_position < $tokens[3]; $current_position++) {
            print $open . "\n";
        }
        $current_position = $tokens[3];
        
    } else { ##Sequence region already open
        if ($tokens[3] > $current_position) { ##If the new read's position is greater than the current position
                                                ##cycle through to catch up to the current position
            for (my $i = $current_position; $i < $tokens[3]; $i++) {
                if (defined($close{$i})) {
                    $open = $open - $close{$i};
                    delete $close{$i};
                }
                print $open . "\n";
            }
            $current_position = $tokens[3];
        }
    }
    $open++; ##Increment the number of open reads
    
    if (($tokens[1] & 0x0080 || $tokens[1] & 0x0040) && $tokens[1] & 0x0010 && $tokens[1] & 0x0002) { ##if second read of mate pair, add close condition
        $open--;
        my $parsed_cig = &parseCigar($tokens[5]);
        my $seq_region_end = $tokens[3] + $parsed_cig->{'M'} + $parsed_cig->{'D'} - 1;
        if (!defined($close{$seq_region_end + 1})) { $close{$seq_region_end + 1} = 0; }
        $close{$seq_region_end + 1} = $close{$seq_region_end + 1} + 1;
    } elsif (!($tokens[1] & 0x0001) || !($tokens[1] & 0x0002)) { ##if unpaired, add close condition
        my $parsed_cig = &parseCigar($tokens[5]);
        my $seq_region_end = $tokens[3] + $parsed_cig->{'M'} + $parsed_cig->{'D'} - 1;
        if (!defined($close{$seq_region_end + 1})) { $close{$seq_region_end + 1} = 0; }
        $close{$seq_region_end + 1} = $close{$seq_region_end + 1} + 1;
    } else {
        #do nothing
    }
}
for (my $i = $current_position; $i <= $current_size; $i++) {  ##Finish up the last sequence region
    if (defined($close{$i})) {
        $open = $open - $close{$i};
        delete $close{$i};
    }
    print $open . "\n";
}
print "\n";
close(SAM);
exit(0);

##reads and tokenizes simple cigarline
sub parseCigar() {
    my $cigar_line = shift;
    $cigar_line =~ s/([0-9]*[A-Z]{1})/$1\t/g;
    my @cigar_tokens = split /\t/, $cigar_line;
    my %parsed = ('M' => 0,
                  'I' => 0,
                  'D' => 0);
    my @events = ();
    for(my $i = 0; $i < scalar(@cigar_tokens); $i++) {
        if ($cigar_tokens[$i] =~ /([0-9]+)([A-Z]{1})/g) {
            if (!defined($parsed{$2})) { $parsed{$2} = 0; }
            my $nt = $2;
            if ($nt ne "M" && $nt ne "D"  && $nt ne "I") { $nt = "M"; }
            $parsed{$nt} += $1;
            my %event_el = ("t" => $nt,
                            "n" => $1);
            push @events, \%event_el;
        }
    }
    $parsed{'events'} = \@events;
    return \%parsed;
}
