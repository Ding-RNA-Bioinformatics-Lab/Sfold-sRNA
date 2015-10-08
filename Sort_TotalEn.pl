#!/usr/bin/perl

# Sort_TotalEn.pl
# sort and finalize output format for sRNA protocol
# Copyright (C) 2015 Ding RNA Bioinformatics Lab
# Author William Rennie

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not see <http://www.gnu.org/licenses/>.
#
# email: william.rennie@health.ny.gov

#
# arguments
#     -i - name of input file, output of total energy script
#     -o - name of output file
#     -t - number of nucleitides to consider, starting at 5' end -- -1 means consider the whole sequence
#     -l - binding site length limit -- optional


use warnings;
use strict;

use Cwd qw/ abs_path /; #use absolute directory
use Getopt::Std;

# get the directory in which all scripts are located, including the current one.
 my $script_dir = `dirname $0`;
 chomp $script_dir;
 $script_dir = abs_path($script_dir);

my $DEBUG = 0; # change to one for debugging information

# Globals
sub usage;

# for input arguments

my $inFile;
my $outFile;
my $stopPos;
my $maxLen;

# constants

my $SORTFIELD = 6;
my $STARTFIELD = 2;

print STDERR "\nEntering $0\n" if $DEBUG;
print STDERR "Arguments:\n\t" if $DEBUG;
print STDERR join("\n\t", @ARGV) . "\n" if $DEBUG;

# process arguments

my %opts;

getopts('i:o:t:l:', \%opts);

die( "You must provide an input filename (-i)\n usage()\n") if (! defined $opts{'i'});
$inFile = $opts{'i'};

die( "You must provide an output filename (-o)\n usage()\n") if (! defined $opts{'o'});
$outFile = $opts{'o'};

# optional arguments
if(! defined $opts{'t'} or 0 > $opts{'t'}) {
    $stopPos = -1;
}
else {
    $stopPos = $opts{'t'}
}

if(! defined $opts{'l'} or 0 > $opts{'l'}) {
    $maxLen = -1;
}
else {
    $maxLen = $opts{'l'}
}


# Now we read the input file into an array.  At this point, before sorting, we filter
# out sites not in the region we are interested in.

open(INFILE, $inFile) or die("Could not open $inFile for input\nusage()\n");

my @outArr=();

# discard header
<INFILE>;
while(<INFILE>) {
    chomp;
    my @record = split(":");
    (my $startP, my $endP) = split('-' , $record[$STARTFIELD]);
    if (($stopPos < 0) or ($startP <= $stopPos)) {
        if (( $maxLen < 0) or (($endP - $startP) < $maxLen)) { 
            # process this record
            push @outArr, [@record];
        }
    }
    # else ignore
} # end while lines in file

close INFILE;

my @sortArr = sort {$a->[$SORTFIELD] <=> $b->[$SORTFIELD]} @outArr;

# print sorted filtered array
open(OUTFILE, ">$outFile") or die("Could not open $outFile for output\nusage()\n");

my $header = "sRNA\tTarget\tPosition\tTotal Energy\n";

print OUTFILE $header;

for (my $i = 0; $i <  scalar @sortArr; $i++) {
    print OUTFILE join("\t", @{$sortArr[$i]}[0,1,2,6]);
    print OUTFILE "\n";
} # end for each record

close OUTFILE;

exit 0;

sub usage() 
{
    return "SortTotalEn.pl [-t <truncate sites at> ] -i <input filename> -o <output filename> [-l  <binding site max length>]\n";
}

