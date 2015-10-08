#!/usr/bin/perl -w

# Calc_disruptEn_sRNA.pl
# Run Sfold disen program to calculate disruption Energy 
# Copyright (C) 2015 Ding RNA Bioinformatics Lab
# Authors William Rennie, Adam Wolnec

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
# 2014-2-6 (war)
#

# usage > perl Calc_disruptEn_sRNA.pl <mRNA_fname> <bsite_file> <output_file> <SFold_outdir>


use warnings;
use strict;

use Cwd qw/ abs_path /;

my $folding_seq_path; 
my $bposf;  
my $outfname;
my $SFold_outdir;  
my $SFold_bin_dir;

my $DEBUG = 0; # set to one for debugging information

print STDERR "\nEntering $0\n" if $DEBUG;
print STDERR "Arguments:\n\t" if $DEBUG;
print STDERR join("\n\t", @ARGV) . "\n" if $DEBUG;

# get the directory in which this script is located
my($base_dir) = `dirname $0`;
chomp $base_dir;
$base_dir = abs_path($base_dir);


if ($#ARGV == 4) {
    $SFold_bin_dir = $ARGV[0];
    $folding_seq_path = $ARGV[1] ; 
    $bposf = $ARGV[2] ;  
    $outfname = $ARGV[3] ;
    $SFold_outdir  = $ARGV[4] ;  
} else {
    die("usage: > perl Calc_disruptEn_sRNA.pl <Sfold bin directory> <RNA_fname> <bsite_file> <output_file> <SFold_outdir>\n") ;
}

if (-e "$SFold_outdir/output") {$SFold_outdir = "$SFold_outdir/output";}

#check the existence of input files for disruptEn_prog.
die " Error: unable to locate $bposf" if (!-e $bposf) ;

my $structf = "$SFold_outdir/bp.out" ;
die " Error: unable to locate $structf" if (!-e $structf) ;

my $febf = "$SFold_outdir/fe.out" ;
die " Error: unable to locate $febf" if (!-e $febf) ;

die " Error: unable to locate $folding_seq_path" if (!-e $folding_seq_path) ;

my $disruptEn_prog = "$SFold_bin_dir/disruptEn/disruptEn";

system("SFOLDBIN=$SFold_bin_dir/bin/ $disruptEn_prog -m $folding_seq_path -p $bposf -s $structf -c 0 -f $febf -t 1 > $outfname") ;

print STDERR "Exiting $0\n" if $DEBUG;

exit 0;

