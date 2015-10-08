#!/usr/bin/perl 

# sRNA_protocol.pl
# used for running current sRNA protocol.
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
# arguments
#     sRNA file in FASTA format
#     mRNA file in FASTA format
#     output filename
#     sfold binary directory, for distruption energy calculation
#     sfold output directory, containing output of Sfold run 
#     temp directory -- optional

# 2015-20-7 (war)
# removed depedency on parameter file
# assumption is RNAhybrid is installed to be globally
# available, and sfold has been run when the protocol is launched.

use warnings;
use strict;

use Cwd qw/ abs_path /;  # one way of getting the absolute path of a file.
my $DEBUG = 0;

# get the current directory.
# assumption is that all required scripts and executables are in the 
# same directory as the launch script

my $script_dir = `dirname $0`;
chomp $script_dir;
$script_dir = abs_path($script_dir);

print STDERR "\nEntering $0\n" if $DEBUG;
print STDERR "Arguments:\n\t" if $DEBUG;
print STDERR join("\n\t", @ARGV) . "\n" if $DEBUG;

die(" Usage: sRNA_protocol.pl <sRNA_file> <mRNA_file> <output filename> <Sfold_bin_dir> <Sfold_outdir> [<temp directory>]\n") unless scalar @ARGV >= 4;

# test the existence of various necessary files.

my ($queryFname, $targetFname, $outFName, $SFold_bin_dir, $SFold_outdir) = @ARGV[0..4];

my $tmpDir= '.'; # if no temporary is defined put intermediate files in current directory

if (defined $ARGV[5]) {
    # use a temporary directory for intermediate files
    $tmpDir = $ARGV[5];
    if (! -e $tmpDir) {
	die("directory $tmpDir not exist!"); # created by launch script
    }
}

-e $SFold_outdir or die("sfold output directory does not exist!");

if (! -e $queryFname) {
    die("File $queryFname does not exist!");
}

if (! -e $targetFname) {
    die("File $targetFname does not exist!");
}


# All required elements are present

# 1: run RNAhybrid
#    constant cutoff of -10 in current incarnation

my $HybFname = "$tmpDir/Hyb-tmp.out"; # output name

# 1a: Run RNAhybrid
#     input -- target sequence filename
#              query sequence filename
system("perl run.RNAhybrid.pl $targetFname $queryFname $HybFname 3utr_human");

# 2: run ParseHybFile.pl
#     input -- $tmpDir/Hyb-tmp.out, output of run.RNAHybrid
my $ResHybFname = "$tmpDir/ResHyb-tmp.out"; # output name
system("perl $script_dir/ParseHybFileSRNA.pl $HybFname $ResHybFname");

# 3: write out the possible binding sites
#    input -- $tmpDir/ResHyb-tmp.out, ParseHybFile.pl output
#             position adjustmentis always zero for this protocol
#    output $tmp_dir/tmp-bsites.txt
my $BsitesFname = "$tmpDir/tmp-bsites.txt"; # output file
system("perl $script_dir/Bsites.pl $ResHybFname 0 $BsitesFname");

# 4: calculate disruption energy
#    input - fe.out, bp.out, from SFold result
#    output - $tmpDir/disruptEn_tmp.out
# Note: $targetfname is the file containing the sequence we gave to Sfold.
my $disruptEnFname = "$tmpDir/disruptEn-tmp.out" ;  #output filename
print $SFold_outdir . "\n";
system("perl $script_dir/Calc_disruptEn_sRNA.pl $SFold_bin_dir $targetFname $BsitesFname $disruptEnFname $SFold_outdir");

# 5: calculate nucleation energies
#    input -- $tmpDir/ResHyb-tmp.out, ParseHybFile.pl output
my $HybFilFname = "$tmpDir/Hyb-Fil-tmp.out"; # output name

system("perl $script_dir/CalcNE.4.2.pl $ResHybFname 0 $HybFilFname $SFold_outdir");


# 6: assemble the data to calculate the total interaction energy
#    input -- $tmp_dir/$mRNA_name-bsites.txt, from Bsites.pl
#             $tmp_dir/disruptEn_$mRNA_name.out, from Calc_disruptEn.pl
my $TotalEnFname = "$tmpDir/EnAU-tmp.out"; # output filename

# position adjustment is always zero in this incarnation
system("perl $script_dir/En_AU_Calc.pl $HybFilFname $targetFname 0 $BsitesFname $disruptEnFname $TotalEnFname 1");


exit 0;
