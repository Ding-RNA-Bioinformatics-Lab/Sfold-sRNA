#!/usr/bin/perl

# main.pl

# main dispatch script for launching programs used in SRNA protocol
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
# --srna - srna sequence file in FASTA format
# --target - mRNA sequence file in FASTA format
# --output - output file name
# --sfold - path to sfold executables including disen
# --tmpPath - path to directory for temporary files
# --retain - do not delete temporary files (to avoid repeating Sfold runs)
# --truncate - do not look for binding sites beyond this position
# --maxLen -- maximum length for binding site
#
#
use warnings;
use strict;

use Cwd qw/ abs_path /; #use absolute directory
use Getopt::Long;

# get the directory in which all scripts are located, including the current one.
my $script_dir = `dirname $0`;
chomp $script_dir;
$script_dir = abs_path($script_dir);

# Globals
sub usage;

# for input arguments

my $srnaFile;
my $targetFile;
my $outFile;
my $sfoldPath = '';
my $stopPos = -1;
my $maxLen = -1;
my $tmpPath = '';
my $retain = '';


# process input arguments

GetOptions(
    "srna=s"     => \$srnaFile,
    "target=s"   => \$targetFile,
    "output=s"   => \$outFile,
    "truncate=i" => \$stopPos,
    "length=i"   => \$maxLen,
    "sfold=s"    => \$sfoldPath,
    "temp=s"     => \$tmpPath,
    "retain!"    => \$retain
    ) or die( usage() . "\n");

# check for existence of input files
if (! -e  $srnaFile) {
    die("You must give an sRNA sequence file in FASTA format.\n usage()");
}

if (! -e  $targetFile) {
    die("You must give an target sequence file in FASTA format.\n usage()");
}

if ($sfoldPath eq '') {
    $sfoldPath = "./";
}

if ($tmpPath eq '') {
    $tmpPath = 'tmp';

}
# if the temp directory does not exist create it
if (! -e $tmpPath) {
    system("mkdir $tmpPath")
}

-e $tmpPath or die("Could not create temp directory $tmpPath\n");

# confirm or create directory for SFold output
my $SFold_outdir = "$tmpPath/sfoldOutput";

# if the output directory does not exist create it
if (! -e $SFold_outdir) {
    system("mkdir $SFold_outdir")
}

-e $SFold_outdir or die("Could not create sfold output directory $SFold_outdir\n");

$outFile = abs_path($outFile);

# my $Pos_Adj = 0;

# check twice to determine if sfold is successful or not.
# If SFold has already run, the script does not rerun it.
# following code contributed to by Chaochun Liu

if (!(-e "$SFold_outdir/clusters/ch.index.out" 
      && -e "$SFold_outdir/bp.out" 
      && -e "$SFold_outdir/fe.out" 
      && -e "$SFold_outdir/sstrand.out")) {
    # it seems likely we have not previously run sfold
  system("$sfoldPath/bin/sfold -o $SFold_outdir -a 1 $targetFile");
}

if (!(-e "$SFold_outdir/clusters/ch.index.out" 
      && -e "$SFold_outdir/bp.out" 
      && -e "$SFold_outdir/fe.out" 
      && -e "$SFold_outdir/sstrand.out")) {
    # to get here first run failed
  system("$sfoldPath/bin/sfold -o $SFold_outdir -a 1 $targetFile");
}

if (!(-e "$SFold_outdir/clusters/ch.index.out" 
      && -e "$SFold_outdir/bp.out" 
      && -e "$SFold_outdir/fe.out" 
      && -e "$SFold_outdir/sstrand.out")) {
    # to get here both attempts to run sfold have failed
    die("Could not run SFold. Terminating.\n");
}

print STDERR "SFold output found.\n";


# run SRNA protocol
system("perl $script_dir/sRNA_protocol.pl $srnaFile $targetFile $outFile $sfoldPath $SFold_outdir $tmpPath");

# Now sort the output
my $TotalEnFname = "$tmpPath/EnAU-tmp.out";

system("perl $script_dir/Sort_TotalEn.pl -t $stopPos -l $maxLen -i $TotalEnFname -o $outFile");

# create the basepairing file
system("perl $script_dir/format_sitelisting.pl $outFile $tmpPath/Hyb-Fil-tmp.out");

if (! $retain) {
    system("rm -r $tmpPath");
}

print STDERR "Done!\n";
exit 0;


sub usage() 
{
    return "main.pl --srna <sRNA filename> --target <target filename> --output <output filename> --sfold <path to sfold executable> --truncate <end of region> [ --maxLen <max binding site length>] [--tmpPath <path directory for temporary files]";
}
