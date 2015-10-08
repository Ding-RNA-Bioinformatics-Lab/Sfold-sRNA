#!/usr/bin/perl -w

# format_sitelisting.pl
# Creates base pairing diagrams
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

use warnings;
use strict;

use Cwd qw/ abs_path /;
use File::Basename;

my $DEBUG = 0; # change to one for debugging information

# get the directory in which this script is located
$0=~/^(.+[\\\/])[^\\\/]+[\\\/]*$/;
my $base_dir = abs_path($1);

# This script uses Hyb-Fil-tmp.out 
# for base pairing information  and produces s text output file
# containing the structural conformation for the selected binding sites

die(" Usage: $0 <sRNA outfile> <hyb fil file> \n") if $#ARGV != 1;

print STDERR "\nEntering $0\n" if $DEBUG;
print STDERR "Arguments:\n\t" if $DEBUG;
print STDERR join("\n\t", @ARGV) . "\n" if $DEBUG;

my $sRNAFile = $ARGV[0];
my $hybFile = $ARGV[1];
my $plotFile = dirname($sRNAFile) . "/" . "plot-" . basename($sRNAFile);

# open files
open(SRNAFILE, $sRNAFile) or die("Can not open $sRNAFile for input\n");
open(HYBFILE, $hybFile) or die("Can not open $hybFile for input\n");
open(PLOTFILE, ">$plotFile") or die("Can not open $plotFile for output\n");

my $target_name = "";
my $sRNA_name = "";

my @dG_total = ();
my @dG_disrupt = ();
my @dG_hybrid = ();
my @dG_N = ();
my @sRNA_len = ();

my @tar_startpos = ();
my @tar_endpos = ();

my @tarconf1 = ();
my @tarconf2 = ();
my @tarconf3 = ();
my @tarconf4 = ();

my @srnaid = ();
my $nsites = 0;

# Create structure information hash
my %structHash = ();

while(<HYBFILE>) {
    chomp;
    my @record = split ':';
    $structHash{($record[0] . '#' . $record[2] . '#' . $record[5])} = [@record];
}
close HYBFILE;

# for each of the selected binding sites, look up the structure information and plot it.
# throw away header
<SRNAFILE>;

while (<SRNAFILE>) {
  chomp;
  next if /^\s*$/;

  my(@temp) = split("\t");

  $target_name = $temp[0];
  $sRNA_name = $temp[1];
  my $targPos =  $temp[2];

  my $e = $structHash{($target_name . '#' . $sRNA_name . '#' .$targPos)};
  die ( "Structure not found\n" ) if not $e;

  # push @dG_total, $e->[0];
  # push @dG_disrupt, $e->[1];
  push @dG_hybrid, $e->[4];
  # push @dG_N, $e->[9];

  push @sRNA_len, $e->[3];

  # break out position
  $e->[5] =~ /(\d+)-(\d+)/;

  push @tar_startpos, $1;
  push @tar_endpos, $2;
  push @tarconf1, $e->[6];
  push @tarconf2, $e->[7];
  push @tarconf3, $e->[8];
  push @tarconf4, $e->[9];
  push @srnaid, $e->[2];

  $nsites++;
} # end while
close SRNAFILE;

die("No target site was identified based on the current energy thresholds.\n") if !$nsites;

print PLOTFILE
	"Target name   = $target_name\n",
	"\n";
print PLOTFILE
	"sRNA name   = $sRNA_name\n",
	"\n";

for (my $i=1; $i<=$nsites; $i++) {
  my $bp = $tarconf2[$i-1];
  $bp =~ s/[^\s]/|/g;
  my $my_spos = sprintf("%5d", $tar_startpos[$i-1]);
  my $my_len = sprintf("%5d", $sRNA_len[$i-1]); # mod here for sRNA start position

  print PLOTFILE <<EndOfOutput;
Site $i -- 

    5'->3'        $tarconf1[$i-1]
    Target $my_spos  $tarconf2[$i-1]  $tar_endpos[$i-1]
                  $bp
    sRNA   $my_len  $tarconf3[$i-1]  1
    3'->5'        $tarconf4[$i-1]
 

EndOfOutput

} # end for each site

close PLOTFILE;

exit 0;
