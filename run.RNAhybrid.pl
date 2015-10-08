#!/usr/bin/perl -w

# run.RNAhybrid.pl
# Run RNAhybrid to identify candidate sites.
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
#

# usage > perl run.RNAhybrid.pl <target file_name> <query file name> <Out_file_name> <species>
# -----------------------------------------------------

use strict;
use warnings;

my $DEBUG = 0; # set to one when debugging.

print STDERR "\nEntering $0\n" if $DEBUG;
print STDERR "Arguments:\n\t" if $DEBUG;
print STDERR join("\n\t", @ARGV) . "\n" if $DEBUG;

# Globals

my $HYB_TH = -10; 

# for arguments
my $DBFname; 
my $MFname; 
my $OutFname; 
my $Hyb_species = 'human'; # this is the default

if ($#ARGV == 3 or $#ARGV == 4) {
  $DBFname =  $ARGV[0];
  $MFname = $ARGV[1];
  $OutFname = $ARGV[2];
  $Hyb_species = $ARGV[3] if defined $ARGV[3]; # we are assuming value is valid
} 
else {
  die("usage: > perl run.RNAhybrid.pl <target file> <query file> <Out_file_name> <Training species>\n");
}


system("RNAhybrid -c -n 200 -m 5000 -s 3utr_human        -e $HYB_TH -t $DBFname -q $MFname > $OutFname");

print STDERR "Exiting $0\n" if $DEBUG;

exit 0;
