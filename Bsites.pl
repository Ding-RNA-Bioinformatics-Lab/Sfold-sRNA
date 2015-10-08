#!/usr/bin/perl -w

# Bsites.pl
# outputs listing of binding sites from processed RNAhybrid data.
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

# 2015-7-20 (war) Adapted for sRNA protocol
# 2014-2-6 (war)
# 10/18/2011 (war) Added standard tracking
# 2010-12-30 Adam Wolenc for Wadsworth Center.

use strict;
use warnings;

use File::Basename;

my $DEBUG = 0;

print STDERR "\nEntering $0\n" if $DEBUG;
print STDERR "Arguments:\n\t" if $DEBUG;
print STDERR join("\n\t", @ARGV) . "\n" if $DEBUG;

my $ResHybFname;
my $CDlen; 
my $BsiteOutFname;

if (scalar @ARGV == 3) {
  $ResHybFname = $ARGV[0];
  $CDlen = $ARGV[1];
  $BsiteOutFname = $ARGV[2];
} else {
  die("usage: ./Bsites.pl <ResHyb_File> <length_offset> <Output_file>\n");
}

open(IN,"<$ResHybFname") || die "Unable to open the $ResHybFname file to read";
my %site;
while (my $line = <IN>) {
  chomp $line;
  my @rec = split(/:/, $line);
  my $genename = $rec[0];
  $genename =~ s/-3pUTR//;
  my $pos = $rec[5];
  my ($st, $en) = split(/-/, $pos);
  $st += $CDlen;
  $en += $CDlen;
  my $key = sprintf "%08d-%08d", $st, $en;
  $site{$key} = "$st\t$en\n";  # prevent duplicates.
}
close IN;

#output in order with no duplicates.
open(OUT,">$BsiteOutFname") || die "Unable to open the $BsiteOutFname file to write";
for my $key (sort keys %site) {
  print OUT $site{$key};
}
close OUT;

print STDERR "Exiting $0\n" if $DEBUG;

exit 0;
