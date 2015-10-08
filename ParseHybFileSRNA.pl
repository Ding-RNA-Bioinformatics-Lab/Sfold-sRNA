#!/usr/bin/perl -w

# ParseHybFileSRNA.pl 
# This script parses RNAhybrid output and produces an encoding of the
# component helixes.
# Copyright (C) 2015 Ding RNA Bioinformatics Lab
# Author William Rennie, Adam Wolnec

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

# usage  >perl ParseHybFileSRNA.pl <Input_File> <Output_filename> 


use strict;
use warnings;

use File::Basename;

my $DEBUG = 0;

print STDERR "\nEntering $0\n" if $DEBUG;
print STDERR "Arguments:\n\t" if $DEBUG;
print STDERR join("\n\t", @ARGV) . "\n" if $DEBUG;

# Globals
# for arguments
my $InFname;
my $OutFname;

# parse arguments
if ($#ARGV == 1) {
    $InFname = $ARGV[0] ;
    $OutFname = $ARGV[1] ;
} else {
    my $bsname = basename($0);
   die("usage: >perl $bsname <Input_File> <Output_file>\n"); 
}

open(IN,"<$InFname") || die "Unable to open the $InFname file to read" ;
open(OUT,">$OutFname") || die "Unable to open the $OutFname file to write" ;

# processing the input file.  We need to isolate helices
# the output files will be used to calculate disruption energy
# code in this section by Adam Wolenc
while (<IN>) {
  chomp;
  my @rec = split(/:/);
  my $sbj_name = $rec[0];
  my $sbj_len = $rec[1];
  if ($sbj_name =~ /\|/) {
    $sbj_name = (split(/\|/, $sbj_name))[1];
  }
  my $q_name = $rec[2];
  my $q_len = $rec[3];
  my $hyb_ener = $rec[4];
  my $rnahyb_be_tar = $rec[6];

  my $tar_misses = uc($rec[7]);
  my $tar_matches = uc($rec[8]);
  my $mir_matches = uc($rec[9]);
  # (mir misses are not relevant to this script.)
  my $pairing = join(":",$rec[7],$rec[8],$rec[9],$rec[10]);

  my @tar_misses_c = split(//, $tar_misses);
  my @tar_matches_c = split(//, $tar_matches);
  my @mir_matches_c = split(//, $mir_matches);

  my $pos_on_tar = $rnahyb_be_tar - 1;

  my $cur_helix;  # this is updated as we accumulate a longer unbroken helix
  my @helices;  # all helixes of the current diagram
  my $be_tar  = 0;  # undef
  my $end_tar = 0;  # undef
  for (my $p = 0; $p < scalar @tar_misses_c; ++$p) {
    my $tar_mis_ch   = $tar_misses_c[$p];
    my $tar_match_ch = $tar_matches_c[$p];
    my $mir_match_ch = $mir_matches_c[$p];
    #there are 3 possibilities for each position in the diagram
    #tar_miss
    # mir_match
    #  mir_miss
    #010 = tar position increments, begin or continue a helix
    #10? = tar position increments, end a helix if one is underway.
    #001 = end a helix if one is underway.
    if ($tar_match_ch ne ' ') {
      # case 010
      ++$pos_on_tar;
      #maintain site begin and end.
      $be_tar = $pos_on_tar if ($be_tar == 0);
      $end_tar = $pos_on_tar;
      #update current helix or start a new one.
      if (!$cur_helix->{'len'}) {
        $cur_helix->{'tar_start'} = $pos_on_tar;
        $cur_helix->{'len'} = 0;  # about to be incremented & appended.
        $cur_helix->{'mir_seq'} = "";
        $cur_helix->{'tar_seq'} = "";
      }
      $cur_helix->{'mir_seq'} .= $mir_match_ch;
      $cur_helix->{'tar_seq'} .= $tar_match_ch;
      $cur_helix->{'tar_end'} = $pos_on_tar;
      ++$cur_helix->{'len'};
    } else {
      # case 10? or 001
      # end a helix if one is underway
      if ($cur_helix->{'len'}) {
        push @helices, $cur_helix;
        undef $cur_helix;
      }
      #update positions only in case of 10?
      ++$pos_on_tar if ($tar_mis_ch ne ' ');
    }
  }
  # one last output in case the diagram ends on a helix.
  if ($cur_helix->{'len'}) {
    push @helices, $cur_helix;
    undef $cur_helix;
  }

  #print
  print OUT $sbj_name,":",$sbj_len,":",$q_name,":",$q_len,":",$hyb_ener,":";
  print OUT $be_tar,"-",$end_tar,":";
  for my $helix (@helices) {
    print OUT $helix->{'tar_seq'},",",$helix->{'mir_seq'},"_";
    print OUT $helix->{'tar_start'},",",$helix->{'tar_end'},":";
  }
  #2011-02-02 Adam: Different between this script and prediction steps, append diagram to output.
  print OUT $pairing,"\n" ;
}
close IN;
close OUT;

print STDERR "Exiting $0\n" if $DEBUG;

exit 0;
