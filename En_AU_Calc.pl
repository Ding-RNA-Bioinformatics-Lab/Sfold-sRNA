#!/usr/bin/perl -w

# En_AU_Calc.pl
# Used for the final energy calculation and the AU calculation for the sRNA project
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

# Adapted from Adam's Get_totalEn script from the starmir protocol
# 2014-2-6 (war)
#
use strict;

use Cwd qw/ abs_path /;
use Bio::SeqIO;
use Bio::Seq;

my $DEBUG = 1;

# globals
# for input arguments
my $HybFname;
my $mRNAFname;
my $mRNABsiteFname;
my $mRNADisEnFname;
my $mRNA_Pos_Adj;
my $sRNAFname;
my $sRNABsiteFname;
my $sRNADisEnFname;
my $sRNA_Pos_Adj;
my $OutFname;

print STDERR "\nEntering $0\n" if $DEBUG;
print STDERR "Arguments:\n\t" if $DEBUG;
print STDERR join("\n\t", @ARGV) . "\n" if $DEBUG;

# get the directory in which this script is located
my $base_dir = `dirname $0`;
chomp $base_dir;
$base_dir = abs_path($base_dir);

#-------------------------------------------


if (scalar @ARGV >= 6) {
    $HybFname = $ARGV[0];
    $mRNAFname = $ARGV[1];
    $mRNA_Pos_Adj = $ARGV[2];
    $mRNABsiteFname = $ARGV[3];
    $mRNADisEnFname = $ARGV[4];
    $OutFname = $ARGV[5];
} else {
  die(
		"args: " . (scalar @ARGV) . ": [@ARGV]\n" .
		"usage: >perl En_AU_Calc.pl <Hyb_Nuc_Ener_File> <mRNA seq file> <Pos_Adj> <bsites_File> <disen_File> <output_file>\n"
	); 
}

# load disEn data. Key is position range.
my %disEnDBTarg; # hash for lookup of disruption energy by site.
# my %disEnDBQuery; # hash for lookup of disruption energy by site. # not used in this version

# we need to associate the binding sites in the Bsite file, with the disruption energy
# the files must have the same sites in the same order for this to work.

my $seq_io =  Bio::SeqIO->new(
			     -file => $mRNAFname,
			     -format => 'fasta'
			     );
			     
open(MRNABSITES, "<$mRNABsiteFname") or die("Could not open $mRNABsiteFname for read.");
open(MRNADISENS, "<$mRNADisEnFname") or die("Could not open $mRNADisEnFname for read.");

# key has to be the position of the binding site on the mRNA  Not sure where that gets translated

while (my $bsite = <MRNABSITES>) {
  my $disen = <MRNADISENS>;
  $disen =~ s/^\s+//;
  die("Disen and Bsites files should have the same number of records.") if ($disen eq "");

  my $disruptEn = (split(/\s+/, $disen))[1];
  my ($beg, $end) = split(/\s+/, $bsite);

  $beg -= $mRNA_Pos_Adj;
  $end -= $mRNA_Pos_Adj;
  my $key = sprintf "%08d-%08d", $beg, $end;
  $disEnDBTarg{$key} = $disruptEn;
}
close MRNABSITES;
close MRNADISENS;

# associate disen data with each record. in the Hyb file, produced by associating binding sites with
# RNAhybrid data.

my %recordsBySiteStart;
my $range_ubound = 0;
open(HYBFIL,"<$HybFname") || die "Unable to open the $HybFname file to read";
while (<HYBFIL>) {
    # 0 dhaK  -- transcript name
    # 1 1371  -- transcript length
    # 2 ryhB  -- sRNA name
    # 3 15    -- sRNA length
    # 4 -13.5 -- hybridization energy
    # 5 410-420   -- site position
    # 6 C    U    C C       -- target mismatch
    # 7 GAGC GAUG C         -- target match
    # 8 CUCG UUAC G         -- sRNA match
    # 9 UACA         A C    -- sRNA mismatch
    # 10 4-410,420          -- no clue
    # 11 -0.041             -- nucleation energy

  chomp;
  tr/\r//d;

  my $site = {};
  $site->{'full_record'} = $_;
  my @rec = split(/:/);

  my $range = $rec[5];
  my ($beg, $end) = split(/-/, $range);
  $range_ubound = $end if ($end > $range_ubound);
  $site->{'beg'} = $beg;
  $site->{'end'} = $end;

  $site->{'hybEn'} = $rec[4];
  $site->{'nuclEn'} = $rec[11];

  my $key = sprintf "%08d-%08d", $beg, $end;
  # if (defined $disEnDBTarg{$key} & defined $disEnDBQuery{$key}) {
  if (defined $disEnDBTarg{$key}) {

    $site->{'disEnTarg'} = $disEnDBTarg{$key};
#    $site->{'disEnQuery'} = $disEnDBQuery{$key};
     $site->{'totalEn'} = $site->{'hybEn'} + $site->{'disEnTarg'};   # output from disruptEn program is negated
    # (equivalent to the more logically clear expression (($hybEn) - (-$disEn)) )

    push @{$recordsBySiteStart{$site->{'beg'}}}, $site;

  } else {
    die("Could not find $key in the disEn results.");
  }
}
close HYBFIL;

# for sRNA all sites are output, and we sort by site start.

open(OUT,">$OutFname") || die "Unable to open the $OutFname file to write";

# print header
print OUT "mRNA:sRNA:Site_position:dG_Hybrid:dG_Nucl + dG_Init:dG_Disrupt_T:dG_Total:5p_AU_pos:3p_AU_pos:Site probable\n";

# first fetch and store the sequence with we need for the AU calculation.
my $seq = $seq_io->next_seq();
my $seq_string = $seq->seq();

for my $sitesKey (sort { $a <=> $b } keys %recordsBySiteStart) {
    for my $site (@{$recordsBySiteStart{$sitesKey}}) {

	my @rec = split(/:/, $site->{'full_record'});
	# now for this site.  Locate an AU region if it exists.
	# get the beginning and end of each site for AU calculation
	
	my ($beg, $end) = split(/-/, $rec[5]);
	my $downOffset;
	my $upOffset;
	
	# figure upstream AU
	my $subseq = reverse(substr($seq_string, ($beg - 101), 100));
	if ( $subseq =~ m/([AaUuTt]{4,})/) {
	    my $AUseq = $1;
	    # print 'DownStream:' . $AUseq . "\n";
	    $downOffset = -(index($subseq, $AUseq));
	}
	else {
	    $downOffset = 'none';
	}

	# figure upstream AU
	$subseq = substr($seq_string, ($end + 100), 100);
	if ( $subseq =~ m/([AaUuTt]{4,})/) {
	    my $AUseq = $1;
	    $upOffset = index($subseq, $AUseq);
	}
	else {
	    $upOffset = 'none';
	}
	
	# create output line
	print OUT $rec[0] . ':'; # transcript name
	print OUT $rec[2] . ':'; # sRNA name
	print OUT $rec[5] . ':'; # binding site
	print OUT $rec[4] . ':'; # hybridization energy
	printf OUT "%.3f:", $site->{'nuclEn'} + 4.09; # Nucleation energy + init
	printf OUT "%.3f:%.3f:", ($site->{'disEnTarg'}), $site->{'totalEn'}; # disrupt and total energy
	print OUT $downOffset . ':'; # downstream AU offset, should be negative
	print OUT $upOffset . ':'; # upstream AU offset, should be positive.
	# now we print what I will call the recomendation.  1 if site is likely 0 if not.
 	if ( ((($site->{'nuclEn'}) + 4.09) < 0)
 	     and ($site->{'totalEn'} < 0)
 	     and (($downOffset ne 'none') or ($upOffset ne 'none'))) {
 	    print OUT '1' ;
 	}
 	else {
 	    print OUT '0' ;
 	}
 	print OUT "\n";
#	print  OUT join(':', @rec[6..9]) . ":\n"; # hybrid conformation
    
    } # for each key
} # for each record matching key

close OUT;

print STDERR "Exiting $0\n" if $DEBUG;

exit 0;
