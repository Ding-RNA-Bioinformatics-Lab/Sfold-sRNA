# sRNA protocol

This project contains scripts that can be used for predicting and ranking binding sites for an sRNA on a target mRNA.

Requirements
The Sfold executables required for this system of scripts can be obtained from http://sfold.wadsworth.org
They are free for academic and non commercial use.  Comercial use requires a license.
Currently Sfold only runs on Linux and Solaris platforms.  The executable used for calculating disruption energy is
part of the Sfold distribution.

RNAhybrid must be installed and in the execution path. RNAhybrid may be obtained from, 
http://bibiserv.techfak.uni-bielefeld.de/rnahybrid/.

The Bioperl library must be installed and available. Instructions for installing Bioperl can be found at http://bioperl.org

Notes:
Sfold is very memory and computationally expensive.  Sequences over 1,000 nts in length can take hours to run.  Sequences
larger than 5,000 nts can take days to run, depending on the resources available.

The Sfold disruption energy executable is very memory intensive.  Anything more than a trivial sequence, (more than a couple of 
hundred nts), requires at least 16 GB of memory.  Most of our work is done on machines with 64 GB of memory.

Running
The program is run by running the "main.pl" script, which calls all the remaining scripts.  The output is a file ranking
binding sites by energy, and a second file prefixed by "plot" that contains base pair diagrams for the same set of binding sites.

All perl scripts should be in the same directory as the main.pl script.  

The main script takes the following options, prefixed by the following flags
  --srna - The sRNA sequence in FASTA format
  --target - The target mRNA sequence in FASTA format
  --output - The name for the output file
  --sfold - The absolute path to the sfold executable
  --tmpPath - The name of the directory for temporary file, this can be created in advance of the run
                 if it is missing the directory will be created.  If this option is not specified
                 the directory created will be a subdirectory of the current directory called "tmp".
  --retain - do not delete the temporary directories and files.  If this option is omitted, the temporary 
                 directories and files will be deleted.
  --truncate - do not look for binding sites after this position on the mRNA.  This was used in our research to limit
                 the binding sites we considered to the 5'UTR of the mRNA.
  --maxLen - do not consider binding sites longer than this length.  We found RNAhybrid will predict some extremely long
                 and very unlikely binding sites.  This allows you to ignore binding sites longer than the argument.  The 
                 default is 70, if this option is not used.
 
An example command line is;
perl main.pl --srna fred.fa --target wilma.fa --output FlintSites.txt --sfold /home/me/development/Sfold/bin/sfold --tmpPath FlintTmp --retain --truncate 172 --maxLen 60

The sRNA sequence is in fred.fa
The target sequence is in wilma.fa
The output will be produced in FlintSites.txt, the basepair diagrams in plot_FlintSites.txt
The path to sfold is /home/me/development/Sfold/bin/sfold.  Its best to use an absolute path here.
Temporary files will be created in a subdirectory of this directory "FlintTmp"
The temporary files will be retained.
We are only looking for binding sites completely contained in the region from the 5' end to nt 172
We are not considering binding sites longer than 60 nts
