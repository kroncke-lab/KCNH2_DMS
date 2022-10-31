# Deep Mutational Scan of KCNH2

Here we develop an experimental method to quantitate variant-induced mistrafficking of missense variants in KCNH2. 
All code and data used and referenced in the manuscript (PMID: 35688148 PMCID: PMC9300756 DOI: 10.1016/j.ajhg.2022.05.003) are included here.

Usage and Summary:

This suit of python scripts are useful for deep mutational scanning experiments, specifically, those which use a barcode linked to a mutation in the protein 
coding sequence. 
First, the script identifies the barcodes from the consensus of the input reads. It then identifies the mutation from the reads by comparing it to the wild-type
sequence and then links the barcode to the mutation to create a barcode-mutation index. 
The script then proceeds to the .fastq files from sort-seq experiments, containing the barcodes from enriched pools of cells. The script quantitates barcode reads in 
each pool to calculate levels of enrichment. 
This version is stringent on filtering reads based on the base pairs flanking the barcode region and it is imperative all barcodes have identical flanking regions.

Set Up and Install: 

First, install any dependencies. This script requires python3 and a virtual environment with the following packages:

$ # if necessary, update or call python3

$ module load python3 

$ virtualenv /path/we/want/

$ source /path/we/want/bin/activate

$ pip install matplotlib==3.4.2

$ pip install biopython==1.77

$ #############the biopython must be version 1.77 to allow for use of the Alphabet.IUPAC module!!! ########

$ ##pandas and numpy are not critical which version

$ pip install pandas

$ pip install numpy

Create a directory that will hold the .fastq reads and the reference sequence. The reference sequence must be formatted as a fasta and include two entries:
the first entry will be the dna sequence of the entire plasmid. The second entry in this fasta file must be the dna only of the protein coding sequence.
The only files in this directory should be the fastq files and reference fasta file, which should be named plasmid.fasta 
The reference fasta file should read like so:

$ >dna_of_full_plasmid

$ atcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcg

$ >dna_of_protein_coding_region

$ atcgatcgatcg

modify the first lines of the main.py script to search in the correct directory 

$ folder = "/Users/directory/with/files/"

This script can be run in the custom virtual environment by calling

$ /path/we/want/bin/python3 main.py

NOTE: If you only wish to analyze the barcode pairings with the variants and generate an index of barcode-variant pairs, comment out the second half of main.py after line 15. 
The default setting is barcode = 18 base pairs. If needed, this can be changed manually by changing 18 to your custom barcode length at consensus.py - line 65, line 81, line 105, and line 158. 



Further note for users new to python: directory structure is important when working with the python virtual environment. Following directory structure recommended: 

$ cd /working/directory/

$ move the entire KCNH2_DMS folder into this directory

create the dna files directory

$ mkdir ./DNA/

modify main.py to point to the /working/directory/DNA/

$ virtualenv /working/directory/KCNH2_DMS

$ source ./KCNH2_DMS/bin/activate

$ ./KCNH2_DMS/bin/python3 -m pip install install packages

$ ./KCNH2_DMS/bin/python3 ./KCNH2_DMS/main.py

This directory structure should allow for correct recognition of all dependencies and subscripts needed in the virtual environment

