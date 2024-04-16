# flank_grabber

Author: Murat Buyukyoruk

        flank_grabber help:

This script is developed to fetch flank sequences of gene/CRISPR loci of interest by using the fasta file with provided position and strand information. 

SeqIO and Seq packages from Bio is required to fetch flank sequences. Additionally, tqdm is required to provide a progress bar since some multifasta files can contain long and many sequences.

Syntax:

        python flank_grabber.py -i demo.fasta -o demo_gene_flanks.fasta -d flank_info_dataframe

        OR

        python flank_grabber.py -i demo.fasta -o demo_gene_flanks.fasta -d demo_flank_info_dataframe.txt -f 200 -c Y -p Y

Example Dataframe (tab separated excel file is required):

        Accession   Protein_accession/array_no  Start   Stop    Strand
        NZ_CP006019 NZ_CP006019_1756            1875203 1877050 -1
        CP000472.1  CP000472.1_235              123     975     1

flank_grabber dependencies:

Bio module, SeqIO and Seq available in this package     refer to https://biopython.org/wiki/Download

tqdm                                                    refer to https://pypi.org/project/tqdm/

Input Paramaters (REQUIRED):
----------------------------
	-i/--input		FASTA			Specify a fasta file. FASTA file requires headers starting with accession number. (i.e. >NZ_CP006019 [fullname])

	-o/--output		Output file	    Specify a output file name that should contain fetched sequences.

	-d/--data		Dataframe		Specify a list of accession (Accession only). Each accession should be included in a new line (i.e. generated with Excel spreadsheet). Script works with or without '>' symbol before the accession.

Parameters [optional]:
----------------------
	-f/--flank		200			This is the default length of flanks that is fetched.

	-c/--circular	N			This is the default option that is assuming the sequence is not circular. Type "Y" insted if you know it is cirgulat genome or plasmid.

	-p/--promoter	N			Option "N" allows to get both flanks from each end of the gene. Type "Y" instead if you just want to fetch promoter region.

Basic Options:
--------------
	-h/--help	HELP			Shows this help text and exits the run.

Output header will contain protein accession/array no, original accession number, positions and fullname (if included in original fasta), observed stand and the 
information about flank in terms of which side of the gene is fetched.

