#!/Users/muratbuyukyorukmsu/PycharmProjects/python/venv/bin/python

import argparse
import os
import sys
import subprocess
import re
import textwrap

try:
    from Bio import SeqIO
except:
    print("SeqIO module is not installed! Please install SeqIO and try again.")
    sys.exit()

try:
    from Bio.Seq import Seq
except:
    print("Seq module is not installed! Please install Seq and try again.")
    sys.exit()

try:
    import tqdm
except:
    print("tqdm module is not installed! Please install tqdm and try again.")
    sys.exit()

parser = argparse.ArgumentParser(prog='python flank_grabber.py',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=textwrap.dedent('''\

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
	-f/--flank		200			    This is the default length of flanks that is fetched.

	-c/--circular	N			    This is the default option that is assuming the sequence is not circular. Type "Y" insted if you know it is cirgulat genome or plasmid.

	-p/--promoter	N			    Option "N" allows to get both flanks from each end of the gene. Type "Y" instead if you just want to fetch promoter region.

Basic Options:
--------------
	-h/--help		HELP			Shows this help text and exits the run.

Output header will contain protein accession/array no, original accession number, positions and fullname (if included in original fasta), observed stand and the information about flank in terms of which side of the gene is fetched.

      	'''))
parser.add_argument('-i', '--input', required=True, type=str, dest='filename',
                    help='Specify a fastafile to fetch regions from.\n')
parser.add_argument('-o', '--output', required=True, dest='out',
                    help='Specify a output file name.\n')
parser.add_argument('-d', '--data', required=True, dest='data',
                    help='Specify a dataframe with accession, prot_accession/array_no, start, stop, strand info in that order. Leave NA if not known.\n')
parser.add_argument('-f', '--flank', type=int, required=False, default=200, dest='flank',
                    help='Specify length of flanking sequence to fetch.\n')
parser.add_argument('-c', '--circular', type=str, required=False, default="N", dest='circular',
                    help='Circular genome? Y/N (Default: N). Split circluar and linear sequences if you have a mixture.\n')
parser.add_argument('-p', '--promoter', type=str, required=False, default="N", dest='promoter',
                    help='Get promoter only genome? Y/N (Default: N).\n')

results = parser.parse_args()
filename = results.filename
out = results.out
data = results.data
flank = results.flank
circ = results.circular
promoter = results.promoter

if circ in ("Y", "N"):
    pass
else:
    print("""\nWrong argument! Type "Y" for Yes(cicular) or "N" for No (linear).\n""")

    os.system('ORF_flank_fetch.py -h')
    sys.exit()

if promoter in ("Y", "N"):
    pass
else:
    print("""\nWrong argument! Type "Y" for Yes(cicular) or "N" for No (linear).\n""")

    os.system('ORF_flank_fetch.py -h')
    sys.exit()

os.system('> ' + out)

seq_id_list = []
seq_list = []
seq_description_list = []

proc = subprocess.Popen("grep -c '>' " + filename, shell=True, stdout=subprocess.PIPE, text=True)
length = int(proc.communicate()[0].split('\n')[0])

with tqdm.tqdm(range(length)) as pbar:
    pbar.set_description('Reading...')
    for record in SeqIO.parse(filename, "fasta"):
        pbar.update()
        seq_id_list.append(record.id)
        seq_list.append(record.seq)
        seq_description_list.append(record.description)

proc = subprocess.Popen("wc -l < " + data, shell=True, stdout=subprocess.PIPE, text=True)
length = int(proc.communicate()[0].split('\n')[0])

with tqdm.tqdm(range(length), desc='Fetching...') as pbar:
    with open(data, 'r') as file:
        f = open(out, 'a')
        sys.stdout = f
        for line in file:
            pbar.update()
            if "Accession" not in line:
                acc = line.split('\t')[0]
                prot_acc = line.split('\t')[1]
                start = line.split('\t')[2]
                stop = line.split('\t')[3]
                strand = line.split('\t')[4].split('\n')[0]
                if strand == "1":
                    ind = seq_id_list.index(acc)
                    fullname = seq_description_list[ind].replace(seq_id_list[ind] + ' ', '')
                    if circ == 'Y':
                        if int(start) - flank > 0:
                            print(
                                '>' + prot_acc + ' | ' + seq_id_list[ind] + ' | ' + str(int(start) - flank) + '-' + str(
                                    int(start) - 1) + ' | ' + prot_acc + ' | ' + fullname + ' | ' + strand + ' | Left_flank')
                            seq = Seq(str(seq_list[ind][int(start) - (flank + 1):int(start) - 1]))
                            print(re.sub("(.{60})", "\\1\n", str(seq), 0, re.DOTALL))
                        else:
                            print('>' + prot_acc + ' | ' + seq_id_list[ind] + ' | ' + str(
                                len(seq_list[ind]) - (flank - int(start))) + '-' + str(
                                len(seq_list[ind])) + ' ~ 1-' + str(
                                int(start) - 1) + ' | ' + prot_acc + ' | ' + fullname + ' | ' + strand + ' | Left_flank')
                            seq = Seq(str(seq_list[ind][0:int(start) - 1]))
                            seq_round = Seq(
                                str(seq_list[ind][len(seq_list[ind]) - ((flank + 1) - int(start)):len(seq_list[ind])]))
                            seq_cat = seq_round + seq
                            print(re.sub("(.{60})", "\\1\n", str(seq_cat), 0, re.DOTALL))
                        if promoter == 'N':
                            if int(stop) + flank < len(seq_list[ind]):
                                print(
                                    '>' + prot_acc + ' | ' + seq_id_list[ind] + ' | ' + str(int(stop) + 1) + '-' + str(
                                        int(stop) + flank) + ' | ' + prot_acc + ' | ' + fullname + ' | ' + strand + ' | Right_flank')
                                seq = Seq(str(seq_list[ind][int(stop):int(stop) + flank]))
                                print(re.sub("(.{60})", "\\1\n", str(seq), 0, re.DOTALL))
                            else:
                                print(
                                    '>' + prot_acc + ' | ' + seq_id_list[ind] + ' | ' + str(int(stop) + 1) + '-' + str(
                                        len(seq_list[ind])) + ' ~ 1-' + str(
                                        flank - (len(seq_list[ind]) - int(
                                            stop))) + ' | ' + prot_acc + ' | ' + fullname + ' | ' + strand + ' | Right_flank')
                                seq = Seq(str(seq_list[ind][0:(flank - (len(seq_list[ind]) - int(stop)))]))
                                seq_round = Seq(str(seq_list[ind][int(stop):]))
                                seq_cat = seq_round + seq
                                print(re.sub("(.{60})", "\\1\n", str(seq_cat), 0, re.DOTALL))

                    elif circ == 'N':
                        if int(start) - flank > 0:
                            print(
                                '>' + prot_acc + ' | ' + seq_id_list[ind] + ' | ' + str(int(start) - flank) + '-' + str(
                                    int(start) - 1) + ' | ' + prot_acc + ' | ' + fullname + ' | ' + strand + ' | Left_flank')
                            seq = Seq(str(seq_list[ind][int(start) - (flank + 1):int(start) - 1]))
                            print(re.sub("(.{60})", "\\1\n", str(seq), 0, re.DOTALL))
                        else:
                            print('>' + prot_acc + ' | ' + seq_id_list[ind] + ' | ' + '1-' + str(
                                int(start) - 1) + ' | ' + prot_acc + ' | ' + fullname + ' | ' + strand + ' | Left_flank')
                            seq = Seq(str(seq_list[ind][0:int(start) - 1]))
                            print(re.sub("(.{60})", "\\1\n", str(seq), 0, re.DOTALL))

                        if promoter == 'N':
                            if int(stop) + flank < len(seq_list[ind]):
                                print(
                                    '>' + prot_acc + ' | ' + seq_id_list[ind] + ' | ' + str(int(stop) + 1) + '-' + str(
                                        int(stop) + flank) + ' | ' + prot_acc + ' | ' + fullname + ' | ' + strand + ' | Right_flank')
                                seq = Seq(str(seq_list[ind][int(stop):int(stop) + flank]))
                                print(re.sub("(.{60})", "\\1\n", str(seq), 0, re.DOTALL))
                            else:
                                print(
                                    '>' + prot_acc + ' | ' + seq_id_list[ind] + ' | ' + str(int(stop) + 1) + '-' + str(
                                        len(seq_list[
                                                ind])) + ' | ' + prot_acc + ' | ' + fullname + ' | ' + strand + ' | Right_flank')
                                seq_round = Seq(str(seq_list[ind][int(stop):]))
                                print(re.sub("(.{60})", "\\1\n", str(seq_round), 0, re.DOTALL))

                elif strand == '-1':

                    ind = seq_id_list.index(acc)
                    fullname = seq_description_list[ind].replace(seq_id_list[ind] + ' ', '')

                    if circ == 'Y':
                        if int(stop) + flank < len(seq_list[ind]):
                            print('>' + prot_acc + ' | ' + seq_id_list[ind] + ' | ' + str(int(stop) + 1) + '-' + str(
                                int(stop) + flank) + ' | ' + prot_acc + ' | ' + fullname + ' | ' + strand + ' | Left_flank')
                            seq = Seq(str(seq_list[ind][int(stop):int(stop) + flank])).reverse_complement()
                            print(re.sub("(.{60})", "\\1\n", str(seq), 0, re.DOTALL))
                        else:
                            print('>' + prot_acc + ' | ' + seq_id_list[ind] + ' | ' + str(int(stop) + 1) + '-' + str(
                                len(seq_list[ind])) + ' ~ 1-' + str(int(stop) + flank - len(seq_list[
                                                                                                ind])) + ' | ' + prot_acc + ' | ' + fullname + ' | ' + strand + ' | Left_flank')
                            seq = Seq(str(seq_list[ind][int(stop):int(stop) + flank])).reverse_complement()
                            seq_round = Seq(
                                str(seq_list[ind][0:int(stop) + flank - len(seq_list[ind])])).reverse_complement()
                            seq_cat = seq_round + seq
                            print(re.sub("(.{60})", "\\1\n", str(seq_cat), 0, re.DOTALL))

                        if promoter == 'N':
                            if int(start) - flank > 0:
                                print('>' + prot_acc + ' | ' + seq_id_list[ind] + ' | ' + str(
                                    int(start) - flank) + '-' + str(
                                    int(start) - 1) + ' | ' + prot_acc + ' | ' + fullname + ' | ' + strand + ' | Right_flank')
                                seq = Seq(
                                    str(seq_list[ind][int(start) - flank - 1:int(start) - 1])).reverse_complement()
                                print(re.sub("(.{60})", "\\1\n", str(seq), 0, re.DOTALL))
                            else:
                                print('>' + prot_acc + ' | ' + seq_id_list[ind] + ' | ' + str(
                                    len(seq_list[ind]) - (flank - int(start))) + '-' + str(
                                    len(seq_list[ind])) + ' ~ 1-' + str(
                                    int(start) - 1) + ' | ' + prot_acc + ' | ' + fullname + ' | ' + strand + ' | Right_flank')
                                seq = Seq(str(
                                    seq_list[ind][len(seq_list[ind]) - (flank - int(start)) - 1:])).reverse_complement()
                                seq_round = Seq(
                                    str(seq_list[ind][0:int(start) - 1])).reverse_complement()
                                seq_cat = seq_round + seq
                                print(re.sub("(.{60})", "\\1\n", str(seq_cat), 0, re.DOTALL))

                    elif circ == 'N':
                        if int(stop) + flank < len(seq_list[ind]):
                            print('>' + prot_acc + ' | ' + seq_id_list[ind] + ' | ' + str(int(stop) + 1) + '-' + str(
                                int(stop) + flank) + ' | ' + prot_acc + ' | ' + fullname + ' | ' + strand + ' | Left_flank')
                            seq = Seq(str(seq_list[ind][int(stop):int(stop) + flank])).reverse_complement()
                            print(re.sub("(.{60})", "\\1\n", str(seq), 0, re.DOTALL))
                        else:
                            print('>' + seq_id_list[ind] + ' | ' + str(int(stop) + 1) + '-' + str(
                                len(seq_list[
                                        ind])) + ' | ' + prot_acc + ' | ' + fullname + ' | ' + strand + ' | Left_flank')
                            seq = Seq(str(seq_list[ind][int(stop):int(stop) + flank])).reverse_complement()
                            print(re.sub("(.{60})", "\\1\n", str(seq), 0, re.DOTALL))

                        if promoter == 'N':
                            if int(start) - flank > 0:
                                print('>' + prot_acc + ' | ' + seq_id_list[ind] + ' | ' + str(
                                    int(start) - flank) + '-' + str(
                                    int(start) - 1) + ' | ' + prot_acc + ' | ' + fullname + ' | ' + strand + ' | Right_flank')
                                seq = Seq(
                                    str(seq_list[ind][int(start) - flank - 1:int(start) - 1])).reverse_complement()
                                print(re.sub("(.{60})", "\\1\n", str(seq), 0, re.DOTALL))
                            else:
                                print('>' + prot_acc + ' | ' + seq_id_list[ind] + ' | 1-' + str(
                                    int(start) - 1) + ' | ' + prot_acc + ' | ' + fullname + ' | ' + strand + ' | Right_flank')
                                seq = Seq(str(seq_list[ind][0:(int(start) - 1)])).reverse_complement()
                                print(re.sub("(.{60})", "\\1\n", str(seq), 0, re.DOTALL))



