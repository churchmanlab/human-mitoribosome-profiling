HELP_STRING = """
get_gene_sequence_densities.py

This takes in chromosome regions, e.g. gene exons, and pulls out the sequences for each exon
and ribosome density data. It puts all exons together and outputs the sequence data and 
density data. 

The chromosome regions have to be listed as:

MT	MT-ND1	+	956	3307-4262
MT	MT-ND2	+	1042	4470-5511
MT	MT-CO1	+	1542	5904-7445
MT	MT-CO2	+	684	7586-8269
MT	MT-ATP8	+	159	8366-8524
"hMitoExonRanges.txt

THe genome sequence should just be of the mito genome. 
The ribosome profiling data should be Asite shifted bed file.

track type=bedGraph color=100,0,100
NC_012920	0	22	0
NC_012920	22	23	4
NC_012920	23	24	2
NC_012920	24	1586	0
NC_012920	1586	1587	81



Author: Stirling Churchman
Date: January 22, 2014
Modified by Mary Couvillion to use different input file formats, April 2016, and for hMitoRP October 2018

     -h     print this help message
     -f     chromosome regions (required)
     -o     output file base (required)
     -g     genome sequence (required)
     -p     plus bedGraph (required)*
     -m     minus bedGraph (required)

 
Run from within CodonFrequency folder
     
example: python get_gene_sequence_densities_MTC_hMitoRP.py -f hMitoExonRanges.txt -l libName -g hMitoGenome.fasta -p Mito_mRNA.noDups.Asite_31to33_P.bedGraph -m Mito_mRNA.noDups.Asite_31to33_M.bedGraph
"""
import sys
from getopt import getopt
import math
import numpy as np
# import matplotlib.pyplot as plt
import itertools
from Bio.Seq import Seq

def main(argv=None):
    if argv is None:
        argv = sys.argv

    inputFile = ""
    libName = ""
    genome = ""
    riboFilePlus = ""
    riboFileMinus = ""
    try:
        optlist, args = getopt(argv[1:], "hf:l:g:p:m:")
    except:
        print("")
        print(HELP_STRING)
        sys.exit(1)
       
    if len(optlist) == 0:
        print("")
        print(HELP_STRING)
        sys.exit(1)
       
    for (opt, opt_arg) in optlist:
        #print opt
        #print opt_arg
        if opt == '-h':
            print("")
            print(HELP_STRING)
            sys.exit(1)
        elif opt == '-f':
            inputFile = opt_arg
        elif opt == '-g':
            genome = opt_arg
        elif opt == '-l':
            libName = opt_arg
        elif opt == '-p':
            riboFilePlus = opt_arg
        elif opt == '-m':
            riboFileMinus = opt_arg
    if inputFile == "" or libName == "" or genome == "" or riboFilePlus == "" or riboFileMinus == "":
        print(HELP_STRING)
        sys.exit(1)


# get ribosome profiling data from a specific region along a genome. Needs a bed coverage file. 
    def getRiboProfilePlus(riboFilePlus, start, stop):
        profile = []
        for line in open(riboFilePlus):
            # This part by Brad
            if line[0] != 't':
                theline = line.split()
                firstPos = int(theline[1])
                secondPos = int(theline[2])
                value = float(theline[3])
                for i in range(firstPos+1,secondPos+1):
                    if i>=start and i<=stop:
                        profile.append(value)
                if secondPos>=stop:
                    break
        return profile
    def getRiboProfileMinus(riboFileMinus, start, stop):
        profile = []
        for line in open(riboFileMinus):
            # This part by Brad
            if line[0] != 't':
                theline = line.split()
                firstPos = int(theline[1])
                secondPos = int(theline[2])
                value = float(theline[3])
                for i in range(firstPos+1,secondPos+1):
                    if i>=start and i<=stop:
                        profile.append(value)
                if secondPos>=stop:
                    break
        return profile


# get the sequence of a region along a genome
    def grabSequence(genome, start, stop):
        index = 1
        InFile = open(genome, 'r')
        seq = ''
        seq_length = stop-start+1
        # find the start of the sequence
        while index<start:
            line = InFile.readline().replace('\n','')
            #usually 50, but just in case...
            line_length = len(line)
            # don't analyze the first line
            if line[0]=='>':continue
            # increase the index
            index+= line_length
            # ask whether we will need more than one line
            if seq_length > (index-start):
                seq += line[(line_length-(index-start)):]
            else:
                seq += line[(line_length-(index-start)):line_length-(index-stop)+1]
        # then read in the lines that will be necessary
        while len(seq)<seq_length:
            line = InFile.readline().replace('\n','')
            index+= len(line)
            if (seq_length-len(seq))>line_length:
                seq+=line
            else:
                seq+= line[:(seq_length-len(seq))]
        InFile.close()
        return seq
        
    gene_seq = {}
    gene_profile = {}
    for line in open(inputFile):
        line = line.replace('\n','').split('\t')
        num_exons = int(len(line[4].replace(';','-').split('-'))/2)
# if this is another exon, then we add it to the other exons
        for i in range(1, num_exons+1):
            if i == 1 and line[2] == '+':
                start = int(line[4].split('-')[i*2-2])
                stop = int(line[4].replace(';','-').split('-')[i*2-1])
                seq = grabSequence(genome,start,stop)
                profile = getRiboProfilePlus(riboFilePlus, start, stop)
                gene_seq[line[1]] = seq
                gene_profile[line[1]] = profile
            if i == 1 and line[2] == '-':
                start = int(line[4].split('-')[i*2-2])
                stop = int(line[4].replace(';','-').split('-')[i*2-1])
                seq = Seq(grabSequence(genome,start,stop))
                profile = getRiboProfileMinus(riboFileMinus, start, stop)
                profile.reverse()
                gene_seq[line[1]] = str(seq.reverse_complement()) #.split(',')[0].tostring() #Seq() is an object from the library Bio.Seq and tostring is a method that returns the first thing from the list
                gene_profile[line[1]] = profile

            if i > 1:
                start = int(line[4].replace(';','-').split('-')[i*2-2])
                stop = int(line[4].replace(';','-').split('-')[i*2-1])
                seq = grabSequence(genome,start,stop)
                profile = getRiboProfile(riboFile, start, stop)
                gene_seq[line[1]] += seq
                gene_profile[line[1]] += profile


    size_range = riboFilePlus.split('Asite_')[1][:6]
    outFile1 = open(libName + '_' + size_range + '_MitoGenes_seq.txt', 'w')
    outFile2 = open(libName + '_' + size_range + '_AsiteFP.txt', 'w')

    for g in gene_seq:
        outFile1.write('>%s\n'% g)
        outFile1.write(gene_seq[g])
        outFile1.write('\n')
        outFile2.write('>%s\n' % g)
        [outFile2.write('%s,' % x) for x in gene_profile[g]]
        outFile2.write('\n')
    outFile1.close()
    outFile2.close()
                 

##############################################
if __name__ == "__main__":
    sys.exit(main())
