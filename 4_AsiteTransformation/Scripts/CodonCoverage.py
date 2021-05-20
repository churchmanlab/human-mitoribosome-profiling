# Count reads in each frame for exons


import sys
infilename_plus = sys.argv[1]
infilename_minus = sys.argv[2]
LibName = sys.argv[3]
CodonFile = open('./Annotations/BEDfiles/cutCodons_hMito_noOverlap.bed','r')  


#Open the bedGraph file to read
InputFilePlus = open(infilename_plus, 'r')
InputFileMinus = open(infilename_minus, 'r')

#Make dictionary for positions and read counts
PlusCounts = {}
MinusCounts = {}

for line in InputFilePlus:
    if line[0:2] != 'tr':
        clist = line.replace('\n','').split('\t')
        PlusCounts[clist[2]] = [clist[3]]
for line in InputFileMinus:
    if line[0:2] != 'tr':
        clist = line.replace('\n','').split('\t')
        MinusCounts[clist[2]] = [clist[3]]
# InputFilePlus.seek(0)
# InputFileMinus.seek(0)

# Go through codon file and make list entry for counts at each codon
PlusCodonCounts = []
MinusCodonCounts = []

for line in CodonFile:
    cols = line.replace('\n','').split('\t')
    POSs = range(int(cols[1]),int(cols[2])+1,1)
    POSs = [str(i) for i in POSs]
    if cols[3] == '+' and int(cols[4].split('.')[-1]) > 6:
        codonCount = int(PlusCounts[POSs[0]][0]) + int(PlusCounts[POSs[1]][0]) + int(PlusCounts[POSs[2]][0])
        PlusCodonCounts.append(codonCount)
    elif cols[3] == '-' and int(cols[4].split('.')[-1]) > 6:
        codonCount = int(MinusCounts[POSs[0]][0]) + int(MinusCounts[POSs[1]][0]) + int(MinusCounts[POSs[2]][0])
        MinusCodonCounts.append(codonCount)

# Get total number of codons and number with >0
CodonCounts = MinusCodonCounts + PlusCodonCounts
Tot = len(CodonCounts)
Covered = Tot - CodonCounts.count(0)
PerCovered = round(float(Covered)/float(Tot) * 100, 1)

def mean(lst):
    return sum(lst) / len(lst)

MeanCoverage = mean(CodonCounts)


#Create a file to write
OutputFile = open((LibName + '_CodonCoverage_ignore1st6.txt'), 'w')

OutputFile.write('Percent_codons_covered' + '\t' + str(PerCovered) + '\n' + 'Average_codon_coverage' + '\t' + str(MeanCoverage) + '\n')


InputFilePlus.close()
InputFileMinus.close()
OutputFile.close()

