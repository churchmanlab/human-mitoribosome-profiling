HELP_STRING = """
FillMissingPositionsBedGraph.py


Author: Mary Couvillion
Date: 6/16/20 updated from MakeMaxSubcodonProfiles_flat.py to split into fillMissingPositions and maxFrame


     -h     print this help message
     -p     plus strand bedGraph (required)
     -m     minus strand bedGraph (required)
     -s     size range (e.g. 31to33)
     
command line, e.g.: python /Users/Mary/Desktop/Data/hMitoRP/PanAnalysis/ScriptsForReanalysis/FillMissingPositionsBedGraph.py -p *_Mito_mRNA.noDups.Asite_31to33_P.bedGraph -m *_Mito_mRNA.noDups.Asite_31to33_M.bedGraph -s 31to33

"""
import sys
import itertools
from getopt import getopt



def main(argv=None):
    if argv is None:
        argv = sys.argv

    plusFile = ""
    minusFile = ""
    sizeRange = ""
    
    try:
        optlist, args = getopt(argv[1:], "hp:m:s:")
    except:
        print ""
        print HELP_STRING
        sys.exit(1)
       
    if len(optlist) == 0:
        print ""
        print HELP_STRING
        sys.exit(1)
       
    for (opt, opt_arg) in optlist:
        if opt == '-h':
            print ""
            print HELP_STRING
            sys.exit(1)
        elif opt == '-p':
            plusFile = opt_arg
        elif opt == '-m':
            minusFile = opt_arg
        elif opt == '-s':
            sizeRange = opt_arg

    if plusFile == '' or minusFile == '':
        print HELP_STRING
        sys.exit(1)

    PlusFile = open(plusFile,'r')
    MinusFile = open(minusFile, 'r')
  
    tempPlusFile = open(plusFile.replace('.bedGraph','All.txt'), 'w')
    tempMinusFile = open(minusFile.replace('.bedGraph', 'All.txt'), 'w')



    # First fill in missing positions in bedGraph file
    with PlusFile as f:
        for line in f:
            first = line[0]
            if first=='c':
                value = int(line.split()[3])
                for position in range(int(line.split()[1])+1,int(line.split()[2])+1):
                    tempPlusFile.write(line.split()[0] + '\t' + str(position - 1) + '\t' + str(position) + '\t' + str(value) + '\n')
            else:
                tempPlusFile.write(line)

    with MinusFile as f:
        for line in f:
            first = line[0]
            if first=='c':
                value = int(line.split()[3])
                for position in range(int(line.split()[1])+1,int(line.split()[2])+1):
                    tempMinusFile.write(line.split()[0] + '\t' + str(position - 1) + '\t' + str(position) + '\t' + str(value) + '\n')
            else:
                tempMinusFile.write(line)


##############################################
if __name__ == "__main__":
    sys.exit(main())
