# Count reads in each frame for exons


#On command line write python FrameCountBed.py <plus bedGraph> <minus bedGraph> <outfile base name>
import sys
infilename_plus = sys.argv[1]
infilename_minus = sys.argv[2]
outfilename = sys.argv[3]
CodonFile = open('./3_CountFramePerLength/Annotations/BEDfiles/cutCodons_hMito_noOverlap.bed','r') 

#Make lists of positions for each frame
Frame1Positions_plus = []
Frame2Positions_plus = []
Frame3Positions_plus = []
Frame1Positions_minus = []
Frame2Positions_minus = []
Frame3Positions_minus = []
#chromosome = CodonFile.readlines()[0].split('\t')[0]
for line in CodonFile.readlines():
    clist = line.split('\t')
    #chr = clist[0]
    if clist[3] == '+':
        Frame1Positions_plus.append('chrM' + clist[1])
        Frame2Positions_plus.append('chrM' + str(int(clist[1])+1))
        Frame3Positions_plus.append('chrM' + clist[2])
    elif clist[3] == '-':
        Frame1Positions_minus.append('chrM' + clist[2])
        Frame2Positions_minus.append('chrM' + str(int(clist[2])-1))
        Frame3Positions_minus.append('chrM' + clist[1])

#Open the bedGraph file to read
InputFilePlus = open(infilename_plus, 'r')
InputFileMinus = open(infilename_minus, 'r')

#Create a file to write
OutputFile = open((outfilename + '_FrameCount.txt'), 'w')

#Go through lines and add up read counts
Frame1Count = 0.0
Frame2Count = 0.0
Frame3Count = 0.0

count = 0
for line in InputFilePlus.readlines():
    count = count +1
    first = line[0]
    if first=='c':
        column = line.split()
        if count % 2000 == 0:
            print 'Plus ', count, Frame1Count, Frame2Count, Frame3Count
        if column[0] + column[2] in Frame2Positions_plus: #Because most match the second frame maybe this will make it go faster
            Frame2Count = Frame2Count + float(column[3])
            Frame2Positions_plus.remove(column[0] + column[2])
        elif column[0] + column[2] in Frame1Positions_plus:
            Frame1Count = Frame1Count + float(column[3])
            Frame1Positions_plus.remove(column[0] + column[2])
        elif column[0] + column[2] in Frame3Positions_plus:
            Frame3Count = Frame3Count + float(column[3])
            Frame3Positions_plus.remove(column[0] + column[2])
        else:
            continue
print 'Plus file done'
print 'Plus Frame 1 count:', Frame1Count
print 'Plus Frame 2 count:', Frame2Count
print 'Plus Frame 3 count:', Frame3Count

count = 0
for line in InputFileMinus.readlines():
    count = count +1
    first = line[0]
    if first=='c':
        column = line.split()
        if count % 2000 == 0:
            print 'Minus ', count, Frame1Count, Frame2Count, Frame3Count
        if column[0] + column[2] in Frame2Positions_minus:
            Frame2Count = Frame2Count + float(column[3])
            Frame2Positions_minus.remove(column[0] + column[2])
        elif column[0] + column[2] in Frame1Positions_minus:
            Frame1Count = Frame1Count + float(column[3])
            Frame1Positions_minus.remove(column[0] + column[2])
        elif column[0] + column[2] in Frame3Positions_minus:
            Frame3Count = Frame3Count + float(column[3])
            Frame3Positions_minus.remove(column[0] + column[2])
        else:
            continue
print 'Minus file done'
print 'Minus Frame 1 count:', Frame1Count
print 'Minus Frame 2 count:', Frame2Count
print 'Minus Frame 3 count:', Frame3Count


OutputFile.write('Frame 1:' + '\t' + str(Frame1Count) + '\n' + 'Frame 2:' + '\t' + str(Frame2Count) + '\n' + 'Frame 3:' + '\t' + str(Frame3Count) + '\n')
    

InputFilePlus.close()
InputFileMinus.close()
OutputFile.close()

