#Sort a .sam file by length of reads

#On command line write python SortSAMbyLength.py infilename
import sys 
infilename=sys.argv[1]

#Open the fasta file to read
InputFile = open(infilename, 'r')

#Create a file to write
OutputFile17 = open(infilename.replace('.sam','_17.sam'), 'w')
OutputFile18 = open(infilename.replace('.sam','_18.sam'), 'w')
OutputFile19 = open(infilename.replace('.sam','_19.sam'), 'w')
OutputFile20 = open(infilename.replace('.sam','_20.sam'), 'w')
OutputFile21 = open(infilename.replace('.sam','_21.sam'), 'w')
OutputFile22 = open(infilename.replace('.sam','_22.sam'), 'w')
OutputFile23 = open(infilename.replace('.sam','_23.sam'), 'w')
OutputFile24 = open(infilename.replace('.sam','_24.sam'), 'w')
OutputFile25 = open(infilename.replace('.sam','_25.sam'), 'w')
OutputFile26 = open(infilename.replace('.sam','_26.sam'), 'w')
OutputFile27 = open(infilename.replace('.sam','_27.sam'), 'w')
OutputFile28 = open(infilename.replace('.sam','_28.sam'), 'w')
OutputFile29 = open(infilename.replace('.sam','_29.sam'), 'w')
OutputFile30 = open(infilename.replace('.sam','_30.sam'), 'w')
OutputFile31 = open(infilename.replace('.sam','_31.sam'), 'w')
OutputFile32 = open(infilename.replace('.sam','_32.sam'), 'w')
OutputFile33 = open(infilename.replace('.sam','_33.sam'), 'w')
OutputFile34 = open(infilename.replace('.sam','_34.sam'), 'w')
OutputFile35 = open(infilename.replace('.sam','_35.sam'), 'w')
OutputFile36 = open(infilename.replace('.sam','_36.sam'), 'w')
OutputFile37 = open(infilename.replace('.sam','_37.sam'), 'w')
OutputFile38 = open(infilename.replace('.sam','_38.sam'), 'w')

#Also could do:
#outfilename=infilename.split('.')[0]+'whatever'+infilename.split('.')[1]

count = 0   

for line in InputFile.readlines():
    
    
    if line[0] == '@':
        OutputFile17.write(line)
        OutputFile18.write(line)
        OutputFile19.write(line)
        OutputFile20.write(line)
        OutputFile21.write(line)
        OutputFile22.write(line)
        OutputFile23.write(line)
        OutputFile24.write(line)
        OutputFile25.write(line)
        OutputFile26.write(line)
        OutputFile27.write(line)
        OutputFile28.write(line)
        OutputFile29.write(line)
        OutputFile30.write(line)
        OutputFile31.write(line)
        OutputFile32.write(line)
        OutputFile33.write(line)
        OutputFile34.write(line)
        OutputFile35.write(line)
        OutputFile36.write(line)
        OutputFile37.write(line)
        OutputFile38.write(line)

    if line[0] != '@':
        column = line.split()
        CIGAR = column[5]
        if CIGAR != '*' and 'D' not in CIGAR and 'I' not in CIGAR and 'N' not in CIGAR:
            if 'S' not in CIGAR or (CIGAR[1] != 'S' and CIGAR[2] != 'S'):
                seq_length = int(CIGAR.split('M')[0])
            elif CIGAR[1] == 'S' or CIGAR[2] == 'S':
                seq_length = int(CIGAR.split('S')[1].split('M')[0])

            if seq_length == 17:
               OutputFile17.write(line)
            elif seq_length ==18:
               OutputFile18.write(line)
            elif seq_length ==19:
               OutputFile19.write(line)
            elif seq_length ==20:
               OutputFile20.write(line)
            elif seq_length == 21:
               OutputFile21.write(line)
            elif seq_length ==22:
               OutputFile22.write(line)
            elif seq_length ==23:
               OutputFile23.write(line)
            elif seq_length ==24:
               OutputFile24.write(line)
            elif seq_length ==25:
               OutputFile25.write(line)
            elif seq_length ==26:
               OutputFile26.write(line)
            elif seq_length ==27:
               OutputFile27.write(line)
            elif seq_length ==28:
               OutputFile28.write(line)
            elif seq_length ==29:
               OutputFile29.write(line)
            elif seq_length ==30:
                OutputFile30.write(line)
            elif seq_length ==31:
                OutputFile31.write(line)
            elif seq_length == 32:
                OutputFile32.write(line)
            elif seq_length ==33:
                OutputFile33.write(line)
            elif seq_length ==34:
                 OutputFile34.write(line)
            elif seq_length ==35:
                 OutputFile35.write(line)
            elif seq_length ==36:
                 OutputFile36.write(line)
            elif seq_length ==37:
                 OutputFile37.write(line)
            elif seq_length ==38:
                 OutputFile38.write(line)

            else:
                count = count + 1
            
            
print count, 'records are out of length bounds'
   
        

InputFile.close()
OutputFile17.close()
OutputFile18.close()
OutputFile19.close()
OutputFile20.close()
OutputFile21.close()
OutputFile22.close()
OutputFile23.close()
OutputFile24.close()
OutputFile25.close()
OutputFile26.close()
OutputFile27.close()
OutputFile28.close()
OutputFile29.close()
OutputFile30.close()
OutputFile31.close()
OutputFile32.close()
OutputFile33.close()
OutputFile34.close()
OutputFile35.close()
OutputFile36.close()
OutputFile37.close()
OutputFile38.close()

#OutputFile28to29.close()
#OutputFile32to34.close()
