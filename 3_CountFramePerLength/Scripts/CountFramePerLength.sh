#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-06:00
#SBATCH -p short
#SBATCH --mem=30G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<youremailaddress>

### WRITTEN:		05/1/2020 by Mary Couvillion

### USE: 			This script will separate a mapped bam(sam) file into lengths  
###					specified and calculate frame frequencies for each  
###					size. Outputs should be checked before running "AsiteAndCountFrame"
###					LibName 


### REQUIREMENTS:	./logs/ directory 
###					*_Aligned.Mito_mRNA.noDups.bam file
###					
###					./3_CountFramePerLength/Scripts/SortSAMbyLength.py
###					./3_CountFramePerLength/Scripts/FrameCountBed_chrM_hMito.py
###					./3_CountFramePerLength/Annotations/BEDfiles/ 
###					cutCodons_hMito_noOverlap.bed (used for frame count)

LibName=$1

# Define sequence lengths being quarried
# sizes='18 19 20 21 29 30 31 32 33 34 35'
sizes='25 26 27 28 29 30 31 32'



echo ${LibName}

# BAM to SAM
echo Converting to sam
samtools view -@ 3 -h ${LibName}_Aligned.Mito_mRNA.noDups.bam -o ${LibName}_Aligned.Mito_mRNA.noDups_temp.sam

# Split SAM according to length
echo Sorting sam by length
python ./3_CountFramePerLength/Scripts/SortSAMbyLength.py ${LibName}_Aligned.Mito_mRNA.noDups_temp.sam
echo Done
rm ${LibName}_Aligned.Mito_mRNA.noDups_temp.sam



# SAM to BAM
for x in $sizes
do
samtools view -@ 3 -b ${LibName}_Aligned.Mito_mRNA.noDups_temp_${x}.sam -o ${LibName}_Aligned.Mito_mRNA.noDups_${x}.bam
rm ${LibName}_Aligned.Mito_mRNA.noDups_temp_${x}.sam
done

# Now remove the rest that aren't in the lengths to be queried list
rm ${LibName}_Aligned.Mito_mRNA.noDups_temp_1*.sam
rm ${LibName}_Aligned.Mito_mRNA.noDups_temp_2*.sam
rm ${LibName}_Aligned.Mito_mRNA.noDups_temp_3*.sam
rm ${LibName}_Aligned.Mito_mRNA.noDups_temp_4*.sam

# Sort BAM
for x in $sizes
do
samtools sort -@ 3 ${LibName}_Aligned.Mito_mRNA.noDups_${x}.bam -o ${LibName}_Aligned.Mito_mRNA.noDups_${x}.sort.bam
done

# Index BAM
for x in $sizes
do
samtools index -@ 3 ${LibName}_Aligned.Mito_mRNA.noDups_${x}.sort.bam
done

# Make bedGraphs for frame count
# 5'
for x in $sizes
do
echo Making ${x}.5p bedGraphs
genomeCoverageBed -bga -5 -strand + -trackline -trackopts 'color=100,0,0' -ibam ${LibName}_Aligned.Mito_mRNA.noDups_${x}.sort.bam > ${LibName}_Mito_mRNA.noDups_${x}.5pP.bedGraph
 
genomeCoverageBed -bga -5 -strand - -trackline -trackopts 'color=0,100,0' -ibam ${LibName}_Aligned.Mito_mRNA.noDups_${x}.sort.bam > ${LibName}_Mito_mRNA.noDups_${x}.5pM.bedGraph
done

# 3'
for x in $sizes
do
echo Making ${x}.3p bedGraphs
genomeCoverageBed -bga -3 -strand + -trackline -trackopts 'color=100,0,0' -ibam ${LibName}_Aligned.Mito_mRNA.noDups_${x}.sort.bam > ${LibName}_Mito_mRNA.noDups_${x}.3pP.bedGraph

genomeCoverageBed -bga -3 -strand - -trackline -trackopts 'color=0,100,0' -ibam ${LibName}_Aligned.Mito_mRNA.noDups_${x}.sort.bam > ${LibName}_Mito_mRNA.noDups_${x}.3pM.bedGraph
done

# Frame Count from 5'
for x in $sizes
do
echo Doing frame count for ${x}.5p
python ./3_CountFramePerLength/Scripts/FrameCountBed_chrM_hMito.py ${LibName}_Mito_mRNA.noDups_${x}.5pP.bedGraph ${LibName}_Mito_mRNA.noDups_${x}.5pM.bedGraph ${LibName}_${x}_5p
done

# # Frame Count from 3'
for x in $sizes
do
echo Doing frame count for ${x}.3p
python ./3_CountFramePerLength/Scripts/FrameCountBed_chrM_hMito.py ${LibName}_Mito_mRNA.noDups_${x}.3pP.bedGraph ${LibName}_Mito_mRNA.noDups_${x}.3pM.bedGraph ${LibName}_${x}_3p
done

for x in $sizes
do
echo Removing ${x}nt bedGraph files
rm ${LibName}_Mito_mRNA.noDups_${x}.5pP.bedGraph
rm ${LibName}_Mito_mRNA.noDups_${x}.5pM.bedGraph
rm ${LibName}_Mito_mRNA.noDups_${x}.3pP.bedGraph
rm ${LibName}_Mito_mRNA.noDups_${x}.3pM.bedGraph
done

echo Done calculating frame counts

# remove rest of temp files
rm ${LibName}_Aligned.Mito_mRNA.noDups_*bam*

# Make a first column for combined frame count files
echo Making combined output file
echo size > ${LibName}_read_sizes.txt
for x in $sizes
do
for i in `seq 3`; do echo ${x} >> ${LibName}_read_sizes.txt ;done
done

# Combine Frame Count files
echo -e Frame:'\t'${LibName} > ${LibName}_FrameCounts_5p.txt
for x in $sizes
do
cat ${LibName}_${x}_5p_FrameCount.txt >> ${LibName}_FrameCounts_5p.txt
done

echo -e Frame:'\t'${LibName} > ${LibName}_FrameCounts_3p.txt
for x in $sizes
do
cat ${LibName}_${x}_3p_FrameCount.txt >> ${LibName}_FrameCounts_3p.txt
done

paste <(awk '{print $0}' ${LibName}_read_sizes.txt) <(awk -F':\t' '{print $2}' ${LibName}_FrameCounts_5p.txt) > ${LibName}_5p_FrameCounts.txt
paste <(awk '{print $0}' ${LibName}_read_sizes.txt) <(awk -F':\t' '{print $2}' ${LibName}_FrameCounts_3p.txt) > ${LibName}_3p_FrameCounts.txt

# Remove individual size frame count files
rm ${LibName}_*_FrameCount.txt
# Remove initial concatenated files without read sizes
rm ${LibName}_FrameCounts_5p.txt
rm ${LibName}_FrameCounts_3p.txt


################################################
### After running for all libs, combine all
Experiment="Li" # OXPHOSinh GalShift MyoDiff1 HelaIP Hela OXPHOSinh2 OXPHOSinh456 HEK HelaFracs Rooijers HEKpearce
# # # Libs=("Glu_1" "Glu_2" "Gal_30m_1" "Gal_30m_2" "Gal_90m_1" "Gal_90m_2" "Gal_4h_1" "Gal_4h_2" "Gal_72h_1" "Gal_72h_2")
# # # Libs=("DMSO_1" "ROT_1nM" "ROT_5nM" "ROT_10nM" "ROT_100nM" "AA_1nM" "AA_5nM" "AA_10nM" "AA_100nM" )
# Libs=("DMSO_1" "DMSO_2" "AA_10nM_1" "AA_10nM_2" "AA_100nM_1" "AA_100nM_2" "AA_500nM_1" "AA_500nM_2" "AA_1uM_1" "AA_1uM_2" "AA_10uM_1" "AA_10uM_2")
# # # Libs=("NonDiff_1" "NonDiff_2" "Diff_4h_1" "Diff_4h_2" "Diff_2d_1" "Diff_2d_2" "Diff_4d_1" "Diff_4d_2" "Diff_7d_1" "Diff_7d_2" "Diff_10d_1" "Diff_10d_2")
# # # Libs=("Hela_IP_4U_1" "Hela_IP_4U_2" "Hela_IP_8U_1" "Hela_IP_8U_2")
# # # Libs=("Hela_noDrug" "Hela_plDrug")
# Libs=("NT4_1" "DMSO4_1" "AA4_100nM" "DMSO5_1" "DMSO5_2" "AA5_100nM_1" "AA5_100nM_2" "DMSO6_1" "DMSO6_2" "AA6_100nM_1" "AA6_100nM_2")
# Libs=("WT1_2" "WT2_2")
# Libs=("H_1_6" "H_2_6" "HM_1_6" "HM_2_6")
# Libs=("BJ_1" "BJ_2")
Libs=("HCT_WT1" "HCT_WT2")

ends='5 3'
for end in $ends
do
paste <(awk '{print $1"\t"$2}' ${Libs[0]}_${end}p_FrameCounts.txt) <(awk '{print $2}' ${Libs[1]}_${end}p_FrameCounts.txt) > ${Experiment}_${end}p_allFrameCounts.txt #    <(awk '{print $2}' 
done