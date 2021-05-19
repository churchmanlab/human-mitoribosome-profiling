#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-4:00
#SBATCH -p short
#SBATCH --mem=5G
#SBATCH --mail-type=END
#SBATCH --mail-user=<youremailaddress>

### WRITTEN:		05/14/2020 by Mary Couvillion

### USE: 			This script will make 5' and 3' bed files with hMito position and 
###					read length
###					sbatch MakeLengthBeds.sh LibName
###					

LibName=$1

# BAM to SAM
samtools view -@ 3 -h -o ${LibName}_Aligned.Mito_mRNA.noDups_forLengthBed.sam  ${LibName}_Aligned.Mito_mRNA.noDups.bam


# 5' bed
python ./2_MakeLengthBeds/Scripts/SAM2length5pBED_softClip.py -f ${LibName}_Aligned.Mito_mRNA.noDups_forLengthBed.sam  -p ${LibName}_Mito_mRNA_lengths5p_P.bed -m ${LibName}_Mito_mRNA_lengths5p_M.bed

# 3' bed
python ./2_MakeLengthBeds/Scripts/SAM2length3pBED_softClip.py -f ${LibName}_Aligned.Mito_mRNA.noDups_forLengthBed.sam  -p ${LibName}_Mito_mRNA_lengths3p_P.bed -m ${LibName}_Mito_mRNA_lengths3p_M.bed

rm ${LibName}_Aligned.Mito_mRNA.noDups_forLengthBed.sam  
