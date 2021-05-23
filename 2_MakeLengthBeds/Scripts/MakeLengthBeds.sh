#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-4:00
#SBATCH -p short
#SBATCH --mem=5G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<youremailaddress>

### WRITTEN:		05/14/2020 by Mary Couvillion

### USE: 			This script will make 5' and 3' bed files with hMito position and 
###					read length
###					sbatch MakeLengthBeds.sh LibName
###
###					./2_MakeLengthBeds/Scripts/
###					SAM2length5pBED_softClip.py
###					SAM2length3pBED_softClip.py
###

# Load required modules
module load gcc/6.2.0
module load python/2.7.12
module load samtools/1.9


LibName=$1

# Go to output directory
cd 2_MakeLengthBeds

# BAM to SAM
samtools view -@ 3 -h -o ${LibName}_Aligned.Mito_mRNA.noDups_forLengthBed.sam  ../1_AlignData/${LibName}_Aligned.Mito_mRNA.noDups.bam


# 5' bed
python ./Scripts/SAM2length5pBED_softClip.py -f ${LibName}_Aligned.Mito_mRNA.noDups_forLengthBed.sam  -p ${LibName}_Mito_mRNA_lengths5p_P.bed -m ${LibName}_Mito_mRNA_lengths5p_M.bed

# 3' bed
python ./Scripts/SAM2length3pBED_softClip.py -f ${LibName}_Aligned.Mito_mRNA.noDups_forLengthBed.sam  -p ${LibName}_Mito_mRNA_lengths3p_P.bed -m ${LibName}_Mito_mRNA_lengths3p_M.bed

rm ${LibName}_Aligned.Mito_mRNA.noDups_forLengthBed.sam  
