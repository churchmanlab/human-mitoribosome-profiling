#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=1G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<youremailaddress>

### WRITTEN:		05/1/2020 by Mary Couvillion

### USE: 			This script will transform .bam human mito ribosome profiling output 
###					to A site bedGraphs for viewing on IGV. 
###					And will count frame preferences across mito mRNAs
###					
###
### REQUIREMENTS:	./logs directory 
###					*_Aligned.Mito_mRNA.noDups.bam file
###					
###					./4_AsiteTransformation/Scripts/SAM2hMitoFBED_softClip.py
###					./4_AsiteTransformation/Scripts/FillMissingPositionsBedGraph.py
###					./4_AsiteTransformation/Scripts/FrameCountBed_chrM_ignore1st3.py
###					./Annotations/chrM.chrom.sizes 
###					./Annotations/BEDfiles/cutCodons_hMito_noOverlap.bed

# load required modules
module load gcc/6.2.0
module load python/2.7.12
module load samtools/1.9
module load bedtools/2.27.1



# Input from command line. LibName is a short descriptive name you want as prefix to all files
LibName=$1
sizeRange=$2 # e.g. 31to33 
offsets=$3 # e.g. 12,13,14

echo $sizeRange
echo $offsets

Ascript="./Scripts/SAM2hMitoFBED_softClip.py"

cd 4_AsiteTransformation

### A site bedGraph
# bam to sam
samtools view -h -o ${LibName}_Aligned.Mito_mRNA.noDups_for${sizeRange}.sam  ../1_AlignData/${LibName}_Aligned.Mito_mRNA.noDups.bam
# ## Custom script to make bed file
python $Ascript -i ${LibName}_Aligned.Mito_mRNA.noDups_for${sizeRange}.sam -s ${sizeRange} -t ${offsets}
# 
# # # convert to bedGraph
# # # plus
genomeCoverageBed -bga -trackline -i ${LibName}_Aligned.Mito_mRNA.noDups.Asite_${sizeRange}_P.bed -g ../Annotations/chrM.chrom.sizes > ${LibName}_Mito_mRNA.noDups.Asite_${sizeRange}_P.bedGraph
# # # minus
genomeCoverageBed -bga -trackline -i ${LibName}_Aligned.Mito_mRNA.noDups.Asite_${sizeRange}_M.bed -g ../Annotations/chrM.chrom.sizes > ${LibName}_Mito_mRNA.noDups.Asite_${sizeRange}_M.bedGraph
# 
rm ${LibName}_Aligned.Mito_mRNA.noDups_for${sizeRange}.sam 
rm ${LibName}_Aligned.Mito_mRNA.noDups.Asite_${sizeRange}_P.bed
rm ${LibName}_Aligned.Mito_mRNA.noDups.Asite_${sizeRange}_M.bed

############ Fill all positions in bedGraph files ###########
# Needed for some downstream scripts (e.g. CodonCoverage.py)

python ./Scripts/FillMissingPositionsBedGraph.py -p ${LibName}_Mito_mRNA.noDups.Asite_${sizeRange}_P.bedGraph -m ${LibName}_Mito_mRNA.noDups.Asite_${sizeRange}_M.bedGraph -s ${sizeRange} 

# # ######## Frame Count ###########
# 
# # First 6 codons are ignored here because in some samples a disproportionate amount of signal is on the 6th codon (and for most genes there is no signal on the first 5)
python ./Scripts/FrameCountBed_chrM_ignore1st6.py ${LibName}_Mito_mRNA.noDups.Asite_${sizeRange}_P.bedGraph ${LibName}_Mito_mRNA.noDups.Asite_${sizeRange}_M.bedGraph ${LibName}_Asite_${sizeRange}
# # Add header
sed -i '1s/^/Frame Frame\t'$LibName'\n/' ${LibName}_Asite_${sizeRange}_FrameCount_ignore1st6.txt

# # ######## Coverage ###########
python ./Scripts/CodonCoverage.py ${LibName}_Mito_mRNA.noDups.Asite_${sizeRange}_PAll.txt ${LibName}_Mito_mRNA.noDups.Asite_${sizeRange}_MAll.txt ${LibName}_Asite_${sizeRange}
# # Add header
sed -i '1s/^/output\t'$LibName'\n/' ${LibName}_Asite_${sizeRange}_CodonCoverage_ignore1st6.txt


