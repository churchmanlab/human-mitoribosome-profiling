#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-06:00
#SBATCH -p short
#SBATCH --mem=80G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<youremailaddress>

### WRITTEN:				03/16/2020 by Mary Couvillion

### USE: 				This script will trim reads, map to the genome, remove PCR duplicates, 
###					and produce files for visualizing on IGV
###					sbatch -e logs/LibName.err -o logs/LibName.log ProcessFASTQ_hMitoRP.sh
###					LibName FASTQfile UMItype
### e.g.
###
### Libs=("Glu_1" "Glu_2" "Gal_30m_1" "Gal_30m_2")
### fastqs=("S1.fastq.gz" "S2.fastq.gz" "S3.fastq.gz" "S4.fastq.gz")
### scriptpath="/n/groups/churchman/user123/Scripts/ProcessFASTQ_hMitoRP.sh"
### UMI="3p6"
### for i in {0..4}
### do
### sbatch -e logs/ProcessFASTQ_${Libs[i]}.err -o logs/ProcessFASTQ_${Libs[i]}.log $scriptpath ${Libs[i]} ${fastqs[i]} $UMI
### done

### REQUIREMENTS:			./logs/ directory
###					*fastq.gz files (gzipped)
###					
###					./Scripts/
###					extractMolecularBarcodeFrom<UMItype>.py
###					removePCRdupsFromBAM_MC20180913.py
###
###					./AnnotationFiles/bowtieindex/
###					rRNA_SSU
###					rRNA_LSU
###					hMito_rRNA
###					mMito_mRNAs
###					tRNAs
###					hMito_tRNA
###					oGAB
###
###					./AnnotationFiles/STARindex/
###					GRCh38_ncRNAs_mM17_merge_riboSeq
###					--> Create this first using instructions in 0_CreateSTARindex
###					or if on O2, can use: /n/groups/churchman/hMitoRP/SequenceFiles/GRCh38_ncRNAs_mM17_merge_riboSeq
###					
###					
###					./AnnotationFiles/BEDfiles/
###					hMito_mRNAsANDsurrounding.bed
###					GENCODE_hg38_proteincoding.sorted.bed
###

# load required modules
module load gcc/6.2.0
module load python/2.7.12
module load samtools/1.9
module load bedtools/2.27.1
module load cutadapt/1.14
module load bowtie/1.2.1.1
module load bamtools/2.4.1
module load fastqc/0.11.3
module load fastx/0.0.13
module load igvtools/2.3.88
module load java/jdk-1.8u112
module load star/2.7.3a



# Input from command line. LibName is a short descriptive name you want as prefix to all files
LibName=$1
InputFASTQ=$2
UMI=$3 # 3p6, 3p6_5p4, 3p10_5p4

# ########## CLEAN READS #################
# ########################################
if [ "${UMI}" = "3p10_5p4" ]
then
	# Cutadapt (minimum is 25-10-4=11 after trimming --untrimmed-output writes reads with no adaptor found INSTEAD of writing to final output file. Do this here so that we keep only reads that have adapter, which means we know where the barcode is
	cutadapt -e 0.2 --untrimmed-output ${LibName}_noAdaptor.fastq -a CTGTAGGCACCATCAAT -m 25 ${InputFASTQ} -o ${LibName}_trimmed.fastq
	rm ${LibName}_noAdaptor.fastq
	###### To extract barcode ######
	python ./Scripts/extractMolecularBarcodeFrom3pr10AND5pr4.py ${LibName}_trimmed.fastq ${LibName}_Cleaned.fastq
	rm ${LibName}_trimmed.fastq

elif [ "${UMI}" = "3p6_5p4" ]
then
	# Cutadapt (minimum is 21-6-4=11 after trimming --untrimmed-output writes reads with no adaptor found INSTEAD of writing to final output file. Do this here so that we keep only reads that have adapter, which means we know where the barcode is
	cutadapt -e 0.2 --untrimmed-output ${LibName}_noAdaptor.fastq -a CTGTAGGCACCATCAAT -m 21 ${InputFASTQ} -o ${LibName}_trimmed.fastq
	rm ${LibName}_noAdaptor.fastq
	###### To extract barcode ######
	python ./Scripts/extractMolecularBarcodeFrom3pr6AND5pr.py ${LibName}_trimmed.fastq ${LibName}_Cleaned.fastq
	rm ${LibName}_trimmed.fastq

elif [ "${UMI}" = "3p6" ]
then 
# 	Cutadapt (minimum is 17-6=11 after trimming --untrimmed-output writes reads with no adaptor found INSTEAD of writing to final output file. Do this here so that we keep only reads that have adapter, which means we know where the barcode is
	cutadapt -e 0.2 --untrimmed-output ${LibName}_noAdaptor.fastq -a CTGTAGGCACCATCAAT -m 17 ${InputFASTQ} -o ${LibName}_trimmed.fastq
	rm ${LibName}_noAdaptor.fastq
	##### To extract barcode ######
	python ./Scripts/extractMolecularBarcodeFrom3pr.py ${LibName}_trimmed.fastq ${LibName}_Cleaned.fastq ${LibName}_Barcodes.txt ${LibName}_Ligation.txt
	rm ${LibName}_trimmed.fastq
	rm ${LibName}_Barcodes.txt
	rm ${LibName}_Ligation.txt
fi
 
############## FASTQC ################
# Make an html report of quality and characteristics of trimmed reads
fastqc ${LibName}_Cleaned.fastq


############## Mapping strategy ################
# Just to get numbers:
# 1. mouse mito mRNA bowtie 5' trim 0mm (Cleaned reads)
# 1. cyto SSU rRNA bowtie 5' trim 1mm
# 2. cyto LSU rRNA bowtie 5' trim 1mm
# 3. mito rRNA bowtie 5' trim 1mm --> output of this to be used for genome mapping with STAR
# 4. cyto tRNA bowtie 5' and 3' trim 1mm
# 5. mito tRNA bowtie 5' and 3' trim 1mm
# 4. STAR concatenated genome soft clip 3mm <-- mito rRNA un
# 5. oGAB bowtie 5' trim 1mm

# Cyto SSU rRNA
bowtie -v 1 -5 1 -S ./AnnotationFiles/bowtieindex/rRNA_SSU ${LibName}_Cleaned.fastq ${LibName}_Nuc_SSUrRNA_almts.sam --un ${LibName}_Nuc_SSUrRNA_un.fastq
rm ${LibName}_Nuc_SSUrRNA_almts.sam

# Cyto LSU rRNA
bowtie -v 1 -5 1 -S ./AnnotationFiles/bowtieindex/rRNA_LSU ${LibName}_Nuc_SSUrRNA_un.fastq ${LibName}_Nuc_LSUrRNA_almts.sam --un ${LibName}_Nuc_LSUrRNA_un.fastq
rm ${LibName}_Nuc_SSUrRNA_un.fastq
rm ${LibName}_Nuc_LSUrRNA_almts.sam

# Mito rRNA
bowtie -v 1 -5 1 -S ./AnnotationFiles/bowtieindex/hMito_rRNA ${LibName}_Nuc_LSUrRNA_un.fastq ${LibName}_hMito_rRNA_almts.sam --al ${LibName}_hMito_rRNA_al.fastq --un ${LibName}_hMito_rRNA_un.fastq
rm ${LibName}_Nuc_LSUrRNA_un.fastq
rm ${LibName}_hMito_rRNA_almts.sam

# Mouse mito mRNA
# First filter input fastq so that only read lengths >= 20 are kept for this part
cutadapt -m 20 ${LibName}_hMito_rRNA_un.fastq -o ${LibName}_hMito_rRNA_20up_un.fastq
bowtie -v 0 -5 1 -S ./AnnotationFiles/bowtieindex/mMito_mRNAs ${LibName}_hMito_rRNA_20up_un.fastq ${LibName}_mMito_mRNAs_0mm.sam --al ${LibName}_mMito_mRNAs_0mm_al.fastq
rm ${LibName}_hMito_rRNA_20up_un.fastq
# Cyto tRNA
bowtie -v 1 -3 3 -5 1 -S ./AnnotationFiles/bowtieindex/tRNAs ${LibName}_hMito_rRNA_un.fastq ${LibName}_Nuc_tRNA_almts.sam --un ${LibName}_Nuc_tRNA_un.fastq
rm ${LibName}_Nuc_tRNA_almts.sam

# Mito tRNA
bowtie -v 1 -3 3 -5 1 -S ./AnnotationFiles/bowtieindex/hMito_tRNA ${LibName}_Nuc_tRNA_un.fastq ${LibName}_hMito_tRNA_almts.sam --al ${LibName}_hMito_tRNA_al.fastq
rm ${LibName}_hMito_tRNA_almts.sam


################ MAP to genome with STAR #################

STAR --runThreadN 4 --genomeDir ./AnnotationFiles/STARindex/GRCh38_ncRNAs_mM17_merge_riboSeq --readFilesIn ${LibName}_hMito_rRNA_un.fastq --outFileNamePrefix ${LibName}_ --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM NM MD --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMultimapNmax 1000 --outFilterMismatchNoverReadLmax 0.07 --outFilterMismatchNmax 3 --outReadsUnmapped Fastx
# 
rm -r ${LibName}__STARtmp
rm ${LibName}_SJ.out.tab

# Map to oGAB
bowtie -v 1 -5 1 -S ./AnnotationFiles/bowtieindex/oGAB ${LibName}_Unmapped.out.mate1 ${LibName}_oGAB_almts.sam
rm ${LibName}_oGAB_almts.sam
# 
# # ############ Get counts from mapping ###########
# # ## Make list of mapping references, in order
echo Input > ${LibName}_MapList.txt
echo Cyto SSU >> ${LibName}_MapList.txt
echo Cyto LSU >> ${LibName}_MapList.txt
echo Mito rRNA >> ${LibName}_MapList.txt
echo Mouse 0mm >> ${LibName}_MapList.txt
echo Cyto tRNA >> ${LibName}_MapList.txt
echo Mito tRNA >> ${LibName}_MapList.txt
echo oGAB >> ${LibName}_MapList.txt
# Grab mapped numbers from bowtie error file
grep -m 1 'reads processed:' logs/Map_${LibName}.err > ${LibName}_bowtieCounts.txt
grep 'reads with at least one reported alignment:' logs/Map_${LibName}.err >> ${LibName}_bowtieCounts.txt
paste <(awk '{print $0}' ${LibName}_MapList.txt) <(awk -F ': ' '{print $2}' ${LibName}_bowtieCounts.txt) > ${LibName}_Counts.txt
# Get mapped reads from STAR log file
echo Mapped from STAR >> ${LibName}_Counts.txt
egrep 'Uniquely mapped reads number | Number of reads mapped to multiple loci' ${LibName}_Log.final.out | awk -F '|\t' '{print $2}' | awk '{sum += $1} END {print sum}' >> ${LibName}_Counts.txt

# Get unmapped reads from STAR log file
echo Unmapped from STAR \(oGAB counts should be subtracted from this\): >> ${LibName}_Counts.txt
egrep 'Number of reads mapped to too many loci | Number of reads unmapped: too short | Number of reads unmapped: other' ${LibName}_Log.final.out | awk -F '|\t' '{print $2}' | awk '{sum += $1} END {print sum}' >> ${LibName}_Counts.txt
 
rm ${LibName}_MapList.txt
rm ${LibName}_bowtieCounts.txt

# # ############ REMOVE PCR DUPLICATES ###########
# # 
# # ## Convert mMito mRNA sam to bam
samtools view -@ 3 -b ${LibName}_mMito_mRNAs_0mm.sam -o ${LibName}_mMito_mRNAs_0mm.bam 
## Sort
samtools sort -@ 3 ${LibName}_mMito_mRNAs_0mm.bam -o ${LibName}_mMito_mRNAs_0mm.sort.bam

## Make BAM index (samtools)
samtools index -@ 3 ${LibName}_Aligned.sortedByCoord.out.bam
samtools index -@ 3 ${LibName}_mMito_mRNAs_0mm.sort.bam

## Removal
python ./Scripts/removePCRdupsFromBAM_MC20180913.py ${LibName}_Aligned.sortedByCoord.out.bam ${LibName}_Aligned_noDups.bam
python ./Scripts/removePCRdupsFromBAM_MC20180913.py ${LibName}_mMito_mRNAs_0mm.sort.bam ${LibName}_mMito_mRNA_noDups.bam



############ MAKE MITO mRNA ONLY BAM FILE ###########

## These should be the footprints
## Make sure they are sorted by coordinate
samtools sort -@ 3 -o ${LibName}_Aligned_noDups.sort.bam ${LibName}_Aligned_noDups.bam

## Choose coordinates by specifying mito bed file
samtools view -@ 3 -b -o ${LibName}_Aligned.Mito_mRNA.noDups.bam -L ./AnnotationFiles/BEDfiles/hMito_mRNAsANDsurrounding.bed ${LibName}_Aligned_noDups.sort.bam
# And hNuc mRNA only file
samtools view -@ 3 -b -o ${LibName}_Aligned.Nuc_mRNA.noDups.bam -L ./AnnotationFiles/BEDfiles/GENCODE_hg38_proteincoding.sorted.bed ${LibName}_Aligned_noDups.sort.bam

# Also for non-deduplicated file, to get counts
samtools view -@ 3 -b -o ${LibName}_hMito_mRNA_Aligned.sortedByCoord.out.bam -L ./AnnotationFiles/BEDfiles/hMito_mRNAsANDsurrounding.bed ${LibName}_Aligned.sortedByCoord.out.bam

# ############ Get counts from noDups ###########
echo '' >> ${LibName}_Counts.txt
echo STAR mapping locations with duplicates: >> ${LibName}_Counts.txt
samtools view -@ 3 -F 0x4 -c ${LibName}_Aligned.sortedByCoord.out.bam >> ${LibName}_Counts.txt
echo STAR mapping locations noDups: >> ${LibName}_Counts.txt
samtools view -@ 3 -F 0x4 -c ${LibName}_Aligned_noDups.bam >> ${LibName}_Counts.txt
echo hMito mRNA with duplicates: >> ${LibName}_Counts.txt
samtools view -@ 3 -F 0x4 -c ${LibName}_hMito_mRNA_Aligned.sortedByCoord.out.bam >> ${LibName}_Counts.txt
echo hMito mRNA noDups \(this is with local alignment and mismatches allowed\): >> ${LibName}_Counts.txt
samtools view -@ 3 -F 0x4 -c ${LibName}_Aligned.Mito_mRNA.noDups.bam >> ${LibName}_Counts.txt
echo mMito mRNA noDups \(this is with end-to-end alignment and 0mm allowed\): >> ${LibName}_Counts.txt
samtools view -@ 3 -F 0x4 -c ${LibName}_mMito_mRNA_noDups.bam >> ${LibName}_Counts.txt

# Add header to Counts file
sed -i '1i\ \n'${LibName} ${LibName}_Counts.txt

# # ############## FASTQC ################
# # Make an html report of quality and characteristics of mito ribo footprints
# 
# # First convert bam to fastq
bedtools bamtofastq -i ${LibName}_Aligned.Mito_mRNA.noDups.bam -fq ${LibName}_Aligned.Mito_mRNA.noDups.fastq
bedtools bamtofastq -i ${LibName}_Aligned.Nuc_mRNA.noDups.bam -fq ${LibName}_Aligned.Nuc_mRNA.noDups.fastq
# Keep only unique lines in fastq (see https://edwards.sdsu.edu/research/sorting-fastq-files-by-their-sequence-identifiers/)
cat ${LibName}_Aligned.Mito_mRNA.noDups.fastq | paste - - - - | sort -k1,1 -t " " | uniq | tr "\t" "\n" > ${LibName}_Aligned.Mito_mRNA.noDups.uniq.fastq 
cat ${LibName}_Aligned.Nuc_mRNA.noDups.fastq | paste - - - - | sort -k1,1 -t " " | uniq | tr "\t" "\n" > ${LibName}_Aligned.Nuc_mRNA.noDups.uniq.fastq 
rm ${LibName}_Aligned.Nuc_mRNA.noDups.bam
rm ${LibName}_Aligned.Mito_mRNA.noDups.fastq
rm ${LibName}_Aligned.Nuc_mRNA.noDups.fastq


# Fastq report
fastqc ${LibName}_Aligned.Mito_mRNA.noDups.uniq.fastq
rm ${LibName}_Aligned.Mito_mRNA.noDups.uniq.fastq
fastqc ${LibName}_Aligned.Nuc_mRNA.noDups.uniq.fastq
rm ${LibName}_Aligned.Nuc_mRNA.noDups.uniq.fastq
fastqc ${LibName}_hMito_tRNA_al.fastq
rm ${LibName}_hMito_tRNA_al.fastq
fastqc ${LibName}_hMito_rRNA_al.fastq
rm ${LibName}_hMito_rRNA_al.fastq
fastqc ${LibName}_mMito_mRNAs_0mm_al.fastq
rm ${LibName}_mMito_mRNAs_0mm_al.fastq

# Get read lengths from fastq report (this will count soft-clipped bases)
list='Aligned.Mito_mRNA.noDups.uniq Aligned.Nuc_mRNA.noDups.uniq hMito_tRNA_al hMito_rRNA_al mMito_mRNAs_0mm_al'
for type in $list
do
	unzip ${LibName}_${type}_fastqc.zip
	sed -n -e '/#Length/,/>>/ p' ${LibName}_${type}_fastqc/fastqc_data.txt | sed -e '1d;$d' > ${LibName}_${type}_LengthDistribution.txt
	# Add header (variable has to be outside of quotes for sed)
	sed '1i\Length\t'${LibName} ${LibName}_${type}_LengthDistribution.txt > ${LibName}_${type}_LengthDist.txt
	rm ${LibName}_${type}_LengthDistribution.txt
	rm -r ${LibName}_${type}_fastqc
done


# Get read length distributions after removing soft clipped bases
list='Aligned.Mito_mRNA.noDups'
for type in $list
do

# convert from .bam to .sam
samtools view -@ 3 -h -O SAM -o ${LibName}_${type}.sam ${LibName}_${type}.bam

# remove soft clipped bases from reads
java -jar /n/groups/churchman/bms36/programs/jvarkit/dist/biostar84452.jar -o ${LibName}_${type}.noSoft.sam ${LibName}_${type}.sam
rm ${LibName}_${type}.sam

# convert back to .bam
samtools view -@ 3 -b -h -o ${LibName}_${type}.noSoft.bam ${LibName}_${type}.noSoft.sam
rm ${LibName}_${type}.noSoft.sam

# sort by coord
samtools sort -@ 3 -o ${LibName}_${type}.noSoft.sort.bam   ${LibName}_${type}.noSoft.bam  
rm ${LibName}_${type}.noSoft.bam

samtools stats ${LibName}_${type}.noSoft.sort.bam > ${LibName}_${type}_samstats.txt
grep ^RL ${LibName}_${type}_samstats.txt | cut -f 2- > ${LibName}_${type}_samstatsReadlength.noSoft.txt
# Add header (variable has to be outside of quotes for sed)
sed -i '1i\Length\t'${LibName} ${LibName}_${type}_samstatsReadlength.noSoft.txt
rm ${LibName}_${type}_samstats.txt

done

# ############ MAKE BEDGRAPH FOR VIEWING ON IGV ###########

### Make bam index
samtools index ${LibName}_Aligned.Mito_mRNA.noDups.bam

## 5' bedGraph (bedtools genomeCoverageBed)
# plus
genomeCoverageBed -bga -5 -strand + -trackline -ibam ${LibName}_Aligned.Mito_mRNA.noDups.bam > ${LibName}_Mito_mRNA.noDups.5p.plus.bedGraph
# minus
genomeCoverageBed -bga -5 -strand - -trackline -ibam ${LibName}_Aligned.Mito_mRNA.noDups.bam > ${LibName}_Mito_mRNA.noDups.5p.minus.bedGraph







### After running for all libs, combine files for all
##################################################
##################################################
##################################################
##########  Do not run as part of script #########
##################################################
##################################################

# Experiment="GalShift" 
# Libs=("Glu_1" "Glu_2" "Gal_30m_1" "Gal_30m_2")

### Read length distributions from fastQC (this will count soft-clipped nts as part of the length)
# list='Aligned.Mito_mRNA.noDups.uniq Aligned.Nuc_mRNA.noDups.uniq hMito_tRNA_al hMito_rRNA_al mMito_mRNAs_0mm_al'
# for type in $list
# do
# paste <(awk '{print $1"\t"$2}' ${Libs[0]}_${type}_LengthDist.txt) <(awk '{print $2}' ${Libs[1]}_${type}_LengthDist.txt) <(awk '{print $2}' ${Libs[2]}_${type}_LengthDist.txt) <(awk '{print $2}' ${Libs[3]}_${type}_LengthDist.txt) > ${Experiment}_${type}_LengthDist.txt 
# done


### Read length distributions from samstats (this does not count soft-clipped nts)
### First make file with all possible lengths to join together
# echo 'Length' > ${Experiment}_lengths.txt
# for x in {13..46}; do
#     echo $x >> ${Experiment}_lengths.txt
# done 

# list='Aligned.Mito_mRNA.noDups'
# # 
# for type in $list
# do
# 
# join -j 1 -a1 -e 0 -o 1.1,2.2 ${Experiment}_lengths.txt ${Libs[0]}_${type}_samstatsReadlength.noSoft.txt > ${Libs[0]}_${type}_samstatsReadlengthAll.noSoft.txt
# join -j 1 -a1 -e 0 -o 1.1,2.2 ${Experiment}_lengths.txt ${Libs[1]}_${type}_samstatsReadlength.noSoft.txt > ${Libs[1]}_${type}_samstatsReadlengthAll.noSoft.txt
# join -j 1 -a1 -e 0 -o 1.1,2.2 ${Experiment}_lengths.txt ${Libs[2]}_${type}_samstatsReadlength.noSoft.txt > ${Libs[2]}_${type}_samstatsReadlengthAll.noSoft.txt
# join -j 1 -a1 -e 0 -o 1.1,2.2 ${Experiment}_lengths.txt ${Libs[3]}_${type}_samstatsReadlength.noSoft.txt > ${Libs[3]}_${type}_samstatsReadlengthAll.noSoft.txt
# done

#### -j 1: Match on column 1
#### -a1: Print all lines in file1 even if it doesn't have match in file2
#### -e 0: Fill missing with '0'

# for type in $list
# do
# paste <(awk '{print $1"\t"$2}' ${Libs[0]}_${type}_samstatsReadlengthAll.noSoft.txt) <(awk '{print $2}' ${Libs[1]}_${type}_samstatsReadlengthAll.noSoft.txt) <(awk '{print $2}' ${Libs[2]}_${type}_samstatsReadlengthAll.noSoft.txt) <(awk '{print $2}' ${Libs[3]}_${type}_samstatsReadlengthAll.noSoft.txt) > ${Experiment}_${type}_noSoftLengthDist.txt 
# done

### Counts files with read mapping stats
# Samps="Glu_1 Glu_2 Gal_30m_1 Gal_30m_2"
# 
# echo All Counts > ${Experiment}_Counts.txt
# for Samp in $Samps
# do
# cat ${Samp}_Counts.txt >> ${Experiment}_Counts.txt
# done



