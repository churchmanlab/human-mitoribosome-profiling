#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-06:00
#SBATCH -p short
#SBATCH --mem=80G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<youremailaddress>

### WRITTEN:		12/09/2020 by Mary Couvillion

### USE: 				This script will trim reads, map to the genome, remove PCR      ###					duplicates, 
###					and produce files for visualizing on IGV
###					sbatch -e logs/1_AlignData_CytoRP_LibName.err -o logs/1_AlignData_CytoRP_LibName.log ProcessFASTQ_hMitoRP.sh
###					LibName FASTQfile UMItype
###
### REQUIREMENTS:			./logs/ directory
###					*fastq.gz files (gzipped)
###					
###					./1_AlignData/Scripts/
###					extractMolecularBarcodeFrom<UMItype>.py
###					removePCRdupsFromBAM_MC20180913.py
###
###					./Annotations/bowtieindex/
###					rRNA_SSU
###					rRNA_LSU
###					hMito_rRNA
###					tRNAs
###					hMito_tRNA
###					oGAB
###
###					./Annotations/STARindex/
###					GRCh38_ncRNAs_riboSeq
###					--> Create this first using instructions in 0_CreateSTARindex
###					(leave out mouse genome here though)
###					or if on O2, can use: 
###					/n/groups/churchman/hMitoRP/SequenceFiles/GRCh38_ncRNAs_riboSeq
###					
###					
###					./Annotations/BEDfiles/
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

STARindexPATH="../Annotations/STARindex/GRCh38_ncRNAs_riboSeq"
# STARindexPATH="/n/groups/churchman/hMitoRP/SequenceFiles/GRCh38_ncRNAs_riboSeq"


# Input from command line. LibName is a short descriptive name you want as prefix to all files
LibName=$1
InputFASTQ=$2
UMI=$3 # 3p6, 3p6_5p4, 3p10_5p4

cd ./1_AlignData_CytoRP

# ########## CLEAN READS #################
# ########################################
if [ "${UMI}" = "3p10_5p4" ]
then
	# Cutadapt (minimum is 25-10-4=11 after trimming --untrimmed-output writes reads with     		  no adaptor found INSTEAD of writing to final output file. Do this here so that we keep only reads that have adapter, which means we know where the barcode is
	cutadapt -e 0.2 --untrimmed-output ${LibName}_noAdaptor.fastq -a CTGTAGGCACCATCAAT -m 25 ../${InputFASTQ} -o ${LibName}_trimmed.fastq
	rm ${LibName}_noAdaptor.fastq
	###### To extract barcode ######
	python ../1_AlignData/Scripts/extractMolecularBarcodeFrom3pr10AND5pr4.py ${LibName}_trimmed.fastq ${LibName}_Cleaned.fastq
	rm ${LibName}_trimmed.fastq
elif [ "${UMI}" = "3p6_5p4" ]
then
	# Cutadapt (minimum is 21-6-4=11 after trimming --untrimmed-output writes reads with no adaptor found INSTEAD of writing to final output file. Do this here so that we keep only reads that have adapter, which means we know where the barcode is
	cutadapt -e 0.2 --untrimmed-output ${LibName}_noAdaptor.fastq -a CTGTAGGCACCATCAAT -m 21 ../${InputFASTQ} -o ${LibName}_trimmed.fastq
	rm ${LibName}_noAdaptor.fastq
	###### To extract barcode ######
	python ../1_AlignData/Scripts/extractMolecularBarcodeFrom3pr6AND5pr.py ${LibName}_trimmed.fastq ${LibName}_Cleaned.fastq
	rm ${LibName}_trimmed.fastq

elif [ "${UMI}" = "3p6" ]
then 
# 	Cutadapt (minimum is 17-6=11 after trimming --untrimmed-output writes reads with no adaptor found INSTEAD of writing to final output file. Do this here so that we keep only reads that have adapter, which means we know where the barcode is
	cutadapt -e 0.2 --untrimmed-output ${LibName}_noAdaptor.fastq -a CTGTAGGCACCATCAAT -m 17 ../${InputFASTQ} -o ${LibName}_trimmed.fastq
	rm ${LibName}_noAdaptor.fastq
	##### To extract barcode ######
	python ../1_AlignData/Scripts/extractMolecularBarcodeFrom3pr.py ${LibName}_trimmed.fastq ${LibName}_Cleaned.fastq ${LibName}_Barcodes.txt ${LibName}_Ligation.txt
	rm ${LibName}_trimmed.fastq
	rm ${LibName}_Barcodes.txt
	rm ${LibName}_Ligation.txt
	
fi
#  
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
bowtie -v 1 -5 1 -S ../Annotations/bowtieindex/rRNA_SSU ${LibName}_Cleaned.fastq ${LibName}_Nuc_SSUrRNA_almts.sam --un ${LibName}_Nuc_SSUrRNA_un.fastq
rm ${LibName}_Nuc_SSUrRNA_almts.sam

# Cyto LSU rRNA
bowtie -v 1 -5 1 -S ../Annotations/bowtieindex/rRNA_LSU ${LibName}_Nuc_SSUrRNA_un.fastq ${LibName}_Nuc_LSUrRNA_almts.sam --un ${LibName}_Nuc_LSUrRNA_un.fastq --al ${LibName}_Nuc_LSUrRNA_al.fastq
rm ${LibName}_Nuc_SSUrRNA_un.fastq
rm ${LibName}_Nuc_LSUrRNA_almts.sam

# Mito rRNA
bowtie -v 1 -5 1 -S ../Annotations/bowtieindex/hMito_rRNA ${LibName}_Nuc_LSUrRNA_un.fastq ${LibName}_hMito_rRNA_almts.sam --al ${LibName}_hMito_rRNA_al.fastq --un ${LibName}_hMito_rRNA_un.fastq
rm ${LibName}_hMito_rRNA_almts.sam

# Cyto tRNA
bowtie -v 1 -3 3 -5 1 -S ../Annotations/bowtieindex/tRNAs ${LibName}_hMito_rRNA_un.fastq ${LibName}_Nuc_tRNA_almts.sam --un ${LibName}_Nuc_tRNA_un.fastq --al ${LibName}_Nuc_tRNA_al.fastq
rm ${LibName}_Nuc_tRNA_almts.sam

# Mito tRNA
bowtie -v 1 -3 3 -5 1 -S ../Annotations/bowtieindex/hMito_tRNA ${LibName}_Nuc_tRNA_un.fastq ${LibName}_hMito_tRNA_almts.sam --al ${LibName}_hMito_tRNA_al.fastq
rm ${LibName}_hMito_tRNA_almts.sam

###################### OPTIONAL for later running RiboQC #####################

############ Make STAR index for RiboseQC
############ Don't include ncRNAs because gtf formatting isn't right for those
############ mkdir GRCh38_riboSeq
############ sbatch -p short -t 0-3:00 -c 6 --mem=80G --wrap="STAR --runThreadN 6 --runMode genomeGenerate --genomeDir GRCh38_riboSeq --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile gencode.v30.primary_assembly.annotation.gtf --sjdbOverhang 50"

############ Also need to mkdir ./1_AlignData_CytoRP/STAR_all

## RiboQC_STARindexPATH="../Annotations/STARindex/GRCh38_riboSeq"

################ MAP to genome with STAR #################

STAR --runThreadN 4 --genomeDir $STARindexPATH --readFilesIn ${LibName}_hMito_rRNA_un.fastq --outFileNamePrefix ${LibName}_ --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM NM MD --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMultimapNmax 1000 --outFilterMismatchNoverReadLmax 0.07 --outFilterMismatchNmax 3 --outReadsUnmapped Fastx

rm -r ${LibName}__STARtmp
rm ${LibName}_SJ.out.tab

# Map to oGAB
bowtie -v 1 -5 1 -S ../Annotations/bowtieindex/oGAB ${LibName}_Unmapped.out.mate1 ${LibName}_oGAB_almts.sam
rm ${LibName}_oGAB_almts.sam

###### OPTIONAL for later running RiboQC #####################
################ MAP FOR RIBOseqQC #################

# Use rRNA-filtered reads to not take so much memory/space

# Multi
# STAR --runThreadN 4 --genomeDir $RiboQC_STARindexPATH --readFilesIn ${LibName}_Nuc_LSUrRNA_un.fastq --outFileNamePrefix STAR_all/${LibName}_all_multi_ --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM NM MD --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMultimapNmax 1000 --outFilterMismatchNoverReadLmax 0.07 --outFilterMismatchNmax 3 --outReadsUnmapped Fastx
# Unique
# STAR --runThreadN 4 --genomeDir $RiboQC_STARindexPATH --readFilesIn ${LibName}_Nuc_LSUrRNA_un.fastq --outFileNamePrefix STAR_all/${LibName}_all_unique_ --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM NM MD --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMultimapNmax 1 --outFilterMismatchNoverReadLmax 0.07 --outFilterMismatchNmax 3 # --outReadsUnmapped Fastx



############ Get counts from mapping ###########
## Make list of mapping references, in order
echo Input > ${LibName}_MapList.txt
echo Cyto SSU >> ${LibName}_MapList.txt
echo Cyto LSU >> ${LibName}_MapList.txt
echo Mito rRNA >> ${LibName}_MapList.txt
echo Cyto tRNA >> ${LibName}_MapList.txt
echo Mito tRNA >> ${LibName}_MapList.txt
echo oGAB >> ${LibName}_MapList.txt
# Grab mapped numbers from bowtie error file
grep -m 1 'reads processed:' ../logs/1_AlignData_CytoRP_${LibName}.err > ${LibName}_bowtieCounts.txt
grep 'reads with at least one reported alignment:' ../logs/1_AlignData_CytoRP_${LibName}.err >> ${LibName}_bowtieCounts.txt
paste <(awk '{print $0}' ${LibName}_MapList.txt) <(awk -F ': ' '{print $2}' ${LibName}_bowtieCounts.txt) > ${LibName}_Counts.txt
# Get mapped reads from STAR log file
echo Mapped from STAR >> ${LibName}_Counts.txt
egrep 'Uniquely mapped reads number | Number of reads mapped to multiple loci' ${LibName}_Log.final.out | awk -F '|\t' '{print $2}' | awk '{sum += $1} END {print sum}' >> ${LibName}_Counts.txt

# Get unmapped reads from STAR log file
echo Unmapped from STAR \(oGAB counts should be subtracted from this\): >> ${LibName}_Counts.txt
egrep 'Number of reads mapped to too many loci | Number of reads unmapped: too short | Number of reads unmapped: other' ${LibName}_Log.final.out | awk -F '|\t' '{print $2}' | awk '{sum += $1} END {print sum}' >> ${LibName}_Counts.txt
 
rm ${LibName}_MapList.txt
rm ${LibName}_bowtieCounts.txt

############ REMOVE PCR DUPLICATES ###########

## Make BAM index (samtools)
samtools index -@ 3 ${LibName}_Aligned.sortedByCoord.out.bam
## Uncomment below for running RiboSeQC later
# samtools index -@ 3 STAR_all/${LibName}_all_multi_Aligned.sortedByCoord.out.bam
# samtools index -@ 3 STAR_all/${LibName}_all_unique_Aligned.sortedByCoord.out.bam

## Removal
python ../1_AlignData/Scripts/removePCRdupsFromBAM_MC20180913.py ${LibName}_Aligned.sortedByCoord.out.bam ${LibName}_Aligned_noDups.bam
## Uncomment below for running RiboSeQC later
# python ../1_AlignData/Scripts/removePCRdupsFromBAM_MC20180913.py STAR_all/${LibName}_all_multi_Aligned.sortedByCoord.out.bam STAR_all/${LibName}_all_multi_Aligned_noDups.bam
# python ../1_AlignData/Scripts/removePCRdupsFromBAM_MC20180913.py STAR_all/${LibName}_all_unique_Aligned.sortedByCoord.out.bam STAR_all/${LibName}_all_unique_Aligned_noDups.bam


############ MAKE MITO mRNA ONLY BAM FILE ###########

## These should be the footprints
## Make sure they are sorted by coordinate
samtools sort -@ 3 -o ${LibName}_Aligned_noDups.sort.bam ${LibName}_Aligned_noDups.bam

## Choose coordinates by specifying mito bed file
samtools view -@ 3 -b -o ${LibName}_Aligned.Mito_mRNA.noDups.bam -L ../Annotations/BEDfiles/hMito_mRNAsANDsurrounding.bed ${LibName}_Aligned_noDups.sort.bam
# And hNuc mRNA only file
samtools view -@ 3 -b -o ${LibName}_Aligned.Nuc_mRNA.noDups.bam -L ../Annotations/BEDfiles/GENCODE_hg38_proteincoding.sorted.bed ${LibName}_Aligned_noDups.sort.bam

# Also for non-deduplicated file, to get counts
samtools view -@ 3 -b -o ${LibName}_Nuc_mRNA_Aligned.sortedByCoord.out.bam -L ../Annotations/BEDfiles/GENCODE_hg38_proteincoding.sorted.bed  ${LibName}_Aligned.sortedByCoord.out.bam

# ############ Get counts from noDups ###########
echo '' >> ${LibName}_Counts.txt
echo STAR mapping locations with duplicates: >> ${LibName}_Counts.txt
samtools view -@ 3 -F 0x4 -c ${LibName}_Aligned.sortedByCoord.out.bam >> ${LibName}_Counts.txt
echo STAR mapping locations noDups: >> ${LibName}_Counts.txt
samtools view -@ 3 -F 0x4 -c ${LibName}_Aligned_noDups.bam >> ${LibName}_Counts.txt
echo Nuc mRNA with duplicates: >> ${LibName}_Counts.txt
samtools view -@ 3 -F 0x4 -c ${LibName}_Nuc_mRNA_Aligned.sortedByCoord.out.bam >> ${LibName}_Counts.txt
echo hNuc mRNA noDups \(this is with local alignment and mismatches allowed\): >> ${LibName}_Counts.txt
samtools view -@ 3 -F 0x4 -c ${LibName}_Aligned.Nuc_mRNA.noDups.bam >> ${LibName}_Counts.txt
echo mMito mRNA noDups \(this is with end-to-end alignment and 0mm allowed\): >> ${LibName}_Counts.txt
samtools view -@ 3 -F 0x4 -c ${LibName}_mMito_mRNA_noDups.bam >> ${LibName}_Counts.txt

# Add header to Counts file
sed -i '1i\ \n'${LibName} ${LibName}_Counts.txt

# ############## FASTQC ################
# Make an html report of quality and characteristics of mito ribo footprints

# First convert bam to fastq
bedtools bamtofastq -i ${LibName}_Aligned.Mito_mRNA.noDups.bam -fq ${LibName}_Aligned.Mito_mRNA.noDups.fastq
bedtools bamtofastq -i ${LibName}_Aligned.Nuc_mRNA.noDups.bam -fq ${LibName}_Aligned.Nuc_mRNA.noDups.fastq
## Keep only unique lines in fastq (see https://edwards.sdsu.edu/research/sorting-fastq-files-by-their-sequence-identifiers/)
cat ${LibName}_Aligned.Mito_mRNA.noDups.fastq | paste - - - - | sort -k1,1 -t " " | uniq | tr "\t" "\n" > ${LibName}_Aligned.Mito_mRNA.noDups.uniq.fastq 
cat ${LibName}_Aligned.Nuc_mRNA.noDups.fastq | paste - - - - | sort -k1,1 -t " " | uniq | tr "\t" "\n" > ${LibName}_Aligned.Nuc_mRNA.noDups.uniq.fastq 
rm ${LibName}_Aligned.Mito_mRNA.noDups.fastq
rm ${LibName}_Aligned.Nuc_mRNA.noDups.fastq


# Fastq report
fastqc ${LibName}_Aligned.Mito_mRNA.noDups.uniq.fastq
rm ${LibName}_Aligned.Mito_mRNA.noDups.uniq.fastq
fastqc ${LibName}_Aligned.Nuc_mRNA.noDups.uniq.fastq
rm ${LibName}_Aligned.Nuc_mRNA.noDups.uniq.fastq
fastqc ${LibName}_Nuc_LSUrRNA_al.fastq
rm ${LibName}_Nuc_LSUrRNA_al.fastq
fastqc ${LibName}_Nuc_tRNA_al.fastq
rm ${LibName}_Nuc_tRNA_al.fastq

# Get read lengths from fastq report (this will count soft-clipped bases, not ideal for footprints)
list='Aligned.Mito_mRNA.noDups.uniq Aligned.Nuc_mRNA.noDups.uniq Nuc_LSUrRNA_al Nuc_tRNA_al'
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

list='Aligned.Mito_mRNA.noDups Aligned.Nuc_mRNA.noDups'
for type in $list
do

# convert from .bam to .sam
samtools view -@ 3 -h -O SAM -o ${LibName}_${type}.sam ${LibName}_${type}.bam

# remove soft clipped bases from reads
java -jar ../1_AlignData/Scripts/biostar84452.jar -o ${LibName}_${type}.noSoft.sam ${LibName}_${type}.sam
rm ${LibName}_${type}.sam
# /n/groups/churchman/bms36/programs/jvarkit/dist/biostar84452.jar

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
samtools index ${LibName}_Aligned_noDups.bam

## 5' bedGraph (bedtools genomeCoverageBed)
# plus
genomeCoverageBed -bga -5 -strand + -trackline -ibam ${LibName}_Aligned_noDups.bam > ${LibName}_noDups.5p.plus.bedGraph
# minus
genomeCoverageBed -bga -5 -strand - -trackline -ibam ${LibName}_Aligned_noDups.bam > ${LibName}_noDups.5p.minus.bedGraph

# index files
igvtools toTDF ${LibName}_noDups.5p.plus.bedGraph ${LibName}_noDups.5p.plus.tdf ../Annotations/hg38.chrom.sizes
igvtools toTDF ${LibName}_noDups.5p.minus.bedGraph ${LibName}_noDups.5p.minus.tdf ../Annotations/hg38.chrom.sizes

### for unique bedgraphs ### 
samtools view -@ 3 -q 255 -o ${LibName}_Aligned_noDups.unique.bam ${LibName}_Aligned_noDups.bam 
## sort by coordinate
samtools sort -@ 3 -o ${LibName}_Aligned_noDups.unique.sort.bam ${LibName}_Aligned_noDups.unique.bam
# index
samtools index ${LibName}_Aligned_noDups.unique.sort.bam 
## 5'
genomeCoverageBed -bga -5 -strand + -trackline -ibam ${LibName}_Aligned_noDups.unique.sort.bam > ${LibName}_noDups.unique.5p.plus.bedGraph
genomeCoverageBed -bga -5 -strand - -trackline -ibam ${LibName}_Aligned_noDups.unique.sort.bam > ${LibName}_noDups.unique.5p.minus.bedGraph
## Indexed files
igvtools toTDF ${LibName}_noDups.unique.5p.plus.bedGraph ${LibName}_noDups.unique.5p.plus.tdf ../Annotations/hg38.chrom.sizes
igvtools toTDF ${LibName}_noDups.unique.5p.minus.bedGraph ${LibName}_noDups.unique.5p.minus.tdf ../Annotations/hg38.chrom.sizes

