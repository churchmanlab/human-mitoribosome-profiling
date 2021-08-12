#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-10:00
#SBATCH -p short
#SBATCH --mem=100G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<youremailaddress>

### WRITTEN:		11/19/2020 by Mary Couvillion

### USE: 				This script will trim reads, map to the genome, remove PCR      ###					duplicates, 
###					and produce files for visualizing on IGV
###					sbatch -e logs/1_AlignData_RNAseq_LibName.err -o logs/1_AlignData_RNAseq_LibName.log ProcessFASTQ_hMitoRP.sh
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

STARindexPATH="../Annotations/STARindex/GRCh38_ncRNAs_mM17_merge_RNAseq"
# STARindexPATH="/n/groups/churchman/hMitoRP/SequenceFiles/GRCh38_ncRNAs_mM17_merge_RNASeq"

# Input from command line. LibName is a short descriptive name you want as prefix to all files
LibName=$1
InputFASTQ=$2
UMI=$3 # 3p6, 3p6_5p4, 3p10_5p4 Truseq (for Truseq use read2 (R2) if available so it will be forward strand)

cd ./1_AlignData_RNAseq

########## CLEAN READS #################
########################################
if [ "${UMI}" = "3p10_5p4" ]
then
	# Cutadapt (minimum is 35-10-4-1=20 after trimming --untrimmed-output writes reads with no adaptor found INSTEAD of writing to final output file. Do this here so that we keep only reads that have adapter, which means we know where the barcode is
	cutadapt -e 0.2 --untrimmed-output ${LibName}_noAdaptor.fastq -a CTGTAGGCACCATCAAT -m 35 ${InputFASTQ} -o ${LibName}_trimmed.fastq
	rm ${LibName}_noAdaptor.fastq
	###### To extract barcode ######
	python ../1_AlignData/Scripts/extractMolecularBarcodeFrom3pr10AND5pr4.py ${LibName}_trimmed.fastq ${LibName}_Cleaned.fastq
	rm ${LibName}_trimmed.fastq
elif [ "${UMI}" = "3p6_5p4" ]
then
	# Cutadapt (minimum is 31-6-4=21 after trimming --untrimmed-output writes reads with no adaptor found INSTEAD of writing to final output file. Do this here so that we keep only reads that have adapter, which means we know where the barcode is
	cutadapt -e 0.2 --untrimmed-output ${LibName}_noAdaptor.fastq -a CTGTAGGCACCATCAAT -m 31 ${InputFASTQ} -o ${LibName}_trimmed.fastq
	rm ${LibName}_noAdaptor.fastq
	###### To extract barcode ######
	python ../1_AlignData/Scripts/extractMolecularBarcodeFrom3pr6AND5pr.py ${LibName}_trimmed.fastq ${LibName}_Cleaned.fastq
	rm ${LibName}_trimmed.fastq

elif [ "${UMI}" = "3p6" ]
then 
# 	Cutadapt (minimum is 26-6=20 after trimming --untrimmed-output writes reads with no adaptor found INSTEAD of writing to final output file. Do this here so that we keep only reads that have adapter, which means we know where the barcode is
	cutadapt -e 0.2 --untrimmed-output ${LibName}_noAdaptor.fastq -a CTGTAGGCACCATCAAT -m 26 ${InputFASTQ} -o ${LibName}_trimmed.fastq
	rm ${LibName}_noAdaptor.fastq
	##### To extract barcode ######
	python ../1_AlignData/Scripts/extractMolecularBarcodeFrom3pr.py ${LibName}_trimmed.fastq ${LibName}_Cleaned.fastq ${LibName}_Barcodes.txt ${LibName}_Ligation.txt
	rm ${LibName}_trimmed.fastq
	rm ${LibName}_Barcodes.txt
	rm ${LibName}_Ligation.txt

elif [ "${UMI}" = "Truseq" ]
then 
	# Need newer cutadapt to handle paired end reads, and need python3 for that
	module load python/3.7.4
	module load cutadapt/2.5
	cutadapt -j 4 -m 30 -o ${LibName}_R1_Cleaned.fastq -p ${LibName}_R2_Cleaned.fastq ${InputFASTQ}_R2_001.fastq.gz ${InputFASTQ}_R1_001.fastq.gz  # Reverse order of R1 and R2 so orientation will be forward
	/n/groups/churchman/mc348/TimelapseSeq/Scripts/NGmerge-master/NGmerge -1 ${LibName}_R1_Cleaned.fastq -2 ${LibName}_R2_Cleaned.fastq -n 4 -o ${LibName}_Cleaned.fastq -p 0.1 -m 10 -v -g
	# Reload python2
	module load python/2.7.12

fi

############## FASTQC ################
# Make an html report of quality and characteristics of trimmed reads
fastqc ${LibName}_Cleaned.fastq


############# MAP ###########

# Cyto SSU rRNA
bowtie -v 1 -5 1 -S ../Annotations/bowtieindex/rRNA_SSU ${LibName}_Cleaned.fastq ${LibName}_Nuc_SSUrRNA_almts.sam --un ${LibName}_Nuc_SSUrRNA_un.fastq
rm ${LibName}_Nuc_SSUrRNA_almts.sam

# Cyto LSU rRNA
bowtie -v 1 -5 1 -S ../Annotations/bowtieindex/rRNA_LSU ${LibName}_Nuc_SSUrRNA_un.fastq ${LibName}_Nuc_LSUrRNA_almts.sam --un ${LibName}_Nuc_LSUrRNA_un.fastq
rm ${LibName}_Nuc_SSUrRNA_un.fastq
rm ${LibName}_Nuc_LSUrRNA_almts.sam

# Mito rRNA
bowtie -v 1 -5 1 -S ../Annotations/bowtieindex/hMito_rRNA ${LibName}_Nuc_LSUrRNA_un.fastq ${LibName}_hMito_rRNA_almts.sam  --un ${LibName}_hMito_rRNA_un.fastq
rm ${LibName}_Nuc_LSUrRNA_un.fastq
rm ${LibName}_hMito_rRNA_almts.sam

# Mouse mito mRNA
# Here reads are already filtered for read lengths >= 20
bowtie -v 0 -5 1 -S ../Annotations/bowtieindex/mMito_mRNAs ${LibName}_hMito_rRNA_un.fastq ${LibName}_mMito_mRNAs_0mm.sam --al ${LibName}_mMito_mRNAs_0mm_al.fastq

# Cyto tRNA
bowtie -v 1 -3 3 -5 1 -S../Annotations/bowtieindex/tRNAs ${LibName}_hMito_rRNA_un.fastq ${LibName}_Nuc_tRNA_almts.sam --un ${LibName}_Nuc_tRNA_un.fastq
rm ${LibName}_Nuc_tRNA_almts.sam

# Mito tRNA
bowtie -v 1 -3 3 -5 1 -S ../Annotations/bowtieindex/hMito_tRNA ${LibName}_Nuc_tRNA_un.fastq ${LibName}_hMito_tRNA_almts.sam --al ${LibName}_hMito_tRNA_al.fastq
rm ${LibName}_hMito_tRNA_almts.sam


STAR --runThreadN 4 --genomeDir $STARindexPATH --readFilesIn ${LibName}_hMito_rRNA_un.fastq --outFileNamePrefix ${LibName}_ --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outFilterMultimapNmax 100
rm -r ${LibName}__STARtmp
rm ${LibName}_SJ.out.tab

# Map to oGAB
bowtie -v 1 -5 1 -S ../Annotations/bowtieindex/oGAB ${LibName}_Unmapped.out.mate1 ${LibName}_oGAB_almts.sam
rm ${LibName}_oGAB_almts.sam

# ############ Get counts from mapping ###########
# ## Make list of mapping references, in order
echo Input > ${LibName}_MapList.txt
echo Cyto SSU >> ${LibName}_MapList.txt
echo Cyto LSU >> ${LibName}_MapList.txt
echo Mito rRNA >> ${LibName}_MapList.txt
echo Mouse 0mm >> ${LibName}_MapList.txt
echo Cyto tRNA >> ${LibName}_MapList.txt
echo Mito tRNA >> ${LibName}_MapList.txt
echo oGAB >> ${LibName}_MapList.txt
# Grab mapped numbers from bowtie error file
grep -m 1 'reads processed:' ../logs/1_AlignData_RNAseq_${LibName}.err > ${LibName}_bowtieCounts.txt
grep 'reads with at least one reported alignment:' ../logs/1_AlignData_RNAseq_${LibName}.err >> ${LibName}_bowtieCounts.txt
paste <(awk '{print $0}' ${LibName}_MapList.txt) <(awk -F ': ' '{print $2}' ${LibName}_bowtieCounts.txt) > ${LibName}_Counts.txt
# Get mapped reads from STAR log file
echo Mapped from STAR >> ${LibName}_Counts.txt
egrep 'Uniquely mapped reads number | Number of reads mapped to multiple loci' ${LibName}_Log.final.out | awk -F '|\t' '{print $2}' | awk '{sum += $1} END {print sum}' >> ${LibName}_Counts.txt

# Get unmapped reads from STAR log file
echo Unmapped from STAR \(oGAB counts should be subtracted from this\): >> ${LibName}_Counts.txt
egrep 'Number of reads mapped to too many loci | Number of reads unmapped: too short | Number of reads unmapped: other' ${LibName}_Log.final.out | awk -F '|\t' '{print $2}' | awk '{sum += $1} END {print sum}' >> ${LibName}_Counts.txt
 
rm ${LibName}_MapList.txt
rm ${LibName}_bowtieCounts.txt

# ############ REMOVE PCR DUPLICATES ###########

## Convert mMito mRNA sam to bam
samtools view -@ 3 -b ${LibName}_mMito_mRNAs_0mm.sam -o ${LibName}_mMito_mRNAs_0mm.bam 
## Sort
samtools sort -@ 3 ${LibName}_mMito_mRNAs_0mm.bam -o ${LibName}_mMito_mRNAs_0mm.sort.bam

## Make BAM index (samtools)
samtools index -@ 3 ${LibName}_Aligned.sortedByCoord.out.bam
samtools index -@ 3 ${LibName}_mMito_mRNAs_0mm.sort.bam

## Removal
if [ "${UMI}" = "Truseq" ]
then
	# Just rename to noDups so it can go through rest of script
	cp ${LibName}_Aligned.sortedByCoord.out.bam ${LibName}_Aligned_noDups.bam
	cp ${LibName}_mMito_mRNAs_0mm.sort.bam ${LibName}_mMito_mRNA_noDups.bam
else
	python ../1_AlignData/Scripts/removePCRdupsFromBAM_MC20180913.py ${LibName}_Aligned.sortedByCoord.out.bam ${LibName}_Aligned_noDups.bam
	python ../1_AlignData/Scripts/removePCRdupsFromBAM_MC20180913.py ${LibName}_mMito_mRNAs_0mm.sort.bam ${LibName}_mMito_mRNA_noDups.bam
fi


############ MAKE MITO ONLY BAM FILE ###########

## Make sure they are sorted by coordinate
samtools sort -@ 3 -o ${LibName}_Aligned_noDups.sort.bam ${LibName}_Aligned_noDups.bam

## Make index
samtools index -@ 3 ${LibName}_Aligned_noDups.sort.bam

## chrM only
samtools view -@ 3 -b ${LibName}_Aligned_noDups.sort.bam chrM -o ${LibName}_Aligned.Mito.noDups.bam 
## chrM mRNA only
samtools view -@ 3 -b -o ${LibName}_Aligned.Mito_mRNA.noDups.bam -L ../Annotations/BEDfiles/hMito_mRNAsANDsurrounding.bed ${LibName}_Aligned_noDups.sort.bam
# And hNuc mRNA only file
samtools view -@ 3 -b -o ${LibName}_Aligned.Nuc_mRNA.noDups.bam -L ../Annotations/BEDfiles/GENCODE_hg38_proteincoding.sorted.bed ${LibName}_Aligned_noDups.sort.bam

# Also for non-deduplicated file, to get counts
samtools view -@ 3 -b -o ${LibName}_hMito_mRNA_Aligned.sortedByCoord.out.bam -L ../Annotations/BEDfiles/hMito_mRNAsANDsurrounding.bed ${LibName}_Aligned.sortedByCoord.out.bam


 ############ Get counts from noDups ###########
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

# ############## FASTQC ################

# First convert bam to fastq
bedtools bamtofastq -i ${LibName}_Aligned.Mito_mRNA.noDups.bam -fq ${LibName}_Aligned.Mito_mRNA.noDups.fastq
bedtools bamtofastq -i ${LibName}_Aligned.Nuc_mRNA.noDups.bam -fq ${LibName}_Aligned.Nuc_mRNA.noDups.fastq
## Keep only unique lines in fastq (see https://edwards.sdsu.edu/research/sorting-fastq-files-by-their-sequence-identifiers/)
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



########## MAKE BEDGRAPHS FOR VIEWING ON IGV ###########

### Make bam index
samtools index ${LibName}_Aligned_noDups.bam

## 5' bedGraph (bedtools genomeCoverageBed)
# plus
genomeCoverageBed -bga -5 -strand + -trackline -ibam ${LibName}_Aligned.Mito_mRNA.noDups.bam > ${LibName}_Mito_mRNA.noDups.5p.plus.bedGraph
# minus
genomeCoverageBed -bga -5 -strand - -trackline -ibam ${LibName}_Aligned.Mito_mRNA.noDups.bam > ${LibName}_Mito_mRNA.noDups.5p.minus.bedGraph


## Coverage
genomeCoverageBed -bga -strand + -split -trackline -trackopts 'color=0,100,50' -ibam ${LibName}_Aligned_noDups.sort.bam > ${LibName}_Aligned_noDups_CovPlus.bedGraph
genomeCoverageBed -bga -strand - -split -trackline -trackopts 'color=0,100,50' -ibam ${LibName}_Aligned_noDups.sort.bam > ${LibName}_Aligned_noDups_CovMinus.bedGraph
## Indexed files
igvtools toTDF ${LibName}_Aligned_noDups_CovPlus.bedGraph ${LibName}_Aligned_noDups_CovPlus.tdf ../Annotations/hg38.chrom.sizes
igvtools toTDF ${LibName}_Aligned_noDups_CovMinus.bedGraph ${LibName}_Aligned_noDups_CovMinus.tdf ../Annotations/hg38.chrom.sizes

### for unique bedgraphs ### 
samtools view -@ 3 -q 255 -o ${LibName}_Aligned_noDups.unique.bam ${LibName}_Aligned_noDups.bam 
## sort by coordinate
samtools sort -@ 3 -o ${LibName}_Aligned_noDups.unique.sort.bam ${LibName}_Aligned_noDups.unique.bam
# index
samtools index ${LibName}_Aligned_noDups.unique.sort.bam 
## Coverage
genomeCoverageBed -bga -strand + -split -trackline -trackopts 'color=0,100,50' -ibam ${LibName}_Aligned_noDups.unique.sort.bam > ${LibName}_uniqueAligned_noDups_CovPlus.bedGraph
genomeCoverageBed -bga -strand - -split -trackline -trackopts 'color=0,100,50' -ibam ${LibName}_Aligned_noDups.unique.sort.bam > ${LibName}_uniqueAligned_noDups_CovMinus.bedGraph
## Indexed files
igvtools toTDF ${LibName}_uniqueAligned_noDups_CovPlus.bedGraph ${LibName}    ../Annotations/hg38.chrom.sizes
igvtools toTDF ${LibName}_uniqueAligned_noDups_CovMinus.bedGraph ${LibName}_uniqueAligned_noDups_CovMinus.tdf ../Annotations/hg38.chrom.sizes



##################################################
### After running for all libs, combine all

# Experiment="FibroNT" # OXPHOSinh GalShift MyoDiff1 HEK_RNAseq2 MyoDiff2
# # # # 
# # # # # Libs=("Glu_1" "Glu_2" "Gal_30m_1" "Gal_30m_2" "Gal_90m_1" "Gal_90m_2" "Gal_4h_1" "Gal_4h_2" "Gal_72h_1" "Gal_72h_2")
# # # # Libs=("DMSO_1" "ROT_1nM" "ROT_5nM" "ROT_10nM" "ROT_100nM" "AA_1nM" "AA_5nM" "AA_10nM" "AA_100nM" )
# # # Libs=("NonDiff_1" "NonDiff_2" "Diff_4h_1" "Diff_4h_2" "Diff_2d_1" "Diff_2d_2")
# # Libs=("WT1_2" "WT2_2" "KO1_2" "KO2_2" "KOplWT1_2" "KOplWT2_2")
# Libs=("FibroNT")
# list='multi unique'
# # 
# for type in $list
# do
# 
# paste <(awk '{print $1"\t"$2}' ${Libs[0]}_featureCounts_${type}.txt)  <(awk '{print $2}' ${Libs[0]}_featureCounts_Length.txt) > ${Experiment}_featureCounts_${type}_noDups.txt #<(awk '{print $2}' ${Libs[1]}_featureCounts_${type}.txt) <(awk '{print $2}' ${Libs[2]}_featureCounts_${type}.txt) <(awk '{print $2}' ${Libs[3]}_featureCounts_${type}.txt) <(awk '{print $2}' ${Libs[4]}_featureCounts_${type}.txt) <(awk '{print $2}' ${Libs[5]}_featureCounts_${type}.txt) <(awk '{print $2}' ${Libs[9]}_featureCounts_allSizes5p_ignore1stLast_multi_noDups.txt)
# # # 
# # # rm *1_featureCounts_allSizes5p_ignore1stLast_multi_noDups.txt
# # # rm *2_featureCounts_allSizes5p_ignore1stLast_multi_noDups.txt
# # # 
# 
# done






# Experiment="HEK_RNAseq2"
# Samps="WT1_2 WT2_2 KO1_2 KO2_2 KOplWT1_2 KOplWT2_2"
# 
# echo All Counts > ${Experiment}_Counts.txt
# for Samp in $Samps
# do
# cat ${Samp}_Counts.txt >> ${Experiment}_Counts.txt
# done


# # ########## For quick batch submission ##########
# Libs=("WT1_2" "WT2_2" "KO1_2" "KO2_2" "KOplWT1_2" "KOplWT2_2")
# fastqs=("WT-R1_S1_R1_001.fastq.gz" "WT-R2_S2_R1_001.fastq.gz" "KO-R1_S3_R1_001.fastq.gz" "KO-R2_S4_R1_001.fastq.gz" "KOplusWT-R1_S5_R1_001.fastq.gz" "KOplusWT-R2_S6_R1_001.fastq.gz")
# Libs=("NT_1" "NT_2" "CAP4h_1" "CAP4h_2" "CAP72h_1" "CAP72h_2")
# fastqs=("RNA_seq_Cyto1_S3_R1_001.fastq.gz" "RNA_seq_Cyto2_S4_R1_001.fastq.gz" "RNA_seq_Cyto3_S5_R1_001.fastq.gz" " RNA_seq_Cyto4_S6_R1_001.fastq.gz" "RNA_seq_Cyto5_S7_R1_001.fastq.gz" "RNA_seq_Cyto6_S8_R1_001.fastq.gz")

Libs=("Tot_1" "Tot_2")
fastqs=("Total-1_S3_L001" "Total-2_S4_L001")
for i in {0..1}
do
sbatch -e logs/Map_${Libs[i]}.err -o logs/Map_${Libs[i]}.log /n/groups/churchman/mc348/hMitoRP/PanAnalysis/Scripts/ProcessFASTQraw_hRNAseq.sh ${Libs[i]} ${fastqs[i]} Truseq
done


# Libs=("FibroNT")
# fastqs=("Exp4_NT.fastq.gz")
# for i in {0..0}
# do
# sbatch -e logs/Map2_${Libs[i]}.err -o logs/Map2_${Libs[i]}.log /n/groups/churchman/mc348/hMitoRP/PanAnalysis/Scripts/ProcessFASTQraw_hRNAseq.sh ${Libs[i]} ${fastqs[i]} 3p10_5p4
# done
