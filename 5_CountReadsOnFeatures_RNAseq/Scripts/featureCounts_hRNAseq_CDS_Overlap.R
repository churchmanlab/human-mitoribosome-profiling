#!/n/app/R/3.5.1/bin/Rscript
					


library(Rsubread)


# Get arguments
args <- commandArgs(trailingOnly = TRUE)
LibName <- args[1]
PWD <- args[2]


file = paste0(PWD,'/../1_AlignData_RNAseq/',LibName,'_Aligned_noDups.bam')

# See https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts

gtffile=paste0(PWD,'/../Annotations/GRCh38_ncRNAs_mM17_merge_noMitoOverlap.gtf')

# Using CDS
DT_uniqueCDS = featureCounts(file, annot.ext=gtffile, isGTFAnnotationFile = TRUE, useMetaFeatures = TRUE, allowMultiOverlap = TRUE, minOverlap=22, countMultiMappingReads = FALSE, GTF.featureType='CDS', strandSpecific = 1, isPairedEnd = FALSE, nthreads = 4)

DT_multiCDS = featureCounts(file, annot.ext=gtffile, isGTFAnnotationFile = TRUE, useMetaFeatures = TRUE, allowMultiOverlap = TRUE, minOverlap=22, countMultiMappingReads = TRUE, GTF.featureType='CDS', strandSpecific = 1, isPairedEnd = FALSE, nthreads = 4)


write.table(DT_uniqueCDS$counts, file=paste0(PWD, '/',LibName,'_featureCounts_CDS_unique_noDups.txt'), col.names=FALSE, sep=("\t"), quote=FALSE)
# 
write.table(DT_multiCDS$counts, file=paste0(PWD, '/',LibName,'_featureCounts_CDS_multi_noDups.txt'), col.names=FALSE, sep=("\t"), quote=FALSE)

# To get feature lengths, for RPK. All of these will be the same
write.table(DT_uniqueCDS$annotation$Length, file=paste0(PWD,'/',LibName,'_featureCounts_CDS_Length.txt'), col.names=FALSE, sep=("\t"), quote=FALSE)




