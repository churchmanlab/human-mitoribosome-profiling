#!/n/app/R/4.0.1/bin/Rscript
					


library(Rsubread)


# Get arguments
args <- commandArgs(trailingOnly = TRUE)
LibName <- args[1]
PWD <- args[2]
offset <- args[3]

file = paste0(PWD,'/../1_AlignData_CytoRP/',LibName,'_Aligned_noDups.bam')
# See https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts

# Want readExtension5 = -17 to get A site, but negative numbers not allowed. Instead first shift read downstream then use 5' readpos

DT_unique = featureCounts(file, annot.ext=paste0(PWD,'/../Annotations/GRCh38_ncRNAs.gtf'), readShiftType  = "downstream", readShiftSize = offset, read2pos = 5, isGTFAnnotationFile = TRUE, allowMultiOverlap = TRUE, useMetaFeatures = TRUE, countMultiMappingReads = FALSE, strandSpecific = 1, isPairedEnd = FALSE, nthreads = 4, GTF.featureType='CDS')

DT_multi = featureCounts(file, annot.ext=paste0(PWD,'/../Annotations/GRCh38_ncRNAs.gtf'), readShiftType  = "downstream", readShiftSize = offset, read2pos = 5, isGTFAnnotationFile = TRUE, allowMultiOverlap = TRUE, useMetaFeatures = TRUE, countMultiMappingReads = TRUE, strandSpecific = 1, isPairedEnd = FALSE, nthreads = 4, GTF.featureType='CDS')



write.table(DT_unique$counts, file=paste0(PWD, '/',LibName,'_featureCounts_allSizesAsite_CDS_unique_noDups.txt'), col.names=FALSE, sep=("\t"), quote=FALSE)
 
write.table(DT_multi$counts, file=paste0(PWD, '/',LibName,'_featureCounts_allSizesAsite_CDS_multi_noDups.txt'), col.names=FALSE, sep=("\t"), quote=FALSE)


# To get feature lengths, for RPK. All of these will be the same
write.table(DT_unique$annotation$Length, file=paste0(PWD,'/',LibName,'_featureCounts_CDS_Length.txt'), col.names=FALSE, sep=("\t"), quote=FALSE)




