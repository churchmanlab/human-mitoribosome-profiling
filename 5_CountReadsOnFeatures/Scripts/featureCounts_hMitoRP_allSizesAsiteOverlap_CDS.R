#!/n/app/R/4.0.1/bin/Rscript

library(Rsubread)


# Get arguments
args <- commandArgs(trailingOnly = TRUE)
LibName <- args[1]
PWD <- args[2]
offset <- args[3]

file = paste0(PWD,'/',LibName,'_Aligned.Mito_mRNA.noDups.bam')
# See https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts

# For mito genes, remove first 6 codons and last 5 codons, as well as all overlapping regions from gtf
# Want readExtension3 = -15 to get A site, but negative numbers not allowed. Instead shift reads upstream by the offset then use 3' readpos

DT_unique = featureCounts(file, annot.ext=paste0(PWD,'/Annotations/GRCh38_mitoCDSnoOverlapPlus6minus5.gtf', readShiftType  = "upstream", readShiftSize = offset, read2pos = 3, isGTFAnnotationFile = TRUE, allowMultiOverlap = TRUE, useMetaFeatures = TRUE, countMultiMappingReads = FALSE, strandSpecific = 1, isPairedEnd = FALSE, nthreads = 4, GTF.featureType='CDS')

DT_multi = featureCounts(file, annot.ext=paste0(PWD,'/Annotations/GRCh38_mitoCDSnoOverlapPlus6minus5.gtf', readShiftType  = "upstream", readShiftSize = offset, read2pos = 3, isGTFAnnotationFile = TRUE, allowMultiOverlap = TRUE, useMetaFeatures = TRUE, countMultiMappingReads = TRUE, strandSpecific = 1, isPairedEnd = FALSE, nthreads = 4, GTF.featureType='CDS')



write.table(DT_unique$counts, file=paste0(PWD, '/',LibName,'_mito_featureCounts_allSizesAsite_CDSignore1stlast_unique_noDups.txt'), col.names=FALSE, sep=("\t"), quote=FALSE)
 
write.table(DT_multi$counts, file=paste0(PWD, '/',LibName,'_mito_featureCounts_allSizesAsite_CDSignore1stlast_multi_noDups.txt'), col.names=FALSE, sep=("\t"), quote=FALSE)


# To get feature lengths, for RPK. All of these will be the same
write.table(DT_unique$annotation$Length, file=paste0(PWD,'/',LibName,'_mito_featureCounts_CDSignore1stlast_Length.txt'), col.names=FALSE, sep=("\t"), quote=FALSE)




