#!/n/app/R/4.0.1/bin/Rscript

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
Exp <- args[1]
type <- args[2] # multi unique
geneSet <- args[3] # all mito


library(data.table)
library(rlist)

# The following parameters should not be changed for running on mito-specific reads. 
pseudo = 'exclude' # include exclude
ncRNA = 'exclude'

fileName = paste0(Exp, '_featureCounts_CDS_',type,'_noDups.txt')
path = paste0(getwd(),'/')

DT <- data.table(read.table(paste0(path, fileName), sep = '\t', quote = '', header=TRUE, stringsAsFactors = FALSE))

SampNum=ncol(DT)-2

samplenames = colnames(DT)[2:(SampNum+1)]

# Convert GeneID_stable to GeneID by removing the .xx
DT[, GeneID := lapply(strsplit(GeneID, '[.]'),'[[',1)]
DT[, GeneID := as.character(GeneID)]


# GeneID	GeneID_stable	TxptID	TxptID_stable	Desc	GeneName	maxLength
GeneConvTable <- data.table(read.table(paste0(path, 'Annotations/hENStoGeneName_maxLength.txt'), header = T, sep = '\t', quote = '', stringsAsFactors = FALSE))

GeneConvTable <- GeneConvTable[!grep('antisense', Desc)] 
GeneConvTableNoPseudo <- GeneConvTable[!grep('pseudogene', Desc)] 
GeneConvTableNoPseudo <- GeneConvTableNoPseudo[!grep('novel transcript', Desc)] 


# Merge by GeneID(stable), only adding the GeneName column from the conversion table
DT_geneNames <- merge(DT, GeneConvTable[, c('GeneID', 'GeneName')], by='GeneID', all = FALSE)
DT_geneNames_noPseudo <- merge(DT, GeneConvTableNoPseudo[, c('GeneID', 'GeneName')], by='GeneID', all = FALSE)

if (pseudo == 'include'){
hDT <- DT_geneNames_noPseudo
}
if (pseudo == 'exclude'){
hDT <- DT_geneNames_noPseudo
 }
 
if (ncRNA == 'exclude') {
# Remove genes that I know to be or overlap with ncRNA: MT-R.., MT-T.., RNU.., FBN1, MALAT1
hDT <- hDT[GeneName != 'FBN1' & GeneName != 'MALAT1' & GeneName != 'CTC1' & GeneName != 'U1' & GeneName != 'U2' & GeneName != 'U3' & GeneName != 'U4' & GeneName != 'U5' & GeneName != 'U6' & GeneName != 'U7' & GeneName != 'U8' & GeneName != 'AL928646.1' & GeneName != 'NEAT1']
hDT <- hDT[!grep('^MT-T', GeneName)]
hDT <- hDT[!grep('^MT-R', GeneName)]
hDT <- hDT[!grep('^RNY', GeneName)]
hDT <- hDT[!grep('^RNU', GeneName)]
hDT <- hDT[!grep('^Y_RNA', GeneName)]
hDT <- hDT[!grep('^Metazoa_SRP', GeneName)]
hDT <- hDT[!grep('^SNO', GeneName)]
hDT <- hDT[!grep('^MIR', GeneName)]
hDT <- hDT[!grep('^RNA5', GeneName)]
hDT <- hDT[!grep('^RN7', GeneName)]
hDT <- hDT[!grep('^RNV', GeneName)]
hDT <- hDT[!grep('^VTRNA', GeneName)]
hDT <- hDT[!grep('^LINC', GeneName)]
hDT <- hDT[!grep('^RNA28S$', GeneName)]
hDT <- hDT[!grep('^RNA18S$', GeneName)]
hDT <- hDT[!grep('^TRNAALL$', GeneName)]


}

if (ncRNA == 'include') {
hDT <- hDT
}


# Keep all genes
if (geneSet == 'all') {
hDT <- hDT
}

# Keep only mito genes
if (geneSet == 'mito') {
hDT <- hDT[grep('^MT-', GeneName)]
}
 



# ######## Calculate RPK, RPKM, and make tables ######### 

pernumber=1000
DT_RPK <- copy(hDT)
RPKcols <- list()
for (i in c(1:SampNum)) {
	RPKcols <- list.append(RPKcols, round(hDT[[samplenames[i]]]/hDT$Length*pernumber, digits=2))
	}
DT_RPK[, (samplenames) := RPKcols]


# Change order of columns so GeneName is first
colorder = c('GeneName', samplenames, 'Length', 'GeneID')
setcolorder(hDT, colorder)
setcolorder(DT_RPK, colorder)



# Write tables to file
if (pseudo == 'include' & ncRNA == 'include') {
write.table(hDT, file=paste0(path, fileName, '_', geneSet,'_readCount_all.txt'), row.names=FALSE, sep=("\t"), quote=FALSE)
write.table(DT_RPK, file=paste0(path, fileName, '_', geneSet, '_RPK_all.txt'), row.names=FALSE, sep=("\t"), quote=FALSE)
}

if (pseudo == 'exclude' & ncRNA == 'include') {
write.table(hDT, file=paste0(path, fileName, '_', geneSet, '_readCount_noPseudo.txt'), row.names=FALSE, sep=("\t"), quote=FALSE)
write.table(DT_RPK, file=paste0(path, fileName, '_', geneSet, '_RPK_noPseudo.txt'), row.names=FALSE,  sep=("\t"), quote=FALSE)
}

if (pseudo == 'exclude' & ncRNA == 'exclude') {
write.table(hDT, file=paste0(path, fileName, '_', geneSet, '_readCount.txt'), row.names=FALSE, sep=("\t"), quote=FALSE)
write.table(DT_RPK, file=paste0(path, fileName, '_', geneSet, '_RPK.txt'), row.names=FALSE,  sep=("\t"), quote=FALSE)
}

# source('./5_CountReadsOnFeatures_RNAseq/Scripts/AddGeneName_RPK_hRNAseq.R')

