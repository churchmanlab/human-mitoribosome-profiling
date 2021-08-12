#!/n/app/R/4.0.1/bin/Rscript


library(data.table)
library('plotrix')
library(rlist)
library(scales)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
data <- args[1] # RP RNAseq
mitoExperiment <- args[2]
cytoExperiment <- args[3]
RNAseqExperiment <- args[4]
soi <- unlist(strsplit(args[5],","))
soiName <- args[6]

# data='RP' # RP RNAseq
# mitoExperiment = 'hMitoRP_HEK'
# cytoExperiment = 'CytoRP_HEK'
# RNAseqExperiment = 'RNAseq_HEK'
# soi = c('WT1', 'WT2')
# soiName = 'WT'

path = paste0(getwd(),'/')




############## First, get mito data from mitoRP or RNAseq ##########################
######################### (mito goes through line 211) ########################


if (data == 'RP') {
Experiment = mitoExperiment
modifier = 'featureCounts_allSizesAsite_CDSignore1stlast_multi_noDups_mito_RPK_noPseudo_ncRNA'
}
if (data == 'RNAseq') {
Experiment = RNAseqExperiment
modifier = 'featureCounts_CDS_multi_noDups_mito_RPK_noPseudo_ncRNA'
}


file = paste0(path, Experiment, '_', modifier,'.txt')
DT <- data.table(read.table(file, header=TRUE, stringsAsFactors=FALSE))


# Normalize to library size (tp10k here)
pernumber = 10000

#Remove unnecessary columns
DT[, c('Length', 'GeneID') := NULL]


# Put rows in order by complexes
to_ord = c('MT-ND1','MT-ND2','^MT-ND4L$','^MT-ND4$','MT-ND6','MT-ND5','MT-ND3','MT-CYB','MT-CO1','MT-CO2','MT-CO3','MT-ATP8','MT-ATP6')

ord_index = c()
for (gene in to_ord) {
ord_index = append(ord_index, which(grepl(gene, DT$GeneName)))
}
setorder(DT[, .r := order(ord_index)], .r)[, .r := NULL]

SampNum = ncol(DT) - 1
sampCols = c(2:(SampNum + 1))

# Get x axis labels
genes=DT$GeneName

# Calculate normalized values
newCols <- list()
for (i in c(1:SampNum)) {
	newCols <- list.append(newCols, DT[[sampCols[i]]]/sum(DT[[sampCols[i]]])*pernumber)
	}
RPKX_DT <- copy(DT)
RPKX_DT[, (sampCols) := newCols]

# Get sample names for legend
samples=colnames(RPKX_DT[,sampCols, with=FALSE])

genes = c('ND1', 'ND2','ND4L','ND4','ND6','ND5','ND3','CYB','CO1','CO2','CO3','ATP8','ATP6')

Color_1 = 'darkgoldenrod1'
Color_2 = 'dodgerblue'
Color_3 = 'indianred3'
Color_4 = 'forestgreen'

soiNum = length(soi)

# Make barplot
pdf(paste0(path, '7_OXPHOSstoichiometry/OXPHOSbarplots/', Experiment, '_', soiName,'_', modifier,'_barplot_mito.pdf'), width = 5, height = 4, pointsize=11)

# Set axis labels
if (data == 'RP') {label = 'Relative synthesis'}
if (data == 'RNAseq') {label = 'Relative RNA abundance'}

barplot(t(RPKX_DT[, ..soi]),beside=TRUE,names.arg = genes, col = c(rep(Color_1,(soiNum)*7), rep(Color_2,(soiNum)*1),  rep(Color_3,(soiNum)*3), rep(Color_4, (soiNum)*2)), border=NA, ylab = label, cex.axis=.9, cex.names = 1, las=2) 

dev.off()

# Write values from barplots

printDT <- cbind(RPKX_DT$GeneName, round(RPKX_DT[, ..soi], 1))
colnames(printDT)[1] <- 'GeneName'
write.table(printDT, file=paste0(path, '7_OXPHOSstoichiometry/OXPHOSbarplots/', Experiment, '_', soiName,'_', modifier,'_barplotData_mito.txt'), row.names=FALSE, sep=("\t"), quote=FALSE)



####################################################################################
############################## Complex averages ####################################
CIave = mean(as.matrix(RPKX_DT[1:7][,..soi]))
CIIIave = mean(as.matrix(RPKX_DT[8][,..soi]))
CIVave = mean(as.matrix(RPKX_DT[9:11][,..soi]))
CVave = mean(as.matrix(RPKX_DT[12:13][,..soi]))

# Get ranges
CIaves = c()
CIIIaves = c()
CIVaves = c()
CVaves = c()
for (i in c(1:length(soi))) {
CIaves = c(CIaves, mean(as.matrix(RPKX_DT[1:7][,..soi][,..i])))
CIIIaves = c(CIIIaves, mean(as.matrix(RPKX_DT[8][,..soi][,..i])))
CIVaves = c(CIVaves, mean(as.matrix(RPKX_DT[9:11][,..soi][,..i])))
CVaves = c(CVaves, mean(as.matrix(RPKX_DT[12:13][,..soi][,..i])))
}
CIrange = max(CIaves)-min(CIaves)
CIIIrange = max(CIIIaves)-min(CIIIaves)
CIVrange = max(CIVaves)-min(CIVaves)
CVrange = max(CVaves)-min(CVaves)


# Store values
Mitosoi = soi
Mitodata = data
MitoExpt = Experiment
MitoAves = c(CIave, CIIIave, CIVave, CVave)
MitoRanges = c(CIrange, CIIIrange, CIVrange, CVrange)


####################################################################################
####################################################################################



####################################################################################
############################ Cumulative fraction plots #############################

# Set data table to RPKX_DT above
dt = copy(RPKX_DT)
# Take means across soi
MeanRPKM <- rowMeans(dt[,..soi])
# Add means
dt[, soiMean := MeanRPKM]

# Get stoichiometries
CImitoStoich = rep(1, 7)
CISmitoStoich = rep(1, 7)
CIIImitoStoich = rep(1, 1)
CIVmitoStoich = rep(1, 3)
CVmitoStoich = rep(1, 2)

# Add stoichiometry
dt[, Stoich := c(CImitoStoich, CIIImitoStoich, CIVmitoStoich, CVmitoStoich)]

# Make columns of interest numeric
dt$soiMean <- as.numeric(dt$soiMean)

# Normalize synth rate to stoichiometry
dt[, NormSynthRate := soiMean/Stoich]

# Define each complex
CImito = c('MT-ND1','MT-ND2','MT-ND4L','MT-ND4','MT-ND6','MT-ND5','MT-ND3')
CIIImito = c('MT-CYB')
CIVmito = c('MT-CO1','MT-CO2','MT-CO3')
CVmito = c('MT-ATP8','MT-ATP6')

# Find median value for each complex
CImed <- median(dt[GeneName %in% CImito]$NormSynthRate, na.rm=TRUE)
CIIImed <- median(dt[GeneName %in% CIIImito]$NormSynthRate, na.rm=TRUE)
CIVmed <- median(dt[GeneName %in% CIVmito]$NormSynthRate, na.rm=TRUE)
CVmed <- median(dt[GeneName %in% CVmito]$NormSynthRate, na.rm=TRUE)

# Find rate limiting value for each complex

CIrl <- min(dt[GeneName %in% CImito]$NormSynthRate, na.rm=TRUE)
CIIIrl <- min(dt[GeneName %in% CIIImito]$NormSynthRate, na.rm=TRUE)
CIVrl <- min(dt[GeneName %in% CIVmito]$NormSynthRate, na.rm=TRUE)
CVrl <- min(dt[GeneName %in% CVmito]$NormSynthRate, na.rm=TRUE)

# Add medians to DT
dt[, Median := c(rep(CImed, 7), rep(CIIImed, 1), rep(CIVmed, 3), rep(CVmed, 2))]

# Add limiting rate to DT
dt[, rl := c(rep(CIrl, 7), rep(CIIIrl, 1), rep(CIVrl, 3), rep(CVrl, 2))]

# Fold diff from median
dt[, log2med := log(NormSynthRate/Median, 2)]

# Fold diff from rate-limiting
dt[, log2rl := log(NormSynthRate/rl, 2)]


pdf(paste0(path, '7_OXPHOSstoichiometry/CumFracPlots/', Experiment, '_', soiName,'_', modifier,'_CumFractionSubunitStoich_mito.pdf'), width=3, height=2.5, pointsize=11)
      
mitocomplexes=c(list(CImito), list(CIIImito), list(CIVmito), list(CVmito))  
colors=c(Color_1, Color_2, 'red', Color_4)


# Plot for mito
compartment='Mito'
complexes=mitocomplexes

# Spread around rl
plot(1, xlim=c(0,5.2), ylim=c(0,1), pch=16, type='n', xlab='Synth rate relative to rate limiting (log2)', ylab='Cumulative fraction of subunits', col=colors[i], lwd=2, main=paste0(Experiment, ' ', compartment, soiName), yaxt='n')
# Make shaded box at 2-fold
polygon(c(0,0,1,1), c(0,1,1,0), col = alpha('black', .1), border = NA)
# Draw lines for CI, CIV, CV (not CIII since there's only 1
for(i in c(1,3,4)) { 
x=dt[GeneName %in% complexes[[i]]]$log2rl
x=sort(x)
x=c(x[1], x) # To make the stair start at the 0 line
y=c(0:(length(x)-1))
y=y/max(y)
points(x,y, pch=16, type='s', col=colors[i], lwd=2)
}
axis(2, at=c(0,.5, 1), labels=c(0, .5, 1))
legend('bottomright', legend=c('ComplexI', 'ComplexIV', 'ComplexV'), text.col=colors[c(1,3,4)], bty='n', cex=.7)

dev.off()




############## Now, get cyto data from cytoRP or RNAseq ##########################
######################### And plot correlations ########################

if (data == 'RP') {
Experiment = cytoExperiment
RPKmodifier = 'featureCounts_allSizesAsite_CDS_multi_noDups_all_RPK_noPseudo_ncRNA'
RCmodifier = 'featureCounts_allSizesAsite_CDS_multi_noDups_all_readCount_noPseudo_ncRNA'
}
if (data == 'RNAseq') {
Experiment = RNAseqExperiment
RPKmodifier = 'featureCounts_CDS_multi_noDups_all_RPK_noPseudo_ncRNA'
RCmodifier = 'featureCounts_CDS_multi_noDups_all_readCount_noPseudo_ncRNA'
}

file = paste0(path, Experiment, '_', RPKmodifier,'.txt')
rcfile = paste0(path, Experiment, '_', RCmodifier,'.txt')
DT <- data.table(read.table(file, header=TRUE, stringsAsFactors=FALSE))
rcDT <- data.table(read.table(rcfile, header=TRUE, stringsAsFactors=FALSE))


# Delete unnecessary columns
DT[, c('Length', 'GeneID') := NULL]
rcDT[, c('Length', 'GeneID') := NULL]

SampNum = ncol(DT)-1
samples=colnames(DT[,c(2:(SampNum+1)), with=FALSE])

# Normalize to library size (here RPKM)
pernumber = 1000000
newCols <- list()
for (i in c(1:SampNum)) {
	newCols <- list.append(newCols, DT[[samples[i]]]/sum(rcDT[[samples[i]]])*pernumber)
	}
RPKX_DT <- copy(DT)
RPKX_DT[, (samples) := newCols]

# Get OXPHOS complex info
CI_DT <- data.table(read.table(paste0(path, 'Annotations/OXPHOScomplexes/hComplexI_core.txt'), sep = '\t', header = TRUE))
CI_nuc = CI_DT[!grep('MT-',CI_DT[[2]])][[2]]
CI_mito = CI_DT[grep('MT-',CI_DT[[2]])][[2]]

CIsupp_DT <- data.table(read.table(paste0(path, 'Annotations/OXPHOScomplexes/hComplexI_supernumery.txt'), sep = '\t', header = TRUE))
CIs_nuc = CIsupp_DT[!grep('MT-',CIsupp_DT[[2]])][[2]]
CIs_mito = CIsupp_DT[grep('MT-',CIsupp_DT[[2]])][[2]]

CII_DT <- data.table(read.table(paste0(path, 'Annotations/OXPHOScomplexes/hComplexII.txt'), sep = '\t', header = TRUE))
CII_nuc = CII_DT[!grep('MT-',CII_DT[[2]])][[2]]
CII_mito = CII_DT[grep('MT-',CII_DT[[2]])][[2]]

CIII_DT <- data.table(read.table(paste0(path, 'Annotations/OXPHOScomplexes/hComplexIII.txt'), sep = '\t', header = TRUE))
CIII_nuc = CIII_DT[!grep('MT-',CIII_DT[[2]])][[2]]
CIII_mito = CIII_DT[grep('MT-',CIII_DT[[2]])][[2]]

CIV_DT <- data.table(read.table(paste0(path, 'Annotations/OXPHOScomplexes/hComplexIV.txt'), sep = '\t', header = TRUE))
CIV_nuc = CIV_DT[!grep('MT-',CIV_DT[[2]])][[2]]
CIV_mito = CIV_DT[grep('MT-',CIV_DT[[2]])][[2]]

CV_DT <- data.table(read.table(paste0(path, 'Annotations/OXPHOScomplexes/hComplexV.txt'), sep = '\t', header = TRUE))
CV_nuc_approved = CV_DT[!grep('MT-',CV_DT[[2]])][[2]]
CV_mito = CV_DT[grep('MT-',CV_DT[[2]])][[2]]
# In case RPK file has old gene names
CV_nuc = c('ATP5A1','ATP5B','ATP5C1','ATP5D','ATP5E','ATPIF1','ATP5G1','ATP5G2','ATP5G3','USMG5','ATP5I','ATP5J2','ATP5L','C14orf2','ATP5F1','ATP5H','ATP5J','ATP5O')

CV_nuc_final = CV_nuc_approved
CV_DTconv = data.table(new = CV_DT$Approved.Symbol, old = CV_DT$Previous.Symbols)


OXPHOSmito = c(as.vector(CI_mito), as.vector(CIII_mito), as.vector(CIV_mito), as.vector(CV_mito))
OXPHOSnuc = c(as.vector(CI_nuc), as.vector(CIII_nuc), as.vector(CIV_nuc), as.vector(CV_nuc_final))

# Remove nuclear OXHPHOS genes not expressed in these cell types
OXPHOSnuc = OXPHOSnuc[OXPHOSnuc != 'COX8C'] # Germ cell specific?

###### Nuc barplot

nucDT <- RPKX_DT[GeneName %in% OXPHOSnuc]

# For isoforms of the same subunit, make a row with combined synth (e.g. ATP5MC1-2-3 (3 different loci of subunit c))
COX4I <- nucDT[GeneName=='COX4I1'][,..samples] + nucDT[GeneName=='COX4I2'][,..samples]
COX6A <- nucDT[GeneName=='COX6A1'][,..samples] + nucDT[GeneName=='COX6A2'][,..samples]
COX6B <- nucDT[GeneName=='COX6B1'][,..samples] + nucDT[GeneName=='COX6B2'][,..samples]
COX7A <- nucDT[GeneName=='COX7A1'][,..samples] + nucDT[GeneName=='COX7A2'][,..samples]
COX7B <- nucDT[GeneName=='COX7B'][,..samples] + nucDT[GeneName=='COX7B2'][,..samples]
ATP5MC123 <- nucDT[GeneName=='ATP5MC1'][,..samples] + nucDT[GeneName=='ATP5MC2'][,..samples] + nucDT[GeneName=='ATP5MC3'][,..samples]

nucDT <- rbind(nucDT, c('COX4I12', COX4I), use.names=FALSE)
nucDT <- rbind(nucDT, c('COX6A12', COX6A), use.names=FALSE)
nucDT <- rbind(nucDT, c('COX6B12', COX6B), use.names=FALSE)
nucDT <- rbind(nucDT, c('COX7A12', COX7A), use.names=FALSE)
nucDT <- rbind(nucDT, c('COX7B12', COX7B), use.names=FALSE)
nucDT <- rbind(nucDT, c('ATP5MC123', ATP5MC123), use.names=FALSE)

# Put them in order
# remove individual MCs
forOrder = OXPHOSnuc[OXPHOSnuc != 'ATP5MC1' & OXPHOSnuc != 'ATP5MC2' & OXPHOSnuc != 'ATP5MC3' & OXPHOSnuc != 'COX4I1' & OXPHOSnuc != 'COX4I2' & OXPHOSnuc != 'COX6A1' & OXPHOSnuc != 'COX6A2' & OXPHOSnuc != 'COX6B1' & OXPHOSnuc != 'COX6B2' & OXPHOSnuc != 'COX7A1' & OXPHOSnuc != 'COX7A2'  & OXPHOSnuc != 'COX7B'  & OXPHOSnuc != 'COX7B2']

nucDT <- nucDT[GeneName != 'ATP5MC1' & GeneName != 'ATP5MC2' & GeneName != 'ATP5MC3' & GeneName != 'COX4I1' & GeneName != 'COX4I2' & GeneName != 'COX6A1' & GeneName != 'COX6A2' & GeneName != 'COX6B1' & GeneName != 'COX6B2' & GeneName != 'COX7A1' & GeneName != 'COX7A2'  & GeneName != 'COX7B'  & GeneName != 'COX7B2']

# insert MC123 in their place
forOrder <- append(forOrder, 'COX4I12', after=16)
forOrder <- append(forOrder, 'COX6A12', after=19)
forOrder <- append(forOrder, 'COX6B12', after=20)
forOrder <- append(forOrder, 'COX7A12', after=21)
forOrder <- append(forOrder, 'COX7B12', after=22)
forOrder <- append(forOrder, 'ATP5MC123', after=32)
to_ord = paste0('^',forOrder,'$')

# Make new complex vectors
CI_nucUpdated <- forOrder[1:7]
CIII_nucUpdated <- forOrder[8:16]
CIV_nucUpdated <- forOrder[17:26]
CV_nucUpdated <- forOrder[27:42]

numCI <- length(CI_nucUpdated)
numCIII <- length(CIII_nucUpdated)
numCIV <- length(CIV_nucUpdated)
numCV <- length(CV_nucUpdated)


ord_index = c()
for (gene in to_ord) {
ord_index = append(ord_index, which(grepl(gene, nucDT$GeneName)))
}
setorder(nucDT[, .r := order(ord_index)], .r)[, .r := NULL]

num = nrow(nucDT)

Color_1 = 'darkgoldenrod1'
Color_2 = 'dodgerblue'
Color_3 = 'indianred3'
Color_4 = 'forestgreen'

n=length(soi)
colors = c(rep(Color_1,numCI*n), rep(Color_2,numCIII*n),  rep(Color_3,numCIV*n), rep(Color_4,numCV*n))

# Make barplot for cyto subunits
pdf(paste0(path, '7_OXPHOSstoichiometry/OXPHOSbarplots/', Experiment, '_', soiName,'_', modifier,'_barplot_cyto.pdf'),
     width=6.5, 
     height=4, 
     pointsize=11)

genes = OXPHOSnuc
plotDT = nucDT[,..soi]
if (n > 1) {
barplot(t(plotDT), beside=TRUE, names = forOrder, col = colors, width=c(rep(1.2,num)), border=NA, ylab = 'RPKM', cex.axis=.7, cex.names = .8, las=2) 
} else {
barplot(t(plotDT), beside=TRUE, names = forOrder, col = colors, width=c(rep(1.2,num)), border=NA, ylab = 'RPKM', cex.axis=.7, cex.names = .8, las=2, space=c(rep(0.1,num)))
}

dev.off()

# Write values from bar plot
printDT <- cbind(nucDT$GeneName, round(plotDT, 1))
colnames(printDT)[1] <- 'GeneName'
write.table(printDT, file=paste0(path, '7_OXPHOSstoichiometry/OXPHOSbarplots/', Experiment, '_', soiName,'_', modifier,'_barplotData_cyto.txt'), row.names=FALSE, sep=("\t"), quote=FALSE)


####################################################################################
######################### Complex averages (cyto) ####################################

CIave = mean(as.matrix(nucDT[GeneName %in% CI_nucUpdated][,..soi]), na.rm=TRUE)
CIIIave = mean(as.matrix(nucDT[GeneName %in% CIII_nucUpdated][,..soi]), na.rm=TRUE)
CIVave = mean(as.matrix(nucDT[GeneName %in% CIV_nucUpdated][,..soi]), na.rm=TRUE)
CVave = mean(as.matrix(nucDT[GeneName %in% CV_nucUpdated][,..soi]), na.rm=TRUE)


# Get ranges
CIaves = c()
CIIIaves = c()
CIVaves = c()
CVaves = c()
for (i in c(1:length(soi))) {
CIaves = c(CIaves, mean(as.matrix(nucDT[1:numCI][,..soi][,..i])))
CIIIaves = c(CIIIaves, mean(as.matrix(nucDT[(numCI+1):(numCI+numCIII)][,..soi][,..i])))
CIVaves = c(CIVaves, mean(as.matrix(nucDT[(numCI+numCIII+1):(numCI+numCIII+numCIV)][,..soi][,..i])))
CVaves = c(CVaves, mean(as.matrix(nucDT[(numCI+numCIII+numCIV+1):(numCI+numCIII+numCIV+numCV)][,..soi][,..i])))
}
CIrange = max(CIaves)-min(CIaves)
CIIIrange = max(CIIIaves)-min(CIIIaves)
CIVrange = max(CIVaves)-min(CIVaves)
CVrange = max(CVaves)-min(CVaves)

NucRanges = c(CIrange, CIIIrange, CIVrange, CVrange)

NucAves = c(CIave, CIIIave, CIVave, CVave)

MitoAves = MitoAves # From mito section of script, see above

# Plot correlations
pdf(paste0(path, '7_OXPHOSstoichiometry/ComplexAveCorrelations/', Experiment, '_vs_',MitoExpt, '_', soiName,'_Stoich_Corr_rangebars.pdf'), 
     width=3, 
     height=3.5,
     pointsize=11)

pearson_r = cor.test(~NucAves+MitoAves, method = c('pearson')) 
mylabel = bquote(italic(r) == .(format(pearson_r$estimate, digits = 2)))

# Set axis labels
if (data == 'RP') {label = 'Average Synthesis'}
if (data == 'RNAseq') {label = 'Average Abundance'}

x=MitoAves
y=NucAves

xmax = max(x)*1.1
ymax = max(y)*1.1
xlimits = c(0, xmax)
ylimits = c(0, ymax)
rlabx=xmax*1/4
rlaby=ymax*3/4

plot(x,y, main = paste0(data, ' ',soiName), cex.main = 0.8, cex.axis = 1, col = c(Color_1,Color_2,Color_3,Color_4), pch = 16, cex = 1.5, ylab = paste0(label,' (cyto)'), xlab = paste0(label, ' (mito)'), ylim=ylimits, xlim = xlimits) 

# x replicate ranges
arrows(x-(MitoRanges/2), y, x+(MitoRanges/2), y, length=0.05, angle=90, code=3, lwd=.3)
# y replicate ranges
arrows(x, y-(NucRanges/2), x, y+(NucRanges/2), length=0.05, angle=90, code=3, lwd=.3)


text(x = rlabx, y = rlaby, labels = mylabel, cex = 1)
legend('bottomright', legend=c('ComplexI', 'ComplexIII', 'ComplexIV', 'ComplexV'), cex= .8, bty="n", text.col =  c(Color_1, Color_2, Color_3, Color_4), ncol=1)

dev.off()

####################################################################################
####################################################################################



####################################################################################
####################### Cumulative fraction plots (cyto) #############################

# Set data table to RPKX_DT above
dt = copy(nucDT)
# Take means across soi
MeanRPKM <- rowMeans(dt[,..soi])
# Add means
dt[, soiMean := MeanRPKM]

# Get stoichiometries
CIcytoStoich = rep(1, numCI)
CIIIcytoStoich = rep(1, numCIII)
CIVcytoStoich = rep(1, numCIV)
CVcytoStoich = c(3,3,rep(1, 4), 8, rep(1,9))

# Add stoichiometry
dt[, Stoich := c(CIcytoStoich, CIIIcytoStoich, CIVcytoStoich, CVcytoStoich)]

# Make columns of interest numeric
dt$soiMean <- as.numeric(dt$soiMean)

# Normalize synth rate to stoichiometry
dt[, NormSynthRate := soiMean/Stoich]


# Find median value for each complex
CImed <- median(dt[GeneName %in% CI_nucUpdated]$NormSynthRate, na.rm=TRUE)
CIIImed <- median(dt[GeneName %in% CIII_nucUpdated]$NormSynthRate, na.rm=TRUE)
CIVmed <- median(dt[GeneName %in% CIV_nucUpdated]$NormSynthRate, na.rm=TRUE)
CVmed <- median(dt[GeneName %in% CV_nucUpdated]$NormSynthRate, na.rm=TRUE)

# Find rate limiting value for each complex

CIrl <- min(dt[GeneName %in% CI_nucUpdated]$NormSynthRate, na.rm=TRUE)
CIIIrl <- min(dt[GeneName %in% CIII_nucUpdated]$NormSynthRate, na.rm=TRUE)
CIVrl <- min(dt[GeneName %in% CIV_nucUpdated]$NormSynthRate, na.rm=TRUE)
CVrl <- min(dt[GeneName %in% CV_nucUpdated]$NormSynthRate, na.rm=TRUE)

# Add medians to DT
dt[, Median := c(rep(CImed, numCI), rep(CIIImed, numCIII), rep(CIVmed, numCIV), rep(CVmed, numCV))]

# Add limiting rate to DT
dt[, rl := c(rep(CIrl, numCI), rep(CIIIrl, numCIII), rep(CIVrl, numCIV), rep(CVrl, numCV))]

CytoRL = unique(dt$rl)


# Fold diff from median
dt[, log2med := log(NormSynthRate/Median, 2)]

# Fold diff from rate-limiting
dt[, log2rl := log(NormSynthRate/rl, 2)]


pdf(paste0(path, '7_OXPHOSstoichiometry/CumFracPlots/', Experiment, '_', soiName,'_', modifier,'_CumFractionSubunitStoich_cyto.pdf'),     		  
	 width=3, 
     height=2.5,
     pointsize=11)
      
cytocomplexes=c(list(CI_nucUpdated), list(CIII_nucUpdated), list(CIV_nucUpdated), list(CV_nucUpdated))  
colors=c(Color_1, Color_2, 'red', Color_4)


# Then plot for mito
compartment='Cyto'
complexes=cytocomplexes

# Spread around rl
plot(1, xlim=c(0,5.2), ylim=c(0,1), pch=16, type='n', xlab='Synth rate relative to rate limiting (log2)', ylab='Cumulative fraction of subunits', col=colors[i], lwd=2, main=paste0(Experiment, ' ', compartment, soiName), yaxt='n')
# Make shaded box at 2-fold
polygon(c(0,0,1,1), c(0,1,1,0), col = alpha('black', .1), border = NA)
# Draw lines for CI, CIV, CV (not CIII since there's only 1
for (i in c(1,2,3,4)) { 
x=dt[GeneName %in% complexes[[i]]]$log2rl
x=sort(x)
x=c(x[1], x) # To make the stair start at the 0 line
y=c(0:(length(x)-1))
y=y/max(y)
points(x,y, pch=16, type='s', col=colors[i], lwd=2)
}
axis(2, at=c(0,.5, 1), labels=c(0, .5, 1))
legend('bottomright', legend=c('ComplexI','ComplexIII', 'ComplexIV', 'ComplexV'), text.col=colors, bty='n', cex=.7)

dev.off()
