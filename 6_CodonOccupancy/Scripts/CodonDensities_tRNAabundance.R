#!/n/app/R/4.0.1/bin/Rscript


library(data.table)
library('scales')
library(Rfast) # to use rowMins


# Get arguments
args <- commandArgs(trailingOnly = TRUE)
Experiment <- args[1]
sizeRange <- args[2]


path = paste0(getwd(),'/')

cdDT <- data.table(read.table(paste0(path, Experiment, '_', sizeRange,'_CodonDensities.txt'),header=TRUE))

ylimits = c(0,8)


# For tRNA abundances and more codon info
taDT <- data.table(read.table(paste('Annotations/hMitoGeneticCode_BattersbyMyoblast_tRNAabundance.txt'), header=TRUE, , sep = ','))

# Merge
DT = data.table(merge(cdDT,taDT, by='codon'))

setorder(DT, cols = 'AminoAcid')

SampNum = (ncol(DT) - 7) / 2

# Get samplenames from column names
samplenames = colnames(DT)[seq(2,(SampNum*2), by=2)]
samples = c()
for (i in c(1:SampNum)) {
samples = c(samples, strsplit(samplenames,'_ave')[[i]][1])
}

# Get densities and std deviation
dens = paste0(samples, '_aveDensity')
sd = paste0(samples, '_SD')

# Get means and ranges across all samples in experiment
means = rowMeans(DT[, ..dens])
ranges = rowMaxs(as.matrix(DT[, ..dens]), value=TRUE) - rowMins(as.matrix(DT[, ..dens]), value=TRUE)

# plotDT = DT[, ..dens]
# sdDT = DT[, ..sd]

# names = paste0(DT$CodonFreq,' ',DT$codon,' ',DT$AminoAcid)
names = paste0(DT$codon,' ',DT$AminoAcid)

colors = c('red', 'darkorange','darkgoldenrod1', 'forestgreen','dodgerblue', 'blue','darkblue','purple','violet','grey50','black','brown','red', 'darkorange','darkgoldenrod1', 'forestgreen','dodgerblue', 'darkblue','purple','violet','grey50','black')

colRep = c(4,4,2,2,2,2,2,4,2,2,6,2,2,2,4,6,2,2,4,2,2,4)
repTimes = c(4,4,2,2,2,2,2,4,2,2,6,2,2,2,4,6,2,2,4,2,2,4) * SampNum # Different color for each amino acid


################ Barplot of codon densities ################
pdf(paste0(path, Experiment,'_', sizeRange,'_CodonDens_barplot.pdf'),
     width=6.8, # width=4.5, width=2, width=2.5
     height=4, #4
     pointsize=11)

# Plot a bar for each sample
xx = barplot(t(as.matrix(DT[, ..dens])), beside=TRUE, names = names, col = rep(colors,repTimes), border=NA, ylab = 'Mean codon occupancy', cex.axis=.7, cex.names = .65, ylim=ylimits, las=2) # t(plotDT)

# Means with range bars
# xx = barplot(means, beside=TRUE, names = names, col = rep(colors,repTimes), border=NA, ylab = 'Mean codon occupancy', cex.axis=.7, cex.names = .65, ylim=ylimits, las=2) # t(plotDT)
# 
# if (length(samples) > 1) {
# # range bars
# arrows(xx, means-(ranges/2), xx, means+(ranges/2), length=0.03, angle=90, code=3, lwd=.2)
# }

# standard deviation
# arrows(xx, t(plotDT)-t(sdDT), xx, t(plotDT)+t(sdDT), length=0.01, angle=90, code=3, lwd=.2)

abline(h=1, lwd = .2, lty = 2, col = alpha('grey60', .5))

# Plot SD

# x=c(1.5:((nrow(DT)*(SampNum+1))+.5))
# # Make 0 column at end of data table to make error bars line up
# fill_plotDT <- copy(plotDT)
# fill_plotDT[, fill := 0]
# y=as.vector(t(fill_plotDT))
# # Get sd
# sdIdx = seq(3,2*SampNum+2,2)
# sdDT <- DT[,sdIdx, with=FALSE]
# # Make 0 column at end of data table to make error bars line up
# fill_sdDT <- copy(sdDT)
# fill_sdDT[, fill := 0]
# sd=as.vector(t(fill_sdDT))
# 
# arrows(x, y-sd, x, y+sd, length=0.005, angle=90, code=3, lwd = 0.1)

dev.off()





################ Individual codon densities ################
pdf(paste0(path,Experiment, '_', sizeRange,'_indCodonDens_boxplot.pdf'),
     width=3, # width=4.5, width=2, width=2.5
     height=3*SampNum, #4
     pointsize=11)

par(mfrow=c(SampNum,1))

for (i in c(1:SampNum)) {

libName = samples[i]
indDT <- data.table(read.table(paste0(path, '6_CodonOccupancy/',libName, '_', sizeRange,'_indCodonDensities.txt'),header=FALSE, sep='['))

## Asp
codon = 'GAT'
occ = as.character(indDT[V1==paste0(codon,',')]$V2)
GAToccs = as.double(strsplit(strsplit(occ,']')[[1]],', ')[[1]]) # Split by ',' and ']'
codon = 'GAC'
occ = as.character(indDT[V1==paste0(codon,',')]$V2)
GACoccs = as.double(strsplit(strsplit(occ,']')[[1]],', ')[[1]]) # Split by ',' and ']'
Asp_occs = c(GAToccs, GACoccs)

# Leu
codon = 'TTG'
occ = as.character(indDT[V1==paste0(codon,',')]$V2)
TTGoccs = as.double(strsplit(strsplit(occ,']')[[1]],', ')[[1]]) # Split by ',' and ']'
codon = 'TTA'
occ = as.character(indDT[V1==paste0(codon,',')]$V2)
TTAoccs = as.double(strsplit(strsplit(occ,']')[[1]],', ')[[1]]) # Split by ',' and ']'
Leu_occs = c(TTGoccs, TTAoccs)

# Phe
codon = 'TTT'
occ = as.character(indDT[V1==paste0(codon,',')]$V2)
TTToccs = as.double(strsplit(strsplit(occ,']')[[1]],', ')[[1]]) # Split by ',' and ']'
codon = 'TTC'
occ = as.character(indDT[V1==paste0(codon,',')]$V2)
TTCoccs = as.double(strsplit(strsplit(occ,']')[[1]],', ')[[1]]) # Split by ',' and ']'
Phe_occs = c(TTToccs, TTCoccs)

# Arg
# codon = 'CGT'
# occ = as.character(indDT[V1==paste0(codon,',')]$V2)
# CGToccs = as.double(strsplit(strsplit(occ,']')[[1]],', ')[[1]]) # Split by ',' and ']'
# Arg_occs = CGToccs

# His
# codon = 'CAC'
# occ = as.character(indDT[V1==paste0(codon,',')]$V2)
# CACoccs = as.double(strsplit(strsplit(occ,']')[[1]],', ')[[1]]) # Split by ',' and ']'
# His_occs = CACoccs


     
boxplot(Asp_occs, Leu_occs, Phe_occs,names = c('Asp', 'Leu', 'Phe'), outline = FALSE, ylab = 'Codon Occupancy', ylim = c(0,16), main=libName) # ylim = c(0,9)
points(jitter(rep(1,length(Asp_occs)),12), Asp_occs, pch = 16,cex = .7, col = alpha('forestgreen',.6))
points(jitter(rep(2,length(Leu_occs)),5), Leu_occs, pch = 16,cex = .7, col = alpha('black', .6))
points(jitter(rep(3,length(Phe_occs)),3), Phe_occs, pch = 16,cex = .7, col = alpha('darkorange', .6))
# points(jitter(rep(4,length(Arg_occs)),3), Arg_occs, pch = 16,cex = .7, col = alpha('darkorange', .6))
# points(jitter(rep(5,length(His_occs)),3), His_occs, pch = 16,cex = .7, col = alpha('violet', .6))

abline(h=1, lty = 2, lwd = 0.2)
}

dev.off()

# Scatterplot of codon occupancies vs tRNA abundance
pdf(paste0(path, Experiment,'_', sizeRange,'_CodonDens_vs_tRNAabundance.pdf'),
     width=4, # width=4.5, width=2, width=2.5
     height=3*SampNum, #4
     pointsize=11)

par(mfrow=c(SampNum,1))

for (i in c(1:SampNum)) {

tRNAabund = DT$meanAbundance
# For plotting each sample separately
codonOcc = DT[[dens[i]]]
# For plotting mean of all samples
# codonOcc = means

pearson_r = cor.test(~tRNAabund[tRNAabund>0]+codonOcc[tRNAabund>0], method = c('pearson')) 
mylabel = bquote(italic(r) == .(format(pearson_r$estimate, digits = 2)))

plot(tRNAabund, codonOcc, main = samples[i], cex.main = 1, cex.axis = 1, col = rep(colors,colRep), pch = 16, cex = 1, xlab = 'tRNA abundance (tpm) (Richter et al. 2018)', ylab = 'Mean codon occupancy', ylim=ylimits) 


# tRNA abundance SD (3 replicates)
arrows(tRNAabund-DT$SD, codonOcc, tRNAabund+DT$SD, codonOcc, length=0.03, angle=90, code=3, lwd=.3)
# codon occupancy ranges
# arrows(tRNAabund, codonOcc-(ranges/2), tRNAabund, codonOcc+(ranges/2), length=0.03, angle=90, code=3, lwd=.3)

# plot points again to make them on top of error bars
points(tRNAabund, codonOcc,col = rep(colors,colRep), pch = 16, cex = 1.5)

text(x = .002, y = max(ylimits)-.5, labels = mylabel, cex = 1)
}
dev.off()

# source('./6_CodonOccupancy/Scripts/CodonDensities_tRNAabundance.R')

