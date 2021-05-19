

library(data.table)
library(pheatmap)
library(inlmisc)
library(RColorBrewer)
library(stringr)
# display.brewer.all()


Folder = 'hMitoRP1' 
Experiment = 'hMitoRP1' 
Species = 'Aligned.Mito_mRNA.noDups_noSoft' # hMito_tRNA_al_ hMito_rRNA_al_ Aligned.Nuc_mRNA.noDups.uniq_ Aligned.Mito_mRNA.noDups.uniq_ Aligned.Mito_mRNA.noDups_noSoft Cleaned_
path = 'pathToLengthDistributionsFile/'
colors=c('dodgerblue', 'forestgreen', 'gold1', 'red')


# Length distributions 
DT <- data.table(read.table(paste0(path, Experiment,'_',Species, 'LengthDist.txt'), sep = '\t', header=TRUE, stringsAsFactors = FALSE, fill=TRUE))

# Replace NA values with 0
DT[is.na(DT)] <- 0

SampNum=ncol(DT)-1
samples = colnames(DT)[2:(SampNum+1)]

n = SampNum


# Normalize
normDT <- copy(DT)
for (i in c(1:n)) {
normDT[, c(samples[i]) := DT[[samples[i]]]/sum(DT[[samples[i]]])*100]
}

# Remove very large values, which are all 0
lower <- normDT[Length <=42]
sampNormDT <- copy(lower)
# Remove length column
sampNormDT[, c('Length') := NULL]
# Reverse order of columns
revSampNames = rev(samples)


# setcolorder(sampNormDT, revSampNames)

# Heat Plot
pdf(paste0(path, Experiment,'_',Species, 'LengthDist_heatplot.pdf'),
     width=4,
     height=2.3,
     )

   
pheatmap(t(as.matrix(sampNormDT, rownames=lower$Length)), cluster_rows = FALSE, cluster_cols = FALSE, border_color = NA, col = GetColors(256, "YlOrBr"), na_col = "white",legend=TRUE, scale = "none", cex=.9) # cm.colors(256)  GetColors(256, "YlOrBr") brewer.pal(n=9,name='Blues')

dev.off()


# Line Plot
xlimits=c(11,40)
ylimits=c(0,15)
pdf(paste0(path, Experiment,'_',Species, 'LengthDist_lineplot.pdf', sep = ''),
     width=3,
     height=3.5,
     pointsize=10.5
     )

plot(spline(normDT$Length,normDT[[samples[1]]], method = 'n', n=150), type = 'l', lwd = 3, col = colors[1], ylim = ylimits, xlab = 'Length', ylab = '%', xlim = xlimits)

for (i in c(2:n)) {
lines(spline(normDT$Length,normDT[[sampNames[i]]], method = 'n', n=150), lwd = 3, col = colors[i], lty=1)
}

legend('topleft', legend=c(samples), col=c(colors), lty=1, lwd = 3, cex=0.7, bty='n') 


dev.off()

# source('./1_AlignData/Scripts/ReadLengthPlot.R')
