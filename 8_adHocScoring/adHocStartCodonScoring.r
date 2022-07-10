
### Written by M. Couvillion 10/2020
### This script will take a bedGraph file of A site RPF density as input, 
### Use A site density on putative start codons, enrichment of in-frame periodicity, and overall RPF density to assign to putative start codons a likelihood score of being used as an initiation codon
### Plot density in each of three frames on plus strand across the mito genome, label known starts and ends, and label scores of putative start codons
### Plot distribution of all scores

### Paths will need to be changed throughout


library(data.table)
library(Biostrings)
library(scales)

# Define experiment
Folder = 'HEK' # Fibroblasts GalShift HelaIP HelaFracs HEK
Experiment = 'HEK_WTall' # GalShift GalShift_Glu MyoDiff1 OXPHOSinh (HelaIP) Gal_90m_1_reseq Fibroblast GalShiftReseq HelaFracs_F6 HelaIP_8U HelaIP_8Uplreseq HEK_WT
sizeRange = '32to34' # 33to35 (HelaIP) 31to35 (OXPHOSinh) 31to33or34 (Fibroblast) 30to33 (GalShift) 31to33 (GalShiftReseq) 32to34 (HEK) 31to33(MyoDiff1) 31to33(MyoblastsExp1-5)
modifier = 'Mito_mRNA' # Mito chrMall
strand = 'Pall' # Pall (GalShift HelaIP) P (Fibroblasts)
region = 'all' # ND5_UTR humanin MOTS-c ND5_UTRall
ylimits = NULL

# Get footprint data * Make sure there is a header line (e.g. 'track type...') in bedGraph
bG <- data.table(read.table(paste0('/Users/Mary/Desktop/Data/hMitoRP/PanAnalysis/', Folder, '/hMitoRP/bedGraphs/', Experiment, '_', modifier, '.noDups.Asite_', sizeRange,'_', strand,'.bedGraph'),skip = 1))


# Get sequence info for region of interest
dna <- readDNAStringSet('/Users/Mary/Desktop/Data/hMitoRP/SequenceFiles/hMitoGenome.fasta')

if (region == 'all') {
begin = 1
end = length(dna[[1]])-3}
if (region == 'ND5_UTRall') {
begin = 14134
end = 14748
ylimits = c(1,600)
}



# To get sequence out
ref1 <- c(begin:end)
ref2 <- c((begin+1):(end+1))
ref3 <- c((begin+2):(end+2))
seq <- dna[[1]][begin:(end+2)]

seq1 <- dna[[1]][begin:end]
seq2 <- dna[[1]][(begin+1):(end+1)]
seq3 <- dna[[1]][(begin+2):(end+2)]

# Get codons for each ref
codons = c()
for (i in c(1:(length(seq)-2))){ # -3 so that they'll all be the same length
codons = c(codons, as.character(seq[i:(i+2)]))
}
codon = 1:length(codons)
codon1 = codons[(codon+2)%%3==0]
codon2 = codons[(codon+1)%%3==0]
codon3 = codons[(codon)%%3==0]

# Get 2nd position for each ref (this is where A site density is assigned)
ref1pos2 = ref1[seq(2, length(ref1), 3)]
ref2pos2 = ref2[seq(2, length(ref2), 3)]
ref3pos2 = ref3[seq(2, length(ref3), 3)]
val1 = bG$V4[ref1pos2]
val2 = bG$V4[ref2pos2]
val3 = bG$V4[ref3pos2]

# add pseudocount of 1 for calculating mean/median later
val1 = val1 + 1
val2 = val2 + 1
val3 = val3 + 1



# List starts and ends for genes
starts = c(648, 1671, 3307, 4470, 5904, 7586, 8366, 8527, 9207, 10059, 10470, 10760, 12337, 14196, 14255, 14297, 14673, 14747)
ends = c(1601, 3229, 4262, 5511, 7445, 8269, 8572, 9207, 9990, 10404, 10766, 12137, 14148, 14258, 14296, 14500, 14149, 15887)
putends = c(14210)
genes = c('RNR1','RNR2','ND1','ND2','CO1','CO2','ATP8','ATP6','CO3','ND3','ND4L','ND4', 'ND5','dORF04','dORF06','dORF065','ND6','CYB')


# Plot read counts
pdf(paste0('/Users/Mary/Desktop/Data/hMitoRP/PanAnalysis/', Folder, '/hMitoRP/Frame/', Experiment, '_', sizeRange,'_', modifier,'_',strand,'_StackedFrames_assignScores.pdf'), width = length(seq)/40, height = 6, pointsize=11)

# par(mfrow=c(3,1))

plot(ref1pos2, val1, log='y',type = 'n', pch=16, cex = .5,  xlab = 'chrM position', ylab = 'pos 2 read count', main = paste0('Stacked frames ',region), ylim = ylimits)

linewidth=1.5
col1 = 'orange'
col2 = 'limegreen'
col3 = 'red'
# col4 = 'limegreen'
# col5 = 'dodgerblue'
# col6 = 'forestgreen'

lines(ref1pos2, val1, type = 'l', lwd = linewidth, col = col1)
lines(ref2pos2, val2, type = 'l', lwd = linewidth, col = col2)
lines(ref3pos2, val3, type = 'l', lwd = linewidth, col = col3)


# Mark start codons in each frame
if (strand == 'P' | strand == 'Pall') {
starts1 = which(codon1 == 'ATG' | codon1 == 'ATT' | codon1 == 'ATA' | codon1 == 'GTG')*3-1
starts2 = which(codon2 == 'ATG' | codon2 == 'ATT' | codon2 == 'ATA' | codon2 == 'GTG')*3
starts3 = which(codon3 == 'ATG' | codon3 == 'ATT' | codon3 == 'ATA' | codon3 == 'GTG')*3+1
ends1 = which(codon1 == 'TAA' | codon1 == 'TAG' | codon1 == 'AGG' | codon1 == 'AGA')*3-1
ends2 = which(codon2 == 'TAA' | codon2 == 'TAG' | codon2 == 'AGG' | codon2 == 'AGA')*3
ends3 = which(codon3 == 'TAA' | codon3 == 'TAG' | codon3 == 'AGG' | codon3 == 'AGA')*3+1
}
if (strand == 'M' | strand == 'Mall') {
starts1 = which(codon1 == 'CAT' | codon1 == 'TAA' | codon1 == 'TAT' | codon1 == 'CAC')*3-1
starts2 = which(codon2 == 'CAT' | codon2 == 'TAA' | codon2 == 'TAT' | codon2 == 'CAC')*3
starts3 = which(codon3 == 'CAT' | codon3 == 'TAA' | codon3 == 'TAT' | codon3 == 'CAC')*3+1
ends1 = which(codon1 == 'ATT' | codon1 == 'ATC' | codon1 == 'TCC' | codon1 == 'TCT')*3-1
ends2 = which(codon2 == 'ATT' | codon2 == 'ATC' | codon2 == 'TCC' | codon2 == 'TCT')*3
ends3 = which(codon3 == 'ATT' | codon3 == 'ATC' | codon3 == 'TCC' | codon3 == 'TCT')*3+1
}
# Remove any starts that are within the window from beginning of known gene start (or will cause problems later)
window = 16 # if this is any longer will get into the ND5 frame from dORF
forwindow = 20
starts1 = starts1[starts1 > (forwindow*3)]
starts2 = starts2[starts2 > (forwindow*3)]
starts3 = starts3[starts3 > (forwindow*3)]

# abline(v=starts1+3, col=col1, lwd=.5)
# abline(v=starts2+3, col=col2, lwd=.5)
# abline(v=starts3+3, col=col3, lwd=.5)
abline(v=ends1, col=col1, lwd=.5, lty=2)
abline(v=ends2, col=col2, lwd=.5, lty=2)
abline(v=ends3, col=col3, lwd=.5, lty=2)

text(starts+3, max(val1), label=genes) 

# Find frame of each known gene
frame1genes = genes[which((starts+1) %in% starts1)]
frame2genes = genes[which((starts+1) %in% starts2)]
frame3genes = genes[which((starts+1) %in% starts3)]
# And get the starts of those
frame1starts = starts[which(genes %in% frame1genes)]
frame2starts = starts[which(genes %in% frame2genes)]
frame3starts = starts[which(genes %in% frame3genes)]


###################### Find putative ORFs ######################
# 1. For each reference, find peaks at start codons
# Define fold-change cutoffs
fivech = 0.8
threech = 0.6
# log transform
logfivech = log(fivech,2)
logthreech = log(threech,2)

vals = c(deparse(substitute(val1)),deparse(substitute(val2)),deparse(substitute(val3)))
cols = c('darkorange', 'forestgreen', 'darkred')
refs = c(deparse(substitute(ref1pos2)), deparse(substitute(ref2pos2)), deparse(substitute(ref3pos2)))
startss = c(deparse(substitute(starts1)), deparse(substitute(starts2)), deparse(substitute(starts3)))
# Now for each frame
for (i in c(1:3)) {
# Get differences 5' and 3' of each position (ND5dORF in ref3, index 4733)
diffval_5p = c(0, diff(get(vals[i]))) # add 0 to get value at 5' change
diffval_3p = c(diff(get(vals[i])), 0)
# Get fold-change 5' and 3' of each position
foldval_5p = log(c(0,tail(get(vals[i]), length(get(vals[i]))-1)/head(get(vals[i]), length(get(vals[i]))-1)),2)
foldval_3p = log(c(head(get(vals[i]), length(get(vals[i]))-1)/tail(get(vals[i]), length(get(vals[i]))-1),0),2)
# Get index of start codons in ref
refstarts = match(get(startss[i])+3, get(refs[i]))
# Get locations where fold-change cutoffs are met
peaks=intersect(which(foldval_5p >= logfivech & foldval_5p != Inf), which(foldval_3p >= logthreech & foldval_3p != Inf))
# Now which of these are at starts (+3)?
peaksANDstarts=intersect(peaks, refstarts)
assign(paste0('peaksANDstarts', i), peaksANDstarts)
# And get the positions
locRef = get(refs[i])[peaksANDstarts]
# abline(v=locRef, col=cols[i], lwd=2)
# Score for fold-change at start (fold-change on 5' end + fold-change on 3' end
stScore = (foldval_5p[peaksANDstarts] + foldval_3p[peaksANDstarts])
assign(paste0('stScore', i), stScore)
}

# 164


# 2. For each reference, get percentage in frame
# For ref1, compare to previous ref3 and same ref2
# For ref2, compare to same ref1 and same ref3
# For ref3, compare to same ref2 and next ref1
val1perc = c(val1,0)/(c(val1,0) + c(val2,0) + c(0,val3))*100 # offset vectors to compare in the correct frame
val1perc = head(val1perc, length(val1perc)-1) # Remove the added element
val2perc = val2/(val1+val2+val3)*100 
val3perc = c(0,val3)/(c(val1,0) + c(0,val2) + c(0,val3))*100 # offset vectors to compare in the correct frame
val3perc = tail(val3perc, length(val3perc)-1) # Remove the added element

# For each location of peaks and starts calculate average % before and after start

FrameIncrease = .1
w = window # windowToAverage in amino acids
forw = forwindow
measure = mean

AveBefPeaksANDstarts1=c()
AveAftPeaksANDstarts1=c()
# Use [peaksANDstarts>w] to avoid getting negative indices
for (peak in peaksANDstarts1) {
AveBefPeaksANDstarts1 = c(AveBefPeaksANDstarts1, measure(val1perc[(peak-w):(peak-1)]))
AveAftPeaksANDstarts1 = c(AveAftPeaksANDstarts1, measure(val1perc[(peak+1):(peak+forw)]))
}
AveBefPeaksANDstarts2=c()
AveAftPeaksANDstarts2=c()
for (peak in peaksANDstarts2) {
AveBefPeaksANDstarts2 = c(AveBefPeaksANDstarts2, measure(val2perc[(peak-w):(peak-1)]))
AveAftPeaksANDstarts2 = c(AveAftPeaksANDstarts2, measure(val2perc[(peak+1):(peak+forw)]))
}
AveBefPeaksANDstarts3=c()
AveAftPeaksANDstarts3=c()
for (peak in peaksANDstarts3) {
AveBefPeaksANDstarts3 = c(AveBefPeaksANDstarts3, measure(val3perc[(peak-w):(peak-1)]))
AveAftPeaksANDstarts3 = c(AveAftPeaksANDstarts3, measure(val3perc[(peak+1):(peak+forw)]))
}
# Get change in frame before and after start (Frame score)
Diff1 = (AveAftPeaksANDstarts1 - AveBefPeaksANDstarts1) 
Diff2 = (AveAftPeaksANDstarts2 - AveBefPeaksANDstarts2) 
Diff3 = (AveAftPeaksANDstarts3 - AveBefPeaksANDstarts3) 



# 3. Add to score overall change in density
measure = mean
wind = 5
TotDensBefPeaksANDstarts1=c()
TotDensAftPeaksANDstarts1=c()
# Use [peaksANDstarts>w] to avoid getting negative indices
for (peak in peaksANDstarts1) {
TotDensBefPeaksANDstarts1 = c(TotDensBefPeaksANDstarts1, measure(val1[(peak-wind):(peak-1)] + val2[(peak-wind):(peak-1)] + val3[(peak-wind-1):(peak-2)])) # to go back one full codon have to substract 2 from frame 3
TotDensAftPeaksANDstarts1 = c(TotDensAftPeaksANDstarts1, measure(val1[(peak+1):(peak+wind)] + val2[(peak+1):(peak+wind)] + val3[(peak):(peak+wind-1)]))
}
TotDensBefPeaksANDstarts2=c()
TotDensAftPeaksANDstarts2=c()
for (peak in peaksANDstarts2) {
TotDensBefPeaksANDstarts2 = c(TotDensBefPeaksANDstarts2, measure(val2[(peak-wind):(peak-1)] + val1[(peak-wind):(peak-1)] + val3[(peak-wind):(peak-1)]))
TotDensAftPeaksANDstarts2 = c(TotDensAftPeaksANDstarts2, measure(val2[(peak+1):(peak+wind)] + val1[(peak+1):(peak+wind)] + val3[(peak+1):(peak+wind)]))
}
TotDensBefPeaksANDstarts3=c()
TotDensAftPeaksANDstarts3=c()
for (peak in peaksANDstarts3) {
TotDensBefPeaksANDstarts3 = c(TotDensBefPeaksANDstarts3, measure(val3[(peak-wind):(peak-1)] + val1[(peak-wind+1):(peak)] + val2[(peak-wind):(peak-1)]))
TotDensAftPeaksANDstarts3 = c(TotDensAftPeaksANDstarts3, measure(val3[(peak+1):(peak+wind)] + val1[(peak+2):(peak+wind+1)] + val2[(peak+1):(peak+wind)]))
}
# Raw values are too large and take over the scoring
densDiff1 = log(TotDensAftPeaksANDstarts1/TotDensBefPeaksANDstarts1)*2
densDiff2 = log(TotDensAftPeaksANDstarts2/TotDensBefPeaksANDstarts2)*2
densDiff3 = log(TotDensAftPeaksANDstarts3/TotDensBefPeaksANDstarts3)*2


# Combine scores
knownGeneFr=35 # Expect known frames to average over 30%

stScore1 = stScore1[which(AveBefPeaksANDstarts1 < knownGeneFr)]
stScore2 = stScore2[which(AveBefPeaksANDstarts2 < knownGeneFr)]
stScore3 = stScore3[which(AveBefPeaksANDstarts3 < knownGeneFr)]
Diff1 = Diff1[which(AveBefPeaksANDstarts1 < knownGeneFr)]
Diff2 = Diff2[which(AveBefPeaksANDstarts2 < knownGeneFr)]
Diff3 = Diff3[which(AveBefPeaksANDstarts3 < knownGeneFr)]
densDiff1 = densDiff1[which(AveBefPeaksANDstarts1 < knownGeneFr)]
densDiff2 = densDiff2[which(AveBefPeaksANDstarts2 < knownGeneFr)]
densDiff3 = densDiff3[which(AveBefPeaksANDstarts3 < knownGeneFr)]

TotScore1 = round(stScore1 + Diff1 + densDiff1, 1)
TotScore2 = round(stScore2 + Diff2 + densDiff2, 1)
TotScore3 = round(stScore3 + Diff3 + densDiff3, 1)

# Get positions & remove all where densDiff=0 (no signal at all)
peaksANDstarts1 = peaksANDstarts1[which(AveBefPeaksANDstarts1 < knownGeneFr)]
peaksANDstarts2 = peaksANDstarts2[which(AveBefPeaksANDstarts2 < knownGeneFr)]
peaksANDstarts3 = peaksANDstarts3[which(AveBefPeaksANDstarts3 < knownGeneFr)]

ref1nt = ref1pos2[peaksANDstarts1]
ref2nt = ref2pos2[peaksANDstarts2]
ref3nt = ref3pos2[peaksANDstarts3]

ref1nt = ref1nt[which(densDiff1 != 0)]
ref2nt = ref2nt[which(densDiff2 != 0)]
ref3nt = ref3nt[which(densDiff3 != 0)]

TotScore1 = TotScore1[which(densDiff1 != 0)]
TotScore2 = TotScore2[which(densDiff2 != 0)]
TotScore3 = TotScore3[which(densDiff3 != 0)]

Diff1 = Diff1[which(densDiff1 != 0)]
Diff2 = Diff2[which(densDiff2 != 0)]
Diff3 = Diff3[which(densDiff3 != 0)]
stScore1 = stScore1[which(densDiff1 != 0)]
stScore2 = stScore2[which(densDiff2 != 0)]
stScore3 = stScore3[which(densDiff3 != 0)]
densDiff1 = densDiff1[which(densDiff1 != 0)]
densDiff2 = densDiff2[which(densDiff2 != 0)]
densDiff3 = densDiff3[which(densDiff3 != 0)]

# Now remove scores that are in known gene frame but too close to beginning to be filtered out with 35% frame cutoff
`%notin%` <- Negate(`%in%`)

# Make lists of values to be removed
nostarts1start = c(frame1starts+5)
nostarts1end = c(frame1starts+(window*3))
removestarts1 = c()
for (i in c(1:length(nostarts1start))) {removestarts1 = c(removestarts1, c(nostarts1start[i]: nostarts1end[i]))}

nostarts2start = c(frame2starts+5)
nostarts2end = c(frame2starts+(window*3))
removestarts2 = c()
for (i in c(1:length(nostarts2start))) {removestarts2 = c(removestarts2, c(nostarts2start[i]: nostarts2end[i]))}

nostarts3start = c(frame3starts+5)
nostarts3end = c(frame3starts+(window*3))
removestarts3 = c()
for (i in c(1:length(nostarts3start))) {removestarts3 = c(removestarts3, c(nostarts3start[i]: nostarts3end[i]))}


TotScore1 = TotScore1[which(ref1nt %notin% removestarts1)]
TotScore2 = TotScore2[which(ref2nt %notin% removestarts2)]
TotScore3 = TotScore3[which(ref3nt %notin% removestarts3)]
Diff1 = Diff1[which(ref1nt %notin% removestarts1)]
Diff2 = Diff2[which(ref2nt %notin% removestarts2)]
Diff3 = Diff3[which(ref3nt %notin% removestarts3)]
stScore1 = stScore1[which(ref1nt %notin% removestarts1)]
stScore2 = stScore2[which(ref2nt %notin% removestarts2)]
stScore3 = stScore3[which(ref3nt %notin% removestarts3)]
densDiff1 = densDiff1[which(ref1nt %notin% removestarts1)]
densDiff2 = densDiff2[which(ref2nt %notin% removestarts2)]
densDiff3 = densDiff3[which(ref3nt %notin% removestarts3)]

ref1nt = ref1nt[ref1nt %notin% removestarts1]
ref2nt = ref2nt[ref2nt %notin% removestarts2]
ref3nt = ref3nt[ref3nt %notin% removestarts3]


# Label putative start codons on plot
abline(v=ref1nt, col=alpha(col1, .3), lwd=10)
text(ref1nt,  .7*max(val1), label=TotScore1, cex = .6)
abline(v=ref2nt, col=alpha(col2, .3), lwd=10)
text(ref2nt,  .7*max(val1), label=TotScore2, cex = .6) 
abline(v=ref3nt, col=alpha(col3, .3), lwd=10)
text(ref3nt,  .7*max(val1), label=TotScore3, cex = .6) 

# Label high-scoring positions that overlap between experiments
# GalshiftReseq & HelaIP
# NewPut = c(10254, 15482, 12925, 7964, 15344)
# abline(v=NewPut, col='black', lwd=1)
# abline(v=NewPut, col=alpha('yellow', .2), lwd=15)

legend('topleft', c('Frame ref 1 start with in-frame peak', 'Frame ref 2 start with in-frame peak', 'Frame ref 3 start with in-frame peak', 'stop codon'), col = c(alpha(c(col1, col2, col3), .3), 'black') ,lty = c(1,1,1,2), box.lty = 0, lwd = c(10,10,10,1))


dev.off()


# GOIs = c('ATP6', 'COX3', 'ND4', 'ND5\ndORF04')
GOIs = c('ND1','ND2','CO1','CO2','ATP8','ATP6','CO3','ND3','ND4L','ND4', 'ND5','ND6','CYB','ND5-dORF')
# StartsOI = c(8531, 9211, 10764, 14200)
StartsOI = c(3307, 4470, 5904, 7586, 8366, 8527, 9207, 10059, 10470, 10760, 12337, 14673, 14747, 14196)+4
cols = c('blue','blue','blue','blue','blue','red','red','blue','blue','red', 'blue','blue','blue','red')
# find scores for GOI
ScoresOI = c()
for (i in c(1:length(GOIs))) {
test1 = 0
test2 = 0
test3 = 0
test1 = which(ref1nt==StartsOI[i])
test2 = which(ref2nt==StartsOI[i])
test3 = which(ref3nt==StartsOI[i])
if (length(test1) == 1 & length(test2) == 0 & length(test3) == 0) {GOIsRef = TotScore1; index = test1}
if (length(test1) == 0 & length(test2) == 1 & length(test3) == 0) {GOIsRef = TotScore2; index = test2}
if (length(test1) == 0 & length(test2) == 0 & length(test3) == 1) {GOIsRef = TotScore3; index = test3}
if (length(test1) == 0 & length(test2) == 0 & length(test3) == 0) {GOIsRef = NA; index = 1}
ScoresOI = c(ScoresOI, GOIsRef[index])
}


pdf(paste0('/Users/Mary/Desktop/Data/hMitoRP/PanAnalysis/', Folder, '/hMitoRP/Frame/', Experiment, '_', sizeRange,'_', modifier,'_',strand,'_Scores.pdf'), width = 5, height = 4, pointsize=10)

# par(mfrow=c(3,1))
scores <- hist(sort(c(TotScore1, TotScore2, TotScore3)),breaks=50, xlab = 'Score', main = paste0(Experiment, ' ', sizeRange, '\nDistribution of scores'), col='gray95', ylim=c(0,22))

# Label bars with genes
xx = c()
yy = c()
for (i in c(1:length(ScoresOI))) {
xx = c(xx, scores$mids[max(which(scores$breaks<ScoresOI[i]))])
}
yy=xx
for (i in c(1:length(yy))) {
yy[i]=scores$counts[max(which(scores$mids == xx[i]))]} # Take the max (of a single number) just so I will get NA instead of integer(0) which throws an error
# Find entries duplicated in both position vectors that will make gene names overlapping
overlap = which(duplicated(yy) & duplicated(xx))
for (i in overlap) {yy[i]=yy[i]+1.5}

text(x = xx+((scores$mids[2]-scores$mids[1])/4), y = yy+1, label = GOIs, pos = 3, cex = 0.6, col = cols, srt=90)

dev.off()

# Find locations of other high scores
TotScores = c(TotScore1, TotScore2, TotScore3)
refnts = c(ref1nt, ref2nt, ref3nt)
stScores = c(stScore1, stScore2, stScore3)
Diffs = c(Diff1, Diff2, Diff3)
densDiffs = c(densDiff1, densDiff2, densDiff3)

or = order(TotScores)
ordrefnts = refnts[or]
ordstScores = stScores[or]
ordDiffs = Diffs[or]
orddensDiffs = densDiffs[or]
ordTotScores = TotScores[or]

# scorecutoff = 5
# ordTotScores[ordTotScores > scorecutoff & ordTotScores %notin% ScoresOI]
# assign(paste0(Experiment,'POIs'),ordrefnts[ordTotScores > scorecutoff & ordTotScores %notin% ScoresOI])


# source('/Users/Mary/Desktop/Data/hMitoRP/PanAnalysis/ScriptsForReanalysis/adHocStartCodonScoring.R')

