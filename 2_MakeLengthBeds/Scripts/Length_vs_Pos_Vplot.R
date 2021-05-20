
library('scales')

Experiment='hMitoRP1'
libName <- 'S1'
zoom = 'start' # start none stop

path = paste0(getwd(),'/')
ylimits = c(10,40)


series5pr <- scan(paste0(path, libName, '_Mito_mRNA_lengths5p_P.bed'), list('',0,0,''))
series3pr <- scan(paste0(path, libName, '_Mito_mRNA_lengths3p_P.bed'), list('',0,0,''))
series5prMinus <- scan(paste0(path, libName, '_Mito_mRNA_lengths5p_M.bed'), list('',0,0,''))
series3prMinus <- scan(paste0(path, libName, '_Mito_mRNA_lengths3p_M.bed'), list('',0,0,''))


### Uncomment to mark the following positions
# AspPositions <- scan('/Users/Mary/Desktop/Data/hMitoRP/SequenceFiles/Asp_positions.bedGraph', list('',0,0,0),skip=1)
# LeuTPositions <- scan('/Users/Mary/Desktop/Data/hMitoRP/SequenceFiles/LeuT_positions.bedGraph', list('',0,0,0),skip=1)
# GluPositions <- scan('/Users/Mary/Desktop/Data/hMitoRP/SequenceFiles/Glu_positions.bedGraph', list('',0,0,0),skip=1)
# TransPositions <- scan('/Users/Mary/Desktop/Data/hMitoRP/SequenceFiles/Transmembrane_BaseCoords.bedGraph', list('',0,0,0),skip=1)
# G4Positions <- scan('/Users/Mary/Desktop/Data/hMitoRP/SequenceFiles/GquadruplexPred.bedGraph', list('',0,0,0),skip=1)


r = 20 # distance to plot upstream of start
t = 20 # distance to plot downstream of stop


if (zoom == 'none') {png(paste0(path, libName, '_5pr3pr.png'),units = 'in',width=10,height=26,pointsize=12,res = 300)}
if (zoom == 'start') {png(paste0(path, libName, '_5pr3pr_start.png'),units = 'in',width=3,height=26,pointsize=12,res = 300)}
if (zoom == 'stop') {png(paste0(path, libName, '_5pr3pr_stopPlus.png'),units = 'in',width=3,height=26,pointsize=12,res = 300)}
if (zoom == 'dORF4') {png(paste0(path, libName, '_5pr3pr_dORF4.png'),units = 'in',width=5,height=5,pointsize=12,res = 300)}

if (zoom != 'dORF4') {

# Plot in order of length
par(mfrow=c(14,1),cex.lab = 1.5, mar = c(2,1,1,1), oma = c(1,4,4,1)) # Margins are c(bottom, left, top, right))

genes = c('ND5','CO1','ND4','CYB','ND2','ND1','CO3','CO2','ATP6','ND6','ND3','ND4L','ATP8','antiND6')

for (gene in genes)
{
if (gene == "ND1"){range = c((3307-r):(4262+t))} else if (gene == "ND2"){range = c((4470-r):(5511+t))} else if (gene == "CO1"){range = c((5904-r):(7445+t))} else if (gene == "CO2"){range = c((7586-r):(8269+t))} else if (gene == "ATP8"){range = c((8366-r):(8526+t))} else if (gene == "ATP6"){range = c((8527-r):(9207+t))} else if (gene == "CO3"){range = c((9207-r):(9990+t))} else if (gene == "ND3"){range = c((10059-r):(10404+t))} else if (gene == "ND4L"){range = c((10470-r):(10766+t))} else if (gene == "ND4"){range = c((10760-r):(12137+t))} else if (gene == "ND5"){range = c((12337-r):(14148+t))} else if (gene == "ND6"){range = c((14149--t):(14673+r))} else if (gene == "antiND6"){range = c((14149-r):(14673+t))} else if (gene == "CYB"){range = c((14747-r):(15887+t))} else if (gene == "ND4L-ND4"){range = c((10470-r):(12137+t))} else if (gene == "ATP8-ATP6"){range = c((8366-r):(9207+t))} else {range = c((start-r):(stop+t))}
length = (tail(range, n=1)-t) - (range[1]+r)
ORF = (range[1]+r):(tail(range, n=1)-t)

if (zoom == 'start'){
if (gene != 'ND6'){range = range[1]:range[1+2.5*r]}
if (gene == 'ND6'){range = (tail(range, n=1)-2.5*t):(tail(range,n=1))}
}
if (zoom == 'stop'){
if (gene != 'ND6'){range = (tail(range, n=1)-2.5*t):(tail(range,n=1))}
if (gene == 'ND6'){range = range[1]:range[1+2.5*r]}
}


if (gene == 'ND5'){rangelength = length(range)}

# Plus strand
if (gene != 'ND6')
{
x5pr=series5pr[[2]][series5pr[[2]] %in% range]
y5pr=series5pr[[3]][series5pr[[2]] %in% range]
x3pr=series3pr[[2]][series3pr[[2]] %in% range]
y3pr=series3pr[[3]][series3pr[[2]] %in% range]

### Uncomment to mark the following positions
# Asp=(AspPositions[[3]]-1)[AspPositions[[3]] %in% range]
# Leu=(LeuTPositions[[3]]-1)[LeuTPositions[[3]] %in% range]
# Glu=(GluPositions[[3]]-1)[GluPositions[[3]] %in% range]
# TransStart=(TransPositions[[2]])[TransPositions[[2]] %in% ORF]
# TransEnd=(TransPositions[[3]])[TransPositions[[3]] %in% ORF]


if (zoom == 'none' | zoom == 'start') {start1 = head(range+r,1)} else {start1 = ''}
if (zoom == 'none' | zoom == 'stop') {stop1 = tail(range-t,1)} else {stop1 = ''}
}

# Minus strand
if (gene == 'ND6')
{
x5pr=series5prMinus[[2]][series5prMinus[[2]] %in% range]
y5pr=series5prMinus[[3]][series5prMinus[[2]] %in% range]
x3pr=series3prMinus[[2]][series3prMinus[[2]] %in% range]
y3pr=series3prMinus[[3]][series3prMinus[[2]] %in% range]

### Uncomment to mark the following positions
# Asp=(AspPositions[[3]]-1)[AspPositions[[3]] %in% range]
# Leu=(LeuTPositions[[3]]-1)[LeuTPositions[[3]] %in% range]
# Glu=(GluPositions[[3]]-1)[GluPositions[[3]] %in% range]
# TransStart=(TransPositions[[2]])[TransPositions[[2]] %in% ORF]
# TransEnd=(TransPositions[[3]])[TransPositions[[3]] %in% ORF]
# G4Start=(G4Positions[[2]])[G4Positions[[2]] %in% ORF]
# G4End=(G4Positions[[3]])[G4Positions[[3]] %in% ORF]


if (zoom == 'none' | zoom == 'start') {start1 = tail(range,1)-t} else {start1 = ''}
if (zoom == 'none' | zoom == 'stop') {stop1 = head(range+r,1)} else {stop1 = ''}
}

# Adjust transparency for points depending on density of signal
total = sum(y5pr)/length 
if (total < 100){transparency = .12}
if (total >= 100 & total < 240){transparency = .07}
if (total >= 240 & total < 500){transparency = .05}
if (total >= 500 & total < 700){transparency = .03}
if (total >= 700 & total < 1000){transparency = .02}
if (total >= 1000 & total < 1500){transparency = .01}
if (total >= 1500 & total < 5000){transparency = .008}
if (total >= 5000){transparency = .005}


# Plus strand
if (gene != 'ND6')
{

ptsize = .4

## Plot it, reordering rows so that densest points are plotted on top
plot(x5pr,jitter(y5pr,2.5), col = alpha('black', transparency), ylim = ylimits, xaxt = 'n', xlab = gene, ylab = 'length', pch = 16, cex = ptsize, xlim = c(head(range,1), tail(range,1)+(rangelength-length(range)))) #xlim = c(head(range,1), tail(range,1)+(1811-length)))  c(head(range,1)+(length-200), tail(range,1))
points(x3pr,jitter(y3pr,2.5), col = alpha('firebrick3', transparency), xlab = gene, pch = 16, cex = ptsize)

### Uncomment to mark the following positions
# abline(v=Asp, col = 'black', lwd = .3)
# abline(v=Leu, col = 'black', lwd = .3)
# abline(v=Glu, col = 'black', lwd = .3, lty = 2)
# legend('topright', legend = c('Asp', 'Glu'), col = c('black','black'), lty = c(1,2), bty = 'n')

if (gene == 'antiND6') {
abline(v=14196,col = 'green', lwd = .3)
abline(v=14208,col = 'red', lwd = .3)
abline(v=14256,col = 'red', lwd = .3)
}

}

# Minus strand
if (gene == 'ND6')
{
## Plot it, reordering rows so that densest points are plotted on top
plot(x5pr,jitter(y5pr,2.5), col = alpha('black', transparency), ylim = ylimits, xaxt = 'n', ylab = 'length', pch = 16, cex = ptsize, xlim = rev(c(head(range,1)-(rangelength-length(range)), tail(range,1)))) 
points(x3pr,jitter(y3pr,2.5), col = alpha('firebrick3', transparency), xlab = gene, pch = 16, cex = ptsize)

### Uncomment to mark the following positions
# abline(v=14666, col = 'black', lwd = .3)
# abline(v=14456, col = 'black', lwd = .3)
# segments(14666,0,14666,15, col='black', lwd = .3)
# segments(14456,0,14456,15, col='black', lwd = .3)

# abline(v=Asp, col = 'black', lwd = .3)
# abline(v=Leu, col = 'black', lwd = .3)
# abline(v=Glu, col = 'black', lwd = .3, lty = 2)
# G4=c()
# for (i in 1:length(G4Start))
# {
# G4=c(G4,G4Start[i]:G4End[i])
# }

}
axis(1, at = c(start1, stop1), labels = FALSE, col= 'black') 


### Uncomment to mark the following positions
# axis(1, at = c(G4), labels = FALSE, col= 'forestgreen') # Make ticks at every position of G4, don't label
# axis(1, at = c(start1, stop1),labels = c('start','stop'), tick = FALSE) # Put labels at +30 to center


title(gene, line = -2) # So that the title overlaps plot, to save space

# Relabel x axis : Choose line to uncomment below based on what you're labeling
axis(1, at = c(start1, stop1), labels = c('start', 'stop'), las = 2)
# axis(1, at = c(start1, Asp, Leu, stop1), labels = c('start', rep('Asp', times = length(Asp)),rep('Leu', times = length(Leu)),'stop'), las = 2)

# trans=c()
# for (i in 1:length(TransStart)) {
# trans=c(trans,TransStart[i]:TransEnd[i])
# }
# axis(1, at = c(start1, trans, stop1), labels = FALSE) # Make ticks at every position of transmembrane domain, don't label
# axis(1, at = c(start1, trans+30, stop1), labels = c('start', rep('Trans', times = length(trans)),'stop'), tick = FALSE) # Put labels at +30 to center


}


mtext("", side = 1, line = 2, outer = TRUE) # line is how close to plot (0-4)
mtext("Length", side = 2, line = 2, outer = TRUE)
mtext(libName, side = 3, line = 2, outer = TRUE)


dev.off()
}



if (zoom == 'dORF4') {
range = c((14196-r):(14200+t))
x5pr=series5pr[[2]][series5pr[[2]] %in% range]
y5pr=series5pr[[3]][series5pr[[2]] %in% range]
x3pr=series3pr[[2]][series3pr[[2]] %in% range]
y3pr=series3pr[[3]][series3pr[[2]] %in% range]
start1 = 14196
stop1 = 14259

plot(x5pr,jitter(y5pr,2.5), col = alpha('black', .2), ylim = ylimits, xaxt = 'n', xlab = 'dORF4', ylab = 'length', pch = 16, cex = .8)
points(x3pr,jitter(y3pr,2.5), col = alpha('firebrick3', .3), xlab = 'dORF4', pch = 16, cex = .8)
axis(1, at = c(start1, stop1), labels = FALSE, col= 'black') 

dev.off()
}


# source('./2_MakeLengthBeds/Scripts/Length_vs_Pos_Vplot.R')