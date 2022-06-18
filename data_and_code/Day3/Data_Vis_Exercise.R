#Read in the input file
bulk<-read.delim(file="tall_allele_counts.txt", head=T, sep="\t")
summary(bulk)

################################################################################
### Look at the distribution of parent 1 allele frequencies with a histogram ###
################################################################################

#Calculate alle frequency rather than count - why??
bulk$p1_freq<-(bulk$p1_count / (bulk$p1_count + bulk$p2_count))
summary(bulk)

#Open an output PDF file
pdf.options(family="Helvetica")
pdf("p1_allele_freq_hist_large.pdf", pointsize=10)

#Set graphical parameters to plot four plots per page
par(mfrow=c(2,2))

#Make a histogram of the allele frequencies
hist(bulk$p1_freq)

#Specify axis labels and color
hist(bulk$p1_freq, xlab="Parent 1 Allele Frequency", ylab="Count", main="Large Plant Pool", col="gray")

#Add lines showing the expected (red) and observed (blue) mean
hist(bulk$p1_freq, xlab="Parent 1 Allele Frequency", ylab="Count", main="Large Plant Pool", col="gray")
lines(c(0.5, 0.5),c(0, 306884), col="red", lty=1, lwd=2)
lines(c(0.4820, 0.4820), c(0, 306884), col="blue", lty=1, lwd=2)

#Modify the axis and breakpoints
bins=seq(0,1,by=0.025)
hist(bulk$p1_freq, xlim=c(0,1), ylim=c(0,200000), xlab="Parent 1 Allele Frequency", ylab="Count", main="Large Plant Pool", col="gray", axes=F, breaks=bins)
axis(1,at=seq(0,1,0.1), labels=T, pos=0)
axis(2,at=seq(0,200000,50000), labels=T, pos=0)
lines(c(0.5, 0.5),c(0, 188903), col="red", lty=1, lwd=2)
lines(c(0.4820, 0.4820), c(0, 188903), col="blue", lty=1, lwd=2)

#Disconnect from the output file
dev.off()

############################################################
### Look at the allele frequencies along each chromosome ###
############################################################

#Make subsets with data for each chromosome
chr1_bulk_sub<-subset(bulk,chr==1)
chr2_bulk_sub<-subset(bulk,chr==2)
chr3_bulk_sub<-subset(bulk,chr==3)
chr4_bulk_sub<-subset(bulk,chr==4)
chr5_bulk_sub<-subset(bulk,chr==5)
chr6_bulk_sub<-subset(bulk,chr==6)
chr7_bulk_sub<-subset(bulk,chr==7)
chr8_bulk_sub<-subset(bulk,chr==8)
chr9_bulk_sub<-subset(bulk,chr==9)
chr10_bulk_sub<-subset(bulk,chr==10)

#Open an output PDF file
pdf.options(family="Helvetica")
pdf("p1_allele_frequency_large_v1.pdf", pointsize=10)

#Set graphical parameters to plot five plots per page
par(mfrow=c(5,1))

#Make a scatter plot for chr1-chr5
plot(chr1_bulk_sub$pos, chr1_bulk_sub$p1_freq)
plot(chr2_bulk_sub$pos, chr2_bulk_sub$p1_freq)
plot(chr3_bulk_sub$pos, chr3_bulk_sub$p1_freq)
plot(chr4_bulk_sub$pos, chr4_bulk_sub$p1_freq)
plot(chr5_bulk_sub$pos, chr5_bulk_sub$p1_freq)

#Disconnect from the output file
dev.off()

#####################################################################

#Now lets make it pretty

#Load necessary libraries and make variables
install.packages("Hmisc")
library(Hmisc)

pos<-c(0,50000000,100000000,150000000,200000000,250000000,300000000)
lab<-c(0,50,100,150,200,250,300)
mark<-c(0,0.5,1)

#Open output PDF file and set file parameters
pdf.options(family="Helvetica")
pdf("p1_allele_frequency_large_v2.pdf", pointsize=10)
par(mfrow=c(5,1), mar=c(5,4,0.2,0.2))

#Plot chromosome 1
plot(chr1_bulk_sub$pos, chr1_bulk_sub$p1_freq, ylim=c(-0.1,1.1), xlim=c(1,301354135), xlab="Chr1 Position (Mb)", ylab="P1 Allele Freq", xaxt="n", yaxt="n", pch=21, col="#009966", bg="#00996625", cex=0.75)

#Add axis to the plot
axis(1,at=pos,labels=lab)
axis(2,at=mark,labels=mark, las=1)
minor.tick(nx=4, ny=0, tick.ratio=0.5)

#Add boxes around selected regions
rect(455861, -0.1, 13854546, 1.1, col="#99330025", border=NA) 
rect(60255361, -0.1, 69849341, 1.1, col="#99330025", border=NA) 
rect(90619709, -0.1, 90959095, 1.1, col="#99330025", border=NA)
rect(170877659, -0.1, 171253827, 1.1, col="#99330025", border=NA)

#Plot chromosome 2-5
plot(chr2_bulk_sub$pos, chr2_bulk_sub$p1_freq, ylim=c(-0.1,1.1), xlim=c(1,301354135), xlab="Chr2 Position (Mb)", ylab="P1 Allele Freq", xaxt="n", yaxt="n", pch=21, col="#009966", bg="#00996625", cex=0.75)
axis(1,at=pos,labels=lab)
axis(2,at=mark,labels=mark, las=1)
minor.tick(nx=4, ny=0, tick.ratio=0.5)
rect(55428942, -0.1, 55961964, 1.1, col="#99330025", border=NA) 
rect(114448932, -0.1, 114957517, 1.1, col="#99330025", border=NA) 
rect(126034984, -0.1, 126475373, 1.1, col="#99330025", border=NA)

plot(chr3_bulk_sub$pos, chr3_bulk_sub$p1_freq, ylim=c(-0.1,1.1), xlim=c(1,301354135), xlab="Chr3 Position (Mb)", ylab="P1 Allele Freq", xaxt="n", yaxt="n", pch=21, col="#009966", bg="#00996625", cex=0.75)
axis(1,at=pos,labels=lab)
axis(2,at=mark,labels=mark, las=1)
minor.tick(nx=4, ny=0, tick.ratio=0.5)
rect(11719757, -0.1, 12002922, 1.1, col="#99330025", border=NA) 
rect(12157474, -0.1, 16907255, 1.1, col="#99330025", border=NA) 
rect(108332482, -0.1, 108658335, 1.1, col="#99330025", border=NA)
rect(144034451, -0.1, 144169093, 1.1, col="#99330025", border=NA) 
rect(149441643, -0.1, 149680694, 1.1, col="#99330025", border=NA) 

plot(chr4_bulk_sub$pos, chr4_bulk_sub$p1_freq, ylim=c(-0.1,1.1), xlim=c(1,301354135), xlab="Chr4 Position (Mb)", ylab="P1 Allele Freq", xaxt="n", yaxt="n", pch=21, col="#009966", bg="#00996625", cex=0.75)
axis(1,at=pos,labels=lab)
axis(2,at=mark,labels=mark, las=1)
minor.tick(nx=4, ny=0, tick.ratio=0.5)
rect(29836106, -0.1, 30360486, 1.1, col="#99330025", border=NA) 
rect(72085352, -0.1, 72456250, 1.1, col="#99330025", border=NA) 
rect(88623650, -0.1, 88964839, 1.1, col="#99330025", border=NA)

plot(chr5_bulk_sub$pos, chr5_bulk_sub$p1_freq, ylim=c(-0.1,1.1), xlim=c(1,301354135), xlab="Chr5 Position (Mb)", ylab="P1 Allele Freq", xaxt="n", yaxt="n", pch=21, col="#009966", bg="#00996625", cex=0.75)
axis(1,at=pos,labels=lab)
axis(2,at=mark,labels=mark, las=1)
minor.tick(nx=4, ny=0, tick.ratio=0.5)
rect(130986673, -0.1, 157960115, 1.1, col="#99330025", border=NA) 

#Close the output file
dev.off()

############################################################################################
### Summarize all significant regions from the small and large plant pools in one figure ### 
############################################################################################

#Open output PDF file and set file parameters
pdf.options(family="Helvetica")
pdf("summary_figure.pdf", pointsize=12)
par(mar=c(0.5,0,0,0))

#Make chromosomes to scale
plot(c(1,1), c(0, 301354135), type="l", lwd=30, col= "gray85", lend=1, xlim=c(0.5,10.5), axes=F, xlab="", ylab="", main="")
lines(c(2,2), c(0, 237068873), type="l", lwd=30, col="gray85", lend=1)
lines(c(3,3), c(0, 232140174), type="l", lwd=30, col="gray85", lend=1)
lines(c(4,4), c(0, 241473504), type="l", lwd=30, col="gray85", lend=1)
lines(c(5,5), c(0, 217872852), type="l", lwd=30, col="gray85", lend=1)
lines(c(6,6), c(0, 169174353), type="l", lwd=30, col="gray85", lend=1)
lines(c(7,7), c(0, 176764762), type="l", lwd=30, col="gray85", lend=1)
lines(c(8,8), c(0, 175793759), type="l", lwd=30, col="gray85", lend=1)
lines(c(9,9), c(0, 156750706), type="l", lwd=30, col="gray85", lend=1)
lines(c(10,10), c(0, 150189435), type="l", lwd=30, col="gray85", lend=1)

#Label the chromosomes
l1 = c(1,2,3,4,5,6,7,8,9,10)
mtext(l1, side=1, line=-0.5, at=l1, col="black", cex=1.25)

#Plot each region
lines(c(0.88,0.88), c(300898274,287499589), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(0.88,0.88), c(130476476,130100308), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(0.88,0.88), c(241098774,231504794), lwd=8.5, type="l", col="#FF0000", lend=1)
lines(c(0.88,0.88), c(210734426,210395040), lwd=8.5, type="l", col="#FF0000", lend=1)
lines(c(1.12,1.12), c(210734426,210495040), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(1.12,1.12), c(123370025,122983135), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(1.12,1.12), c(232004271,231680910), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(1.88,1.88), c(181639931,181106909), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(1.88,1.88), c(122619941,122111356), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(1.88,1.88), c(111033889,110593500), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(2.12,2.12), c(142098450,141868529), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(2.12,2.12), c(231543366,231170530), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(2.12,2.12), c(122619941,122311356), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(2.88,2.88), c(220420417,220137252), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(2.88,2.88), c(219982700,215232919), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(2.88,2.88), c(88105723,87971081), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(2.88,2.88), c(82698531,82459480), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(2.88,2.88), c(123807692,123481839), lwd=8.5, type="l", col="#FF0000", lend=1)
lines(c(3.12,3.12), c(219982700,215072919), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(3.12,3.12), c(195304261,195253553), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(3.12,3.12), c(123807692,123774839), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(3.88,3.88), c(211637398,211113018), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(3.88,3.88), c(152849854,152508665), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(3.88,3.88), c(169388152,169017254), lwd=8.5, type="l", col="#FF0000", lend=1)
lines(c(4.12,4.12), c(203373404,203189599), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(4.12,4.12), c(169388152,169351063), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(4.88,4.88), c(86886179,59912737), lwd=8.5, type="l", col="#FF0000", lend=1)
lines(c(5.12,5.12), c(217645967,217341173), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(5.12,5.12), c(86886179,84212737), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(5.88,5.88), c(131108794,130861870), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(5.88,5.88), c(93657391,93399663), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(5.88,5.88), c(91803114,90906850), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(6.12,6.12), c(61360380,61040004), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(6.12,6.12), c(131108794,130861870), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(6.12,6.12), c(94588543,94241180), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(6.88,6.88), c(158856588,156911203), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(6.88,6.88), c(98351300,98107634), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(6.88,6.88), c(66675936,65687018), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(6.88,6.88), c(43168177,42679906), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(6.88,6.88), c(140298474,139354288), lwd=8.5, type="l", col="#FF0000", lend=1)
lines(c(7.12,7.12), c(158856588,156911203), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(7.12,7.12), c(98351300,98107634), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(7.12,7.12), c(66675936,65687018), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(7.12,7.12), c(140298474,139354288), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(7.88,7.88), c(155838138,155642611), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(7.88,7.88), c(118412484,116507285), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(7.88,7.88), c(64793524,64365266), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(7.88,7.88), c(43400815,43012304), lwd=8.5, type="l", col="#FF0000", lend=1)
lines(c(8.12,8.12), c(118412484,113407285), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(8.12,8.12), c(19693367,19390135), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(8.12,8.12), c(64793524,64765266), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(8.88,8.88), c(138106903,137729896), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(9.12,9.12), c(23013476,22976644), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(9.88,9.88), c(24426311,24033867), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(9.88,9.88), c(44731,2162), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(9.88,9.88), c(128668058,127487002), lwd=8.5, type="l", col="#FF0000", lend=1)
lines(c(10.12,10.12), c(24426311,24023867), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(10.12,10.12), c(44731,22162), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(10.12,10.12), c(128668058,127667002), lwd=8.5, type="l", col="#0000FF", lend=1)

#Make a legend for the figure
legend(7.3, 301354135, c("Large Plant Pool", "Small Plant Pool"), col=c("#FF0000", "#0000FF"), lwd=c(4,4))

#Close the output file
dev.off()

#####################################################################

#Summarize all significant regions from the small and large plant pools in one figure with a scale

#Open output PDF file and set file parameters
pdf.options(family="Helvetica")
pdf("summary_figure_v2.pdf", pointsize=12)
par(mar=c(2,4,0,0))

#Make chromosomes to scale
plot(c(1,1), c(0, 301354135), type="l", lwd=30, col= "gray85", lend=1, xlim=c(0.5,10.5), axes=F, xlab="", ylab="Position (Mb)", main="")
lines(c(2,2), c(0, 237068873), type="l", lwd=30, col="gray85", lend=1)
lines(c(3,3), c(0, 232140174), type="l", lwd=30, col="gray85", lend=1)
lines(c(4,4), c(0, 241473504), type="l", lwd=30, col="gray85", lend=1)
lines(c(5,5), c(0, 217872852), type="l", lwd=30, col="gray85", lend=1)
lines(c(6,6), c(0, 169174353), type="l", lwd=30, col="gray85", lend=1)
lines(c(7,7), c(0, 176764762), type="l", lwd=30, col="gray85", lend=1)
lines(c(8,8), c(0, 175793759), type="l", lwd=30, col="gray85", lend=1)
lines(c(9,9), c(0, 156750706), type="l", lwd=30, col="gray85", lend=1)
lines(c(10,10), c(0, 150189435), type="l", lwd=30, col="gray85", lend=1)

#Label the chromosomes
l1 = c(1,2,3,4,5,6,7,8,9,10)
mtext(l1, side=1, line=-0.5, at=l1, col="black")
mtext("Chromosome", side=1, line=1.0, at=5, col="black")

#Plot each region
lines(c(0.88,0.88), c(300898274,287499589), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(0.88,0.88), c(130476476,130100308), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(0.88,0.88), c(241098774,231504794), lwd=8.5, type="l", col="#FF0000", lend=1)
lines(c(0.88,0.88), c(210734426,210395040), lwd=8.5, type="l", col="#FF0000", lend=1)
lines(c(1.12,1.12), c(210734426,210495040), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(1.12,1.12), c(123370025,122983135), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(1.12,1.12), c(232004271,231680910), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(1.88,1.88), c(181639931,181106909), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(1.88,1.88), c(122619941,122111356), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(1.88,1.88), c(111033889,110593500), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(2.12,2.12), c(142098450,141868529), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(2.12,2.12), c(231543366,231170530), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(2.12,2.12), c(122619941,122311356), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(2.88,2.88), c(220420417,220137252), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(2.88,2.88), c(219982700,215232919), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(2.88,2.88), c(88105723,87971081), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(2.88,2.88), c(82698531,82459480), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(2.88,2.88), c(123807692,123481839), lwd=8.5, type="l", col="#FF0000", lend=1)
lines(c(3.12,3.12), c(219982700,215072919), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(3.12,3.12), c(195304261,195253553), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(3.12,3.12), c(123807692,123774839), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(3.88,3.88), c(211637398,211113018), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(3.88,3.88), c(152849854,152508665), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(3.88,3.88), c(169388152,169017254), lwd=8.5, type="l", col="#FF0000", lend=1)
lines(c(4.12,4.12), c(203373404,203189599), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(4.12,4.12), c(169388152,169351063), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(4.88,4.88), c(86886179,59912737), lwd=8.5, type="l", col="#FF0000", lend=1)
lines(c(5.12,5.12), c(217645967,217341173), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(5.12,5.12), c(86886179,84212737), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(5.88,5.88), c(131108794,130861870), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(5.88,5.88), c(93657391,93399663), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(5.88,5.88), c(91803114,90906850), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(6.12,6.12), c(61360380,61040004), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(6.12,6.12), c(131108794,130861870), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(6.12,6.12), c(94588543,94241180), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(6.88,6.88), c(158856588,156911203), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(6.88,6.88), c(98351300,98107634), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(6.88,6.88), c(66675936,65687018), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(6.88,6.88), c(43168177,42679906), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(6.88,6.88), c(140298474,139354288), lwd=8.5, type="l", col="#FF0000", lend=1)
lines(c(7.12,7.12), c(158856588,156911203), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(7.12,7.12), c(98351300,98107634), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(7.12,7.12), c(66675936,65687018), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(7.12,7.12), c(140298474,139354288), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(7.88,7.88), c(155838138,155642611), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(7.88,7.88), c(118412484,116507285), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(7.88,7.88), c(64793524,64365266), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(7.88,7.88), c(43400815,43012304), lwd=8.5, type="l", col="#FF0000", lend=1)
lines(c(8.12,8.12), c(118412484,113407285), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(8.12,8.12), c(19693367,19390135), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(8.12,8.12), c(64793524,64765266), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(8.88,8.88), c(138106903,137729896), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(9.12,9.12), c(23013476,22976644), lwd=8.5, type="l", col="#0000FF", lend=1)
lines(c(9.88,9.88), c(24426311,24033867), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(9.88,9.88), c(44731,2162), lwd=8.5, type="l", col="#FF000075", lend=1)
lines(c(9.88,9.88), c(128668058,127487002), lwd=8.5, type="l", col="#FF0000", lend=1)
lines(c(10.12,10.12), c(24426311,24023867), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(10.12,10.12), c(44731,22162), lwd=8.5, type="l", col="#0000FF75", lend=1)
lines(c(10.12,10.12), c(128668058,127667002), lwd=8.5, type="l", col="#0000FF", lend=1)

#Make a legend for the figure
legend(7.3, 301354135, c("Large Plant Pool", "Small Plant Pool"), col=c("#FF0000", "#0000FF"), lwd=c(4,4))

#Add y-axis
axis(2, at=c(1354135, 51354135, 101354135, 151354135, 201354135, 251354135, 301354135), labels=c(300, 250, 200, 150, 100, 50, 0), las=1)

#Close the output file
dev.off()

