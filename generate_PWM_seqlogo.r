#list.files()

# Usage: R --no-save --slave --args inputfilename background_perc_A background_perc_C background_perc_G background_perc_T < generate_PWM_seqlogo.R
#====================================================

arg = commandArgs(trailingOnly=T);
input<-arg[1];
pA<-as.numeric(arg[2]);  # background percentage of base A (0-100), equiprobable composition=25
pC<-as.numeric(arg[3]);  # background percentrage of base C (0-100), equiprobable composition=25
pG<-as.numeric(arg[4]);  # background percentage of base G (0-100), equiprobable composition=25
pT<-as.numeric(arg[5]);  # background percentage of base T (0-100), equiprobable composition=25
# e.g.,
#-- XL-----
#pA=20.1
#pC=18.7
#pG=34.9
#pT=26.2

## Step 1: Read in input file:
# For example:
#input="cas9_XL_blunt_3bp_nucleotide_count.txt"
d<-read.table(input, sep="\t", header=T)
#d
# Pos	A	C	G	T
#1	53857	41050	70235	58835
#2	51127	11151	135767	25932
#3	29851	9243	182374	2508
#4	39190	48707	77295	58785
#5	47316	42419	73508	60734
#6	41732	35421	84186	62638
#7	40518	45406	80428	57625
#8	48224	40038	79349	56366
#9	46694	42718	76774	57791
#10	46904	40698	78390	57985

## Step 2: Calculate enrichment score as: log2(observed/background)
d$lA<-log(d$A/((d$A+d$C+d$G+d$T)*(pA/100)), 2)
d$lC<-log(d$C/((d$A+d$C+d$G+d$T)*(pC/100)), 2)
d$lG<-log(d$G/((d$A+d$C+d$G+d$T)*(pG/100)), 2)
d$lT<-log(d$T/((d$A+d$C+d$G+d$T)*(pT/100)), 2)

## Step 3: Plot log ratio of nucleotide frequencies at each position
## Follow the color scheme of weblog:
## G: orange; pch=71
## T: red; pch=84
## C: blue; pch=67
## A: green; pch=65

pdffile<-paste(input, "position_weight_matrix_logo.pdf", sep="_")

pdf(file=pdffile, width=3.4, height=2.8)

par(mar=c(2, 1.5, 0.1,0.1)+1)
plot(d$Pos, d$lA, ylim=c(-7, 3), pch=65, col="green",  xlab="position", ylab="score", cex=1.5, xaxt="n", yaxt="n",ann = FALSE)
points(d$Pos, d$lC, pch=67, col="blue", cex=1.5)
points(d$Pos, d$lT, pch=84,  col="red", cex=1.5)
points(d$Pos, d$lG, pch=71, col="orange", cex=1.5)
abline(h=0, lty=3)

axis(1, at=d$Pos, tck = -0.03,label = rep("", 10), cex.axis=0.8) # first add the tick marks
axis(1, at =d$Pos, labels=d$Pos, line=-0.5, lwd=0, cex.axis=0.8) # then add the labels
axis(2, at=-7:3, tck = -0.03,label = rep("", 11), cex.axis=0.8) # first add the tick marks
axis(2, at =-7:3, labels=seq(-7,3,1), line=-0.5, lwd=0, cex.axis=0.8) # then add the labels
title(ylab = "score", line = 1.5)
title(xlab = "position", line = 1.5)

dev.off()

