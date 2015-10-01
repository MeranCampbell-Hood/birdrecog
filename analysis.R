rm(list=ls())
library(tuneR)
library(seewave)
library(tidyr)
oldpar <- par()

bestguesses <- function(nameddistancematrix){
inorder <- as.data.frame(apply(as.matrix(nameddistancematrix), 1, function(x){(names(x)[order(x)])}))
topref <- data.frame(t(apply(inorder, 2, function(x){(x[grep("reference",x)])[1:3]})))
topref <- topref[grep("test|self",row.names(topref)),]
topref$bird <- row.names(topref)
return(topref)
}


usethese <- read.csv("~/nbu/birding/Sep/usethese.csv")
prepFiles <- list.files(path="xenoprep", pattern=".wav$")
exmplars <- usethese
exmplars$wavs <- paste("~/nbu/birding/Sep/xenoprep/",exmplars$id, "_mod.wav", sep="")
#listOfSoundsE <- lapply(exmplars$wavs, readWave)
#specia <- lapply(listOfSoundsE, function(x){meanspec(x, f=f, plot = FALSE, identity=TRUE)})
#save(specia, file="spectra.RData")
load("spectra.RData")
referenceset <- sort(which(usethese$use == "reference"))
referenceandtestset <- sort(c(which(usethese$use == "reference"),which(usethese$use == "test")))
referenceandmyset <- sort(c(which(usethese$use == "reference"),which(usethese$use == "self")))
SongA <- rep(1:nrow(usethese), each= length(referenceset)) #all songs
SongB <- rep(referenceset, length.out=length(SongA)) #comparison set
pairs <- data.frame(SongA,SongB)

########## get kl.distance closeness of only references
kpairs <- pairs[pairs$SongA %in% referenceset,]
kpairs$kldist <- 99999
for (i in 1:nrow(kpairs)){
  item1 <- (specia[[kpairs[i,"SongA"]]])
  item2 <- (specia[[kpairs[i,"SongB"]]])
  AtoB <- (kl.dist(item1, item2))$D
  BtoA <- (kl.dist(item2, item1))$D
  kpairs$kldist[i] <- (AtoB + BtoA) / 2
}
kpairs$kldist <- kpairs$kldist / max(kpairs$kldist)
klwide <- spread(kpairs, key=SongB, value=kldist)
row.names(klwide) <- paste(exmplars$en[referenceset], exmplars$use[referenceset])
klwide$SongA <- NULL
kd <- dist(as.matrix(klwide))
khc <- hclust(kd) 
png("fig3_refsamples_only_raw.png", width=6, height=4, units="in", res=300)
plot(khc, cex=0.5, xlab="Species")
dev.off()
klwidePCA <- prcomp(klwide, scale=TRUE)
summary(klwidePCA) # first 3 components explain 89.87% of variation
klpd <- dist(as.matrix(klwidePCA$x[,1:3]))
klphc <- hclust(klpd)
png("fig4_refsamples_only_pca.png", width=6, height=4, units="in", res=300)
plot(klphc, cex=0.5, xlab="Species")
dev.off()
firstcomponent <- klwidePCA$x[,1]
secondcomponent <- klwidePCA$x[,2]
xmin <- min(firstcomponent) - 0.5
xmax <- max(firstcomponent) + 2
ymin <- min(secondcomponent) - 0.5
ymax <- max(secondcomponent) + 0.5
pchval <- rep(15, times=nrow(klwide))
gval <- (klwidePCA$x[,3]-min(klwidePCA$x[,3]))/(max(klwidePCA$x[,3])-min(klwidePCA$x[,3]))
png("fig5_refsamples_only_pcascatter.png", width=6, height=4, units="in", res=300)
plot(firstcomponent, secondcomponent, pch=pchval, cex=0.8, col=grey(gval),  frame.plot=F, ylim=c(ymin,ymax), xlim=c(xmin,xmax))
par(mar=c(4,4,1,1))
pchval <- rep(0, times=nrow(klwide))
points(firstcomponent, secondcomponent, pch=pchval, cex=0.8)
text(x = firstcomponent, y = secondcomponent, labels= exmplars$en[referenceset], cex=0.5,pos=4)
text(x=xmax, y=ymax, labels="point shading is 3rd principle component", cex=0.5,pos=2, font=3)
legend(x=xmax-1, y=ymax-0.5, "Species", pch=0, cex=0.5)
dev.off()
par(mar=oldpar$mar)

###########kldist of test samples vs reference
kpairs <- pairs[pairs$SongA %in% referenceandtestset,]
kpairs$kldist <- 99999
for (i in 1:nrow(kpairs)){
  item1 <- (specia[[kpairs[i,"SongA"]]])
  item2 <- (specia[[kpairs[i,"SongB"]]])
  AtoB <- (kl.dist(item1, item2))$D
  BtoA <- (kl.dist(item2, item1))$D
  kpairs$kldist[i] <- (AtoB + BtoA) / 2
}
kpairs$kldist <- kpairs$kldist / max(kpairs$kldist)
klwide <- spread(kpairs, key=SongB, value=kldist)
row.names(klwide) <- paste(exmplars$en[as.integer(row.names(klwide))], exmplars$use[as.integer(row.names(klwide))])
klwide$SongA <- NULL
kd <- dist(as.matrix(klwide))
kdpicks <- bestguesses(kd)
write.csv(kdpicks, file="kdpicks.csv", row.names = FALSE)
khc <- hclust(kd)          
plot(khc, cex=0.5, main="KL distance. Sound similarities of example recordings", xlab="Species")
klwidePCA <- prcomp(klwide, scale=TRUE)
summary(klwidePCA) # first 3 components explain 89.79% of variation
klpd <- dist(as.matrix(klwidePCA$x[,1:3]))
kpcadpicks <- bestguesses(klpd)
write.csv(kpcadpicks, file="kdpcapicks.csv", row.names = FALSE)
klphc <- hclust(klpd)
plot(klphc, cex=0.5, main="KL distance. Sound similarities of example recordings\n based on first 3 Principle Components", xlab="Species")

firstcomponent <- klwidePCA$x[,1]
secondcomponent <- klwidePCA$x[,2]
xmin <- min(firstcomponent) - 0.5
xmax <- max(firstcomponent) + 2
ymin <- min(secondcomponent) - 0.5
ymax <- max(secondcomponent) + 0.5
pchval <- rep(15, times=nrow(klwide))
pchval[grep("test", row.names(klwidePCA$x))] <- 19 
gval <- (klwidePCA$x[,3]-min(klwidePCA$x[,3]))/(max(klwidePCA$x[,3])-min(klwidePCA$x[,3]))
png("fig6_reftest_klpcascatter.png", width=6, height=4, units="in", res=300)
par(mar=c(4,4,1,1))
plot(firstcomponent, secondcomponent, pch=pchval, cex=0.8, col=grey(gval),  frame.plot=F, ylim=c(ymin,ymax), xlim=c(xmin,xmax))
pchval <- rep(0, times=nrow(klwide))
pchval[grep("test", row.names(klwidePCA$x))] <- 1 
points(firstcomponent, secondcomponent, pch=pchval, cex=0.8)
legendtext <- row.names(klwidePCA$x)
legendtext <- gsub(" test", "", legendtext)
legendtext <- gsub(" reference", "", legendtext)
text(x = firstcomponent, y = secondcomponent, labels= legendtext, cex=0.5,pos=4)
text(x=xmax, y=ymax, labels="point shading is 3rd principle component", cex=0.5,pos=2, font=3)
legend(x=xmax-1, y=ymax-0.5, c("Reference","Test"), pch=c(0,1), cex=0.5)
dev.off()
par(mar=oldpar$mar)

#########logspec.dist of test samples vs reference

##############
lpairs <- pairs[pairs$SongA %in% referenceandtestset,]
lpairs$kldist <- 99999
for (i in 1:nrow(lpairs)){
  item1 <- (specia[[lpairs[i,"SongA"]]])
  item2 <- (specia[[lpairs[i,"SongB"]]])
  AtoB <- (logspec.dist(item1, item2))
  BtoA <- (logspec.dist(item2, item1))
  lpairs$kldist[i] <- (AtoB + BtoA) / 2
}
lpairs$kldist <- lpairs$kldist / max(lpairs$kldist)
lgwide <- spread(lpairs, key=SongB, value=kldist)
row.names(lgwide) <- paste(exmplars$en[as.integer(row.names(lgwide))], exmplars$use[as.integer(row.names(lgwide))])
lgwide$SongA <- NULL
lgld <- dist(as.matrix(lgwide))
lgdpicks <- bestguesses(lgld)
write.csv(lgdpicks, file="ldpicks.csv", row.names = FALSE)
lghc <- hclust(lgld)          
plot(lghc, cex=0.5, main="LogSpec distance. Sound similarities of example recordings", xlab="Species")
lgwidePCA <- prcomp(lgwide, scale=TRUE)
summary(lgwidePCA) # first 3 components explain 85.67% var of variation
lgpd <- dist(as.matrix(lgwidePCA$x[,1:3]))
lgpdpicks <- bestguesses(lgpd)
write.csv(lgpdpicks, file="lpcadpicks.csv", row.names = FALSE)
lgphc <- hclust(lgpd)
plot(lgphc, cex=0.5, main="LogSpec distance. Sound similarities of example recordings\n based on first 3 Principle Components", xlab="Species")

firstcomponent <- lgwidePCA$x[,1]
secondcomponent <- lgwidePCA$x[,2]
xmin <- min(firstcomponent) - 0.5
xmax <- max(firstcomponent) + 2
ymin <- min(secondcomponent) - 0.5
ymax <- max(secondcomponent) + 0.5
pchval <- rep(15, times=nrow(lgwidePCA$x))
pchval[grep("test", row.names(lgwidePCA$x))] <- 19 
gval <- (lgwidePCA$x[,3]-min(lgwidePCA$x[,3]))/(max(lgwidePCA$x[,3])-min(lgwidePCA$x[,3]))
png("fig7_reftest_lgpcascatter.png", width=6, height=4, units="in", res=300)
par(mar=c(4,4,1,1))
plot(firstcomponent, secondcomponent, pch=pchval, cex=0.8, col=grey(gval),  frame.plot=F, ylim=c(ymin,ymax), xlim=c(xmin,xmax))
pchval <- rep(0, times=nrow(lgwidePCA$x))
pchval[grep("test", row.names(lgwidePCA$x))] <- 1 
points(firstcomponent, secondcomponent, pch=pchval, cex=0.8)
legendtext <- row.names(lgwidePCA$x)
legendtext <- gsub(" test", "", legendtext)
legendtext <- gsub(" reference", "", legendtext)
text(x = firstcomponent, y = secondcomponent, labels= legendtext, cex=0.5,pos=4)
text(x=xmax, y=ymax, labels="point shading is 3rd principle component", cex=0.5,pos=2, font=3)
legend(x=xmax-1, y=ymax-0.5, c("Reference","Test"), pch=c(0,1), cex=0.5)
dev.off()
par(mar=oldpar$mar)
#diffspec
dpairs <- pairs[pairs$SongA %in% referenceandtestset,]
dpairs$kldist <- 99999
for (i in 1:nrow(dpairs)){
  item1 <- (specia[[dpairs[i,"SongA"]]])
  item2 <- (specia[[dpairs[i,"SongB"]]])
  AtoB <- (diffspec(item1, item2))
  BtoA <- (diffspec(item2, item1))
  dpairs$kldist[i] <- (AtoB + BtoA) / 2
}
dpairs$kldist <- dpairs$kldist / max(dpairs$kldist)
dswide <- spread(dpairs, key=SongB, value=kldist)
row.names(dswide) <- paste(exmplars$en[as.integer(row.names(dswide))], exmplars$use[as.integer(row.names(dswide))])
dswide$SongA <- NULL
dsld <- dist(as.matrix(dswide))
dsldpicks <- bestguesses(dsld)
write.csv(dsldpicks, file="dsldpicks.csv", row.names = FALSE)
dshc <- hclust(dsld)          
plot(dshc, cex=0.5, main="DiffSpec distance. Sound similarities of example recordings", xlab="Species")
dswidePCA <- prcomp(dswide, scale=TRUE)
summary(dswidePCA) # first 3 components explain 77.92% of variation
dspd <- dist(as.matrix(dswidePCA$x[,1:3]))
dspcapicks <- bestguesses(dspd)
write.csv(dspcapicks, file="dspcadpicks.csv", row.names = FALSE)
dsphc <- hclust(dspd)
plot(dsphc, cex=0.5, main="DiffSpec distance. Sound similarities of example recordings\n based on first 3 Principle Components", xlab="Species")

firstcomponent <- dswidePCA$x[,1]
secondcomponent <- dswidePCA$x[,2]
xmin <- min(firstcomponent) - 0.5
xmax <- max(firstcomponent) + 2
ymin <- min(secondcomponent) - 0.5
ymax <- max(secondcomponent) + 0.5
pchval <- rep(15, times=nrow(dswidePCA$x))
pchval[grep("test", row.names(dswidePCA$x))] <- 19 
gval <- (dswidePCA$x[,3]-min(dswidePCA$x[,3]))/(max(dswidePCA$x[,3])-min(dswidePCA$x[,3]))
png("fig8_reftest_dspcascatter.png", width=6, height=4, units="in", res=300)
par(mar=c(4,4,1,1))
plot(firstcomponent, secondcomponent, pch=pchval, cex=0.8, col=grey(gval),  frame.plot=F, ylim=c(ymin,ymax), xlim=c(xmin,xmax))
pchval <- rep(0, times=nrow(dswidePCA$x))
pchval[grep("test", row.names(dswidePCA$x))] <- 1 
points(firstcomponent, secondcomponent, pch=pchval, cex=0.8)
legendtext <- row.names(dswidePCA$x)
legendtext <- gsub(" test", "", legendtext)
legendtext <- gsub(" reference", "", legendtext)
text(x = firstcomponent, y = secondcomponent, labels= legendtext, cex=0.5,pos=4)
text(x=xmax, y=ymax, labels="point shading is 3rd principle component", cex=0.5,pos=2, font=3)
legend(x=xmax-1, y=ymax-0.5, c("Reference","Test"), pch=c(0,1), cex=0.5)
dev.off()
par(mar=oldpar$mar)

################my recordings

kpairs <- pairs[pairs$SongA %in% referenceandmyset,]
kpairs$kldist <- 99999
for (i in 1:nrow(kpairs)){
  item1 <- (specia[[kpairs[i,"SongA"]]])
  item2 <- (specia[[kpairs[i,"SongB"]]])
  AtoB <- (kl.dist(item1, item2))$D
  BtoA <- (kl.dist(item2, item1))$D
  kpairs$kldist[i] <- (AtoB + BtoA) / 2
}
kpairs$kldist <- kpairs$kldist / max(kpairs$kldist)
klwide <- spread(kpairs, key=SongB, value=kldist)
row.names(klwide) <- paste(exmplars$en[referenceandmyset], exmplars$use[referenceandmyset])
klwide$SongA <- NULL
kd <- dist(as.matrix(klwide))
kdmypicks <- bestguesses(kd)
write.csv(kdmypicks, file="kdmypicks.csv", row.names = FALSE)
khc <- hclust(kd)          
plot(khc, cex=0.5, main="KL distance. Sound similarities of example recordings", xlab="Species")
klwidePCA <- prcomp(klwide, scale=TRUE)
summary(klwidePCA) # first 3 components explain 90.16% of variation
klpd <- dist(as.matrix(klwidePCA$x[,1:3]))
kpcadmypicks <- bestguesses(klpd)
write.csv(kpcadmypicks, file="kdmypcapicks.csv", row.names = FALSE)
klphc <- hclust(klpd)
plot(klphc, cex=0.5, main="KL distance. Sound similarities of example recordings\n based on first 3 Principle Components", xlab="Species")

firstcomponent <- klwidePCA$x[,1]
secondcomponent <- klwidePCA$x[,2]
xmin <- min(firstcomponent) - 0.5
xmax <- max(firstcomponent) + 2
ymin <- min(secondcomponent) - 0.5
ymax <- max(secondcomponent) + 0.5
pchval <- rep(15, times=nrow(klwide))
pchval[grep("self", row.names(klwidePCA$x))] <- 19 
gval <- (klwidePCA$x[,3]-min(klwidePCA$x[,3]))/(max(klwidePCA$x[,3])-min(klwidePCA$x[,3]))
png("fig9_refmy_klpcascatter.png", width=6, height=4, units="in", res=300)
par(mar=c(4,4,1,1))
plot(firstcomponent, secondcomponent, pch=pchval, cex=0.8, col=grey(gval),  frame.plot=F, ylim=c(ymin,ymax), xlim=c(xmin,xmax))
pchval <- rep(0, times=nrow(klwidePCA$x))
pchval[grep("self", row.names(klwidePCA$x))] <- 1 
pointcol <- rep(1, times=nrow(klwidePCA$x))
pointcol[grep("phone", row.names(klwidePCA$x))] <- 2
pointcol[grep("ipod", row.names(klwidePCA$x))] <- 3
pointcol[grep("voicer", row.names(klwidePCA$x))] <- 4
pointcol[grep("camera", row.names(klwidePCA$x))] <- 5
points(firstcomponent, secondcomponent, pch=pchval, cex=0.8, col=pointcol)
legendtext <- row.names(klwidePCA$x)
legendtext <- gsub(" self", "", legendtext)
legendtext <- gsub(" reference", "", legendtext)
legendtext <- gsub("ipod", "", legendtext)
legendtext <- gsub("phone", "", legendtext)
legendtext <- gsub("camera", "", legendtext)
legendtext <- gsub("voicer", "", legendtext)
text(x = firstcomponent, y = secondcomponent, labels= legendtext, cex=0.5,pos=4)
text(x=xmax, y=ymax, labels="point shading is 3rd principle component", cex=0.5,pos=2, font=3)
legend(x=xmax-1.5, y=ymax-0.5, c("Reference","My phone","My ipod","My voicerecorder","My camera" ), pch=c(0,1,1,1,1), col=c(1,2,3,4,5), cex=0.5)
dev.off()
par(mar=oldpar$mar)

######preparing comparision
names(kdpicks) <- c("KL1","KL2","KL3","KLBird")
kdpicks <- kdpicks[order(kdpicks$KLBird),]
names(lgdpicks) <- c("LG1","LG2","LG3","LGBird")
lgdpicks <- lgdpicks[order(lgdpicks$LGBird),]
names(dsldpicks) <- c("DS1","DS2","DS3","DSBird")
dsldpicks <- dsldpicks[order(dsldpicks$DSBird),]
names(kpcadpicks) <- c("KLp1","KLp2","KLp3","KLpBird")
kpcadpicks <- kpcadpicks[order(kpcadpicks$KLpBird),]
names(lgpdpicks) <- c("LG1","LG2","LG3","LGpBird")
lgpdpicks <- lgpdpicks[order(lgpdpicks$LGpBird),]
names(dspcapicks) <- c("DS1","DS2","DS3","DSpBird")
dspcapicks <- dspcapicks[order(dspcapicks$DSpBird),]

guesses <- cbind(kdpicks,kpcadpicks, lgdpicks, lgpdpicks, dsldpicks, dspcapicks)
write.csv(guesses, file="guesses.csv", row.names = FALSE)

