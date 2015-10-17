rm(list=ls()) #clearing variables
library(tuneR) #reading in sounds
library(seewave) #analysing sounds
library(tidyr) #reorganising data
library(ggplot2) #making fancy graphs

oldpar <- par() #stores current graph settings

bestguesses <- function(nameddistancematrix){
inorder <- as.data.frame(apply(as.matrix(nameddistancematrix), 1, function(x){(names(x)[order(x)])}))
#rearranges data in order of closeness; does this for each line
topref <- data.frame(t(apply(inorder, 2, function(x){(x[grep("reference",x)])[1:3]})))
#only for comparisons with reference samples rather than phone recordings
topref <- topref[grep("test|self",row.names(topref)),]
#get closests 3 species names 
topref$bird <- row.names(topref)
#stores changes in data + returns
return(topref)
}   #got help with this

# NB not all graphs made got used

###figure1
tui <- readMP3("30 Tui-Song-50.Mp3.mp3") #read in tui file
penguin <- readMP3("35 Yellow-Eyed-Penguin.Mp3.mp3") #read in penguin file
png(filename = "figure1a.png",width = 16, height = 8, units = "cm", res=150) #make image file
spectro(tui, flim=c(0.5,6), main="a) Tui recording") #put spec in image file
dev.off()
png(filename = "figure1b.png",width = 16, height = 8, units = "cm", res=150)
spectro(penguin, flim=c(0.5,6), main="b) Penguin recording")
dev.off() #stop making image
######

usethese <- read.csv("~/nbu/birding/Sep/usethese.csv") #information about recordings
exmplars <- usethese
# prepFiles <- list.files(path="xenoprep", pattern=".wav$")
# exmplars$wavs <- paste("~/nbu/birding/Sep/xenoprep/",exmplars$id, "_mod.wav", sep="")
# listOfSoundsE <- lapply(exmplars$wavs, readWave)
# specia <- lapply(listOfSoundsE, function(x){meanspec(x, f=f, plot = FALSE, identify=TRUE)})
#save(specia, file="spectra.RData")
### having read in .wavfiles, was used as r.data file which was faster
load("spectra.RData") #contains meanspec of all recordings
referenceset <- sort(which(usethese$use == "reference"))
referenceandtestset <- sort(c(which(usethese$use == "reference"),which(usethese$use == "test")))
referenceandmyset <- sort(c(which(usethese$use == "reference"),which(usethese$use == "self"))) #splitting data into sets
SongA <- rep(1:nrow(usethese), each= length(referenceset)) #all songs
SongB <- rep(referenceset, length.out=length(SongA)) #comparison set
pairs <- data.frame(SongA,SongB) #every possible pair of A and B

###figure 2
## silvereye is spectrum 25, fantail is spectrum 15, falcon is spectrum 5
silv <- data.frame(specia[[25]])
fant <- data.frame(specia[[15]])
falc <- data.frame(specia[[5]])
silv$species <- "silvereye"
fant$species <- "fantail"
falc$species <- "falcon"
mnspecs <- rbind(silv, fant, falc)
png(filename = "figure2_meanspecscomp.png",width = 16, height = 8, units = "cm", res=150)
ggplot(data=mnspecs, aes(x=x, y=y, group=species, colour=species)) +
  geom_line() + xlab("frequencies") + ylab("Normalised Amplitude") + 
  ggtitle("Comparison of mean spectrums from 3 different species") +     # Set title
  theme_bw() #making line graph of three bird species
dev.off()
####

########## get kl.distance closeness of only references
kpairs <- pairs[pairs$SongA %in% referenceset,] #only the pairs where song A is in the reference group
kpairs$kldist <- 99999 #needs something there
for (i in 1:nrow(kpairs)){
  item1 <- (specia[[kpairs[i,"SongA"]]]) #get spec of whatever song A is in that row
  item2 <- (specia[[kpairs[i,"SongB"]]]) #same for song B
  AtoB <- (kl.dist(item1, item2))$D #find KLdist between
  BtoA <- (kl.dist(item2, item1))$D #same, might be different
  kpairs$kldist[i] <- (AtoB + BtoA) / 2 #average them
}
kpairs$kldist <- kpairs$kldist / max(kpairs$kldist) #normalise against longest distance
klwide <- spread(kpairs, key=SongB, value=kldist) #changes table from long to wide
row.names(klwide) <- paste(exmplars$en[referenceset], exmplars$use[referenceset]) #fixes row names/ means row name is species
klwide$SongA <- NULL #gets rid of serial numbers of species recording as has name
kd <- dist(as.matrix(klwide)) #builds a distance matrix of data, ie. how distant everything is as one number
khc <- hclust(kd) 
png("fig3_refsamples_only_raw.png", width=6, height=4, units="in", res=300)
plot(khc, cex=0.5, xlab="Species")
dev.off() #builds cluster dendrogram from data
klwidePCA <- prcomp(klwide, scale=TRUE) #making PCA
summary(klwidePCA) # first 3 components explain 89.87% of variation
klpd <- dist(as.matrix(klwidePCA$x[,1:3]))
klphc <- hclust(klpd)
png("fig4_refsamples_only_pca.png", width=6, height=4, units="in", res=300)
plot(klphc, cex=0.5, xlab="Species")
dev.off() #makes cluster dendrogram of PCA
firstcomponent <- klwidePCA$x[,1]
secondcomponent <- klwidePCA$x[,2]
xmin <- min(firstcomponent) - 0.5
xmax <- max(firstcomponent) + 2 #fixing graph axes
ymin <- min(secondcomponent) - 0.5
ymax <- max(secondcomponent) + 0.5
pchval <- rep(15, times=nrow(klwide)) #sets default value for shape
gval <- (klwidePCA$x[,3]-min(klwidePCA$x[,3]))/(max(klwidePCA$x[,3])-min(klwidePCA$x[,3])) #amount of grey = shading based on PCA 3rd value
png("fig5_refsamples_only_pcascatter.png", width=6, height=4, units="in", res=300)
plot(firstcomponent, secondcomponent, pch=pchval, cex=0.8, col=grey(gval),  frame.plot=F, ylim=c(ymin,ymax), xlim=c(xmin,xmax))
par(mar=c(4,4,1,1))
pchval <- rep(0, times=nrow(klwide))
points(firstcomponent, secondcomponent, pch=pchval, cex=0.8) #plot values, then plots shape border
text(x = firstcomponent, y = secondcomponent, labels= exmplars$en[referenceset], cex=0.5,pos=4)
text(x=xmax, y=ymax, labels="point shading is 3rd principle component", cex=0.5,pos=2, font=3) #auxiliary text
legend(x=xmax-1, y=ymax-0.5, "Species", pch=0, cex=0.5) #graph key
dev.off() #making scatter graph
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
} #same again
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

#same as ref except uses different group at the start

firstcomponent <- klwidePCA$x[,1]
secondcomponent <- klwidePCA$x[,2]
xmin <- min(firstcomponent) - 0.5
xmax <- max(firstcomponent) + 2
ymin <- min(secondcomponent) - 0.5
ymax <- max(secondcomponent) + 0.5
pchval <- rep(15, times=nrow(klwide))
pchval[grep("test", row.names(klwidePCA$x))] <- 19 
gval <- (klwidePCA$x[,3]-min(klwidePCA$x[,3]))/(max(klwidePCA$x[,3])-min(klwidePCA$x[,3]))
png("fig7a_reftest_klpcascatter.png", width=6, height=4, units="in", res=300)
par(mar=c(4,4,1,1))
plot(firstcomponent, secondcomponent, pch=pchval, cex=0.8, col=grey(gval),  frame.plot=F, ylim=c(ymin,ymax), xlim=c(xmin,xmax), main="7a) KL distance")
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
png("fig7b_reftest_lgpcascatter.png", width=6, height=4, units="in", res=300)
par(mar=c(4,4,1,1))
plot(firstcomponent, secondcomponent, pch=pchval, cex=0.8, col=grey(gval),  frame.plot=F, ylim=c(ymin,ymax), xlim=c(xmin,xmax), , main="7b) Logspec")
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
par(mar=oldpar$mar) #using logspec
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
png("fig7c_reftest_dspcascatter.png", width=6, height=4, units="in", res=300)
par(mar=c(4,4,1,1))
plot(firstcomponent, secondcomponent, pch=pchval, cex=0.8, col=grey(gval),  frame.plot=F, ylim=c(ymin,ymax), xlim=c(xmin,xmax), main="7c) Diffspec")
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
par(mar=oldpar$mar) #same but using diffspec

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

#only used kl distance, saved results as spreadsheet; used my recordsing vs references

####################
#simlation estimating how likely 3 picks of 2 species and 1 of 3 (or better) is

sim <-function(){
  testset <- 1:14
  test1 <- sample(testset, 4, replace=TRUE)
  test2 <- sample(testset, 4, replace=TRUE)
  test3 <- sample(testset, 4, replace=TRUE)
  l1 <- length(tapply(test1,test1,length))
  l2 <- length(tapply(test2,test2,length))
  l3 <- length(tapply(test3,test3,length))
  lorder <- sort(c(l1,l2,l3))
  if (lorder[1] <= 2 & lorder[2] <= 2 & lorder[3] <= 3){
    return(1)
  } else {
    return(0)
  }
}

isgood <- replicate(100000,sim())

print(sum(isgood)/length(isgood))
#####

# brute force test
