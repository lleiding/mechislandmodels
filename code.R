## Code for producing figures and calculating MDS
## for Review "Biodiversity dynamics on islands: explicitly accounting for causality in mechanistic models"
## by Ludwig Leidinger and Juliano Sarmento Cabral
##
## ludwig.leidinger@uni-wuerzburg.de

##load packages:
library(vegan)
## Note: for package vegan 'tk' and 'gcc-fortran' must be installed on the OS
library(smacof)

## read in data:
mtraw=read.csv("master_table_newnames.csv", header=T)

mt=mtraw
## make factors consistent and if possible binary
## always check 'levels(mt$COLNAME)' (before and after each re-assignment)!
levels(mt$stochastic.deterministic)

colnames(mt)[11]<-"stochastic"
levels(mt$spatially.explicit)

levels(mt$neutral)

levels(mt$evolution)

levels(mt$evolution)[c(1,3)]<-"yes" #speciation, neutral mutation, etc.
mt$evolution<-factor(mt$evolution, levels=c("no","yes"))
levels(mt$agent.level)

levels(mt$static.env.)
names(mt)[names(mt)=="static.env."]<-"static.environment"

names(mt)[names(mt)=="niche.axis"]<-"no.niche.axes"

#######################################################################
## MDS:
mtred=mt[c("study","stochastic","spatially.explicit","neutral",
           "no.niche.axes","evolution","agent.level","static.environment","disturbance")]
mtred[,-1]<-lapply(mtred[,-1],as.numeric) 

distPapers<-vegdist(mtred[,-1], method='bray',na.rm=T)

Papersmds=smacofSym(distPapers)

## Backup:
Papersmds.points.orig=Papersmds$conf

## Jitter:
set.seed(5) #5
Papersmds.points<-Papersmds.points.orig
for(i in 1:nrow(Papersmds.points)){
    Papersmds.points[i,2]<-Papersmds.points[i,2]+((runif(1)-0.5)/1.6) #1.6
}
pdf(file="pcaneu.pdf", bg="white", width=7, height=6)
par(cex=0.7, mar=c(4.2,4.2,0.8,0.2))
plot(Papersmds.points, type="n",xlim=c(-1.1,1.7), xlab="Dimension 1", ylab="Dimension 2", cex.lab=1.5, cex.axis=1.5, font=2, font.lab=2)

ordihull(Papersmds.points, groups=mtred$no.niche.axes>0,draw="polygon",col="lightgreen",show.groups=T, cex=0.5)
ordihull(Papersmds.points,groups=mtred$evolution-1,draw="polygon",col="orange",show.groups=T, cex=0.5)

text(Papersmds.points,labels=mtred[,1], font=1, cex=0.9)

fit<- envfit(Papersmds.points,mtred[,c(2:9)], perm = 999, na.rm=T)

var.names <- rownames(scores(fit, "vectors"))
plot(fit, p.max = 0.1,col = "blue", font=3,cex = 1.2) # plotting axes with p < 0.1
dev.off()


#######################################################################
## MDS for closeness of model characteristics:
library(MASS)
mtred[,c(2,3,4,6,8,9)]<-mtred[,c(2,3,4,6,8,9)]-1 # make categories zero-based

mtredmds<-mtred[,c(2,3,4,6,8)] # reduce data to binary categories
mtredmds$niche.based<-ifelse(mtred$no.niche.axes>0,1,0) # create binary categories...
mtredmds$ind.based<-ifelse(mtred$agent.level==1,1,0) # ... for rest of data
mtredmds$pop.based<-ifelse(mtred$agent.level==2,1,0)
mtredmds$sp.based<-ifelse(mtred$agent.level==3,1,0)
mdsneu<-metaMDS(mtredmds, try=100) # calculate metaMDS for getting species coordinates

## Create Plot:
pdf(file='species.pdf',bg='white',width=4,height=4)
par(cex=0.7)
plot(mdsneu, xlab="Dimension 1", ylab="Dimension 2",display=c('species'), type='n', color='black', font=2, font.lab=2 )
text(mdsneu$species, labels=rownames(mdsneu$species), font=1, cex=0.9)
dev.off()
#######################################################################
## Publication year distribution:

## requires working directory of PCA data
years = read.table("studyyears.txt", sep="\t", header=T)

found=as.data.frame(table(mt$year))
found[,1]<-as.numeric(as.character(found[,1]))


pdf(file="../manuscript/yearplot.pdf", bg="white", width=8, height=6)
par(ps=15,mfrow=c(2,1),
          oma = c(4,4,0,0) + 0.1,
          mar = c(1,0,1,1) + 0.1)
plot(years$records ~ years$Year, xlim=c(1980,2017),
     ylim=c(0,70),type="l", las=1,
     xlab="Year", ylab="Number of published papers", bty="l", col="black", lwd=2,
     font=2, font.lab=2, xaxs="i", yaxs="i",
     xaxt = 'n')
axis(side = 1, labels=F)
mtext("A", side=3, font=2, adj=0.02, line=-1)
plot(found$Var1, found$Freq, xlim=c(1980,2017),
     ylim=c(0,5), type="l", las=1,
     xlab="Year", ylab="Number of published papers", bty="l", col="black", lwd=2,
     font=2, font.lab=2, xaxs="i", yaxs="i")
mtext("B", side=3, font=2, adj=0.02, line=-1)
title(xlab = "Year",
      ylab = "Number of published papers",
      outer = TRUE, line = 3, font.lab=2)
dev.off()

