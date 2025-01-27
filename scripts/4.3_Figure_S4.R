## 4.3 ## Figure S4: Sensitivity analyses

plot.f<-paste0(FigureFolder,"Figure_S4.pdf")
pdf(plot.f,w=9,h=3)

par(mfrow=c(1,3))
par(mar=c(1,2,1,2))
par(oma=c(2,1,2,1))

## Fig 4 A (Beak scan) ####

# Open data Beak Scan:
pc.dat12 <- read.csv(paste0(DataFolder,"data_tmp/BeakScanPC12.csv"))

# Add beak task enrichment:
beakTaskEnrich <- read.csv(paste0(DataFolder,"Data2_BeakTaskEnrichment.csv"))
fdat <- merge(pc.dat12,beakTaskEnrich,by="species")

xcuts<-seq(min(fdat$PC1),max(fdat$PC1)+1.5,length.out=40)
ycuts<-seq(min(fdat$PC2),max(fdat$PC2),length.out=40)

xdiff<-(unique(diff(xcuts),6)/2)[1]
ydiff<-(unique(diff(ycuts),6)/2)[1]

fdat$xCell<-cut(fdat$PC1,xcuts,labels=FALSE,include.lowest = TRUE)
fdat$yCell<-cut(fdat$PC2,ycuts,labels=FALSE,include.lowest = TRUE)
fdat$Cell<-paste(fdat$xCell,fdat$yCell,sep="_")

xcuts<-seq(min(fdat$PC1),max(fdat$PC1)+1.5,length.out=40)
ycuts<-seq(min(fdat$PC2),max(fdat$PC2),length.out=40)

xdiff<-(unique(diff(xcuts),6)/2)[1]
ydiff<-(unique(diff(ycuts),6)/2)[1]

fdat$xCell<-cut(fdat$PC1,xcuts,labels=FALSE,include.lowest = TRUE)
fdat$yCell<-cut(fdat$PC2,ycuts,labels=FALSE,include.lowest = TRUE)
fdat$Cell<-paste(fdat$xCell,fdat$yCell,sep="_")

cellTab<-data.frame(Cell=names(table(fdat$Cell)),btENG=NA,btREA=NA,btCRU=NA)

cellTab$btENG<-as.numeric(sapply(split(fdat$Engulf_task,fdat$Cell),sum))
cellTab$btREA<-as.numeric(sapply(split(fdat$Reach_task,fdat$Cell),sum))
cellTab$btCRU<-as.numeric(sapply(split(fdat$Crush_task,fdat$Cell),sum))

cellTab$x<-as.numeric(unlist(strsplit(cellTab$Cell,"_"))[seq(1,nrow(cellTab)*2,2)])
cellTab$y<-as.numeric(unlist(strsplit(cellTab$Cell,"_"))[seq(2,nrow(cellTab)*2,2)])

cellTab$xMid<-xcuts[cellTab$x]+xdiff
cellTab$yMid<-ycuts[cellTab$y]+ydiff

cellTab$btENGScaled<-cellTab$btENG/sum(fdat$Engulf_task)
cellTab$btREAScaled<-cellTab$btREA/sum(fdat$Reach_task)
cellTab$btCRUScaled<-cellTab$btCRU/sum(fdat$Crush_task)

cellTab$btENGScaled2<-cellTab[,9]/rowSums(cellTab[,9:11])
cellTab$btREAScaled2<-cellTab[,10]/rowSums(cellTab[,9:11])
cellTab$btCRUScaled2<-cellTab[,11]/rowSums(cellTab[,9:11])

ENGCol<-brewer.pal(9, "RdYlBu")[6:9]
REACol<-rev(brewer.pal(9, "BrBG")[1:4])
CRUCol<-brewer.pal(9, "PiYG")[6:9]

cellTab$ENGCol<-ENGCol[cut(cellTab$btENGScaled2,round(seq(0.33,1,length.out=5),2),labels=FALSE,include.lowest = TRUE)]
cellTab$REACol<-REACol[cut(cellTab$btREAScaled2,round(seq(0.33,1,length.out=5),2),labels=FALSE,include.lowest = TRUE)]
cellTab$CRUCol<-CRUCol[cut(cellTab$btCRUScaled2,round(seq(0.33,1,length.out=5),2),labels=FALSE,include.lowest = TRUE)]

cellTabCols <- cellTab[,15:17]
colCol<-apply(cellTab[,9:11],1,which.max)

vals<-rep(NA,nrow(cellTab[,12:14]))
for(i in 1:nrow(cellTab[,12:14])){
  vals[i]<-as.numeric(cellTab[i,12:14][colCol[i]])
}

cellTabCols<-cellTab[,15:17]
colCol<-apply(cellTab[,9:11],1,which.max)

plot(0,0,type="n",xlim=c(min(xcuts)-0.5,max(xcuts)),ylim=c(min(ycuts),max(ycuts)),ylab="",xlab="",yaxt="n",xaxt="n")
axis(1,at=seq(-4,3,2),labels=seq(-4,3,2),cex.axis=0.75,padj=-2.5,tck=-0.0125,lwd=0.75)
axis(2,at=seq(-4,3.5,2),labels=c("-4","-2"," 0"," 2"),cex.axis=0.75,hadj=-0.3,las=2,tck=-0.0125,lwd=0.75)
mtext("PC1",side=1,outer=F,line=0,cex=0.75,at=max(xcuts))
mtext("PC2",side=2,outer=F,line=0,cex=0.75,at=max(ycuts))
mtext("Beak shape (3D Scan)",side=3,outer=F,line=0,cex=0.75)

for(i in 1:nrow(cellTab)){
  polygon(c(cellTab$xMid[i]-xdiff,cellTab$xMid[i]+xdiff,cellTab$xMid[i]+xdiff,cellTab$xMid[i]-xdiff),
          c(cellTab$yMid[i]-ydiff,cellTab$yMid[i]-ydiff,cellTab$yMid[i]+ydiff,cellTab$yMid[i]+ydiff),
          col=cellTabCols[i,colCol[i]],border="white",lwd=0.25)
}

xrange<-range(xcuts)
yrange<-range(ycuts)

xOr<-min(xcuts)+(0.8*diff(xrange))
yOr<-min(ycuts)+(0.15*diff(yrange))

xdisp<-diff(xrange)/10
ydisp<-diff(yrange)/10

xE2<-xOr+xdisp
yE2<-yOr-ydisp

xE3<-xOr-xdisp
yE3<-yOr-ydisp

xE1<-xOr
yE1<-yOr+sqrt(xdisp^2+ydisp^2)

lines(c(xOr,xE1),c(yOr,yE1))
lines(c(xOr,xE2),c(yOr,yE2))
lines(c(xOr,xE3),c(yOr,yE3))

x1<-seq(xOr,xE1,length.out=5)
y1<-seq(yOr,yE1,length.out=5)

x2<-seq(xOr,xE2,length.out=5)
y2<-seq(yOr,yE2,length.out=5)

x3<-seq(xOr,xE3,length.out=5)
y3<-seq(yOr,yE3,length.out=5)

for(i in 2:5){
  points(x1[i],y1[i],pch=21,bg=CRUCol[i-1],cex=1,lwd=0.5)
  points(x2[i],y2[i],pch=21,bg=REACol[i-1],cex=1,lwd=0.5)
  points(x3[i],y3[i],pch=21,bg=ENGCol[i-1],cex=1,lwd=0.5)
}
text(x1[2:5]-0.15,y1[2:5]-0.1,c("50","66","83","100%"),cex=0.4,pos=4)
text(x1[5],y1[5]-0.05,"Crush",cex=0.7,pos=3)
text(x2[5]+0.3,y2[5]+0.05,"Reach",cex=0.7,pos=1)
text(x3[5]-0.3,y3[5]+0.05,"Engulf",cex=0.7,pos=1)

# Plot archetypes points:
load(paste0(OutputFolder,"out4_ArchetypesObsScan_RobustAlg_PC12_Beaks.rda"))
a3Rxy <- archCoords[[3]]
points(a3Rxy,pch=16,cex=1.5)

# Add panel letter:
xr <- c(xrange[2]-xrange[1])
yr <- c(yrange[2]-yrange[1])

# Letter outside corner:
text(min(par("usr")[1])+(xr/30),max(par("usr")[4])-(yr/12), "A", pos = 4, cex = 1.5, font = 2)

## Fig 4 B (PPCA Beak) ####

# Open PPCA beak data:
load(file=paste0(DataFolder,"data_tmp/BeakPPCA.Rdata"))
# Open beak task enrichment data:
beakTaskEnrich <- read.csv(paste0(DataFolder,"Data2_BeakTaskEnrichment.csv"))
fdat <- merge(ppca.spp,beakTaskEnrich,by="species")

xcuts<-seq(min(fdat$PC2),max(fdat$PC2)+1.5,length.out=40)
ycuts<-seq(min(fdat$PC3),max(fdat$PC3),length.out=40)

xdiff<-(unique(diff(xcuts),6)/2)[1]
ydiff<-(unique(diff(ycuts),6)/2)[1]

fdat$xCell<-cut(fdat$PC2,xcuts,labels=FALSE,include.lowest = TRUE)
fdat$yCell<-cut(fdat$PC3,ycuts,labels=FALSE,include.lowest = TRUE)
fdat$Cell<-paste(fdat$xCell,fdat$yCell,sep="_")

xcuts<-seq(min(fdat$PC2),max(fdat$PC2)+1.5,length.out=40)
ycuts<-seq(min(fdat$PC3),max(fdat$PC3),length.out=40)

xdiff<-(unique(diff(xcuts),6)/2)[1]
ydiff<-(unique(diff(ycuts),6)/2)[1]

fdat$xCell<-cut(fdat$PC2,xcuts,labels=FALSE,include.lowest = TRUE)
fdat$yCell<-cut(fdat$PC3,ycuts,labels=FALSE,include.lowest = TRUE)
fdat$Cell<-paste(fdat$xCell,fdat$yCell,sep="_")

cellTab<-data.frame(Cell=names(table(fdat$Cell)),btENG=NA,btREA=NA,btCRU=NA)

cellTab$btENG<-as.numeric(sapply(split(fdat$Engulf_task,fdat$Cell),sum))
cellTab$btREA<-as.numeric(sapply(split(fdat$Reach_task,fdat$Cell),sum))
cellTab$btCRU<-as.numeric(sapply(split(fdat$Crush_task,fdat$Cell),sum))

cellTab$x<-as.numeric(unlist(strsplit(cellTab$Cell,"_"))[seq(1,nrow(cellTab)*2,2)])
cellTab$y<-as.numeric(unlist(strsplit(cellTab$Cell,"_"))[seq(2,nrow(cellTab)*2,2)])

cellTab$xMid<-xcuts[cellTab$x]+xdiff
cellTab$yMid<-ycuts[cellTab$y]+ydiff

cellTab$btENGScaled<-cellTab$btENG/sum(fdat$Engulf_task)
cellTab$btREAScaled<-cellTab$btREA/sum(fdat$Reach_task)
cellTab$btCRUScaled<-cellTab$btCRU/sum(fdat$Crush_task)

cellTab$btENGScaled2<-cellTab[,9]/rowSums(cellTab[,9:11])
cellTab$btREAScaled2<-cellTab[,10]/rowSums(cellTab[,9:11])
cellTab$btCRUScaled2<-cellTab[,11]/rowSums(cellTab[,9:11])

ENGCol<-brewer.pal(9, "RdYlBu")[6:9]
REACol<-rev(brewer.pal(9, "BrBG")[1:4])
CRUCol<-brewer.pal(9, "PiYG")[6:9]

cellTab$ENGCol<-ENGCol[cut(cellTab$btENGScaled2,round(seq(0.33,1,length.out=5),2),labels=FALSE,include.lowest = TRUE)]
cellTab$REACol<-REACol[cut(cellTab$btREAScaled2,round(seq(0.33,1,length.out=5),2),labels=FALSE,include.lowest = TRUE)]
cellTab$CRUCol<-CRUCol[cut(cellTab$btCRUScaled2,round(seq(0.33,1,length.out=5),2),labels=FALSE,include.lowest = TRUE)]

cellTabCols<-cellTab[,15:17]
colCol<-apply(cellTab[,9:11],1,which.max)

vals<-rep(NA,nrow(cellTab[,12:14]))
for(i in 1:nrow(cellTab[,12:14])){
  vals[i]<-as.numeric(cellTab[i,12:14][colCol[i]])
}

cellTabCols<-cellTab[,15:17]
colCol<-apply(cellTab[,9:11],1,which.max)

plot(0,0,type="n",xlim=c(min(xcuts)-0.5,max(xcuts)),ylim=c(min(ycuts),max(ycuts)),ylab="",xlab="",yaxt="n",xaxt="n")
axis(1,at=seq(-4,4,2),labels=seq(-4,4,2),cex.axis=0.75,padj=-2.5,tck=-0.0125,lwd=0.75)
axis(2,at=seq(-4,4,2),labels=seq(-4,4,2),cex.axis=0.75,hadj=-0.3,las=2,tck=-0.0125,lwd=0.75)
mtext("PC2",side=1,outer=F,line=0,cex=0.75,at=max(xcuts))
mtext("PC3",side=2,outer=F,line=0,cex=0.75,at=max(ycuts))
mtext("Beak shape",side=3,outer=F,line=0,cex=0.75)

for(i in 1:nrow(cellTab)){
  polygon(c(cellTab$xMid[i]-xdiff,cellTab$xMid[i]+xdiff,cellTab$xMid[i]+xdiff,cellTab$xMid[i]-xdiff),
          c(cellTab$yMid[i]-ydiff,cellTab$yMid[i]-ydiff,cellTab$yMid[i]+ydiff,cellTab$yMid[i]+ydiff),
          col=cellTabCols[i,colCol[i]],border="white",lwd=0.25)
}

xrange<-range(xcuts)
yrange<-range(ycuts)

xOr<-min(xcuts)+(0.8*diff(xrange))
yOr<-min(ycuts)+(0.15*diff(yrange))

xdisp<-diff(xrange)/10
ydisp<-diff(yrange)/10

xE2<-xOr+xdisp
yE2<-yOr-ydisp

xE3<-xOr-xdisp
yE3<-yOr-ydisp

xE1<-xOr
yE1<-yOr+sqrt(xdisp^2+ydisp^2)

lines(c(xOr,xE1),c(yOr,yE1))
lines(c(xOr,xE2),c(yOr,yE2))
lines(c(xOr,xE3),c(yOr,yE3))

x1<-seq(xOr,xE1,length.out=5)
y1<-seq(yOr,yE1,length.out=5)

x2<-seq(xOr,xE2,length.out=5)
y2<-seq(yOr,yE2,length.out=5)

x3<-seq(xOr,xE3,length.out=5)
y3<-seq(yOr,yE3,length.out=5)

# Plot legend: 

for(i in 2:5){
  points(x1[i],y1[i],pch=21,bg=CRUCol[i-1],cex=1,lwd=0.5)
  points(x2[i],y2[i],pch=21,bg=REACol[i-1],cex=1,lwd=0.5)
  points(x3[i],y3[i],pch=21,bg=ENGCol[i-1],cex=1,lwd=0.5)
}
text(x1[2:5]-0.15,y1[2:5]-0.1,c("50","66","83","100%"),cex=0.4,pos=4)
text(x1[5],y1[5]-0.02,"Crush",cex=0.7,pos=3)
text(x2[5]+0.3,y2[5],"Reach",cex=0.7,pos=1)
text(x3[5]-0.3,y3[5],"Engulf",cex=0.7,pos=1)

# Plot archetype points:
load(paste0(OutputFolder,"out4_ArchetypesObs_RobustAlg_PPC23_Beaks.rda"))
a3Rxy <- archCoords[[3]]
points(a3Rxy,pch=16,cex=1.5)

# Add panel letter:
xr <- c(xrange[2]-xrange[1])
yr <- c(yrange[2]-yrange[1])

# Letter outside corner:
#mtext("B", side = 3, line = 0.1, at = par("usr")[1]-xr/8, adj = 0, cex = 1.25, font=2)
# Letter inside corner:
text(min(par("usr")[1])+(xr/30),max(par("usr")[4])-(yr/12), "B", pos = 4, cex = 1.5, font = 2)

## Body enrichment plot with PPCA ####

# Open PPCA beak data:
load(file=paste0(DataFolder,"data_tmp/BodyPPCA.Rdata"))
# Open beak task enrichment data:
ForTaskEnrich <- read.csv(paste0(DataFolder,"Data3_BodyTaskEnrichment.csv"))
fdat <- merge(ppca.spp,ForTaskEnrich,by="species")
names(fdat)

xcuts<-seq(min(fdat$PC2),max(fdat$PC2)+1.5,length.out=40)
ycuts<-seq(min(fdat$PC3),max(fdat$PC3),length.out=40)

xdiff<-(unique(diff(xcuts),6)/2)[1]
ydiff<-(unique(diff(ycuts),6)/2)[1]

fdat$xCell<-cut(fdat$PC2,xcuts,labels=FALSE,include.lowest = TRUE)
fdat$yCell<-cut(fdat$PC3,ycuts,labels=FALSE,include.lowest = TRUE)
fdat$Cell<-paste(fdat$xCell,fdat$yCell,sep="_")

xcuts<-seq(min(fdat$PC2),max(fdat$PC2)+1.5,length.out=40)
ycuts<-seq(min(fdat$PC3),max(fdat$PC3),length.out=40)

xdiff<-(unique(diff(xcuts),6)/2)[1]
ydiff<-(unique(diff(ycuts),6)/2)[1]

fdat$xCell<-cut(fdat$PC2,xcuts,labels=FALSE,include.lowest = TRUE)
fdat$yCell<-cut(fdat$PC3,ycuts,labels=FALSE,include.lowest = TRUE)
fdat$Cell<-paste(fdat$xCell,fdat$yCell,sep="_")

cellTab<-data.frame(Cell=names(table(fdat$Cell)),ftFLY=NA,ftSWM=NA,ftWLK=NA)

cellTab$ftFLY<-as.numeric(sapply(split(fdat$Fly_task,fdat$Cell),sum))
cellTab$ftSWM<-as.numeric(sapply(split(fdat$Swim_task,fdat$Cell),sum))
cellTab$ftWLK<-as.numeric(sapply(split(fdat$Walk_task,fdat$Cell),sum))

cellTab$x<-as.numeric(unlist(strsplit(cellTab$Cell,"_"))[seq(1,nrow(cellTab)*2,2)])
cellTab$y<-as.numeric(unlist(strsplit(cellTab$Cell,"_"))[seq(2,nrow(cellTab)*2,2)])

cellTab$xMid<-xcuts[cellTab$x]+xdiff
cellTab$yMid<-ycuts[cellTab$y]+ydiff

cellTab$ftFLYScaled<-cellTab$ftFLY/sum(fdat$Fly_task)
cellTab$ftSWMScaled<-cellTab$ftSWM/sum(fdat$Swim_task)
cellTab$ftWLKScaled<-cellTab$ftWLK/sum(fdat$Walk_task)

cellTab$ftFLYScaled2<-cellTab[,9]/rowSums(cellTab[,9:11])
cellTab$ftSWMScaled2<-cellTab[,10]/rowSums(cellTab[,9:11])
cellTab$ftWLKScaled2<-cellTab[,11]/rowSums(cellTab[,9:11])

SWMCol<-brewer.pal(9, "RdYlBu")[6:9]
FLYCol<-brewer.pal(9, "PiYG")[6:9]
WLKCol<-rev(brewer.pal(9, "BrBG")[1:4])

cellTab$FLYCol<-FLYCol[cut(cellTab$ftFLYScaled2,round(seq(0.33,1,length.out=5),2),labels=FALSE,include.lowest = TRUE)]
cellTab$SWMCol<-SWMCol[cut(cellTab$ftSWMScaled2,round(seq(0.33,1,length.out=5),2),labels=FALSE,include.lowest = TRUE)]
cellTab$WLKCol<-WLKCol[cut(cellTab$ftWLKScaled2,round(seq(0.33,1,length.out=5),2),labels=FALSE,include.lowest = TRUE)]

cellTabCols<-cellTab[,15:17]
colCol<-apply(cellTab[,9:11],1,which.max)

vals<-rep(NA,nrow(cellTab[,12:14]))
for(i in 1:nrow(cellTab[,12:14])){
  vals[i]<-as.numeric(cellTab[i,12:14][colCol[i]])
}

cellTabCols<-cellTab[,15:17]
colCol<-apply(cellTab[,9:11],1,which.max)

plot(0,0,type="n",xlim=c(min(xcuts)-0.5,max(xcuts)),ylim=c(min(ycuts),max(ycuts)),ylab="",xlab="",yaxt="n",xaxt="n")
axis(1,at=seq(-4,4,2),labels=seq(-4,4,2),cex.axis=0.75,padj=-2.5,tck=-0.0125,lwd=0.75)
axis(2,at=seq(-4,4,2),labels=c("-4","-2"," 0"," 2"," 4"),cex.axis=0.75,hadj=-0.3,las=2,tck=-0.0125,lwd=0.75)
mtext("PC2",side=1,outer=F,line=0,cex=0.75,at=max(xcuts))
mtext("PC3",side=2,outer=F,line=0,cex=0.75,at=max(ycuts))
mtext("Body shape",side=3,outer=F,line=0,cex=0.75)

for(i in 1:nrow(cellTab)){
  polygon(c(cellTab$xMid[i]-xdiff,cellTab$xMid[i]+xdiff,cellTab$xMid[i]+xdiff,cellTab$xMid[i]-xdiff),
          c(cellTab$yMid[i]-ydiff,cellTab$yMid[i]-ydiff,cellTab$yMid[i]+ydiff,cellTab$yMid[i]+ydiff),
          col=cellTabCols[i,colCol[i]],border="white",lwd=0.25)
}

xrange<-range(xcuts)
yrange<-range(ycuts)

xOr<-min(xcuts)+(0.8*diff(xrange))
yOr<-min(ycuts)+(0.15*diff(yrange))

xdisp<-diff(xrange)/10
ydisp<-diff(yrange)/10

xE2<-xOr+xdisp
yE2<-yOr-ydisp

xE3<-xOr-xdisp
yE3<-yOr-ydisp

xE1<-xOr
yE1<-yOr+sqrt(xdisp^2+ydisp^2)

lines(c(xOr,xE1),c(yOr,yE1))
lines(c(xOr,xE2),c(yOr,yE2))
lines(c(xOr,xE3),c(yOr,yE3))

x1<-seq(xOr,xE1,length.out=5)
y1<-seq(yOr,yE1,length.out=5)

x2<-seq(xOr,xE2,length.out=5)
y2<-seq(yOr,yE2,length.out=5)

x3<-seq(xOr,xE3,length.out=5)
y3<-seq(yOr,yE3,length.out=5)

# Plot legend:

for(i in 2:5){
  points(x1[i],y1[i],pch=21,bg=SWMCol[i-1],cex=1,lwd=0.5)
  points(x2[i],y2[i],pch=21,bg=WLKCol[i-1],cex=1,lwd=0.5)
  points(x3[i],y3[i],pch=21,bg=FLYCol[i-1],cex=1,lwd=0.5)
}
text(x1[2:5]-0.15,y1[2:5]-0.1,c("50","66","83","100%"),cex=0.4,pos=4)
text(x1[5],y1[5]-0.05,"Swim",cex=0.7,pos=3)
text(x2[5]+0.3,y2[5]+0.05,"Walk",cex=0.7,pos=1)
text(x3[5]-0.3,y3[5]+0.05,"Fly",cex=0.7,pos=1)

# Plot archetypes points:
load(paste0(OutputFolder,"out4_ArchetypesObs_RobustAlg_PPC23_Body.rda"))
a3Rxy <- archCoords[[3]]
points(a3Rxy,pch=16,cex=1.5)

# Add panel letter:
xr <- c(xrange[2]-xrange[1])
yr <- c(yrange[2]-yrange[1])

# Letter inside corner:
text(min(par("usr")[1])+(xr/30),max(par("usr")[4])-(yr/12), "C", pos = 4, cex = 1.5, font = 2)

dev.off()

## End script 4.3 ##