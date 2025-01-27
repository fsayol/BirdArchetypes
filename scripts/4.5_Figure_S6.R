## 4.5 ## Figure S6: Sensitivity of task enrichment scoring

outFolder <- paste0(DataFolder,'Data4_Task_sensitivity/')

## Load data for sensitivity plots ####

pc.dat23beak <- read.csv(paste0(DataFolder,"data_tmp/BeakPC23.csv"))
pc.dat23body <- read.csv(paste0(DataFolder,"data_tmp/BodyPC23.csv"))

# Plot beak and body enrichment sensitivity plots (Max) ####

plot.f <- paste0(FigureFolder,"Figure_S6.pdf")
pdf(plot.f,w=5*1.5,h=5)

# Set up multipanel plotting, e.g., 2 rows and 3 columns for 6 plots
par(mfrow=c(2, 3))  # Adjust based on the number of plots
par(oma=c(2,2,2,2))  # Outer margins: can be adjusted based on your needs
par(mar=c(2,2,2,2))  # Inner margins between individual plots

### Beak plots (A,B,C)

load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Beaks.rda"))
a3Rxy<-archCoords[[3]]

arch.list <- c("CRU","ENG","REA")
arch.names <- c("Crush","Engulf","Reach")

for (arch.i in 1:length(arch.list)){
  
  # Add beak task enrichment
  beakTaskEnrich <- read.csv(paste0(outFolder,"Data4_BeakTaskEnrichment_Sensi_Max",arch.list[arch.i],".csv"))
  names(beakTaskEnrich)
  fdat.i <- merge(pc.dat23beak,beakTaskEnrich,by="species")
  
  xcuts<-seq(min(fdat.i$PC2_beak),max(fdat.i$PC2_beak)+1.5,length.out=40)
  ycuts<-seq(min(fdat.i$PC3_beak),max(fdat.i$PC3_beak),length.out=40)
  
  xdiff<-(unique(diff(xcuts),6)/2)[1]
  ydiff<-(unique(diff(ycuts),6)/2)[1]
  
  fdat.i$xCell<-cut(fdat.i$PC2_beak,xcuts,labels=FALSE,include.lowest = TRUE)
  fdat.i$yCell<-cut(fdat.i$PC3_beak,ycuts,labels=FALSE,include.lowest = TRUE)
  fdat.i$Cell<-paste(fdat.i$xCell,fdat.i$yCell,sep="_")
  
  xcuts<-seq(min(fdat.i$PC2_beak),max(fdat.i$PC2_beak)+1.5,length.out=40)
  ycuts<-seq(min(fdat.i$PC3_beak),max(fdat.i$PC3_beak),length.out=40)
  
  xdiff<-(unique(diff(xcuts),6)/2)[1]
  ydiff<-(unique(diff(ycuts),6)/2)[1]
  
  fdat.i$xCell<-cut(fdat.i$PC2_beak,xcuts,labels=FALSE,include.lowest = TRUE)
  fdat.i$yCell<-cut(fdat.i$PC3_beak,ycuts,labels=FALSE,include.lowest = TRUE)
  fdat.i$Cell<-paste(fdat.i$xCell,fdat.i$yCell,sep="_")
  
  cellTab<-data.frame(Cell=names(table(fdat.i$Cell)),btENG=NA,btREA=NA,btCRU=NA)
  
  cellTab$btENG<-as.numeric(sapply(split(fdat.i$Engulf_task,fdat.i$Cell),sum))
  cellTab$btREA<-as.numeric(sapply(split(fdat.i$Reach_task,fdat.i$Cell),sum))
  cellTab$btCRU<-as.numeric(sapply(split(fdat.i$Crush_task,fdat.i$Cell),sum))
  
  cellTab$x<-as.numeric(unlist(strsplit(cellTab$Cell,"_"))[seq(1,nrow(cellTab)*2,2)])
  cellTab$y<-as.numeric(unlist(strsplit(cellTab$Cell,"_"))[seq(2,nrow(cellTab)*2,2)])
  
  cellTab$xMid<-xcuts[cellTab$x]+xdiff
  cellTab$yMid<-ycuts[cellTab$y]+ydiff
  
  cellTab$btENGScaled<-cellTab$btENG/sum(fdat.i$Engulf_task)
  cellTab$btREAScaled<-cellTab$btREA/sum(fdat.i$Reach_task)
  cellTab$btCRUScaled<-cellTab$btCRU/sum(fdat.i$Crush_task)
  
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
  axis(2,at=seq(-4,4,2),labels=c("-4","-2"," 0"," 2"," 4"),cex.axis=0.75,hadj=-0.3,las=2,tck=-0.0125,lwd=0.75)
  mtext("PC2",side=1,outer=F,line=0,cex=0.75,at=max(xcuts))
  mtext("PC3",side=2,outer=F,line=0,cex=0.75,at=max(ycuts)+0.5)
  mtext(paste0("Upweighted ",arch.names[arch.i]),side=3,outer=F,line=0,cex=0.75)
  mtext(LETTERS[arch.i], side=3, line=-1.3, adj=0.02, font=2, cex=1)  
  
  for(i in 1:nrow(cellTab)){
    polygon(c(cellTab$xMid[i]-xdiff,cellTab$xMid[i]+xdiff,cellTab$xMid[i]+xdiff,cellTab$xMid[i]-xdiff),
            c(cellTab$yMid[i]-ydiff,cellTab$yMid[i]-ydiff,cellTab$yMid[i]+ydiff,cellTab$yMid[i]+ydiff),
            col=cellTabCols[i,colCol[i]],border="white",lwd=0.25)
  }
  points(a3Rxy,pch=16,cex=1.5)
  
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
  text(x1[5],y1[5]-0.15,"Crush",cex=0.7,pos=3)
  text(x2[5]+0.3,y2[5]+0.1,"Reach",cex=0.7,pos=1)
  text(x3[5]-0.3,y3[5]+0.1,"Engulf",cex=0.7,pos=1)
  

} # End for sim i

### Body plots (D,E,F)

load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Body.rda"))
a3Rxy<-archCoords[[3]]

arch.list <- c("FLY","WLK","SWM")
arch.names <- c("Fly","Walk","Swim")

for (arch.i in 1:length(arch.list)){
  
  # Add beak task enrichment
  bodyTaskEnrich <- read.csv(paste0(outFolder,"Data4_BodyTaskEnrichment_Sensi_Max",arch.list[arch.i],".csv"))
  names(bodyTaskEnrich)
  fdat.i <- merge(pc.dat23body,bodyTaskEnrich,by="species")

  xcuts<-seq(min(fdat.i$PC2_body),max(fdat.i$PC2_body)+1.5,length.out=40)
  ycuts<-seq(min(fdat.i$PC3_body),max(fdat.i$PC3_body),length.out=40)
  
  xdiff<-(unique(diff(xcuts),6)/2)[1]
  ydiff<-(unique(diff(ycuts),6)/2)[1]
  
  fdat.i$xCell<-cut(fdat.i$PC2_body,xcuts,labels=FALSE,include.lowest = TRUE)
  fdat.i$yCell<-cut(fdat.i$PC3_body,ycuts,labels=FALSE,include.lowest = TRUE)
  fdat.i$Cell<-paste(fdat.i$xCell,fdat.i$yCell,sep="_")
  
  xcuts<-seq(min(fdat.i$PC2_body),max(fdat.i$PC2_body)+1.5,length.out=40)
  ycuts<-seq(min(fdat.i$PC3_body),max(fdat.i$PC3_body),length.out=40)
  
  xdiff<-(unique(diff(xcuts),6)/2)[1]
  ydiff<-(unique(diff(ycuts),6)/2)[1]
  
  fdat.i$xCell<-cut(fdat.i$PC2_body,xcuts,labels=FALSE,include.lowest = TRUE)
  fdat.i$yCell<-cut(fdat.i$PC3_body,ycuts,labels=FALSE,include.lowest = TRUE)
  fdat.i$Cell<-paste(fdat.i$xCell,fdat.i$yCell,sep="_")
  
  cellTab<-data.frame(Cell=names(table(fdat.i$Cell)),ftFLY=NA,ftSWM=NA,ftWLK=NA)
  
  cellTab$ftFLY<-as.numeric(sapply(split(fdat.i$Fly_task,fdat.i$Cell),sum))
  cellTab$ftSWM<-as.numeric(sapply(split(fdat.i$Swim_task,fdat.i$Cell),sum))
  cellTab$ftWLK<-as.numeric(sapply(split(fdat.i$Walk_task,fdat.i$Cell),sum))
  
  cellTab$x<-as.numeric(unlist(strsplit(cellTab$Cell,"_"))[seq(1,nrow(cellTab)*2,2)])
  cellTab$y<-as.numeric(unlist(strsplit(cellTab$Cell,"_"))[seq(2,nrow(cellTab)*2,2)])
  
  cellTab$xMid<-xcuts[cellTab$x]+xdiff
  cellTab$yMid<-ycuts[cellTab$y]+ydiff
  
  cellTab$ftFLYScaled<-cellTab$ftFLY/sum(fdat.i$Fly_task)
  cellTab$ftSWMScaled<-cellTab$ftSWM/sum(fdat.i$Swim_task)
  cellTab$ftWLKScaled<-cellTab$ftWLK/sum(fdat.i$Walk_task)
  
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
  axis(1,at=seq(-4,3,2),labels=seq(-4,3,2),cex.axis=0.75,padj=-2.5,tck=-0.0125,lwd=0.75)
  axis(2,at=seq(-4,4,2),labels=c("-4","-2"," 0"," 2"," 4"),cex.axis=0.75,hadj=-0.3,las=2,tck=-0.0125,lwd=0.75)
  mtext("PC2",side=1,outer=F,line=0,cex=0.75,at=max(xcuts))
  mtext("PC3",side=2,outer=F,line=0,cex=0.75,at=max(ycuts)+0.5)
  mtext(paste0("Upweighted ",arch.names[arch.i]),side=3,outer=F,line=0,cex=0.75)
  mtext(LETTERS[arch.i+3], side=3, line=-1.3, adj=0.02, font=2, cex=1)  
  
  for(i in 1:nrow(cellTab)){
    polygon(c(cellTab$xMid[i]-xdiff,cellTab$xMid[i]+xdiff,cellTab$xMid[i]+xdiff,cellTab$xMid[i]-xdiff),
            c(cellTab$yMid[i]-ydiff,cellTab$yMid[i]-ydiff,cellTab$yMid[i]+ydiff,cellTab$yMid[i]+ydiff),
            col=cellTabCols[i,colCol[i]],border="white",lwd=0.25)
  }
  points(a3Rxy,pch=16,cex=1.5)
  
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
    points(x1[i],y1[i],pch=21,bg=SWMCol[i-1],cex=1,lwd=0.5)
    points(x2[i],y2[i],pch=21,bg=WLKCol[i-1],cex=1,lwd=0.5)
    points(x3[i],y3[i],pch=21,bg=FLYCol[i-1],cex=1,lwd=0.5)
  }
  text(x1[2:5]-0.15,y1[2:5]-0.1,c("50","66","83","100%"),cex=0.4,pos=4)
  text(x1[5],y1[5]-0.15,"Swim",cex=0.7,pos=3)
  text(x2[5]+0.3,y2[5]+0.1,"Walk",cex=0.7,pos=1)
  text(x3[5]-0.3,y3[5]+0.1,"Fly",cex=0.7,pos=1)
  
} # End for sim i

dev.off()

## End script 4.5 ##