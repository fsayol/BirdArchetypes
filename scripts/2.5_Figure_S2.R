## 2.5 ## Plot Figure S2 (Irregular polygons)

plot.f<-paste0(FigureFolder,"Figure_S2.pdf")
pdf(plot.f,w=8,h=6.5)

cols<-brewer.pal(6, "Dark2")
mat<-matrix(ncol=5,nrow=4)
mat[1:2,1:2]<-1
mat[3:4,1:2]<-2
mat[1,3:5]<-3:5
mat[2,3:5]<-6:8
mat[3,3:5]<-9:11
mat[4,3:5]<-12:14
layout(mat)
par(mar=c(2.5,3.7,2,1.5))
for(x in 1:2){
  if(x==1){
    load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Beaks.rda"))
    pc.dat23 <- read.csv(paste0(DataFolder,"data_tmp/BeakPC23.csv"))
  }
  if(x==2){
    load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Body.rda"))
    pc.dat23 <- read.csv(paste0(DataFolder,"data_tmp/BodyPC23.csv"))
  }
  
  rss<-rowMins(rssMat,na.rm=T)[3:7]
  rss<-rowMedians(rssMat,na.rm=T)[3:7]
  rssMin<-rowMins(rssMat,na.rm=T)[3:7]
  rssMax<-rowMaxs(rssMat,na.rm=T)[3:7]
  
  plot(1:5,1:5,type="n",ylim=c(min(rssMin),max(rssMax)),yaxt="n",xaxt="n",ylab="",xlab="")
  polygon(c(1:5,rev(1:5)),c(rssMin,rev(rssMax)),border=FALSE,col="grey90")
  points(1:5,rss,type="l")
  points(1:5,rss,col=cols,pch=16,cex=2)
  
  axis(1,at=1:5,labels=3:7,
       cex.axis=1,padj=-1.5,tck=-0.0125,lwd=0.75)
  axis(2,at=seq(min(rssMin),max(rssMax),length.out=5),labels=round(seq(min(rssMin),max(rssMax),length.out=5),4),
       cex.axis=1,hadj=0.7,las=2,tck=-0.0125,lwd=0.75)
  text(1,max(rssMax),letters[x],font=2,cex=1.5)
  mtext("Number of vertices",side=1,outer=FALSE,line=1.25,cex=0.75)
  mtext("Error",side=2,outer=FALSE,line=2.75,cex=0.75)
  if(x==1){mtext("Beak shape",side=4,outer=FALSE,line=1)}
  if(x==2){mtext("Body shape",side=4,outer=FALSE,line=1)}
  
}

par(mar=c(0,0,0,0))
for(x in 1:2){
  if(x==1){
    load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Beaks.rda"))
    pc.dat23 <- read.csv(paste0(DataFolder,"data_tmp/BeakPC23.csv"))
  }
  if(x==2){
    load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Body.rda"))
    pc.dat23 <- read.csv(paste0(DataFolder,"data_tmp/BodyPC23.csv"))
  }
  rss<-rowMedians(rssMat,na.rm=T)[3:7]
  
  for(ii in 1:5){
    plot(pc.dat23[,2],pc.dat23[,3],col="grey80",pch=".",xlim=c(-4.5,4.5),ylim=c(-4.5,4.5),yaxt="n",xaxt="n",bty="n")
    if(x==1){text(-3,3,letters[ii+2],font=2,cex=1.5)}
    if(x==2){text(-3,3,letters[ii+7],font=2,cex=1.5)}
    
    shapeCoord<-matrix(ncol=2,nrow=ii+2)
    shapeCoord[,1]<-archCoords[[ii+2]][,1]
    shapeCoord[,2]<-archCoords[[ii+2]][,2]
    
    chull<-shapeCoord[chull(shapeCoord),]
    chull.poly <- Polygon(chull, hole=F)
    
    ps <- Polygons(list(chull.poly),1)
    sps <- SpatialPolygons(list(ps))
    plot(sps,add=TRUE,col=FALSE,border=cols[ii],lwd=2)
    text(0,-4.5,paste0("Error = ", round(rss[ii],4),""))
  }
  plot(pc.dat23[,2],pc.dat23[,3],type="n",xlim=c(-4,4),ylim=c(-4,4),yaxt="n",xaxt="n",bty="n")
}
dev.off()

## End script 2.5 ##