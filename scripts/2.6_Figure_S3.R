## 2.6 ## Plot Figure S3 (Regular polygons)

plot.f<-paste0(FigureFolder,"Figure_S3.pdf")
pdf(plot.f,w=8,h=6.5)

cols <- brewer.pal(6, "Dark2")
mat <- matrix(ncol=5,nrow=4)
mat[1:2,1:2]<-1
mat[3:4,1:2]<-2
mat[1,3:5]<-3:5
mat[2,3:5]<-6:8
mat[3,3:5]<-9:11
mat[4,3:5]<-12:14

layout(mat)
par(mar=c(2.5,2.7,2,2))

for(x in 1:2){
  if(x==1){
    load(paste0(OutputFolder,"out2_RegularPolygonFitting_Beak.rda"))
    pc.dat23<-read.csv(paste0(DataFolder,"data_tmp/BeakPC23.csv"))
  }
  if(x==2){
    load(paste0(OutputFolder,"out2_RegularPolygonFitting_Body.rda"))
    pc.dat23<-read.csv(paste0(DataFolder,"data_tmp/BodyPC23.csv"))
  }
  
  plot(borderDist[,1],type="n",ylim=c(min(borderDist),max(borderDist)),yaxt="n",xaxt="n",ylab="",xlab="")
  for(ii in 1:6){
    points(borderDist[,ii],type="l",col=cols[ii])
  }  
  axis(1,at=seq(0,80,10),labels=seq(0,80,10),
       cex.axis=0.75,padj=-2,tck=-0.0125,lwd=0.75)
  axis(2,at=seq(min(borderDist),max(borderDist),length.out=5),labels=round(seq(min(borderDist),max(borderDist),length.out=5),0),
       cex.axis=0.75,hadj=0.3,las=2,tck=-0.0125,lwd=0.75)
  text(1,max(borderDist),letters[x],font=2)
  mtext("Polygon area",side=1,outer=FALSE,line=1.25,cex=0.75)
  mtext("Error",side=2,outer=FALSE,line=1.7,cex=0.75)
  if(x==1){mtext("Beak shape",side=4,outer=FALSE,line=1)}
  if(x==2){mtext("Body shape",side=4,outer=FALSE,line=1)}
  
}

par(mar=c(0,0,0,0))
for(x in 1:2){
  if(x==1){
    load(paste0(OutputFolder,"out2_RegularPolygonFitting_Beak.rda"))
    pc.dat23<-read.csv(paste0(DataFolder,"data_tmp/BeakPC23.csv"))
  }
  if(x==2){
    load(paste0(OutputFolder,"out2_RegularPolygonFitting_Body.rda"))
    pc.dat23<-read.csv(paste0(DataFolder,"data_tmp/BodyPC23.csv"))
  }
  for(ii in 1:6){
    plot(pc.dat23[,2],pc.dat23[,3],col="grey80",pch=".",xlim=c(-5,5),ylim=c(-5,5),yaxt="n",xaxt="n",bty="n")
    if(x==1){text(-3,3,letters[ii+2],font=2)}
    if(x==2){text(-3,3,letters[ii+8],font=2)}
    
    shape<-regular.poly(nSides=OptShape$nVertex[ii], area=OptShape$Area[ii])
    
    shape<-Rotate(x=shape$x, y=shape$y,my=NULL,theta=deg2rad(OptShape$Rotation[ii]))[1:2]
    
    shapeCoord<-matrix(ncol=2,nrow=OptShape$nVertex[ii])
    shapeCoord[,1]<-shape$x
    shapeCoord[,2]<-shape$y
    
    chull<-shapeCoord[chull(shapeCoord),]
    chull.poly <- Polygon(chull, hole=F)
    
    ps <- Polygons(list(chull.poly),1)
    sps <- SpatialPolygons(list(ps))
    #points(datList[[x]],pch=".")
    plot(sps,add=TRUE,col=FALSE,border=cols[ii],lwd=2)
    text(0,-3.5,paste0("Min error = ", round(OptShape$borderDist[ii]),""))
  }
}
dev.off()

## End script 2.6 ##