# functions used along the scripts (loaded with setup script)

## Make transparent function ####
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

## Regular polygons function ####
regular.poly <- function(nSides, area)
{
  # Find the radius of the circumscribed circle
  radius <- sqrt((2*area)/(nSides*sin((2*pi)/nSides)))
  
  # I assume the center is at (0;0) and the first point lies at (0; radius)
  points <- list(x=NULL, y=NULL)
  angles <- (2*pi)/nSides * 1:nSides
  
  points$x <- cos(angles) * radius
  points$y <- sin(angles) * radius
  
  return (points);
}

rsphere <- function(n, r = 1.0, surface_only = FALSE) {
  phi       <- runif(n, 0.0, 2.0 * pi)
  cos_theta <- runif(n, -1.0, 1.0)
  sin_theta <- sqrt((1.0-cos_theta)*(1.0+cos_theta))
  radius <- r
  if (surface_only == FALSE) {
    radius <- r * runif(n, 0.0, 1.0)^(1.0/3.0)
  }
  
  x <- radius * sin_theta * cos(phi)
  y <- radius * sin_theta * sin(phi)
  z <- radius * cos_theta
  
  cbind(x, y, z)
}

## Bubble Plot function ####
BubblePlot<-function(Data,traits,
                     PlotBox=TRUE,PlotAxis=FALSE,
                     PlotAxisLabels=TRUE,bubbleCol="yellow",
                     PlotPoints=FALSE,
                     bubbleSize=0.075,paintWalls=TRUE,
                     GridNBreaks=40,
                     wallColour=colorRampPalette(c("white","black"))(10)){
  
  colors <- wallColour
  
  foccol<-match(traits,names(Data))
  Data<-Data[,foccol]
  
  PC1min<-min(Data[,1])
  PC1max<-max(Data[,1])
  PC2min<-min(Data[,2])
  PC2max<-max(Data[,2])
  PC3min<-min(Data[,3])
  PC3max<-max(Data[,3])
  
  if(PlotAxis==TRUE){
    if(PlotAxisLabels==TRUE){
      plot3d(Data,box=PlotBox,col="grey95",size=0.1,pch=1,xlim=c(PC1min,PC1max),ylim=c(PC2min,PC2max),zlim=c(PC3min,PC3max))
    }else{
      plot3d(Data,box=PlotBox,col="grey95",size=0.1,pch=1,xlim=c(PC1min,PC1max),ylim=c(PC2min,PC2max),zlim=c(PC3min,PC3max),xlab="",ylab="",zlab="")
    }	
  }else{
    plot3d(Data,box=PlotBox,col="grey95",size=0.1,pch=1,xlim=c(PC1min,PC1max),ylim=c(PC2min,PC2max),zlim=c(PC3min,PC3max),xlab="",ylab="",zlab="",xaxt="n",yaxt="n",zaxt="n")
  }
  
  if(paintWalls==TRUE){
    
    nbreaks<-GridNBreaks
    
    brks <- classIntervals(Data[,1], n=nbreaks, style="equal");brks1 <- brks$brks
    Data$PC1int <- findInterval(Data[,1], brks1, all.inside = TRUE)
    
    brks <- classIntervals(Data[,2], n=nbreaks, style="equal");brks2 <- brks$brks
    Data$PC2int <- findInterval(Data[,2], brks2, all.inside = TRUE)
    
    brks <- classIntervals(Data[,3], n=nbreaks, style="equal");brks3 <- brks$brks
    Data$PC3int <- findInterval(Data[,3], brks3, all.inside = TRUE)
    
    Data$cell12<-rep(NA,dim(Data)[1])
    for(i in 1:dim(Data)[1]){Data$cell12[i]<-paste(Data$PC1int[i],Data$PC2int[i],sep="")}
    richness<-table(Data$cell12)
    gridTable<-data.frame(cell=names(richness),richness=as.numeric(richness))
    gridTable$x<-Data$PC1int[match(gridTable$cell,Data$cell12)]
    gridTable$y<-Data$PC2int[match(gridTable$cell,Data$cell12)]
    
    brksr <- classIntervals(gridTable$richness, n=10, style="quantile")
    brks <- brksr$brks
    col1 <- findInterval(gridTable$richness, brks, all.inside = TRUE)
    gridTable$cols<-as.character(findColours(brksr,wallColour)) 
    
    for(i in 1:(length(brks1)-1)){
      xmin<-brks1[i]
      xmax<-brks1[i+1]
      for(it in 1:(length(brks2)-1)){
        ymin<-brks2[it]
        ymax<-brks2[it+1]
        focCell<-which(gridTable$x==i & gridTable$y==it)
        if(length(focCell)>0){
        }else{
          polygon3d(x=c(xmin,xmin,xmax,xmax),y=c(ymin,ymax,ymax,ymin),z=c(PC3min+0.001,PC3min+0.001,PC3min+0.001,PC3min+0.001), col="grey80",fill=FALSE)	
        }
      }
    }
    
    
    for(i in 1:(length(brks1)-1)){
      xmin<-brks1[i]
      xmax<-brks1[i+1]
      for(it in 1:(length(brks2)-1)){
        ymin<-brks2[it]
        ymax<-brks2[it+1]
        focCell<-which(gridTable$x==i & gridTable$y==it)
        if(length(focCell)>0){
          focCol<-gridTable$cols[focCell]
          focRichness<-gridTable$richness[focCell]/max(gridTable$richness)
          focRichness<-focRichness*(PC3max-PC3min)*0.05
          polygon3d(x=c(xmin,xmin,xmax,xmax),y=c(ymin,ymax,ymax,ymin),z=c(PC3min+0.001,PC3min+0.001,PC3min+0.001,PC3min+0.001),col=focCol)
          polygon3d(x=c(xmin,xmin,xmax,xmax),y=c(ymin,ymax,ymax,ymin),z=c(PC3min+0.001,PC3min+0.001,PC3min+0.001,PC3min+0.001),col="white",fill=FALSE)
        }
      }
    }
    
    Data$cell23<-rep(NA,dim(Data)[1])
    for(i in 1:dim(Data)[1]){Data$cell23[i]<-paste(Data$PC2int[i],Data$PC3int[i],sep="")}
    richness<-table(Data$cell23)
    gridTable<-data.frame(cell=names(richness),richness=as.numeric(richness))
    gridTable$x<-Data$PC2int[match(gridTable$cell,Data$cell23)]
    gridTable$y<-Data$PC3int[match(gridTable$cell,Data$cell23)]
    
    brksr <- classIntervals(gridTable$richness, n=10, style="quantile")
    brks <- brksr$brks
    col1 <- findInterval(gridTable$richness, brks, all.inside = TRUE)
    gridTable$cols<-as.character(findColours(brksr,wallColour)) 
    
    for(i in 1:(length(brks2)-1)){
      xmin<-brks2[i]
      xmax<-brks2[i+1]
      for(it in 1:(length(brks3)-1)){
        ymin<-brks3[it]
        ymax<-brks3[it+1]
        focCell<-which(gridTable$x==i & gridTable$y==it)
        if(length(focCell)>0){
        }else{
          polygon3d(x=c(PC1min,PC1min,PC1min+0.001,PC1min+0.001),y=c(xmin,xmax,xmax,xmin),z=c(ymin,ymin,ymax,ymax),col="grey80",fill=FALSE)
        }
      }
    }
    
    
    for(i in 1:(length(brks2)-1)){
      xmin<-brks2[i]
      xmax<-brks2[i+1]
      for(it in 1:(length(brks3)-1)){
        ymin<-brks3[it]
        ymax<-brks3[it+1]
        focCell<-which(gridTable$x==i & gridTable$y==it)
        if(length(focCell)>0){
          focCol<-gridTable$cols[focCell]
          polygon3d(x=c(PC1min,PC1min,PC1min+0.001,PC1min+0.001),y=c(xmin,xmax,xmax,xmin),z=c(ymin,ymin,ymax,ymax),col=focCol)
          polygon3d(x=c(PC1min,PC1min,PC1min+0.001,PC1min+0.001),y=c(xmin,xmax,xmax,xmin),z=c(ymin,ymin,ymax,ymax),col="white",fill=FALSE)
        }
      }
    }
    
    Data$cell13<-rep(NA,dim(Data)[1])
    for(i in 1:dim(Data)[1]){Data$cell13[i]<-paste(Data$PC1int[i],Data$PC3int[i],sep="")}
    richness<-table(Data$cell13)
    gridTable<-data.frame(cell=names(richness),richness=as.numeric(richness))
    gridTable$x<-Data$PC1int[match(gridTable$cell,Data$cell13)]
    gridTable$y<-Data$PC3int[match(gridTable$cell,Data$cell13)]
    
    brksr <- classIntervals(gridTable$richness, n=10, style="quantile")
    brks <- brksr$brks
    col1 <- findInterval(gridTable$richness, brks, all.inside = TRUE)
    gridTable$cols<-as.character(findColours(brksr,wallColour)) 
    
    for(i in 1:(length(brks1)-1)){
      xmin<-brks1[i]
      xmax<-brks1[i+1]
      for(it in 1:(length(brks3)-1)){
        ymin<-brks3[it]
        ymax<-brks3[it+1]
        focCell<-which(gridTable$x==i & gridTable$y==it)
        if(length(focCell)>0){
        }else{
          polygon3d(x=c(xmin,xmax,xmax,xmin),y=c(PC2max,PC2max,PC2max-0.001,PC2max-0.001),z=c(ymax,ymax,ymin,ymin),col="grey80",fill=FALSE)
        }
      }
    }
    
    for(i in 1:(length(brks1)-1)){
      xmin<-brks1[i]
      xmax<-brks1[i+1]
      for(it in 1:(length(brks3)-1)){
        ymin<-brks3[it]
        ymax<-brks3[it+1]
        focCell<-which(gridTable$x==i & gridTable$y==it)
        if(length(focCell)>0){
          focCol<-gridTable$cols[focCell]
          polygon3d(x=c(xmin,xmax,xmax,xmin),y=c(PC2max,PC2max,PC2max-0.001,PC2max-0.001),z=c(ymax,ymax,ymin,ymin),col=focCol)
          polygon3d(x=c(xmin,xmax,xmax,xmin),y=c(PC2max,PC2max,PC2max-0.001,PC2max-0.001),z=c(ymax,ymax,ymin,ymin),col="white",fill=FALSE)
        }
      }
    }
    
  }
  if(PlotPoints==TRUE){
    spheres3d(Data,radius=bubbleSize,col=bubbleCol)
  }
}

## Cell intersect function ####

# find which lineages intersect with each grid cell

cellIntersect <- function(pc.dat23,nCells,toplot=F,printRow=F){
  
  minX<-floor(min(pmin(pc.dat23[,c("pc2.Anc","pc2.Des")])))
  maxX<-ceiling(max(pmax(pc.dat23[,c("pc2.Anc","pc2.Des")])))
  minY<-floor(min(pmin(pc.dat23[,c("pc3.Anc","pc3.Des")])))
  maxY<-ceiling(max(pmax(pc.dat23[,c("pc3.Anc","pc3.Des")])))
  
  gridSize<-(maxX-minX)/nCells
  
  sfc <- st_sfc(st_polygon(list(rbind(c(minX,minY), c(maxX,minY), c(maxX,maxY), c(minX,minY)))))
  poly_dat<-st_make_grid(x=sfc,what = "polygons",cellsize=gridSize)
  poly_datxy<-st_coordinates(poly_dat)
  
  if(toplot==TRUE){plot(poly_dat)}
  
  ## identify which branches evolve through each grid cell
  cellIntersectList<-as.list(rep(NA,nrow(pc.dat23)))
  for(i in 1:nrow(pc.dat23)){
    line1 <- st_linestring(rbind(c(pc.dat23$pc2.Anc[i],pc.dat23$pc3.Anc[i]),c(pc.dat23$pc2.Des[i],pc.dat23$pc3.Des[i])))
    line_dat <- st_sf(linename = c("1"), geometry = st_sfc(line1))
    if(toplot==TRUE){plot(line_dat, add = TRUE, col = "red", lty = 2)}
    cellIntersectList[[i]]<-st_intersects(line_dat, poly_dat)[[1]]
    if(printRow==T){
      if(i%%100==0){print(i)}}
  }
  cellsPerNode<-sapply(cellIntersectList,length)
  cellIntersectDat<-data.frame(Node.Anc=rep(pc.dat23$Node.Anc,cellsPerNode),
                               Node.Des=rep(pc.dat23$Node.Des,cellsPerNode),
                               Cells=unlist(cellIntersectList))
  outlist<-list()
  outlist[[1]]<-cellIntersectDat
  outlist[[2]]<-sfc
  outlist[[3]]<-poly_dat
  outlist[[4]]<-poly_datxy
  return(outlist)
}

## Make cell function ####

makeCellTab<-function(pc.dat23,cellIntersectDat,poly_datxy,printRow=F){
  pc.dat23$UniqueID<-paste(pc.dat23$Node.Anc,pc.dat23$Node.Des,sep="_")
  cellIntersectDat$UniqueID<-paste(cellIntersectDat$Node.Anc,cellIntersectDat$Node.Des,sep="_")
  pc.dat23$DX<-pc.dat23$pc2.Des-pc.dat23$pc2.Anc
  pc.dat23$DY<-pc.dat23$pc3.Des-pc.dat23$pc3.Anc
  pc.dat23$rad2<-atan2(pc.dat23$DY,pc.dat23$DX)
  pc.dat23$deg2<-rad2deg(pc.dat23$rad2)
  pc.dat23$D<-sqrt(pc.dat23$DX^2+pc.dat23$DY^2)
  
  cellsFoc<-sort(unique(cellIntersectDat$Cells))
  keep<-rep(c(1,0,0,0,0),max(poly_datxy[,4]))
  poly_datxy2<-poly_datxy[which(keep==1),]
  
  cellTab<-data.frame(Cells=cellsFoc,X=poly_datxy2[cellsFoc,1],Y=poly_datxy2[cellsFoc,2],xdeg2=NA,nb=NA,xdegV=NA)
  for(i in 1:nrow(cellTab)){
    focRow<-which(cellIntersectDat$Cells==cellTab$Cells[i])
    focAnc<-cellIntersectDat$UniqueID[focRow]
    cellTab$nb[i]<-length(focAnc)
    focB<-match(focAnc,pc.dat23$UniqueID)
    cellTab$xdeg2[i]<-as.numeric(mean(circular(pc.dat23$deg2[focB],units="degrees")))
    cellTab$xdeg2W[i]<-as.numeric(weighted.mean(circular(pc.dat23$deg2[focB],units="degrees"),w=pc.dat23$D[focB]))
    cellTab$xdegV[i]<-as.numeric(angular.variance(circular(pc.dat23$deg2[focB],units="degrees")))
    if(printRow==T){print(i)}
  }
  cellTab$xdegV2<-max(cellTab$xdegV)-cellTab$xdegV
  return(cellTab)
}

## Arrows function ####

# Used to plot arrows (inside GAM plot function)
arrows.0 <- function(x0, y0, length.ar, angle.ar,col,lwd=lwd){
  
  ab <- cos(angle.ar) * length.ar
  bc <- sign(sin(angle.ar)) * sqrt(length.ar^2 - ab^2)
  
  x1 <- x0 + ab
  y1 <- y0 + bc
  
  arrows(x0, y0, x1, y1,length=0.020,angle=30,col=col,lwd=lwd)
}


## Plot GAM function ####

# Used to plot morphospace with GAM model output
# Delete empty cells and add arrows is optional

plot.predictedGAM <- function(pd,space,brks2,emptyCells=F,plotArrows=F,gridWidth=0.25,col.pal="YlOrRd",letter){
  
  if(emptyCells==T){
    pd[pd$spp==0,]$out <- NA}
  
  if(space=="beak"){
    load(paste0(outFolder,"Beak_arrows.rda"))
    pc.dat23<-read.csv(paste0(DataFolder,"data_tmp/BeakPC23.csv"))
    pc.dat23$PC2<-pc.dat23$PC2_beak
    pc.dat23$PC3<-pc.dat23$PC3_beak
    pd$PC2<-pd$PC2_beak
    pd$PC3<-pd$PC3_beak
    a3Rxy<-archCoords[[3]]}
  else{
    load(paste0(outFolder,"Body_arrows.rda"))
    pc.dat23<-read.csv(paste0(DataFolder,"data_tmp/BodyPC23.csv"))
    pc.dat23$PC2<-pc.dat23$PC2_body
    pc.dat23$PC3<-pc.dat23$PC3_body
    pd$PC2<-pd$PC2_body
    pd$PC3<-pd$PC3_body
    a3Rxy<-archCoords[[3]]}
  
  ydiff<-xdiff<-min(diff(sort(unique(pd$PC2))))
  xboundary<-xdiff*4
  pd<-pd[which(pd$PC2>=(min(a3Rxy[,1])-xboundary) & pd$PC2<(xboundary+max(a3Rxy[,1])) &
                 pd$PC3>=(min(a3Rxy[,2])-xboundary) & pd$PC3<(max(a3Rxy[,2])+xboundary)),]
  
  ymax<-max(pd$PC3)+0.1
  ymin<-min(pd$PC3)-0.1
  xmax<-max(pd$PC2)+0.1
  xmin<-min(pd$PC2)-0.1
  
  plot(pc.dat23$PC2,pc.dat23$PC3,type="n",bty="n",
       yaxt="n",xaxt="n",xlim=c(xmin,xmax),ylim=c(ymin,ymax),ylab="",xlab="")
  text(xmin,ymax,letter,font=2,cex=2)
  
  if(length(col.pal)==1){
    cols3 <- brewer.pal(9, col.pal)
  }else{
    cols3 <-col.pal  
  }
  cols4 <- colorRampPalette(cols3)(length(brks2)-1)
  pd$cellCol<-cols4[cut(pd$out,brks2,include.lowest=TRUE)]           
  
  for(i in 1:nrow(pd)){
    polygon(c(pd$PC2[i]-xdiff/2,pd$PC2[i]+xdiff/2,pd$PC2[i]+xdiff/2,pd$PC2[i]-xdiff/2),
            c(pd$PC3[i]-ydiff/2,pd$PC3[i]-ydiff/2,pd$PC3[i]+ydiff/2,pd$PC3[i]+ydiff/2),
            col=pd$cellCol[i],border="grey50",lwd=gridWidth)
  }
  
  points(a3Rxy,col="black",pch=16,cex=1.5)
  #polygon(a3Rxy,lwd=0.5)
  
  ## Arrows
  
  if(plotArrows==T){
    cellTab <- cellTabO
    for(i in 1:nrow(cellTab)){
      if(cellTab$RejectNull[i]==1){
        arrows.0(cellTab$X[i], cellTab$Y[i], 
                 length.ar= cellTab$xdegV2[i]/6, angle.ar=deg2rad(cellTab$xdeg2W[i]),
                 col = "black",lwd=0.6)
      }else{
        points(cellTab$X[i], cellTab$Y[i],pch=16,cex=0.2,col="black")}  
    } # end if inside
  } # end if outside arrows=T
  
} # end function

## End archetype functions script ##