## 1.2 ## Plot Figure S1

load(file=paste0(OutputFolder,"out1_NullModels_1.Rdata"))

## Plot Figure S1 ####

plot.f <- paste0(FigureFolder,"Figure_S1.pdf")
pdf(plot.f,w=3.135,h=3.135)

par(mfrow=c(2,1))
par(mar=c(1.2,2.5,0.3,1))

Metric<-c("CubeNM1.MV","SphereNM2.1.MV","SphereNM2.2.MV")

plot(1:(1+length(Metric)),c(dat$Med[match(Metric,dat$Metric)],EmpiricalVol),xaxt="n",yaxt="n",ylim=c(0,100))

for(i in 1:length(Metric)){
  lines(c(i,i),c(dat$LQ[match(Metric,dat$Metric)][i],dat$UQ[match(Metric,dat$Metric)][i]))
  points(i,dat$Med[match(Metric,dat$Metric)][i],pch=21,bg="white")
}
axis(1,at=seq(1,length(Metric)+1,1),labels=c("NM 1\ncube","NM 2.1\nsphere","NM 2.2\nsphere","Observed\n"),
     cex.axis=0.5,padj=-1,tck=-0.0125,lwd=0.75)
axis(2,at=seq(20,80,20),labels=seq(20,80,20),cex.axis=0.75,hadj=-0.15,las=2,tck=-0.0125,lwd=0.75)
axis(2,at=c(0,100),labels=c("  0","100"),cex.axis=0.75,hadj=0.2,las=2,tck=-0.0125,lwd=0.75)
text(1,100*0.975,"A",font=2)
mtext("Beak morphospace\nvolume",side=2,outer=F,line=1,cex=0.75)

Metric<-c("CubeNM1.BDiff.Cone","SphereNM2.1.BDiff.Cone","SphereNM2.2.BDiff.Cone")

plot(1:5,c(dat$Med[match(Metric,dat$Metric)],EmpiricalBeakVolumeDiffCone,EmpiricalBeakVolumeMaxDiffCone),xaxt="n",yaxt="n",ylim=c(0,7))
for(i in 1:5){
  lines(c(i,i),c(dat$LQ[match(Metric,dat$Metric)][i],dat$UQ[match(Metric,dat$Metric)][i]))
  points(i,dat$Med[match(Metric,dat$Metric)][i],pch=21,bg="white")
}
axis(1,at=seq(1,5,1),
     labels=c("NM 1\ncube","NM 2.1\nsphere","NM 2.2\nsphere","Observed\n","Geometric\nlimit"),
     cex.axis=0.5,padj=-1,tck=-0.0125,lwd=0.75)
axis(2,at=seq(0,7,1),labels=seq(0,7,1),cex.axis=0.75,hadj=-1,las=2,tck=-0.0125,lwd=0.75)
text(1,7*0.975,"B",font=2)
mtext("Variation in\nBeak size",side=2,outer=F,line=1,cex=0.75)

dev.off()

## Plot Figure 1 ####

## Each plot will be displayed in 3D in separate window (i.e. they can be explored/rotated)

# Plot Observed (Empirical):

BubblePlot(Data=as.data.frame(x),traits=c("Length","Width","Depth"),PlotBox=TRUE,PlotAxis=TRUE,
           PlotAxisLabels=TRUE,
           PlotPoints=TRUE,
           bubbleCol="yellow",bubbleSize=0.075,
           paintWalls=TRUE,
           GridNBreaks=10,
           wallColour=colorRampPalette(c("white","black"))(10))

ts.surf <- t(convhulln(x))
xx<-x[ts.surf,]
spheres3d(xx,radius=0.1,col="forestgreen")
triangles3d(x[ts.surf,1],x[ts.surf,2],x[ts.surf,3],col="green",alpha=.4)
segments3d(x[ts.surf,1],x[ts.surf,2],x[ts.surf,3],col="forestgreen")

# Plot Sphere:

BubblePlot(Data=as.data.frame(x2),traits=c("Length","Width","Depth"),PlotBox=TRUE,PlotAxis=TRUE,
           PlotAxisLabels=TRUE,
           PlotPoints=TRUE,
           bubbleCol="yellow",bubbleSize=0.075,
           paintWalls=TRUE,
           GridNBreaks=10,
           wallColour=colorRampPalette(c("white","black"))(10))

ts.surf <- t(convhulln(x2))
xx<-x2[ts.surf,]
spheres3d(xx,radius=0.1,col="forestgreen")
triangles3d(x2[ts.surf,1],x2[ts.surf,2],x2[ts.surf,3],col="green",alpha=.4)
segments3d(x2[ts.surf,1],x2[ts.surf,2],x2[ts.surf,3],col="forestgreen")

# Plot Cube:

BubblePlot(Data=as.data.frame(x3),traits=c("Length","Width","Depth"),PlotBox=TRUE,PlotAxis=TRUE,
           PlotAxisLabels=TRUE,
           PlotPoints=TRUE,
           bubbleCol="yellow",bubbleSize=0.075,
           paintWalls=TRUE,
           GridNBreaks=10,
           wallColour=colorRampPalette(c("white","black"))(10))

ts.surf <- t(convhulln(x3))
xx<-x3[ts.surf,]
spheres3d(xx,radius=0.1,col="forestgreen")
triangles3d(x3[ts.surf,1],x3[ts.surf,2],x3[ts.surf,3],col="green",alpha=.4)
segments3d(x3[ts.surf,1],x3[ts.surf,2],x3[ts.surf,3],col="forestgreen")

## End script 1.2 ##