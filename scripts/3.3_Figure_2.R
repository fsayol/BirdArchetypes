## 3.3 ## Figure 2: Scatterplot of beak and body morphospace, with enrichment analysis.

plot.f<-paste0(FigureFolder,"Figure_2.pdf")
pdf(plot.f,w=3.135*2.5,h=3.135*2)

# Open beak and body shape datasets (PCs):
fdat.beak <- read.csv(paste0(DataFolder,"data_tmp/BeakPC23.csv"))
fdat.body <- read.csv(paste0(DataFolder,"data_tmp/BodyPC23.csv"))

# Define matrix for panels:
mat<-matrix(ncol=6,nrow=5)

mat[1:3,1:3]<-1
mat[1:3,4:6]<-2
mat[4:5,1:2]<-3
mat[4:5,3:4]<-4
mat[4:5,5:6]<-5

layout(mat)
par(mar=c(2,2,1,1))

# Open beak archetypes:
(load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Beaks.rda")))
a3Rxy<-archCoords[[3]]
xmin<-min(fdat.beak$PC2_beak)-1
ymin<-min(fdat.beak$PC3_beak)-1
xmax<-max(fdat.beak$PC2_beak)
ymax<-max(fdat.beak$PC3_beak)

plot(fdat.beak$PC2_beak,fdat.beak$PC3_beak,ylim=c(ymin,ymax),
     xlim=c(xmin,xmax),
     pch=19,cex=0.4,lwd=0.4,col="grey70",ylab="",xlab="",yaxt="n",xaxt="n",bty="n")

arrows(x0=xmin, y0=ymin, x1 = xmin+(xmax-xmin)/10, y1 = ymin, length = 0.05, angle = 30,code = 2)
arrows(x0=xmin, y0=ymin, x1 = xmin, y1 = ymin+(ymax-ymin)/10, length = 0.05, angle = 30,code = 2)
mtext("Beak PC2",side=1,outer=F,line=-0.25,cex=0.75,at=xmin+(xmax-xmin)/10)
mtext("Beak PC3",side=2,outer=F,line=-0.25,cex=0.75,at=ymin+(ymax-ymin)/10)

text(xmin,ymax,"A",font=2,cex=2)

points(a3Rxy,pch=16,col="black",cex=2)
text(a3Rxy[3,1]-0.375,a3Rxy[3,2]+0.35,"A1",font=2,col="black",cex=1)
text(a3Rxy[2,1]+0.45,a3Rxy[2,2]+0.35,"A2",font=2,col="black",cex=1)
text(a3Rxy[1,1],a3Rxy[1,2]-0.5,"A3",font=2,col="black",cex=1)

# Open body archetypes:
load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Body.rda"))
a3Rxy<-archCoords[[3]]
xmin<-min(fdat.body$PC2_body)-1
ymin<-min(fdat.body$PC3_body)
xmax<-max(fdat.body$PC2_body)
ymax<-max(fdat.body$PC3_body)+1.5

plot(fdat.body$PC2_body,fdat.body$PC3_body,ylim=c(ymin,ymax),
     xlim=c(xmin,xmax),
     pch=19,cex=0.4,lwd=0.4,col="grey70",ylab="",xlab="",yaxt="n",xaxt="n",bty="n")

arrows(x0=xmin, y0=ymin, x1 = xmin+(xmax-xmin)/10, y1 = ymin, length = 0.05, angle = 30,code = 2)
arrows(x0=xmin, y0=ymin, x1 = xmin, y1 = ymin+(ymax-ymin)/10, length = 0.05, angle = 30,code = 2)
mtext("Body PC2",side=1,outer=F,line=-0.25,cex=0.75,at=xmin+(xmax-xmin)/10)
mtext("Body PC3",side=2,outer=F,line=-0.25,cex=0.75,at=ymin+(ymax-ymin)/10)

points(a3Rxy,pch=16,col="black",cex=2)
text(a3Rxy[3,1]+0.5,a3Rxy[3,2]-0.4,"A2",font=2,col="black",cex=1)
text(a3Rxy[2,1],a3Rxy[2,2]+0.5,"A1",font=2,col="black",cex=1)
text(a3Rxy[1,1]-0.4,a3Rxy[1,2]-0.5,"A3",font=2,col="black",cex=1)

text(xmin,ymax,"B",font=2,cex=2)

## enrichment plots (panels c / d / e) ####

par(mar=c(2,2,1,1))

load(paste0(OutputFolder,"out3_BeakNicheEnrichment.rda"))

colgr <- brewer.pal(n = 10, name = "Paired")
colCol <- apply(cellTab[,16:25],1,which.max)

plot(0,0,type="n",xlim=c(min(xcuts)-0.5,max(xcuts)),ylim=c(min(ycuts),max(ycuts)),ylab="",xlab="",yaxt="n",xaxt="n")
axis(1,at=seq(-4,5,2),labels=seq(-4,5,2),cex.axis=0.75,padj=-2.5,tck=-0.0125,lwd=0.75)
axis(2,at=seq(-4,4,2),labels=c("-4","-2"," 0"," 2"," 4"),cex.axis=0.75,hadj=-0.3,las=2,tck=-0.0125,lwd=0.75)
mtext("PC2",side=1,outer=F,line=0,cex=0.75,at=max(xcuts))
mtext("PC3",side=2,outer=F,line=0,cex=0.75,at=max(ycuts)+0.5)
mtext("Beak shape & Feeding technique",side=3,outer=F,line=0.2,cex=0.75)

for(i in 1:nrow(cellTab)){
  polygon(c(cellTab$xMid[i]-xdiff,cellTab$xMid[i]+xdiff,cellTab$xMid[i]+xdiff,cellTab$xMid[i]-xdiff),
          c(cellTab$yMid[i]-ydiff,cellTab$yMid[i]-ydiff,cellTab$yMid[i]+ydiff,cellTab$yMid[i]+ydiff),
          col=colgr[colCol[i]],border="white",lwd=0.25)
}
points(a3Rxy,pch=16,cex=1.5)
TaskNames<-c("Tearing","Crushing","Filtering","Sweeping","Engulfing","Grazing","Hammering","Plucking","Probing","Spearing") 
points(rep(4.5,10),seq(-5.5,-2,length.out=10),pch=21,bg=colgr,lwd=0.5)
text(rep(4.5,10),seq(-5.5,-2,length.out=10),TaskNames,cex=0.5,pos=4)

text(min(xcuts)+(diff(range(xcuts))*0.025)-0.5,max(ycuts)-(diff(range(ycuts))*0.025),"C",font=2,cex=2)

load(paste0(OutputFolder,"out3_BeakTaskEnrichment.rda"))

cellTabCols<-cellTab[,15:17]
colCol<-apply(cellTab[,9:11],1,which.max)

plot(0,0,type="n",xlim=c(min(xcuts)-0.5,max(xcuts)),ylim=c(min(ycuts),max(ycuts)),ylab="",xlab="",yaxt="n",xaxt="n")
axis(1,at=seq(-4,5,2),labels=seq(-4,5,2),cex.axis=0.75,padj=-2.5,tck=-0.0125,lwd=0.75)
axis(2,at=seq(-4,4,2),labels=c("-4","-2"," 0"," 2"," 4"),cex.axis=0.75,hadj=-0.3,las=2,tck=-0.0125,lwd=0.75)
mtext("PC2",side=1,outer=F,line=0,cex=0.75,at=max(xcuts))
mtext("PC3",side=2,outer=F,line=0,cex=0.75,at=max(ycuts)+0.5)
mtext("Beak shape & Feeding objective",side=3,outer=F,line=0.2,cex=0.75)

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
text(x1[2:5]-0.2,y1[2:5]-0.1,c("50","66","83","100%"),cex=0.4,pos=4)
text(x1[5],y1[5]-0.2,"Crush",cex=0.7,pos=3)
text(x2[5]+0.3,y2[5]+0.2,"Reach",cex=0.7,pos=1)
text(x3[5]-0.3,y3[5]+0.2,"Engulf",cex=0.7,pos=1)

text(min(xcuts)+(diff(range(xcuts))*0.0125)-0.5,max(ycuts)-(diff(range(ycuts))*0.025),"D",font=2,cex=2)

load(paste0(OutputFolder,"out3_BodyTaskEnrichment.rda"))

cellTabCols<-cellTab[,15:17]
colCol<-apply(cellTab[,9:11],1,which.max)

plot(0,0,type="n",xlim=c(min(xcuts),max(xcuts)),ylim=c(min(ycuts),max(ycuts)),ylab="",xlab="",yaxt="n",xaxt="n")
axis(1,at=seq(-4,8,2),labels=seq(-4,8,2),cex.axis=0.75,padj=-2.5,tck=-0.0125,lwd=0.75)
axis(2,at=seq(-4,4,2),labels=c("-4","-2"," 0"," 2"," 4"),cex.axis=0.75,hadj=-0.3,las=2,tck=-0.0125,lwd=0.75)
mtext("PC2",side=1,outer=F,line=0,cex=0.75,at=max(xcuts))
mtext("PC3",side=2,outer=F,line=0,cex=0.75,at=max(ycuts)+0.5)
mtext("Body shape & Foraging objective",side=3,outer=F,line=0.2,cex=0.75)

for(i in 1:nrow(cellTab)){
  polygon(c(cellTab$xMid[i]-xdiff,cellTab$xMid[i]+xdiff,cellTab$xMid[i]+xdiff,cellTab$xMid[i]-xdiff),
          c(cellTab$yMid[i]-ydiff,cellTab$yMid[i]-ydiff,cellTab$yMid[i]+ydiff,cellTab$yMid[i]+ydiff),
          col=cellTabCols[i,colCol[i]],border="white",lwd=0.25)
}
points(a3Rxy,pch=16,cex=1.5)

xrange<-range(xcuts)
yrange<-range(ycuts)

xOr<-min(xcuts)+(0.75*diff(xrange))
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
text(x1[2:5]-0.3,y1[2:5]-0.1,c("50","66","83","100%"),cex=0.4,pos=4)
text(x1[5],y1[5]-0.2,"Swim",cex=0.7,pos=3)
text(x2[5]+0.3,y2[5]+0.2,"Walk",cex=0.7,pos=1)
text(x3[5]-0.3,y3[5]+0.2,"Fly",cex=0.7,pos=1)

text(min(xcuts)+(diff(range(xcuts))*0.025),max(ycuts)-(diff(range(ycuts))*0.025),"E",font=2,cex=2)

dev.off()

## End script 3.3 ##