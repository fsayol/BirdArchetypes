## 6.3 ## Plot Figure 4 (Spp richness, ER, DR)

# Folder to load Extinction risk and DR models (GAM models)
outFolder <-paste0(OutputFolder,'out6_ER_DR_Models/')

plot.f<-paste0(FigureFolder,"Figure_4.pdf")
pdf(plot.f,w=3.135*2.5,h=3.135*2)

# Load and save beak arrows for plotting:
load(paste0(OutputFolder,"out5_NodeTraits/BeakNodeTraits_Null_cellIntersections.rda"))
load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Beaks.rda"))
save(archCoords,rssMat,rssMatC,cellIntersectDatObs,cellTabO,file=paste0(outFolder,"beak_arrows.rda"))

# Load and save body arrows for plotting:
load(paste0(OutputFolder,"out5_NodeTraits/BodyNodeTraits_Null_cellIntersections.rda"))
load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Body.rda"))
save(archCoords,rssMat,rssMatC,cellIntersectDatObs,cellTabO,file=paste0(outFolder,"body_arrows.rda"))


cols<-c("gold","red", "darkred")

mat<-matrix(ncol=12,nrow=9)
mat[1:4,1:4]<-1
mat[5:8,1:4]<-2
mat[9,1:4]<-3

mat[1:4,5:8]<-4
mat[5:8,5:8]<-5
mat[9,5:8]<-6

mat[1:4,9:12]<-7
mat[5:8,9:12]<-8
mat[9,9:12]<-9

layout(mat)

par(mar=c(0,1.2,1.2,0))
par(oma=c(0,0,0,0))

## Species richness (observed) across beak [with arrows]
pred.er <- read.csv(paste0(outFolder,"GAMBeakExtinctionModelPrediction.csv"))
pred.er$out <- pred.er$spp
print(range(pred.er$spp))
brks <- seq(0,300,20)
plot.predictedGAM(pred.er,"beak",brks,emptyCells=T,plotArrows=T,gridWidth = 0.0001,col.pal = cols,letter="A")
mtext("Species richness",side=3,outer=F,line=0,cex=0.8)
mtext("Beak shape",side=2,outer=F,line=0,cex=0.8)

## Species richness (observed) across body [with arrows]
pred.er <- read.csv(paste0(outFolder,"GAMBodyExtinctionModelPrediction.csv"))
pred.er$out <- pred.er$spp
print(range(pred.er$spp))
brks <- seq(0,300,20)
plot.predictedGAM(pred.er,"body",brks,emptyCells=T,gridWidth = 0.0001,plotArrows = T,col.pal = cols,letter="B")
mtext("Body shape",side=2,outer=F,line=0,cex=0.8)

par(mar=c(0,0,0,0))
brks2 <-  seq(0,300,20)
if(length(cols)==1){
  cols3 <- brewer.pal(9, cols)
}else{
  cols3 <-cols  
}
cols4 <- colorRampPalette(cols3)(length(brks2)-1)
plot(0,0,type="n",yaxt="n",xaxt="n",xlim=c(-0.025,max(brks2)),ylim=c(0,1),bty="n",ylab="",xlab="")
for(i in 1:length(cols4)){polygon(c(brks2[i],brks2[i+1],brks2[i+1],brks2[i]),c(0.4,0.4,0.7,0.7),col=cols4[i],border="white",lwd=0.5)}
for(i in seq(1,length(brks2),by=2)){text(brks2[i],0.5,round(brks2[i],2),pos=1,cex=0.75)}

par(mar=c(0,1.2,1.2,0))

## Extinction risk
pred.er <- read.csv(paste0(outFolder,"GAMBeakExtinctionModelPrediction.csv"))
pred.er2 <- read.csv(paste0(outFolder,"GAMBodyExtinctionModelPrediction.csv"))
pred.er[pred.er$spp==0,]$out <- NA
pred.er2[pred.er2$spp==0,]$out <- NA

load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Beaks.rda"))
a3Rxy<-archCoords[[3]]
ydiff<-xdiff<-min(diff(sort(unique(pred.er$PC2_beak))))
xboundary<-xdiff*4
pred.er<-pred.er[which(pred.er$PC2_beak>=(min(a3Rxy[,1])-xboundary) & pred.er$PC2_beak<(xboundary+max(a3Rxy[,1])) &
                         pred.er$PC3_beak>=(min(a3Rxy[,2])-xboundary) & pred.er$PC3_beak<(max(a3Rxy[,2])+xboundary)),]

load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Body.rda"))
a3Rxy<-archCoords[[3]]
ydiff<-xdiff<-min(diff(sort(unique(pred.er2$PC2_body))))
xboundary<-xdiff*4
pred.er2<-pred.er2[which(pred.er2$PC2_body>=(min(a3Rxy[,1])-xboundary) & pred.er2$PC2_body<(xboundary+max(a3Rxy[,1])) &
                         pred.er2$PC3_body>=(min(a3Rxy[,2])-xboundary) & pred.er2$PC3_body<(max(a3Rxy[,2])+xboundary)),]

brks<-quantile(c(pred.er$out,pred.er2$out),probs=seq(0,1,0.1),na.rm=T)
summary(pred.er$out)
summary(pred.er2$out)
plot.predictedGAM(pred.er,"beak",brks,emptyCells=T,gridWidth = 0.0001,col.pal = cols,letter="C")
mtext("Extinction risk",side=3,outer=F,line=0,cex=0.8)
plot.predictedGAM(pred.er2,"body",brks,emptyCells=T,gridWidth = 0.0001,col.pal = cols,letter="D")

par(mar=c(0,0,0,0))
brks2 <-  seq(min(brks),max(brks),length.out=length(brks))
if(length(cols)==1){
  cols3 <- brewer.pal(9, cols)
}else{
  cols3 <-cols  
}
cols4 <- colorRampPalette(cols3)(length(brks2)-1)
plot(0,0,type="n",yaxt="n",xaxt="n",xlim=c(-0.025,max(brks2)),ylim=c(0,1),bty="n",ylab="",xlab="")
for(i in 1:length(cols4)){polygon(c(brks2[i],brks2[i+1],brks2[i+1],brks2[i]),c(0.4,0.4,0.7,0.7),col=cols4[i],border="white",lwd=0.5)}
for(i in seq(1,length(brks2),by=2)){text(brks2[i],0.5,round(brks[i],2),pos=1,cex=0.75)}

par(mar=c(0,1.2,1.2,0))

## DR
pred.er <- read.csv(paste0(outFolder,"GAM_BeakModelPrediction_DR.csv"))
pred.er2 <- read.csv(paste0(outFolder,"GAM_BodyModelPrediction_DR.csv"))
pred.er[pred.er$spp==0,]$out <- NA
pred.er2[pred.er2$spp==0,]$out <- NA

load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Beaks.rda"))
a3Rxy<-archCoords[[3]]
ydiff<-xdiff<-min(diff(sort(unique(pred.er$PC2_beak))))
xboundary<-xdiff*4
pred.er<-pred.er[which(pred.er$PC2_beak>=(min(a3Rxy[,1])-xboundary) & pred.er$PC2_beak<(xboundary+max(a3Rxy[,1])) &
                         pred.er$PC3_beak>=(min(a3Rxy[,2])-xboundary) & pred.er$PC3_beak<(max(a3Rxy[,2])+xboundary)),]

load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Body.rda"))
a3Rxy<-archCoords[[3]]
ydiff<-xdiff<-min(diff(sort(unique(pred.er2$PC2_body))))
xboundary<-xdiff*4
pred.er2<-pred.er2[which(pred.er2$PC2_body>=(min(a3Rxy[,1])-xboundary) & pred.er2$PC2_body<(xboundary+max(a3Rxy[,1])) &
                           pred.er2$PC3_body>=(min(a3Rxy[,2])-xboundary) & pred.er2$PC3_body<(max(a3Rxy[,2])+xboundary)),]


brks<-quantile(c(pred.er$out,pred.er2$out),probs=seq(0,1,0.1),na.rm=T)
summary(pred.er$out)
summary(pred.er2$out)
plot.predictedGAM(pred.er,"beak",brks,emptyCells=T,gridWidth = 0.0001,col.pal = cols,letter="E")
mtext("Speciation rate",side=3,outer=F,line=0,cex=0.8)
plot.predictedGAM(pred.er2,"body",brks,emptyCells=T,gridWidth = 0.0001,col.pal = cols,letter="F")

par(mar=c(0,0,0,0))
brks2 <-  seq(min(brks),max(brks),length.out=length(brks))
if(length(cols)==1){
  cols3 <- brewer.pal(9, cols)
}else{
  cols3 <-cols  
}
cols4 <- colorRampPalette(cols3)(length(brks2)-1)
plot(0,0,type="n",yaxt="n",xaxt="n",xlim=c(-0.025,max(brks2)),ylim=c(0,1),bty="n",ylab="",xlab="")
for(i in 1:length(cols4)){polygon(c(brks2[i],brks2[i+1],brks2[i+1],brks2[i]),c(0.4,0.4,0.7,0.7),col=cols4[i],border="white",lwd=0.5)}
for(i in seq(1,length(brks2),by=2)){text(brks2[i],0.5,round(brks[i],2),pos=1,cex=0.75)}

dev.off()

## End script 6.3 ##