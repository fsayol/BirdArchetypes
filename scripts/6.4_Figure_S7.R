## 6.4 ## Plot Figure S7 (Insect gleaning and beak/body task)

outFolder <-paste0(OutputFolder,'out6_ER_DR_Models/')

plot.f<-paste0(FigureFolder,"Figure_S7.pdf")
pdf(plot.f,w=3.135*2.5,h=3.135*2)

cols<-c("gold","red", "darkred")

mat<-matrix(ncol=8,nrow=9)
mat[1:4,1:4]<-1
mat[5:8,1:4]<-2
mat[9,1:4]<-3

mat[1:4,5:8]<-4
mat[5:8,5:8]<-5
mat[9,5:8]<-6

layout(mat)

par(mar=c(0,1.2,1.2,0))
par(oma=c(0,0,0,0))

## Invertivore gleaners across beak
pred.iga <- read.csv(paste0(outFolder,"GAM_BeakModelPrediction_InvertGleanSpec.csv"))
pred.iga[pred.iga$spp==0,]$out <- NA
summary(pred.iga$out)

## Invertivore gleaners across body
pred.iga2 <- read.csv(paste0(outFolder,"GAM_BodyModelPrediction_InvertGleanSpec.csv"))
pred.iga2[pred.iga2$spp==0,]$out <- NA
summary(pred.iga2$out)

brks<-quantile(c(pred.iga$out,pred.iga2$out),probs=seq(0,1,0.1),na.rm=T)
brks<-seq(0,0.6,0.05)
plot.predictedGAM(pred.iga,"beak",brks,emptyCells=T,gridWidth = 0.0001,col.pal = cols,letter="A")
mtext("Arboreal gleaning invertivore",side=3,outer=F,line=0,cex=0.8)
mtext("Beak shape",side=2,outer=F,line=0,cex=0.8)
plot.predictedGAM(pred.iga2,"body",brks,emptyCells=T,gridWidth = 0.0001,col.pal = cols,letter="B")
mtext("Body shape",side=2,outer=F,line=0,cex=0.8)

par(mar=c(0,0,0,0))
brks2 <-  seq(min(brks),max(brks),length.out=length(brks))
if(length(cols)==1){
  cols3 <- brewer.pal(9, cols)
}else{
  cols3 <- cols  
}
cols4 <- colorRampPalette(cols3)(length(brks2)-1)
plot(0,0,type="n",yaxt="n",xaxt="n",xlim=c(-0.025,max(brks2)),ylim=c(0,1),bty="n",ylab="",xlab="")
for(i in 1:length(cols4)){polygon(c(brks2[i],brks2[i+1],brks2[i+1],brks2[i]),c(0.4,0.4,0.7,0.7),col=cols4[i],border="white",lwd=0.5)}
for(i in seq(1,length(brks2),by=2)){text(brks2[i],0.5,round(brks[i],3),pos=1,cex=0.75)}

par(mar=c(0,1.2,1.2,0))

## Beak specialization
pred.beakSpec <- read.csv(paste0(outFolder,"GAM_BeakModelPrediction_BeakTask.csv"))
pred.beakSpec[pred.beakSpec$spp==0,]$out <- NA
summary(pred.beakSpec$out)

## Body specialization
pred.bodySpec <- read.csv(paste0(outFolder,"GAM_BodyModelPrediction_BodyTask.csv"))
pred.bodySpec[pred.bodySpec$spp==0,]$out <- NA
summary(pred.bodySpec$out)

brks<-quantile(c(pred.beakSpec$out,pred.bodySpec$out),probs=seq(0,1,0.1),na.rm=T)
plot.predictedGAM(pred.beakSpec,"beak",brks,emptyCells=T,gridWidth = 0.0001,col.pal = cols,letter="C")
mtext("Task specialization",side=3,outer=F,line=0,cex=0.8)
plot.predictedGAM(pred.bodySpec,"body",brks,emptyCells=T,gridWidth = 0.0001,col.pal = cols,letter="D")

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

## End script 6.4 ##