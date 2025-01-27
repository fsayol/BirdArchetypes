## 3.4 ## Figure 3: Enrichment plots (as distance to archetypes)

## Model enrichment as distance to archetype

pc.beak <- read.csv(paste0(DataFolder,"data_tmp/BeakPC23.csv"))
pc.body <- read.csv(paste0(DataFolder,"data_tmp/BodyPC23.csv"))
fdat <- merge(pc.beak,pc.body,by="species")

# Add beak task enrichment
beakTaskEnrich <- read.csv(paste0(DataFolder,"Data2_BeakTaskEnrichment.csv"))
fdat <- merge(fdat,beakTaskEnrich,by="species")

bodyTaskEnrich <- read.csv(paste0(DataFolder,"Data3_BodyTaskEnrichment.csv"))
fdat <- merge(fdat,bodyTaskEnrich,by="species")

dm <- read.csv(paste0(DataFolder,"data_tmp/BeakDistanceMetrics.csv"))
dm <-dm[,c("species","distAR1","distAR2","distAR3")]
names(dm) <-c("species",paste0("Beak_",names(dm)[2:4]))
fdat <- merge(fdat,dm,by="species")

dm <- read.csv(paste0(DataFolder,"data_tmp/BodyDistanceMetrics.csv"))
dm<-dm[,c("species","distAR1","distAR2","distAR3")]
names(dm)<-c("species",paste0("Body_",names(dm)[2:4]))
fdat <- merge(fdat,dm,by="species")

Archetype <- c("Crush","Reach","Engulf","Fly","Walk","Swim")

modelOutput<-data.frame(Archetype=Archetype,Intercept=NA,Slope=NA,Quadratic=NA)

nStarsList<-list()
nStarsList[[1]]<-""
nStarsList[[2]]<-"*"
nStarsList[[3]]<-"**"
nStarsList[[4]]<-"***"

plot.f<-paste0(FigureFolder,"Figure_3.pdf")
pdf(plot.f,w=3.135*1.7,h=3.135)

par(mfrow=c(2,3))
par(mar=c(2,2.25,1,1))

head(fdat)

for(x in 1:6){
  
  if(x==1){
    Col<-brewer.pal(9, "PiYG")[6:9]
    fdat$focVar<-round(fdat$Crush_task/100,1)
    fdat$focD<-fdat$Beak_distAR1
  }
  if(x==2){
    Col<-rev(brewer.pal(9, "BrBG")[1:4])
    fdat$focVar<-round(fdat$Reach_task/100,1)
    fdat$focD<-fdat$Beak_distAR2
  }
  if(x==3){
    Col<-brewer.pal(9, "RdYlBu")[6:9]
    fdat$focVar<-round(fdat$Engulf_task/100,1)
    fdat$focD<-fdat$Beak_distAR3
  }
  
  if(x==4){
    Col<-brewer.pal(9, "PiYG")[6:9]
    fdat$focVar<-round(fdat$Fly_task/100,1)
    fdat$focD<-fdat$Body_distAR3
  }
  if(x==5){
    Col<-rev(brewer.pal(9, "BrBG")[1:4])
    fdat$focVar<-round(fdat$Walk_task/100,1)
    fdat$focD<-fdat$Body_distAR2
  }
  if(x==6){
    Col<-brewer.pal(9, "RdYlBu")[6:9]
    fdat$focVar<-round(fdat$Swim_task/100,1)
    fdat$focD<-fdat$Body_distAR1
  }
 
  fdat$focD<-rescale(fdat$focD,c(0,1))
  vals<-sort(unique(fdat$focVar))
  Col2<-c(Col,rev(Col))
  
  plot(0,0,type="n",ylim=c(-0.05,1.05),xlim=c(0,1),ylab="",xlab="",yaxt="n",xaxt="n")
  points(fdat$focD,fdat$focVar,pch=".",col=Col2[1])
  
  axis(1,at=seq(0,1,0.2),labels=seq(0,1,0.2),cex.axis=0.75,padj=-2.5,tck=-0.0125,lwd=0.75)
  axis(2,at=seq(0,1,0.2),labels=seq(0,1,0.2),cex.axis=0.75,hadj=0.1,las=2,tck=-0.0125,lwd=0.75)
  if(x<4){mtext(paste0("Archetype ",x," distance"),side=1,outer=F,line=0.75,cex=0.75,at=0.5)}
  if(x==4){mtext(paste0("Archetype ",3," distance"),side=1,outer=F,line=0.75,cex=0.75,at=0.5)}
  if(x==5){mtext(paste0("Archetype ",2," distance"),side=1,outer=F,line=0.75,cex=0.75,at=0.5)}
  if(x==6){mtext(paste0("Archetype ",1," distance"),side=1,outer=F,line=0.75,cex=0.75,at=0.5)}
  
  mtext(Archetype[x],side=2,outer=F,line=1.25,cex=0.75,at=0.5)
  text(0.955,1,LETTERS[x],font=2,cex=1.25)
  
  for(i in 1:length(vals)){
    focRow<-which(fdat$focVar==vals[i])
    if(length(focRow)>=10){
      xvals<-quantile(fdat$focD[focRow],probs=seq(0.1,0.9,0.1))
      for(ii in 1:length(Col2)){
        polygon(c(xvals[ii],xvals[ii+1],xvals[ii+1],xvals[ii]),
                c(vals[i]-0.05,vals[i]-0.05,vals[i]+0.05,vals[i]+0.05),col=Col2[ii],border=FALSE)
      }
    }
  }
  fdat$focD<-rescale(fdat$focD,c(0.001,0.999))
  mod<-betareg(focD~focVar+I(focVar^2),data=fdat)
  modSum<-summary(mod)$coef$mean
  
  if(modSum[1,4]>0.05){nStars<-1}
  if(modSum[1,4]<0.05){nStars<-2}
  if(modSum[1,4]<0.01){nStars<-3}
  if(modSum[1,4]<0.001){nStars<-4}
  
  modelOutput$Intercept[x]<-paste0(round(modSum[1,1],2)," ± ",round(modSum[1,2],2),nStarsList[[nStars]],sep="")
  
  if(modSum[2,4]>0.05){nStars<-1}
  if(modSum[2,4]<0.05){nStars<-2}
  if(modSum[2,4]<0.01){nStars<-3}
  if(modSum[2,4]<0.001){nStars<-4}
  
  modelOutput$Slope[x]<-paste0(round(modSum[2,1],2)," ± ",round(modSum[2,2],2),nStarsList[[nStars]],sep="")
  
  if(modSum[3,4]>0.05){nStars<-1}
  if(modSum[3,4]<0.05){nStars<-2}
  if(modSum[3,4]<0.01){nStars<-3}
  if(modSum[3,4]<0.001){nStars<-4}
  
  modelOutput$Quadratic[x]<-paste0(round(modSum[3,1],2)," ± ",round(modSum[3,2],2),nStarsList[[nStars]],sep="")
  
  predDat<-data.frame(focVar=seq(0,1,0.01))
  lines(predict(mod, predDat),predDat$focVar,col="black",lwd=1)
}  

dev.off()

## End script 3.4 ##