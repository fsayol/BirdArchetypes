## 3.2 ## Calculate enrichment data across species

## Calculate enrichment across morphospace ####

# We split morphospace in grid cells and calculate percentages of techniques / tasks

# Open beak shape dataset:
fdat.beak <-read.csv(paste0(DataFolder,"data_tmp/BeakPC23.csv"))

# Add feeding techniques dataset:
FeedingTechniques <- read.csv(paste0(DataFolder,"Data1_Foraging_Feeding.csv"))
names(FeedingTechniques)
fdat <- merge(fdat.beak,FeedingTechniques,by="species")
names(fdat)

# Add beak task dataset:
beakTaskEnrich <- read.csv(paste0(DataFolder,"Data2_BeakTaskEnrichment.csv"))
names(beakTaskEnrich)
fdat <- merge(fdat,beakTaskEnrich,by="species")
names(fdat)

## beak techniques enrichment ####

xcuts<-seq(min(fdat$PC2_beak),max(fdat$PC2_beak)+1.5,length.out=40)
ycuts<-seq(min(fdat$PC3_beak),max(fdat$PC3_beak),length.out=40)

xdiff<-(unique(diff(xcuts),6)/2)[1]
ydiff<-(unique(diff(ycuts),6)/2)[1]

fdat$xCell<-cut(fdat$PC2_beak,xcuts,labels=FALSE,include.lowest = TRUE)
fdat$yCell<-cut(fdat$PC3_beak,ycuts,labels=FALSE,include.lowest = TRUE)
fdat$Cell<-paste(fdat$xCell,fdat$yCell,sep="_")

cellTab<-data.frame(Cell=names(table(fdat$Cell)),Tearing=NA,Crushing=NA,Filtering=NA,
                    Sweeping=NA,Engulfing=NA,Grazing=NA,Hammering=NA,Plucking=NA,Probing=NA,Spearing=NA)

cellTab$Tearing<-as.numeric(sapply(split(fdat$Tearing,fdat$Cell),sum))
cellTab$Crushing<-as.numeric(sapply(split(fdat$Crushing,fdat$Cell),sum))
cellTab$Filtering<-as.numeric(sapply(split(fdat$Filtering,fdat$Cell),sum))
cellTab$Sweeping<-as.numeric(sapply(split(fdat$Sweeping,fdat$Cell),sum))
cellTab$Engulfing<-as.numeric(sapply(split(fdat$Engulfing,fdat$Cell),sum))
cellTab$Grazing<-as.numeric(sapply(split(fdat$Grazing,fdat$Cell),sum))
cellTab$Hammering<-as.numeric(sapply(split(fdat$Hammering,fdat$Cell),sum))
cellTab$Plucking<-as.numeric(sapply(split(fdat$Plucking,fdat$Cell),sum))
cellTab$Probing<-as.numeric(sapply(split(fdat$Probing,fdat$Cell),sum))
cellTab$Spearing<-as.numeric(sapply(split(fdat$Spearing,fdat$Cell),sum))

cellTab$x<-as.numeric(unlist(strsplit(cellTab$Cell,"_"))[seq(1,nrow(cellTab)*2,2)])
cellTab$y<-as.numeric(unlist(strsplit(cellTab$Cell,"_"))[seq(2,nrow(cellTab)*2,2)])

cellTab$xMid<-xcuts[cellTab$x]+xdiff
cellTab$yMid<-ycuts[cellTab$y]+ydiff

cellTab$TearingScaled<-cellTab$Tearing/sum(fdat$Tearing)
cellTab$CrushingScaled<-cellTab$Crushing/sum(fdat$Crushing)
cellTab$FilteringScaled<-cellTab$Filtering/sum(fdat$Filtering)
cellTab$SweepingScaled<-cellTab$Sweeping/sum(fdat$Sweeping)
cellTab$EngulfingScaled<-cellTab$Engulfing/sum(fdat$Engulfing)
cellTab$GrazingScaled<-cellTab$Grazing/sum(fdat$Grazing)
cellTab$HammeringScaled<-cellTab$Hammering/sum(fdat$Hammering)
cellTab$PluckingScaled<-cellTab$Plucking/sum(fdat$Plucking)
cellTab$ProbingScaled<-cellTab$Probing/sum(fdat$Probing)
cellTab$SpearingScaled<-cellTab$Spearing/sum(fdat$Spearing)

cellTab$TearingScaled2<-cellTab[,16]/rowSums(cellTab[,16:25])
cellTab$CrushingScaled2<-cellTab[,17]/rowSums(cellTab[,16:25])
cellTab$FilteringScaled2<-cellTab[,18]/rowSums(cellTab[,16:25])
cellTab$SweepingScaled2<-cellTab[,19]/rowSums(cellTab[,16:25])
cellTab$EngulfingScaled2<-cellTab[,20]/rowSums(cellTab[,16:25])
cellTab$GrazingScaled2<-cellTab[,21]/rowSums(cellTab[,16:25])
cellTab$HammeringScaled2<-cellTab[,22]/rowSums(cellTab[,16:25])
cellTab$PluckingScaled2<-cellTab[,23]/rowSums(cellTab[,16:25])
cellTab$ProbingScaled2<-cellTab[,24]/rowSums(cellTab[,16:25])
cellTab$SpearingScaled2<-cellTab[,25]/rowSums(cellTab[,16:25])

load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Beaks.rda"))
a3Rxy<-archCoords[[3]]

save(xcuts,ycuts,xdiff,ydiff,cellTab,a3Rxy,file=paste0(OutputFolder,"out3_BeakNicheEnrichment.rda"))

## beak task enrichment ####

xcuts<-seq(min(fdat$PC2_beak),max(fdat$PC2_beak)+1.5,length.out=40)
ycuts<-seq(min(fdat$PC3_beak),max(fdat$PC3_beak),length.out=40)

xdiff<-(unique(diff(xcuts),6)/2)[1]
ydiff<-(unique(diff(ycuts),6)/2)[1]

fdat$xCell<-cut(fdat$PC2_beak,xcuts,labels=FALSE,include.lowest = TRUE)
fdat$yCell<-cut(fdat$PC3_beak,ycuts,labels=FALSE,include.lowest = TRUE)
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

load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Beaks.rda"))
a3Rxy<-archCoords[[3]]

save(xcuts,ycuts,xdiff,ydiff,cellTab,ENGCol,REACol,CRUCol,a3Rxy,file=paste0(OutputFolder,"out3_BeakTaskEnrichment.rda"))

## body task enrichment ####

fdat.body <-read.csv(paste0(DataFolder,"data_tmp/BodyPC23.csv"))

# Add beak task dataset:
bodyTaskEnrich <- read.csv(paste0(DataFolder,"Data3_BodyTaskEnrichment.csv"))
names(bodyTaskEnrich)
fdat <- merge(fdat.body,bodyTaskEnrich,by="species")
names(fdat)

xcuts<-seq(min(fdat$PC2_body),max(fdat$PC2_body),length.out=40)
ycuts<-seq(min(fdat$PC3_body),max(fdat$PC3_body),length.out=40)

xdiff<-(unique(diff(xcuts),6)/2)[1]
ydiff<-(unique(diff(ycuts),6)/2)[1]

fdat$xCell<-cut(fdat$PC2_body,xcuts,labels=FALSE,include.lowest = TRUE)
fdat$yCell<-cut(fdat$PC3_body,ycuts,labels=FALSE,include.lowest = TRUE)
fdat$Cell<-paste(fdat$xCell,fdat$yCell,sep="_")

cellTab<-data.frame(Cell=names(table(fdat$Cell)),ftFLY=NA,ftWLK=NA,ftSWM=NA)

cellTab$ftFLY<-as.numeric(sapply(split(fdat$Fly_task,fdat$Cell),sum))
cellTab$ftWLK<-as.numeric(sapply(split(fdat$Walk_task,fdat$Cell),sum))
cellTab$ftSWM<-as.numeric(sapply(split(fdat$Swim_task,fdat$Cell),sum))

cellTab$x<-as.numeric(unlist(strsplit(cellTab$Cell,"_"))[seq(1,nrow(cellTab)*2,2)])
cellTab$y<-as.numeric(unlist(strsplit(cellTab$Cell,"_"))[seq(2,nrow(cellTab)*2,2)])

cellTab$xMid<-xcuts[cellTab$x]+xdiff
cellTab$yMid<-ycuts[cellTab$y]+ydiff

cellTab$ftFLYScaled<-cellTab$ftFLY/sum(fdat$Fly_task)
cellTab$ftWLKScaled<-cellTab$ftWLK/sum(fdat$Walk_task)
cellTab$ftSWMScaled<-cellTab$ftSWM/sum(fdat$Swim_task)

cellTab$ftFLYScaled2<-cellTab[,9]/rowSums(cellTab[,9:11])
cellTab$ftWLKScaled2<-cellTab[,10]/rowSums(cellTab[,9:11])
cellTab$ftSWMScaled2<-cellTab[,11]/rowSums(cellTab[,9:11])

# Specify colors
FLYCol <- brewer.pal(9, "PiYG")[6:9]
WLKCol <- rev(brewer.pal(9, "BrBG")[1:4])
SWMCol <- brewer.pal(9, "RdYlBu")[6:9]
#show_col(FLYCol)

cellTab$FLYCol<-FLYCol[cut(cellTab$ftFLYScaled2,round(seq(0.33,1,length.out=5),2),labels=FALSE,include.lowest = TRUE)]
cellTab$WLKCol<-WLKCol[cut(cellTab$ftWLKScaled2,round(seq(0.33,1,length.out=5),2),labels=FALSE,include.lowest = TRUE)]
cellTab$SWMCol<-SWMCol[cut(cellTab$ftSWMScaled2,round(seq(0.33,1,length.out=5),2),labels=FALSE,include.lowest = TRUE)]

cellTabCols<-cellTab[,15:17]
colCol<-apply(cellTab[,9:11],1,which.max)

load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Body.rda"))
a3Rxy<-archCoords[[3]]

save(xcuts,ycuts,xdiff,ydiff,cellTab,FLYCol,WLKCol,SWMCol,a3Rxy,file=paste0(OutputFolder,"out3_BodyTaskEnrichment.rda"))

## End script 3.2 ##