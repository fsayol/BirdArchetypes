## 1.1 ## Null models of beak morphospace volume and shape

# We calculate volume of convex hull enclosing raw beak length, width and depth
# We compare observed data to null models

## 1a ## Observed data ####

# Open morphology dataset:
fdat <- read.csv(paste0(DataTmpFolder,"Avonet_BirdTree.csv"))

# Remove NA's from beak traits (n=9942 spp)
x <- (na.omit(fdat[,c("species","Beak.Length_Culmen","Beak.Width","Beak.Depth")]))

# Extract length, width and depth of the beak and log-transform.
rownames(x)<-x$species
x <- x[,-1]
names(x)<-c("Length","Width","Depth")
xLog <- log(x)

# empirical beak morphovolume
x2 <- as.matrix(xLog)
ch <- cxhull(x2[!duplicated(x2),])
EmpiricalVol<-ch[[5]]

# empirical beak volume
x2 <- as.matrix(x)
EmpiricalBeakVolume<-rowProds(x2)*1/3

x[which.max(x[,1]),]
x[which.max(x[,2]),]
x[which.max(x[,3]),]

x[which.min(EmpiricalBeakVolume),]
x[which.max(EmpiricalBeakVolume),]

EmpiricalBeakVolumeCone<-EmpiricalBeakVolume*pi
EmpiricalBeakVolumeDiff<-log10(max(EmpiricalBeakVolume))-log10(min(EmpiricalBeakVolume))
EmpiricalBeakVolumeDiffCone<-log10(max(EmpiricalBeakVolumeCone))-log10(min(EmpiricalBeakVolumeCone))
EmpiricalBeakVolumeDiff
EmpiricalBeakVolumeDiffCone

# Maximum possible beak volume
x22<-as.matrix(x)
x22[,1]<-sort(x22[,1])
x22[,2]<-sort(x22[,2])
x22[,3]<-sort(x22[,3])

EmpiricalBeakVolumeMax<-rowProds(x22)*1/3
EmpiricalBeakVolumeMaxCone<-EmpiricalBeakVolumeMax*pi

EmpiricalBeakVolumeMaxDiff<-log10(max(EmpiricalBeakVolumeMax))-log10(min(EmpiricalBeakVolumeMax))
EmpiricalBeakVolumeMaxDiffCone<-log10(max(EmpiricalBeakVolumeMaxCone))-log10(min(EmpiricalBeakVolumeMaxCone))

EmpiricalBeakVolumeDiff/EmpiricalBeakVolumeMaxDiff
EmpiricalBeakVolumeDiffCone/EmpiricalBeakVolumeMaxDiffCone

SimulatedTraits<-data.frame(
  CubeNM1.MV=rep(NA,100),CubeNM1.BMax=NA,CubeNM1.BMin=NA,CubeNM1.BDiff=NA,CubeNM1.BDiff.Cone=NA,
  SphereNM2.1.MV=rep(NA,100),SphereNM2.1.BMax=NA,SphereNM2.1.BMin=NA,SphereNM2.1.BDiff=NA, SphereNM2.1.BDiff.Cone=NA,
  SphereNM2.2.MV=rep(NA,100),SphereNM2.2.BMax=NA,SphereNM2.2.BMin=NA,SphereNM2.2.BDiff=NA,SphereNM2.2.BDiff.Cone=NA)

## 1b ## Null model 1: Cube ####

# cube - uniform trait distribution
for(i in 1:nrow(SimulatedTraits)){
  x3<-as.matrix(x)
  x3[,1]<-sample(seq(min(x[,1]),max(x[,1]),length=nrow(x)))
  x3[,2]<-sample(seq(min(x[,2]),max(x[,2]),length=nrow(x)))
  x3[,3]<-sample(seq(min(x[,3]),max(x[,3]),length=nrow(x)))
  
  SimulatedBeakVolume<-rowProds(x3)*1/3
  
  SimulatedTraits$CubeNM1.BMax[i]<-max(SimulatedBeakVolume)
  SimulatedTraits$CubeNM1.BMin[i]<-min(SimulatedBeakVolume)
  
  x3<-log(x3)
  ch3<-NULL
  try(ch3<-cxhull(x3[!duplicated(x3),]))
  if(is.null(ch3)==FALSE){
    SimulatedTraits$CubeNM1.MV[i]<-ch3[[5]]
  }
}
mean(SimulatedTraits$CubeNM1.MV,na.rm=T)
sd(SimulatedTraits$CubeNM1.MV,na.rm=T)

SimulatedTraits$CubeNM1.BDiff<-log10(SimulatedTraits$CubeNM1.BMax)-log10(SimulatedTraits$CubeNM1.BMin)
SimulatedTraits$CubeNM1.BDiff.Cone<-log10(SimulatedTraits$CubeNM1.BMax*pi)-log10(SimulatedTraits$CubeNM1.BMin*pi)

mean(SimulatedTraits$CubeNM1.BDiff)/EmpiricalBeakVolumeMaxDiff
mean(SimulatedTraits$CubeNM1.BDiff.Cone)/EmpiricalBeakVolumeMaxDiffCone


## 1c ## Null model 2: Sphere ####

# sphere - observed trait clustering along each axis

for(i in 1:nrow(SimulatedTraits)){
  x2<-as.matrix(x)
  x2[,1]<-sample(x2[,1])
  x2[,2]<-sample(x2[,2])
  x2[,3]<-sample(x2[,3])
  
  SimulatedBeakVolume<-rowProds(x2)*1/3
  
  SimulatedTraits$SphereNM2.1.BMax[i]<-max(SimulatedBeakVolume)
  SimulatedTraits$SphereNM2.1.BMin[i]<-min(SimulatedBeakVolume)
  
  x2<-log(x2)
  ch2<-NULL
  try(ch2<-cxhull(x2[!duplicated(x2),]))
  if(is.null(ch2)==FALSE){
    SimulatedTraits$SphereNM2.1.MV[i]<-ch2[[5]]
  }
} 
mean(SimulatedTraits$SphereNM2.1.MV,na.rm=T)
sd(SimulatedTraits$SphereNM2.1.MV,na.rm=T)

SimulatedTraits$SphereNM2.1.BDiff<-log10(SimulatedTraits$SphereNM2.1.BMax)-log10(SimulatedTraits$SphereNM2.1.BMin)
SimulatedTraits$SphereNM2.1.BDiff.Cone<-log10(SimulatedTraits$SphereNM2.1.BMax*pi)-log10(SimulatedTraits$SphereNM2.1.BMin*pi)

mean(SimulatedTraits$SphereNM2.1.BDiff)/EmpiricalBeakVolumeMaxDiff
mean(SimulatedTraits$SphereNM2.1.BDiff.Cone)/EmpiricalBeakVolumeMaxDiffCone

# sphere - uniform trait distribution
radLog<-0.5*mean(c(diff(range(xLog[,1])),diff(range(xLog[,2])),diff(range(xLog[,3]))))

D<-mean(c(diff(range(x[,1])),diff(range(x[,2])),diff(range(x[,3]))))
for(i in 1:nrow(SimulatedTraits)){
  x2uLog <- rsphere(n=nrow(x),r=radLog)
  
  x2u <- exp(x2uLog)
  
  SimulatedBeakVolume<-rowProds(x2u)*1/3
  
  SimulatedTraits$SphereNM2.2.BMax[i]<-max(SimulatedBeakVolume)
  SimulatedTraits$SphereNM2.2.BMin[i]<-min(SimulatedBeakVolume)
  
  ch2u<-NULL
  try(ch2u<-cxhull(x2uLog[!duplicated(x2uLog),]))
  if(is.null(ch2u)==FALSE){
    SimulatedTraits$SphereNM2.2.MV[i]<-ch2u[[5]]
  }
} 
mean(SimulatedTraits$SphereNM2.2.MV,na.rm=T)
sd(SimulatedTraits$SphereNM2.2.MV,na.rm=T)

SimulatedTraits$SphereNM2.2.BDiff<-log10(SimulatedTraits$SphereNM2.2.BMax)-log10(SimulatedTraits$SphereNM2.2.BMin)
SimulatedTraits$SphereNM2.2.BDiff.Cone<-log10(SimulatedTraits$SphereNM2.2.BMax*pi)-log10(SimulatedTraits$SphereNM2.2.BMin*pi)

mean(SimulatedTraits$SphereNM2.2.BDiff)/EmpiricalBeakVolumeMaxDiff
mean(SimulatedTraits$SphereNM2.2.BDiff.Cone)/EmpiricalBeakVolumeMaxDiffCone

## 1d ## Summarize and save data output ####

dat<-data.frame(Metric=names(colMeans(as.matrix(SimulatedTraits),na.rm=T)),
                Med=as.numeric(colMedians(as.matrix(SimulatedTraits),na.rm=T)),
                LQ=as.numeric(colQuantiles(as.matrix(SimulatedTraits),na.rm=T,prob=0.025)),
                UQ=as.numeric(colQuantiles(as.matrix(SimulatedTraits),na.rm=T,prob=0.975)))

save(dat,EmpiricalVol,EmpiricalBeakVolumeDiffCone,EmpiricalBeakVolumeMaxDiffCone,x,x2,x3,
     file=paste0(OutputFolder,"out1_NullModels_1.Rdata"))

## End script 1.1 ##