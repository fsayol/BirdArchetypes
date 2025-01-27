## 2.4 ## Null models for archetypes and fit regular polygons analysis

## We run two types of null models, for both beak and body shape:
# (1) Null models randomizing rows within traits (RND null model)
# (2) Null models using brownian motion simulations (BM null model)

## We first create folders to store outputs

# Save randomized null models (1)
outFolder <-paste0(OutputFolder,'out2_RegPol_Null_Random/')
if(!file.exists(outFolder)) dir.create(outFolder)
# Save BM null models (2)
outFolderBM <-paste0(OutputFolder,'out2_RegPol_Null_BM/')
if(!file.exists(outFolderBM)) dir.create(outFolderBM)

## (1 beak) Null models regular polygons with shifts for Beak ####

## (1a) Reshuffle beak traits, (1b) redo PCA, (1c) run archetype analysis

# Open morphology dataset:
fdat <- read.csv(paste0(DataTmpFolder,"Avonet_BirdTree.csv"))
# Select beak traits:
pc.dat23<-read.csv(paste0(DataFolder,"data_tmp/BeakPC23.csv"))
fdat <- fdat[fdat$species %in% pc.dat23$species,]
x <- log(fdat[,c("Beak.Length_Culmen","Beak.Width","Beak.Depth")]) # # Log transform
x <- scale(x)
rownames(x) <- fdat$species

Nspp <- nrow(fdat)

# (1) Reshufle beak traits 
Nnull <- 100 

for (Ni in 1:Nnull) {
  
  # Step 1: Reshuffle values
  x[,1] <- x[sample(1:Nspp,replace=F),1] # Reshuffle Beak length
  x[,2] <- x[sample(1:Nspp,replace=F),2] # Reshuffle Beak width
  x[,3] <- x[sample(1:Nspp,replace=F),3] # Reshuffle Beak depth
  
  # Step 2: Run PCA
  pca <- prcomp(x)
  pcx <- data.frame(pca$x) #Chose which data to use.
  pcx$species <- rownames(pca$x)
  names(pcx)[1:3] <- c("PC1_beak","PC2_beak","PC3_beak")
  # Scale PCs
  pcx$PC1_beak <- scale(pcx$PC1_beak)
  pcx$PC2_beak <- scale(pcx$PC2_beak)
  pcx$PC3_beak <- scale(pcx$PC3_beak)
  fdat.pc <- merge(fdat,pcx[,c(1:3,ncol(pcx))],by="species")
  pc.dat23 <- fdat.pc[,c("species","PC2_beak","PC3_beak")]
  
  # Step 3: Run archetype analysis

  pc12Foc <- as.matrix(data.frame(PC1=pc.dat23[,2],PC2=pc.dat23[,3]))
  H12Foc <- Hpi(x=pc12Foc)      # optimal bandwidth estimation
  est12Foc<- kde(x=pc12Foc, H=H12Foc, compute.cont=TRUE)     # kernel density estimation
  
  pc12<-pc12Foc
  H12<-H12Foc
  est12<-est12Foc
  
  dat<-expand.grid(x=est12$eval.points[[1]],y=est12$eval.points[[2]])
  
  dat$x2<-rep(1:length(est12$eval.points[[1]]),length(est12$eval.points[[2]]))
  dat$y2<-rep(1:length(est12$eval.points[[1]]),each=length(est12$eval.points[[2]]))
  
  dat$estimate<-as.vector(est12$estimate)
  dat<-dat[which(dat$estimate>est12$cont[which(names(est12$cont)=="10%")]),]
  
  ydiff<-xdiff<-1
  dat$nb<-NA
  for(i in 1:nrow(dat)){
    left<-which(round(dat$x2)==round(dat$x2[i]-xdiff) & round(dat$y2)==round(dat$y2[i]))
    right<-which(round(dat$x2)==round(dat$x2[i]+xdiff) & round(dat$y2)==round(dat$y2[i]))
    
    top<-which(round(dat$y2)==round(dat$y2[i]+ydiff) & round(dat$x2)==round(dat$x2[i]))
    low<-which(round(dat$y2)==round(dat$y2[i]-ydiff) & round(dat$x2)==round(dat$x2[i]))
    
    topleft<-which(round(dat$x2)==round(dat$x2[i]-xdiff) & round(dat$y2)==round(dat$y2[i]+ydiff))
    topright<-which(round(dat$x2)==round(dat$x2[i]+xdiff) & round(dat$y2)==round(dat$y2[i]+ydiff))
    
    lowleft<-which(round(dat$x2)==round(dat$x2[i]-xdiff) & round(dat$y2)==round(dat$y2[i]-ydiff))
    lowright<-which(round(dat$x2)==round(dat$x2[i]+xdiff) & round(dat$y2)==round(dat$y2[i]-ydiff))
    
    dat$nb[i]<-length(c(left,right,top,low,topleft,topright,lowleft,lowright))
  }
  table(dat$nb)
  dat <- dat[dat$nb<8,]
  datList <- dat
  
  MaxArea<-ceiling(diff(range(pc.dat23[,2]))*diff(range(pc.dat23[,3])))/2
  nVertex<-c(3,4,5,6,7,360)
  InputDat<-data.frame()
  for(i in 1:length(nVertex)){
    dat <- expand.grid(nVertex=nVertex[i],Area=seq(1,MaxArea,1),xpos=0,ypos=0,Rotation=seq(0,360/nVertex[i],1),Position=0,borderDist=NA)
    InputDat<-rbind(InputDat,dat)
  }  
  
  for(i in 1:nrow(InputDat)){
    
    shape<-regular.poly(nSides=InputDat$nVertex[i], area=InputDat$Area[i])
    
    shape<-Rotate(x=shape$x, y=shape$y,my=NULL,theta=deg2rad(InputDat$Rotation[i]))[1:2]
    
    shape$x<-shape$x+InputDat$xpos[i]
    shape$y<-shape$y+InputDat$ypos[i]
    
    shapeCoord<-matrix(ncol=2,nrow=InputDat$nVertex[i])
    shapeCoord[,1]<-shape$x
    shapeCoord[,2]<-shape$y
    shapeCoord<-rbind(shapeCoord,shapeCoord[1,])
    shapeCoord[,1]<-shapeCoord[,1]+rnorm(nrow(shapeCoord),0,0.0001)
    shapeCoord[,2]<-shapeCoord[,2]+rnorm(nrow(shapeCoord),0,0.0001)
    shapeCoordListx<-list()
    shapeCoordListy<-list()
    for(ii in 1:(nrow(shapeCoord)-1)){
      xy<-approx(x=c(shapeCoord[ii,1],shapeCoord[ii+1,1]), y=c(shapeCoord[ii,2],shapeCoord[ii+1,2]),method="linear", n=50)
      shapeCoordListx[[ii]]<-xy$x
      shapeCoordListy[[ii]]<-xy$y
    }
    shapeCoord<-data.frame(x=unlist(shapeCoordListx),y=unlist(shapeCoordListy))
    shapeCoord<-unique(shapeCoord)
    
    Morph <- datList[,1:2]
    MS<-rbind(Morph,shapeCoord)
    
    borderDist<-pointDistance(Morph,shapeCoord,F)
    borderDist2<-rowMins(as.matrix(borderDist))
    
    InputDat$borderDist[i]<-sum(borderDist2)
  }
  
  borderDist<-matrix(ncol=length(nVertex),nrow=length(seq(1,MaxArea,1)))
  OptShape<-data.frame()
  for(i in 1:length(nVertex)){
    toplot<-which(InputDat$nVertex==nVertex[i])
    borderDist[,i]<-sapply(split(InputDat$borderDist[toplot],InputDat$Area[toplot]),min)
    mD<-min(borderDist[,i])
    OptShape<-rbind(OptShape,InputDat[which(InputDat$nVertex==nVertex[i] & InputDat$borderDist==mD),])
  }  
  print(OptShape[,c(1,7)])
  save(pc.dat23,OptShape,file=paste0(outFolder,"RegularPolygonFit_Null_",Ni,"_Beak.rda"))
  cat(Ni,"/",Nnull,"\n")
} # End for Ni (Null model i)


## (1 body) Null models regular polygons with shifts for Body ####
## (1a) Reshuffle body traits, (1b) redo PCA, (1c) run archetype analysis

# Open morphology dataset:
fdat <- read.csv(paste0(DataTmpFolder,"Avonet_BirdTree.csv"))
# Select body traits:
pc.dat23<-read.csv(paste0(DataFolder,"data_tmp/BodyPC23.csv"))
fdat <- fdat[fdat$species %in% pc.dat23$species,]
x <- log(fdat[,c("Tarsus.Length","Wing.Length","Kipps.Distance","Tail.Length")]) # # Log transform
x <- scale(x)
rownames(x) <- fdat$species

Nspp <- nrow(fdat)

# Null models to perform  
Nnull <- 100 

for (Ni in 1:15) {
  
  # Step 1: Reshuffle values
  x[,1] <- x[sample(1:Nspp,replace=F),1] # Reshuffle tarsus
  x[,2] <- x[sample(1:Nspp,replace=F),2] # Reshuffle wing
  x[,3] <- x[sample(1:Nspp,replace=F),3] # Reshuffle kipps
  x[,4] <- x[sample(1:Nspp,replace=F),4] # Reshuffle tail
  # Step 2: Run PCA
  pca <- prcomp(x)
  pcx <- data.frame(pca$x) #Chose which data to use.
  pcx$species <- rownames(pca$x)
  names(pcx)[1:3] <- c("PC1_body","PC2_body","PC3_body")
  # Scale PCs
  pcx$PC1_body <- scale(pcx$PC1_body)
  pcx$PC2_body <- scale(pcx$PC2_body)
  pcx$PC3_body <- scale(pcx$PC3_body)
  fdat.pc <- merge(fdat,pcx[,c(1:3,ncol(pcx))],by="species")
  pc.dat23 <- fdat.pc[,c("species","PC2_body","PC3_body")]
  
  # Step 3: Run archetype analysis
  
  pc12Foc <- as.matrix(data.frame(PC1=pc.dat23[,2],PC2=pc.dat23[,3]))
  H12Foc <- Hpi(x=pc12Foc)      # optimal bandwidth estimation
  est12Foc<- kde(x=pc12Foc, H=H12Foc, compute.cont=TRUE)     # kernel density estimation
  
  pc12<-pc12Foc
  H12<-H12Foc
  est12<-est12Foc
  
  dat<-expand.grid(x=est12$eval.points[[1]],y=est12$eval.points[[2]])
  
  dat$x2<-rep(1:length(est12$eval.points[[1]]),length(est12$eval.points[[2]]))
  dat$y2<-rep(1:length(est12$eval.points[[1]]),each=length(est12$eval.points[[2]]))
  
  dat$estimate<-as.vector(est12$estimate)
  dat<-dat[which(dat$estimate>est12$cont[which(names(est12$cont)=="10%")]),]
  
  ydiff<-xdiff<-1
  dat$nb<-NA
  for(i in 1:nrow(dat)){
    left<-which(round(dat$x2)==round(dat$x2[i]-xdiff) & round(dat$y2)==round(dat$y2[i]))
    right<-which(round(dat$x2)==round(dat$x2[i]+xdiff) & round(dat$y2)==round(dat$y2[i]))
    
    top<-which(round(dat$y2)==round(dat$y2[i]+ydiff) & round(dat$x2)==round(dat$x2[i]))
    low<-which(round(dat$y2)==round(dat$y2[i]-ydiff) & round(dat$x2)==round(dat$x2[i]))
    
    topleft<-which(round(dat$x2)==round(dat$x2[i]-xdiff) & round(dat$y2)==round(dat$y2[i]+ydiff))
    topright<-which(round(dat$x2)==round(dat$x2[i]+xdiff) & round(dat$y2)==round(dat$y2[i]+ydiff))
    
    lowleft<-which(round(dat$x2)==round(dat$x2[i]-xdiff) & round(dat$y2)==round(dat$y2[i]-ydiff))
    lowright<-which(round(dat$x2)==round(dat$x2[i]+xdiff) & round(dat$y2)==round(dat$y2[i]-ydiff))
    
    dat$nb[i]<-length(c(left,right,top,low,topleft,topright,lowleft,lowright))
  }
  table(dat$nb)
  dat <- dat[dat$nb<8,]
  datList <- dat
  
  MaxArea<-ceiling(diff(range(pc.dat23[,2]))*diff(range(pc.dat23[,3])))/2
  nVertex<-c(3,4,5,6,7,360)
  InputDat<-data.frame()
  for(i in 1:length(nVertex)){
    dat <- expand.grid(nVertex=nVertex[i],Area=seq(1,MaxArea,1),xpos=0,ypos=0,Rotation=seq(0,360/nVertex[i],1),Position=0,borderDist=NA)
    InputDat<-rbind(InputDat,dat)
  }  
  
  for(i in 1:nrow(InputDat)){
    
    shape<-regular.poly(nSides=InputDat$nVertex[i], area=InputDat$Area[i])
    
    shape<-Rotate(x=shape$x, y=shape$y,my=NULL,theta=deg2rad(InputDat$Rotation[i]))[1:2]
    
    shape$x<-shape$x+InputDat$xpos[i]
    shape$y<-shape$y+InputDat$ypos[i]
    
    shapeCoord<-matrix(ncol=2,nrow=InputDat$nVertex[i])
    shapeCoord[,1]<-shape$x
    shapeCoord[,2]<-shape$y
    shapeCoord<-rbind(shapeCoord,shapeCoord[1,])
    shapeCoord[,1]<-shapeCoord[,1]+rnorm(nrow(shapeCoord),0,0.0001)
    shapeCoord[,2]<-shapeCoord[,2]+rnorm(nrow(shapeCoord),0,0.0001)
    shapeCoordListx<-list()
    shapeCoordListy<-list()
    for(ii in 1:(nrow(shapeCoord)-1)){
      xy<-approx(x=c(shapeCoord[ii,1],shapeCoord[ii+1,1]), y=c(shapeCoord[ii,2],shapeCoord[ii+1,2]),method="linear", n=50)
      shapeCoordListx[[ii]]<-xy$x
      shapeCoordListy[[ii]]<-xy$y
    }
    shapeCoord<-data.frame(x=unlist(shapeCoordListx),y=unlist(shapeCoordListy))
    shapeCoord<-unique(shapeCoord)
    
    Morph <- datList[,1:2]
    MS<-rbind(Morph,shapeCoord)
    
    borderDist<-pointDistance(Morph,shapeCoord,F)
    borderDist2<-rowMins(as.matrix(borderDist))
    InputDat$borderDist[i]<-sum(borderDist2)
  }
  
  borderDist<-matrix(ncol=length(nVertex),nrow=length(seq(1,MaxArea,1)))
  OptShape<-data.frame()
  for(i in 1:length(nVertex)){
    toplot<-which(InputDat$nVertex==nVertex[i])
    borderDist[,i]<-sapply(split(InputDat$borderDist[toplot],InputDat$Area[toplot]),min)
    mD<-min(borderDist[,i])
    OptShape<-rbind(OptShape,InputDat[which(InputDat$nVertex==nVertex[i] & InputDat$borderDist==mD),])
  }  
  print(OptShape[,c(1,7)])
  save(pc.dat23,OptShape,file=paste0(outFolder,"RegularPolygonFit_Null_",Ni,"_Body.rda"))

  cat(Ni,"/",Nnull,"\n")
  
} # End for Ni (Null model i)

## End script.
## (2 beak) Null models regular polygons with simulated BM for Beak ####
## (2a) Simulate beak traits (BM), (2b) redo PCA, (2c) run archetype analysis

# Open morphology dataset:
fdat <- read.csv(paste0(DataTmpFolder,"Avonet_BirdTree.csv"))
# Select beak traits:
fdat <- na.omit(fdat[,c("species","Beak.Length_Culmen","Beak.Width","Beak.Depth")])

tree <- read.nexus(paste0(DataFolder,"phylo/AllBirdsHackett1_summary.tre"))
tree <- drop.tip(tree,tree$tip.label[!tree$tip.label %in% fdat$species])
Nspp <- nrow(fdat)
length(tree$tip.label)

x <- log(fdat[,c("Beak.Length_Culmen","Beak.Width","Beak.Depth")]) # # Log transform
x <- scale(x)
rownames(x) <- fdat$species

# Simulate 100 trees using cvc
vcv.x <- ratematrix(tree,x)
sim.x <- sim.char(tree,vcv.x,100,"BM")

# Null models to perform  
Nnull <- 100 

for (Ni in 1:Nnull) {
  
  # Select simulation Ni
  x <- sim.x[1:9912,1:3,Ni]

  # Step 2: Run PCA
  pca <- prcomp(x)
  pcx <- data.frame(pca$x) #Chose which data to use.
  pcx$species <- rownames(pca$x)
  names(pcx)[1:3] <- c("PC1_beak","PC2_beak","PC3_beak")
  # Scale PCs
  pcx$PC1_beak <- scale(pcx$PC1_beak)
  pcx$PC2_beak <- scale(pcx$PC2_beak)
  pcx$PC3_beak <- scale(pcx$PC3_beak)
  fdat.pc <- merge(fdat,pcx[,c(1:3,ncol(pcx))],by="species")
  pc.dat23 <- fdat.pc[,c("species","PC2_beak","PC3_beak")]
  
  # Step 3: Run archetype analysis
  pc12Foc <- as.matrix(data.frame(PC1=pc.dat23[,2],PC2=pc.dat23[,3]))
  H12Foc <- Hpi(x=pc12Foc)      # optimal bandwidth estimation
  est12Foc<- kde(x=pc12Foc, H=H12Foc, compute.cont=TRUE)     # kernel density estimation
  pc12<-pc12Foc
  H12<-H12Foc
  est12<-est12Foc
  dat<-expand.grid(x=est12$eval.points[[1]],y=est12$eval.points[[2]])
  
  dat$x2<-rep(1:length(est12$eval.points[[1]]),length(est12$eval.points[[2]]))
  dat$y2<-rep(1:length(est12$eval.points[[1]]),each=length(est12$eval.points[[2]]))
  
  dat$estimate<-as.vector(est12$estimate)
  dat<-dat[which(dat$estimate>est12$cont[which(names(est12$cont)=="10%")]),]
  
  ydiff<-xdiff<-1
  dat$nb<-NA
  for(i in 1:nrow(dat)){
    left<-which(round(dat$x2)==round(dat$x2[i]-xdiff) & round(dat$y2)==round(dat$y2[i]))
    right<-which(round(dat$x2)==round(dat$x2[i]+xdiff) & round(dat$y2)==round(dat$y2[i]))
    
    top<-which(round(dat$y2)==round(dat$y2[i]+ydiff) & round(dat$x2)==round(dat$x2[i]))
    low<-which(round(dat$y2)==round(dat$y2[i]-ydiff) & round(dat$x2)==round(dat$x2[i]))
    
    topleft<-which(round(dat$x2)==round(dat$x2[i]-xdiff) & round(dat$y2)==round(dat$y2[i]+ydiff))
    topright<-which(round(dat$x2)==round(dat$x2[i]+xdiff) & round(dat$y2)==round(dat$y2[i]+ydiff))
    
    lowleft<-which(round(dat$x2)==round(dat$x2[i]-xdiff) & round(dat$y2)==round(dat$y2[i]-ydiff))
    lowright<-which(round(dat$x2)==round(dat$x2[i]+xdiff) & round(dat$y2)==round(dat$y2[i]-ydiff))
    
    dat$nb[i]<-length(c(left,right,top,low,topleft,topright,lowleft,lowright))
  }
  table(dat$nb)
  dat <- dat[dat$nb<8,]
  datList <- dat
  
  MaxArea<-ceiling(diff(range(pc.dat23[,2]))*diff(range(pc.dat23[,3])))/2
  nVertex<-c(3,4,5,6,7,360)
  InputDat<-data.frame()
  for(i in 1:length(nVertex)){
    #dat<-expand.grid(nVertex=nVertex[i],Area=seq(1,MaxArea,1),xpos=seq(-0.5,0.5,0.25),ypos=seq(-0.5,0.5,0.25),Rotation=seq(0,360/nVertex[i],1),Position=0,borderDist=NA)
    dat <- expand.grid(nVertex=nVertex[i],Area=seq(1,MaxArea,1),xpos=0,ypos=0,Rotation=seq(0,360/nVertex[i],1),Position=0,borderDist=NA)
    InputDat<-rbind(InputDat,dat)
  }  
  
  for(i in 1:nrow(InputDat)){
    
    shape<-regular.poly(nSides=InputDat$nVertex[i], area=InputDat$Area[i])
    
    shape<-Rotate(x=shape$x, y=shape$y,my=NULL,theta=deg2rad(InputDat$Rotation[i]))[1:2]
    
    shape$x<-shape$x+InputDat$xpos[i]
    shape$y<-shape$y+InputDat$ypos[i]
    
    shapeCoord<-matrix(ncol=2,nrow=InputDat$nVertex[i])
    shapeCoord[,1]<-shape$x
    shapeCoord[,2]<-shape$y
    shapeCoord<-rbind(shapeCoord,shapeCoord[1,])
    shapeCoord[,1]<-shapeCoord[,1]+rnorm(nrow(shapeCoord),0,0.0001)
    shapeCoord[,2]<-shapeCoord[,2]+rnorm(nrow(shapeCoord),0,0.0001)
    shapeCoordListx<-list()
    shapeCoordListy<-list()
    for(ii in 1:(nrow(shapeCoord)-1)){
      xy<-approx(x=c(shapeCoord[ii,1],shapeCoord[ii+1,1]), y=c(shapeCoord[ii,2],shapeCoord[ii+1,2]),method="linear", n=50)
      shapeCoordListx[[ii]]<-xy$x
      shapeCoordListy[[ii]]<-xy$y
    }
    shapeCoord<-data.frame(x=unlist(shapeCoordListx),y=unlist(shapeCoordListy))
    shapeCoord<-unique(shapeCoord)
    
    Morph <- datList[,1:2]
    MS<-rbind(Morph,shapeCoord)
    
    borderDist<-pointDistance(Morph,shapeCoord,F)
    borderDist2<-rowMins(as.matrix(borderDist))
    
    InputDat$borderDist[i]<-sum(borderDist2)
    # print(i)
  }
  
  borderDist<-matrix(ncol=length(nVertex),nrow=length(seq(1,MaxArea,1)))
  OptShape<-data.frame()
  for(i in 1:length(nVertex)){
    toplot<-which(InputDat$nVertex==nVertex[i])
    borderDist[,i]<-sapply(split(InputDat$borderDist[toplot],InputDat$Area[toplot]),min)
    mD<-min(borderDist[,i])
    
    #print(InputDat[which(InputDat$nVertex==nVertex[i] & InputDat$borderDist==mD),])
    OptShape<-rbind(OptShape,InputDat[which(InputDat$nVertex==nVertex[i] & InputDat$borderDist==mD),])
    #print(i)
  }  
  print(OptShape[,c(1,7)])
  save(pc.dat23,OptShape,file=paste0(outFolderBM,"RegPolygonFit_Null_BM_",Ni,"_Beak.rda"))
  
  cat(Ni,"/",Nnull,"\n")
  
} # End for Ni (Null model i)


## (2 body) Null models regular polygons with simulated BM for Body ####
## (2a) Simulate body traits (BM), (2b) redo PCA, (2c) run archetype analysis

# Open morphology dataset:
fdat <- read.csv(paste0(DataTmpFolder,"Avonet_BirdTree.csv"))
# Select body traits:
fdat <- na.omit(fdat[,c("species","Tarsus.Length","Wing.Length","Kipps.Distance","Tail.Length")])
fdat <- fdat[-which(fdat$Wing.Length==0.1),]

# Open tree:
tree <- read.nexus(paste0(DataFolder,"phylo/AllBirdsHackett1_summary.tre"))
tree <- drop.tip(tree,tree$tip.label[!tree$tip.label %in% fdat$species])
Nspp <- nrow(fdat)
length(tree$tip.label)

x <- log(fdat[,c("Tarsus.Length","Wing.Length","Kipps.Distance","Tail.Length")]) # # Log transform
x <- scale(x)
rownames(x) <- fdat$species

# Simulate 100 trees using vcv
vcv.x <- ratematrix(tree,x)
sim.x <- sim.char(tree,vcv.x,100,"BM")

# Null models to perform  
Nnull <- 100 

for (Ni in 1:Nnull) {
  
  # Step 1: Select simulation Ni
  x <- sim.x[1:Nspp,1:3,Ni]

  # Step 2: Run PCA
  pca <- prcomp(x)
  pcx <- data.frame(pca$x) #Chose which data to use.
  pcx$species <- rownames(pca$x)
  names(pcx)[1:3] <- c("PC1_body","PC2_body","PC3_body")
  # Scale PCs
  pcx$PC1_body <- scale(pcx$PC1_body)
  pcx$PC2_body <- scale(pcx$PC2_body)
  pcx$PC3_body <- scale(pcx$PC3_body)
  fdat.pc <- merge(fdat,pcx[,c(1:3,ncol(pcx))],by="species")
  pc.dat23 <- fdat.pc[,c("species","PC2_body","PC3_body")]
  
  # Step 3: Run archetype analysis
  
  pc12Foc <- as.matrix(data.frame(PC1=pc.dat23[,2],PC2=pc.dat23[,3]))
  H12Foc <- Hpi(x=pc12Foc)      # optimal bandwidth estimation
  est12Foc<- kde(x=pc12Foc, H=H12Foc, compute.cont=TRUE)     # kernel density estimation
  
  pc12<-pc12Foc
  H12<-H12Foc
  est12<-est12Foc
  
  dat<-expand.grid(x=est12$eval.points[[1]],y=est12$eval.points[[2]])
  
  dat$x2<-rep(1:length(est12$eval.points[[1]]),length(est12$eval.points[[2]]))
  dat$y2<-rep(1:length(est12$eval.points[[1]]),each=length(est12$eval.points[[2]]))
  
  dat$estimate<-as.vector(est12$estimate)
  dat<-dat[which(dat$estimate>est12$cont[which(names(est12$cont)=="10%")]),]
  
  ydiff<-xdiff<-1
  dat$nb<-NA
  for(i in 1:nrow(dat)){
    left<-which(round(dat$x2)==round(dat$x2[i]-xdiff) & round(dat$y2)==round(dat$y2[i]))
    right<-which(round(dat$x2)==round(dat$x2[i]+xdiff) & round(dat$y2)==round(dat$y2[i]))
    
    top<-which(round(dat$y2)==round(dat$y2[i]+ydiff) & round(dat$x2)==round(dat$x2[i]))
    low<-which(round(dat$y2)==round(dat$y2[i]-ydiff) & round(dat$x2)==round(dat$x2[i]))
    
    topleft<-which(round(dat$x2)==round(dat$x2[i]-xdiff) & round(dat$y2)==round(dat$y2[i]+ydiff))
    topright<-which(round(dat$x2)==round(dat$x2[i]+xdiff) & round(dat$y2)==round(dat$y2[i]+ydiff))
    
    lowleft<-which(round(dat$x2)==round(dat$x2[i]-xdiff) & round(dat$y2)==round(dat$y2[i]-ydiff))
    lowright<-which(round(dat$x2)==round(dat$x2[i]+xdiff) & round(dat$y2)==round(dat$y2[i]-ydiff))
    
    dat$nb[i]<-length(c(left,right,top,low,topleft,topright,lowleft,lowright))
  }
  table(dat$nb)
  dat <- dat[dat$nb<8,]
  datList <- dat
  
  MaxArea<-ceiling(diff(range(pc.dat23[,2]))*diff(range(pc.dat23[,3])))/2
  nVertex<-c(3,4,5,6,7,360)
  InputDat<-data.frame()
  for(i in 1:length(nVertex)){
    #dat<-expand.grid(nVertex=nVertex[i],Area=seq(1,MaxArea,1),xpos=seq(-0.5,0.5,0.25),ypos=seq(-0.5,0.5,0.25),Rotation=seq(0,360/nVertex[i],1),Position=0,borderDist=NA)
    dat <- expand.grid(nVertex=nVertex[i],Area=seq(1,MaxArea,1),xpos=0,ypos=0,Rotation=seq(0,360/nVertex[i],1),Position=0,borderDist=NA)
    InputDat<-rbind(InputDat,dat)
  }  
  
  for(i in 1:nrow(InputDat)){
    
    shape<-regular.poly(nSides=InputDat$nVertex[i], area=InputDat$Area[i])
    
    shape<-Rotate(x=shape$x, y=shape$y,my=NULL,theta=deg2rad(InputDat$Rotation[i]))[1:2]
    
    shape$x<-shape$x+InputDat$xpos[i]
    shape$y<-shape$y+InputDat$ypos[i]
    
    shapeCoord<-matrix(ncol=2,nrow=InputDat$nVertex[i])
    shapeCoord[,1]<-shape$x
    shapeCoord[,2]<-shape$y
    shapeCoord<-rbind(shapeCoord,shapeCoord[1,])
    shapeCoord[,1]<-shapeCoord[,1]+rnorm(nrow(shapeCoord),0,0.0001)
    shapeCoord[,2]<-shapeCoord[,2]+rnorm(nrow(shapeCoord),0,0.0001)
    shapeCoordListx<-list()
    shapeCoordListy<-list()
    for(ii in 1:(nrow(shapeCoord)-1)){
      xy<-approx(x=c(shapeCoord[ii,1],shapeCoord[ii+1,1]), y=c(shapeCoord[ii,2],shapeCoord[ii+1,2]),method="linear", n=50)
      shapeCoordListx[[ii]]<-xy$x
      shapeCoordListy[[ii]]<-xy$y
    }
    shapeCoord<-data.frame(x=unlist(shapeCoordListx),y=unlist(shapeCoordListy))
    shapeCoord<-unique(shapeCoord)
    
    Morph <- datList[,1:2]
    MS<-rbind(Morph,shapeCoord)
    
    borderDist<-pointDistance(Morph,shapeCoord,F)
    borderDist2<-rowMins(as.matrix(borderDist))
    
    InputDat$borderDist[i]<-sum(borderDist2)
    # print(i)
  }
  
  borderDist<-matrix(ncol=length(nVertex),nrow=length(seq(1,MaxArea,1)))
  OptShape<-data.frame()
  for(i in 1:length(nVertex)){
    toplot<-which(InputDat$nVertex==nVertex[i])
    borderDist[,i]<-sapply(split(InputDat$borderDist[toplot],InputDat$Area[toplot]),min)
    mD<-min(borderDist[,i])
    OptShape<-rbind(OptShape,InputDat[which(InputDat$nVertex==nVertex[i] & InputDat$borderDist==mD),])
    #print(i)
  }  
  print(OptShape[,c(1,7)])
  save(pc.dat23,OptShape,file=paste0(outFolderBM,"RegPolygonFit_Null_BM_",Ni,"_Body.rda"))
  
  cat(Ni,"/",Nnull,"\n")
  
} # End for Ni (Null model i)

## Prepare output for Fig S3 (Fitting of regular polygons)

datList<-list()

for(x in 1:2){
  if(x==1){
    pc.dat23<-read.csv(paste0(DataFolder,"data_tmp/BeakPC23.csv"))
  }
  if(x==2){
    pc.dat23<-read.csv(paste0(DataFolder,"data_tmp/BodyPC23.csv"))
  }
  
  pc12Foc<-as.matrix(data.frame(PC1=pc.dat23[,2],PC2=pc.dat23[,3]))
  H12Foc <- Hpi(x=pc12Foc)      # optimal bandwidth estimation
  est12Foc<- kde(x=pc12Foc, H=H12Foc, compute.cont=TRUE)     # kernel density estimation
  
  pc12<-pc12Foc
  H12<-H12Foc
  est12<-est12Foc
  
  dat<-expand.grid(x=est12$eval.points[[1]],y=est12$eval.points[[2]])
  
  dat$x2<-rep(1:length(est12$eval.points[[1]]),length(est12$eval.points[[2]]))
  dat$y2<-rep(1:length(est12$eval.points[[1]]),each=length(est12$eval.points[[2]]))
  
  dat$estimate<-as.vector(est12$estimate)
  dat<-dat[which(dat$estimate>est12$cont[which(names(est12$cont)=="10%")]),]
  
  ydiff<-xdiff<-1
  dat$nb<-NA
  for(i in 1:nrow(dat)){
    left<-which(round(dat$x2)==round(dat$x2[i]-xdiff) & round(dat$y2)==round(dat$y2[i]))
    right<-which(round(dat$x2)==round(dat$x2[i]+xdiff) & round(dat$y2)==round(dat$y2[i]))
    
    top<-which(round(dat$y2)==round(dat$y2[i]+ydiff) & round(dat$x2)==round(dat$x2[i]))
    low<-which(round(dat$y2)==round(dat$y2[i]-ydiff) & round(dat$x2)==round(dat$x2[i]))
    
    topleft<-which(round(dat$x2)==round(dat$x2[i]-xdiff) & round(dat$y2)==round(dat$y2[i]+ydiff))
    topright<-which(round(dat$x2)==round(dat$x2[i]+xdiff) & round(dat$y2)==round(dat$y2[i]+ydiff))
    
    lowleft<-which(round(dat$x2)==round(dat$x2[i]-xdiff) & round(dat$y2)==round(dat$y2[i]-ydiff))
    lowright<-which(round(dat$x2)==round(dat$x2[i]+xdiff) & round(dat$y2)==round(dat$y2[i]-ydiff))
    
    dat$nb[i]<-length(c(left,right,top,low,topleft,topright,lowleft,lowright))
  }
  table(dat$nb)
  dat<-dat[dat$nb<8,]
  datList[[x]]<-dat
}  

for(x in 1:2){
  if(x==1){pc.dat23<-read.csv(paste0(DataFolder,"data_tmp/BeakPC23.csv"))}
  if(x==2){pc.dat23<-read.csv(paste0(DataFolder,"data_tmp/BodyPC23.csv"))}
  
  MaxArea<-ceiling(diff(range(pc.dat23[,2]))*diff(range(pc.dat23[,3])))/2
  nVertex<-c(3,4,5,6,7,360)
  InputDat<-data.frame()
  for(i in 1:length(nVertex)){
    #dat<-expand.grid(nVertex=nVertex[i],Area=seq(1,MaxArea,1),xpos=seq(-0.5,0.5,0.25),ypos=seq(-0.5,0.5,0.25),Rotation=seq(0,360/nVertex[i],1),Position=0,borderDist=NA)
    dat<-expand.grid(nVertex=nVertex[i],Area=seq(1,MaxArea,1),xpos=0,ypos=0,Rotation=seq(0,360/nVertex[i],1),Position=0,borderDist=NA)
    InputDat<-rbind(InputDat,dat)
  }  
  
  for(i in 1:nrow(InputDat)){
    
    shape<-regular.poly(nSides=InputDat$nVertex[i], area=InputDat$Area[i])
    
    shape<-Rotate(x=shape$x, y=shape$y,my=NULL,theta=deg2rad(InputDat$Rotation[i]))[1:2]
    
    shape$x<-shape$x+InputDat$xpos[i]
    shape$y<-shape$y+InputDat$ypos[i]
    
    shapeCoord<-matrix(ncol=2,nrow=InputDat$nVertex[i])
    shapeCoord[,1]<-shape$x
    shapeCoord[,2]<-shape$y
    shapeCoord<-rbind(shapeCoord,shapeCoord[1,])
    shapeCoord[,1]<-shapeCoord[,1]+rnorm(nrow(shapeCoord),0,0.0001)
    shapeCoord[,2]<-shapeCoord[,2]+rnorm(nrow(shapeCoord),0,0.0001)
    shapeCoordListx<-list()
    shapeCoordListy<-list()
    for(ii in 1:(nrow(shapeCoord)-1)){
      xy<-approx(x=c(shapeCoord[ii,1],shapeCoord[ii+1,1]), y=c(shapeCoord[ii,2],shapeCoord[ii+1,2]),method="linear", n=50)
      shapeCoordListx[[ii]]<-xy$x
      shapeCoordListy[[ii]]<-xy$y
    }
    shapeCoord<-data.frame(x=unlist(shapeCoordListx),y=unlist(shapeCoordListy))
    shapeCoord<-unique(shapeCoord)
    
    Morph<-datList[[x]][,1:2]
    MS<-rbind(Morph,shapeCoord)
    
    borderDist<-pointDistance(Morph,shapeCoord,F)
    borderDist2<-rowMins(as.matrix(borderDist))
    
    InputDat$borderDist[i]<-sum(borderDist2)
    print(i)
  }
  
  borderDist<-matrix(ncol=length(nVertex),nrow=length(seq(1,MaxArea,1)))
  OptShape<-data.frame()
  for(i in 1:length(nVertex)){
    toplot<-which(InputDat$nVertex==nVertex[i])
    borderDist[,i]<-sapply(split(InputDat$borderDist[toplot],InputDat$Area[toplot]),min)
    mD<-min(borderDist[,i])
    
    print(InputDat[which(InputDat$nVertex==nVertex[i] & InputDat$borderDist==mD),])
    OptShape<-rbind(OptShape,InputDat[which(InputDat$nVertex==nVertex[i] & InputDat$borderDist==mD),])
    print(i)
  }  
  
  if(x==1){save(InputDat,borderDist,MaxArea,OptShape,file=paste0(OutputFolder,"out2_RegularPolygonFitting_Beak.rda"))}
  if(x==2){save(InputDat,borderDist,MaxArea,OptShape,file=paste0(OutputFolder,"out2_RegularPolygonFitting_Body.rda"))}
}

## End script 2.4 ##