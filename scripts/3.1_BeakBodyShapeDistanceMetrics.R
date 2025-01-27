## 3.1 ## Beak and body shape distance metrics

## 1 ## Beak shape distance metrics ####

# Open morphology dataset:
fdat <- read.csv(paste0(DataTmpFolder,"Avonet_BirdTree.csv"))
fdat <- na.omit(fdat)

# Merge with beak shape data:
pc.dat23 <- read.csv(paste0(DataFolder,"data_tmp/BeakPC23.csv"))
fdat <- merge(fdat,pc.dat23,by="species")

# Open robust archetypes for beak shape:
load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Beaks.rda"))
a3Rxy<-archCoords[[3]]

## Calculate metrics of morphological specialization:

distMat<-dist(fdat[,c("PC2_beak","PC3_beak")])
distMat<-as.matrix(distMat)
diag(distMat)<-NA
fdat$distNN<-log(rowMins(distMat,na.rm=T))

# Metric 1: Distance to morphospace centroid
pCent <- c(mean(fdat$PC2_beak),mean(fdat$PC3_beak))
fdat$distCent <- pointDistance(fdat[c("PC2_beak","PC3_beak")],pCent,lonlat = F, allpairs = T)

# Metric 2: Distance to closest archetype, assuming 3 vertices 
sps.dista <- pointDistance(fdat[,c("PC2_beak","PC3_beak")],a3Rxy,F) # Distance between Sps and Archetypes.

fdat$distAR1 <- sps.dista[,3] # A1 (Cockatoo archetype)
fdat$distAR3 <- sps.dista[,1] # A3 (Frogmouth archetype)
fdat$distAR2 <- sps.dista[,2] # A2 (Jacamar archetype)

# Select minimum for each row.
fdat$distMinAR <- apply(sps.dista,1,min)

distMetrics<-c("distCent",
               "distAR1","distAR2","distAR3",
               "distMinAR")
fdat <- fdat[,c(1,which(names(fdat)%in%distMetrics))]

# Save distance metrics data
write.csv(fdat,file=paste0(DataFolder,"data_tmp/BeakDistanceMetrics.csv"),row.names=FALSE)

## 2 ## Body shape distance metrics ####

# Open morphology dataset:
fdat <- read.csv(paste0(DataTmpFolder,"Avonet_BirdTree.csv"))
fdat <- na.omit(fdat)

# Merge with body shape data:
pc.dat23 <- read.csv(paste0(DataFolder,"data_tmp/BodyPC23.csv"))
fdat <- merge(fdat,pc.dat23,by="species")

# Open robust archetypes for body shape:
load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Body.rda"))
a3Rxy<-archCoords[[3]]

## Calculate metrics of morphological specialization:

# Metric 1: Distance to morphospace centroid
pCent <- c(mean(fdat$PC2_body),mean(fdat$PC3_body))
fdat$distCent <- pointDistance(fdat[c("PC2_body","PC3_body")],pCent,lonlat = F, allpairs = T)

# Metric 2: Distance to closest archetype, assuming 3 vertices 
sps.dista <- pointDistance(fdat[,c("PC2_body","PC3_body")],a3Rxy,F) # Distance between Sps and Archetypes.

fdat$distAR1 <- sps.dista[,2] # A1 (Grebe archetype)
fdat$distAR2 <- sps.dista[,3] # A2 (Ground-Cuckoo archetype)
fdat$distAR3 <- sps.dista[,1] # A3 (Palm swift archetype)

# Select minimum for each row.
fdat$distMinAR <- apply(sps.dista,1,min)

distMetrics<-c("distCent",
               "distAR1","distAR2","distAR3",
               "distMinAR")

fdat <- fdat[,c(1,which(names(fdat)%in%distMetrics))]

# Save body shape distance metrics
write.csv(fdat,file=paste0(DataFolder,"data_tmp/BodyDistanceMetrics.csv"),row.names=FALSE)

## End script 3.1 ##