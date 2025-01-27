## 2.1 ## Calculate Beak and Body shape axes (PCAs)

## 1 ## Beak PCA ####

# Open morphology dataset:
fdat <- read.csv(paste0(DataTmpFolder,"Avonet_BirdTree.csv"))
fdat <- na.omit(fdat[,c("species","Beak.Length_Culmen","Beak.Width","Beak.Depth")])

# PCA of 3 beak traits
x <- log(fdat[,c("Beak.Length_Culmen","Beak.Width","Beak.Depth")]) # # Log transform
x <- scale(x)
rownames(x) <- fdat$species
pca <- prcomp(x)
pca
pca$rotation
summary(pca)
# Invert PCs to facilitate interpretation (i.e., PC1 from small to large)
pca$rotation*(-1)
pcx <- data.frame(pca$x)*(-1)
pcx$species <- rownames(pca$x)
names(pcx)[1:3] <- c("PC1_beak","PC2_beak","PC3_beak")

summary(pcx$PC1_beak)
# Scale PCs
pcx$PC1_beak <- scale(pcx$PC1_beak)
pcx$PC2_beak <- scale(pcx$PC2_beak)
pcx$PC3_beak <- scale(pcx$PC3_beak)

fdat <- merge(fdat,pcx[,c(1:3,ncol(pcx))],by="species")

pc.dat23 <- fdat[,c("species","PC2_beak","PC3_beak")]
write.csv(pc.dat23,file=paste0(DataFolder,"data_tmp/BeakPC23.csv"),row.names=FALSE)

pc.dat123 <- fdat[,c("species","PC1_beak","PC2_beak","PC3_beak")]
write.csv(pc.dat123,file=paste0(DataFolder,"data_tmp/BeakPC123.csv"),row.names=FALSE)

## 2 ## Body PCA ####

# Open morphology dataset:
fdat <- read.csv(paste0(DataTmpFolder,"Avonet_BirdTree.csv"))
fdat <- na.omit(fdat[,c("species","Tarsus.Length","Wing.Length","Kipps.Distance","Tail.Length")])
fdat <- fdat[-which(fdat$Wing.Length==0.1),]

# PCA of 4 beak traits
x <- log(fdat[,c("Tarsus.Length","Wing.Length","Kipps.Distance","Tail.Length")]) # # Log transform

x <- scale(x)
rownames(x) <- fdat$species
pca <- prcomp(x)
pca
summary(pca)
# Invert PCs to facilitate interpretation (i.e., PC1 from small to large)
pca$rotation*(-1)
pcx <- data.frame(pca$x)*(-1)
pcx$species <- rownames(pca$x)
names(pcx)[1:4] <- c("PC1_body","PC2_body","PC3_body","PC4_body")

# Scale PCs
# Invert PCs to facilitate interpretation (i.e., PC1 from small to large)
pcx$PC1_body <- scale(pcx$PC1_body)
pcx$PC2_body <- scale(pcx$PC2_body)
pcx$PC3_body <- scale(pcx$PC3_body)
pcx$PC4_body <- scale(pcx$PC4_body)
fdat <- merge(fdat,pcx[,c(1:4,ncol(pcx))],by="species")

pc.dat23 <- fdat[,c("species","PC2_body","PC3_body")]

write.csv(pc.dat23,file=paste0(DataFolder,"data_tmp/BodyPC23.csv"),row.names=FALSE)

pc.dat123 <- fdat[,c("species","PC1_body","PC2_body","PC3_body")]

write.csv(pc.dat123,file=paste0(DataFolder,"data_tmp/BodyPC123.csv"),row.names=FALSE)

## End script 2.1 ##