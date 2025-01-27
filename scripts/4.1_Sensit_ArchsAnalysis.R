## 4.1 ## Run sensitivity analyses:

# Includes archetype analysis for Beak scan (A), PPCA Beak data (B) and PPCA Body data (C)
# The results will be used to obtain the archetype coordinates for Figure S4.

## 4.1 A: Archetypes with Beak scan PCs ####

# Create folder "PCBeakArchetypes" inside output folder:
outFolder<-paste0(OutputFolder,'out4_PCBeakArchetypesScan/')
if(!file.exists(outFolder)) dir.create(outFolder)

# Open beak shape (PC2, PC3).
pc.dat23 <- read.csv(paste0(DataTmpFolder,"BeakPC23.csv"))
# Open beak shape scan PCs (PC1,PC2 from Cooney et al. 2017)
pc.dat12 <-read.csv(paste0(DataTmpFolder,"BeakScanPC12.csv")) #Morphology data

fdat <- merge(pc.dat23,pc.dat12,by="species")
fdat$PC1<-scale(fdat$PC1)
fdat$PC2<-scale(fdat$PC2)

# Plot correlation of Beak Scan PCs with Beak linear PCs.
plot(fdat$PC2_beak,fdat$PC1)
plot(fdat$PC3_beak,fdat$PC2)
cor(fdat$PC2_beak,fdat$PC1)
cor(fdat$PC3_beak,fdat$PC2)

## Run archetype analysis with PC1 and PC2 (from Beak Scan)

pc.dat12 <- pc.dat12[,-1]

Nrep <- 100 
Nkmin<-1
Nk <- 3
Kvals<-Nkmin:Nk
maxTrials <- 10

for(i in 1:length(Kvals)){
  
  # aim to obtain Nrep estimates
  for(ii in 1:Nrep){
    
    #only run the model if the k by rep combination does not yet exist 
    outfile<-paste0(outFolder,"ArchetypesObsScan_PC12Traits_RobustAlg_k_",
                    Kvals[i],"_rep_",ii,"_Beaks.rda")
    exF<-paste0(outFolder,list.files(outFolder))
    
    if(length(intersect(outfile,exF))==0){
      
      for(iii in 1:maxTrials){
        print(paste0(i,"_",ii,"_",iii))
        #run may fail to converge
        asAR<-NULL
        try(asAR<-robustArchetypes(data=pc.dat12, k=Kvals[i]))
        
        rssV<-rss(asAR)
        arch<-asAR$archetypes
        
        # if it coverges save the result ad move on to the next rep
        if(is.null(arch)==FALSE){
          save(rssV,arch,file=outfile)
          break
        }
      }    
    }
  }  
}

rssMat<-matrix(nrow=length(Kvals),ncol=Nrep)
for(i in 1:length(Kvals)){
  
  # aim to obtain Nrep estimates
  for(ii in 1:Nrep){
    
    outfile<-paste0(outFolder,"ArchetypesObsScan_PC12Traits_RobustAlg_k_",
                    Kvals[i],"_rep_",ii,"_Beaks.rda")
    exF<-paste0(outFolder,list.files(outFolder))
    
    if(length(intersect(outfile,exF))==1){
      (load(outfile))
      rssMat[i,ii]<-rssV
    }  
  } 
  print(i)
}

rssMatC<-rssMat
rssMatC[rssMat>0]<-1
rowSums(rssMatC,na.rm=T)
# K values > 0
Kvals_No0 <- sum(rowSums(rssMatC,na.rm=T)>0)

archCoords<-list()
for(i in 1:Kvals_No0){
  ii<-which.min(rssMat[i,])
  outfile<-paste0(outFolder,"ArchetypesObsScan_PC12Traits_RobustAlg_k_",
                  Kvals[i],"_rep_",ii,"_Beaks.rda")
  load(outfile)
  archCoords[[i]]<-arch
} 

save(rssMat,rssMatC,archCoords,file=paste0(OutputFolder,"out4_ArchetypesObsScan_RobustAlg_PC12_Beaks.rda"))

## 

## 4.1 B: Archetypes with PPCA (Beak) ####

# Create folder "PCBeakArchetypes" inside output folder:
outFolder <- paste0(OutputFolder,'out4_PPCA_Beak/')
if(!file.exists(outFolder)) dir.create(outFolder)

# Open morphological data to calculate beak PPCA:
fdat <- read.csv(paste0(DataTmpFolder,"Avonet_BirdTree.csv"))
fdat <- na.omit(fdat[,c("species","Beak.Length_Culmen","Beak.Width","Beak.Depth")])

# PCA of 3 beak traits
x <- log(fdat[,c("Beak.Length_Culmen","Beak.Width","Beak.Depth")]) # # Log transform
x <- scale(x)
rownames(x) <- fdat$species

# Open tree and prune
tree <- read.nexus(paste0(DataFolder,"phylo/AllBirdsHackett1_summary.tre"))
tree.beak <- drop.tip(tree,tree$tip.label[!tree$tip.label %in% fdat$species])

ppca <- phyl.pca(tree = tree.beak,x)
ppca
ppca.spp <- data.frame(ppca$S)
ppca.spp$species <- rownames(ppca$S)
save(ppca,ppca.spp,file=paste0(DataFolder,"data_tmp/BeakPPCA.Rdata"))

# Load PPC beak shape (PC2, PC3)
load(paste0(DataFolder,"data_tmp/BeakPPCA.Rdata"))

pc.dat23 <- ppca.spp[,c("PC2","PC3")]

Nrep <- 100 
Nkmin<-1
Nk <- 3
Kvals<-Nkmin:Nk
maxTrials<-4

# Run archetype analysis:
for(i in 1:length(Kvals)){
  
  # aim to obtain Nrep estimates
  for(ii in 1:Nrep){
    
    #only run the model if the k by rep combination does not yet exist 
    outfile <- paste0(outFolder,"ArchetypesObs_PC23Traits_RobustAlg_k_",
                      Kvals[i],"_rep_",ii,"_Beaks.rda")
    exF<-paste0(outFolder,list.files(outFolder))
    
    if(length(intersect(outfile,exF))==0){
      
      for(iii in 1:maxTrials){
        print(paste0(i,"_",ii,"_",iii))
        #run may fail to converge
        asAR<-NULL
        try(asAR<-robustArchetypes(data=pc.dat23, k=Kvals[i]))
        
        rssV<-rss(asAR)
        arch<-asAR$archetypes
        
        # if it coverges save the result ad move on to the next rep
        if(is.null(arch)==FALSE){
          save(rssV,arch,file=outfile)
          break
        }
      }    
    }
  }  
}

# extract rss scores and archetypes coords (Beak shape)

rssMat<-matrix(nrow=length(Kvals),ncol=Nrep)
for(i in 1:length(Kvals)){
  
  # aim to obtain Nrep estimates
  for(ii in 1:Nrep){
    
    outfile<-paste0(outFolder,"ArchetypesObs_PC23Traits_RobustAlg_k_",
                    Kvals[i],"_rep_",ii,"_Beaks.rda")
    exF<-paste0(outFolder,list.files(outFolder))
    
    if(length(intersect(outfile,exF))==1){
      (load(outfile))
      rssMat[i,ii]<-rssV
    }  
  } 
  print(i)
}

rssMatC<-rssMat
rssMatC[rssMat>0]<-1
rowSums(rssMatC,na.rm=T)
# K values > 0
Kvals_No0 <- sum(rowSums(rssMatC,na.rm=T)>0)

archCoords<-list()
for(i in 1:Kvals_No0){
  ii<-which.min(rssMat[i,])
  outfile<-paste0(outFolder,"ArchetypesObs_PC23Traits_RobustAlg_k_",
                  Kvals[i],"_rep_",ii,"_Beaks.rda")
  (load(outfile))
  archCoords[[i]]<-arch
} 

save(rssMat,rssMatC,archCoords,file=paste0(OutputFolder,"out4_ArchetypesObs_RobustAlg_PPC23_Beaks.rda"))

## 4.1 C: Archetypes with PPCA (Body) ####

# Create folder "PCBeakArchetypes" inside output folder:
outFolder <- paste0(OutputFolder,'out4_PPCA_Body/')
if(!file.exists(outFolder)) dir.create(outFolder)

# Open morphological data to calculate body PPCA:
fdat <- read.csv(paste0(DataTmpFolder,"Avonet_BirdTree.csv"))
fdat <- na.omit(fdat[,c("species","Tarsus.Length","Wing.Length","Kipps.Distance","Tail.Length")])
fdat <- fdat[-which(fdat$Wing.Length==0.1),]

# PCA of 4 beak traits
x <- log(fdat[,c("Tarsus.Length","Wing.Length","Kipps.Distance","Tail.Length")]) # # Log transform
x <- scale(x)
rownames(x) <- fdat$species

# Prune body tree
tree.body <- drop.tip(tree,tree$tip.label[!tree$tip.label %in% fdat$species])
ppca <- phyl.pca(tree = tree.body,x)
ppca.spp <- data.frame(ppca$S)
ppca.spp$species <- rownames(ppca$S)
save(ppca,ppca.spp,file=paste0(DataFolder,"data_tmp/BodyPPCA.Rdata"))

# Load PPC beak shape (PC2, PC3)
load(paste0(DataFolder,"data_tmp/BodyPPCA.Rdata"))

pc.dat23 <- ppca.spp[,c("PC2","PC3")]

Nrep <- 100 
Nkmin<-1
Nk <- 3
Kvals<-Nkmin:Nk
maxTrials<-4

# Run archetype analysis:
for(i in 1:length(Kvals)){
  
  # aim to obtain Nrep estimates
  for(ii in 1:Nrep){
    
    #only run the model if the k by rep combination does not yet exist 
    outfile <- paste0(outFolder,"ArchetypesObs_PC23Traits_RobustAlg_k_",
                      Kvals[i],"_rep_",ii,"_Body.rda")
    exF<-paste0(outFolder,list.files(outFolder))
    
    if(length(intersect(outfile,exF))==0){
      
      for(iii in 1:maxTrials){
        print(paste0(i,"_",ii,"_",iii))
        #run may fail to converge
        asAR<-NULL
        try(asAR<-robustArchetypes(data=pc.dat23, k=Kvals[i]))
        
        rssV<-rss(asAR)
        arch<-asAR$archetypes
        
        # if it coverges save the result ad move on to the next rep
        if(is.null(arch)==FALSE){
          save(rssV,arch,file=outfile)
          break
        }
      }    
    }
  }  
}

# extract rss scores and archetypes coords (Body shape)

rssMat<-matrix(nrow=length(Kvals),ncol=Nrep)
for(i in 1:length(Kvals)){
  
  # aim to obtain Nrep estimates
  for(ii in 1:Nrep){
    
    outfile<-paste0(outFolder,"ArchetypesObs_PC23Traits_RobustAlg_k_",
                    Kvals[i],"_rep_",ii,"_Body.rda")
    exF<-paste0(outFolder,list.files(outFolder))
    
    if(length(intersect(outfile,exF))==1){
      (load(outfile))
      rssMat[i,ii]<-rssV
    }  
  } 
  print(i)
}

rssMatC<-rssMat
rssMatC[rssMat>0]<-1
rowSums(rssMatC,na.rm=T)
# K values > 0
Kvals_No0 <- sum(rowSums(rssMatC,na.rm=T)>0)

archCoords<-list()
for(i in 1:Kvals_No0){
  ii<-which.min(rssMat[i,])
  outfile<-paste0(outFolder,"ArchetypesObs_PC23Traits_RobustAlg_k_",
                  Kvals[i],"_rep_",ii,"_Body.rda")
  (load(outfile))
  archCoords[[i]]<-arch
} 

save(rssMat,rssMatC,archCoords,file=paste0(OutputFolder,"out4_ArchetypesObs_RobustAlg_PPC23_Body.rda"))

## End script 4.1 ##