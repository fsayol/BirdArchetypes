## 5.1 ## Directionality in trait evolution (Global analysis)

# For beak and body shape (PCs), we will perform an ancestral state reconstruction at each node (ASR).

# For each node, we will compute the distance of the ASR value to the morphospace centroid
# We will compare the change in directionality with null expectations (Evolving under a BM model).

# We create folder to save all ancestral state reconstructions (ASR) outputs
outFolder <-paste0(OutputFolder,'out5_ASR_Empirical_and_Null_RateHet/')
if(!file.exists(outFolder)) dir.create(outFolder)

## Beak shape empirical ####

# 1 # Perform ASR at each node

pc.dat23beak <- read.csv(paste0(DataFolder,"data_tmp/BeakPC23.csv"))

# Open phylogenies (Rate-transformed trees)
ptree2 <- read.tree(paste0(DataFolder,"phylo/tree_rateTrans_PC2beak.tre"))
ptree3 <- read.tree(paste0(DataFolder,"phylo/tree_rateTrans_PC3beak.tre"))

# select only species in tree
pdat <- pc.dat23beak[pc.dat23beak$species %in% ptree2$tip.label,]
# Select PC2 and PC3
pc2 <- pdat$PC2_beak
pc3 <- pdat$PC3_beak
names(pc2) <- names(pc3) <- pdat$species
# Prune tree to species in dataset
ptree2 <- drop.tip(ptree2,ptree2$tip.label[!ptree2$tip.label %in% pdat$species])
ptree3 <- drop.tip(ptree3,ptree3$tip.label[!ptree3$tip.label %in% pdat$species])
pc.obs <- data.frame("species"=pdat$species,"pc2"=pc2,"pc3"=pc3)
pc.obs <- pc.obs[pc.obs$species %in% ptree2$tip.label,] # Order as tips
# ancestral reconstructions PCs
rec.pc2obs <- fastAnc(ptree2,pc.obs$pc2)
rec.pc3obs <- fastAnc(ptree3,pc.obs$pc3)
node.traits <- data.frame("Nodes"=names(rec.pc2obs),cbind(rec.pc2obs,rec.pc3obs))
node.pairs <- data.frame(ptree2$edge,ptree2$edge.length)
names(node.pairs) <- c("Node.Anc","Node.Des","edge.length")
node.traitsA <- node.traits[,c(1:3)]
names(node.traitsA) <- c("Node.Anc","pc2.Anc","pc3.Anc")
node.pairs <- merge(node.pairs,node.traitsA,by="Node.Anc")
node.traitsD <- node.traits[,c(1:3)]
names(node.traitsD) <- c("Node.Des","pc2.Des","pc3.Des")
node.pairs <- merge(node.pairs,node.traitsD,by="Node.Des")

# Save ASR of beak observed traits
save(node.traits,node.pairs,file=paste0(outFolder,"AncestralReconstructionBeakTraits_RateHet_Observed.rda"))

# 2 # Calculate distance to centroid

# Centroid of empirical data
centroid.spp <- data.frame("pc2"=mean(pc.dat23$PC2_beak),"pc3"=mean(pc.dat23$PC3_beak))

# Vector to save results:
ratio.BeakObsCen <- NULL

# Load results ASR Beak Obsserved:
load(file=paste0(outFolder,"AncestralReconstructionBeakTraits_RateHet_Observed.rda"))

# Distance to centroid: 
node.traits$distCen <- pointDistance(node.traits[,c("rec.pc2obs","rec.pc3obs")],centroid.spp,F) # Distance between Sps and centroid
# Merge distance metrics with Node pairs.
node.traitsA <- node.traits[,c("Nodes","distCen")]
names(node.traitsA) <- c("Node.Anc","distCen.Anc")
node.pairs <- merge(node.pairs,node.traitsA,by="Node.Anc")
node.traitsD <- node.traits[,c("Nodes","distCen")]
names(node.traitsD) <- c("Node.Des","distCen.Des")
node.pairs <- merge(node.pairs,node.traitsD,by="Node.Des")

# Distance to centroid change (descendant vs ancestor)
node.pairs$dif.CenDist <- node.pairs$distCen.Des - node.pairs$distCen.Anc
node.pairs$dif.CenSign <- ifelse(node.pairs$dif.CenDist>0,"Reduced (closer Arch)","Increased (dist Arch)")
ratio.BeakObsCen <- table(node.pairs$dif.CenSign)[2]/nrow(node.pairs) 
# Save results of Beak observed traits
save(ratio.BeakObsCen,file=paste0(outFolder,"BeakDirectionEvolved_RateHet_Observed.rda"))

## Beak shape simulated ####

# 1 # Perform ASR at each node

# Open phylogenies (Rate-transformed trees)
ptree2 <- read.tree(paste0(DataFolder,"phylo/tree_rateTrans_PC2beak.tre"))
ptree3 <- read.tree(paste0(DataFolder,"phylo/tree_rateTrans_PC3beak.tre"))

Nsim <- 100

for(sim_i in 1:Nsim){
  
# Simulate BM on the tree:  
  pc2 <- data.frame(fastBM(tree = ptree2))
  pc2$species <- rownames(pc2)
  pc3 <- data.frame(fastBM(tree = ptree3))
  pc3$species <- rownames(pc3)
  pc.sim <- merge(pc2,pc3,by="species")
  names(pc.sim) <- c("species","PC2_beak","PC3_beak")
  
  pc2 <- pc.sim$PC2_beak
  pc3 <- pc.sim$PC3_beak
  names(pc2) <- names(pc3) <- pc.sim$species
  
  # ancestral reconstructions PCs
  rec.pc2sim <- fastAnc(ptree2,pc2)
  rec.pc3sim <- fastAnc(ptree3,pc3)
  node.traits <- data.frame("Nodes"=names(rec.pc2sim),cbind(rec.pc2sim,rec.pc3sim))
  
  node.pairs <- data.frame(ptree2$edge,ptree2$edge.length)
  names(node.pairs) <- c("Node.Anc","Node.Des","edge.length")
  node.traitsA <- node.traits[,c(1:3)]
  names(node.traitsA) <- c("Node.Anc","pc2.Anc","pc3.Anc")
  node.pairs <- merge(node.pairs,node.traitsA,by="Node.Anc")
  node.traitsD <- node.traits[,c(1:3)]
  names(node.traitsD) <- c("Node.Des","pc2.Des","pc3.Des")
  node.pairs <- merge(node.pairs,node.traitsD,by="Node.Des")

  # Save ASR of beak simulated traits
  save(node.traits,node.pairs,pc.sim,file=paste0(outFolder,"AncestralReconstructionBeak_RateHet_",sim_i,".rda"))
  cat("Beak simulated Tree",sim_i,"\n")
} # end focal tree

# 2 # Calculate distance to centroid

Nsim <- 100
ratio.BeakSimCen <- numeric(Nsim)

for(sim_i in 1:Nsim){
  # Load results from simulated reconstructions.
  load(file=paste0(outFolder,"AncestralReconstructionBeak_RateHet_",sim_i,".rda"))
  # Distance to centroid
  centroidSpp <- data.frame("pc2"=mean(pc.sim$PC2_beak),"pc3"=mean(pc.sim$PC3_beak))
  node.traits$distCen <- pointDistance(node.traits[,c("rec.pc2sim","rec.pc3sim")],centroidSpp,F) # Distance between Sps and centroid
  # Merge distance metrics with Node pairs.
  node.traitsA <- node.traits[,c("Nodes","distCen")]
  names(node.traitsA) <- c("Node.Anc","distCen.Anc")
  node.pairs <- merge(node.pairs,node.traitsA,by="Node.Anc")
  node.traitsD <- node.traits[,c("Nodes","distCen")]
  names(node.traitsD) <- c("Node.Des","distCen.Des")
  node.pairs <- merge(node.pairs,node.traitsD,by="Node.Des")
  # Distance change from centroid:
  node.pairs$dif.CenDist <- node.pairs$distCen.Des - node.pairs$distCen.Anc
  node.pairs$dif.CenSign <- ifelse(node.pairs$dif.CenDist>0,"Reduced (closer Arch)","Increased (dist Arch)")
  ratio.BeakSimCen[sim_i] <- table(node.pairs$dif.CenSign)[2]/nrow(node.pairs)
  cat("Tree",sim_i,"ratioCen:",ratio.BeakSimCen[sim_i],"\n")
  
} # End for simulation "i"

# Save results of beak simulated traits
save(ratio.BeakSimCen,file=paste0(outFolder,"BeakDirectionEvolved_RateHet_Simulated.rda"))

## Body shape empirical ####

# 1 # Perform ASR at each node

# Open body shape data
pc.dat23body <- read.csv(paste0(DataFolder,"data_tmp/BodyPC23.csv"))

# Open phylogenies (Rate-transformed trees)
ptree2 <- read.tree(paste0(DataFolder,"phylo/tree_rateTrans_PC2body.tre"))
ptree3 <- read.tree(paste0(DataFolder,"phylo/tree_rateTrans_PC3body.tre"))

# select only species in tree
pdat <- pc.dat23body[pc.dat23body$species %in% ptree2$tip.label,]
pc2 <- pdat$PC2_body
pc3 <- pdat$PC3_body
names(pc2) <- names(pc3) <- pdat$species
  
# Prune tree to species in dataset
ptree2 <- drop.tip(ptree2,ptree2$tip.label[!ptree2$tip.label %in% pdat$species])
ptree3 <- drop.tip(ptree3,ptree3$tip.label[!ptree3$tip.label %in% pdat$species])
pc.obs <- data.frame("species"=pdat$species,"pc2"=pc2,"pc3"=pc3)
pc.obs <- pc.obs[pc.obs$species %in% ptree2$tip.label,] # Order as tips
# ancestral reconstructions PCs
rec.pc2obs <- fastAnc(ptree2,pc.obs$pc2)
rec.pc3obs <- fastAnc(ptree3,pc.obs$pc3)
node.traits <- data.frame("Nodes"=names(rec.pc2obs),cbind(rec.pc2obs,rec.pc3obs))
node.pairs <- data.frame(ptree2$edge,ptree2$edge.length)
names(node.pairs) <- c("Node.Anc","Node.Des","edge.length")
node.traitsA <- node.traits[,c(1:3)]
names(node.traitsA) <- c("Node.Anc","pc2.Anc","pc3.Anc")
node.pairs <- merge(node.pairs,node.traitsA,by="Node.Anc")
node.traitsD <- node.traits[,c(1:3)]
names(node.traitsD) <- c("Node.Des","pc2.Des","pc3.Des")
node.pairs <- merge(node.pairs,node.traitsD,by="Node.Des")

# Save ASR of body simulated traits
save(node.traits,node.pairs,file=paste0(outFolder,"AncestralReconstructionBodyTraits_RateHet_Observed.rda"))

# 2 # Calculate distance to centroid

# Vector to save results:
ratio.BodyObsCen <- NULL
# Centroid of empirical data (Beak):
centroid.spp <- data.frame("pc2"=mean(pc.dat23$PC2_body),"pc3"=mean(pc.dat23$PC3_body))
# Load results ASR Beak Obsserved:
load(file=paste0(outFolder,"AncestralReconstructionBodyTraits_RateHet_Observed.rda"))

# Distance to centroid
node.traits$distCen <- pointDistance(node.traits[,c("rec.pc2obs","rec.pc3obs")],centroid.spp,F) # Distance between Sps and centroid
# Merge distance metrics with Node pairs.
node.traitsA <- node.traits[,c("Nodes","distCen")]
names(node.traitsA) <- c("Node.Anc","distCen.Anc")
node.pairs <- merge(node.pairs,node.traitsA,by="Node.Anc")
node.traitsD <- node.traits[,c("Nodes","distCen")]
names(node.traitsD) <- c("Node.Des","distCen.Des")
node.pairs <- merge(node.pairs,node.traitsD,by="Node.Des")  
# Distance to centroid change (descendant vs ancestor)
node.pairs$dif.CenDist <- node.pairs$distCen.Des - node.pairs$distCen.Anc
node.pairs$dif.CenSign <- ifelse(node.pairs$dif.CenDist>0,"Reduced (closer Arch)","Increased (dist Arch)")
ratio.BodyObsCen <- table(node.pairs$dif.CenSign)[2]/nrow(node.pairs)

# Save directionality of body simulated traits
save(ratio.BodyObsCen,file=paste0(outFolder,"BodyDirectionEvolved_RateHet_Observed.rda"))

## Body shape simulated ####

# 1 # Perform ASR at each node

# Open phylogenies (Rate-transformed trees)
ptree2 <- read.tree(paste0(DataFolder,"phylo/tree_rateTrans_PC2body.tre"))
ptree3 <- read.tree(paste0(DataFolder,"phylo/tree_rateTrans_PC3body.tre"))

Nsim <- 100

for(sim_i in 1:Nsim){

  # Simulate BM on the tree:  
  pc2 <- data.frame(fastBM(tree = ptree2))
  pc2$species <- rownames(pc2)
  pc3 <- data.frame(fastBM(tree = ptree3))
  pc3$species <- rownames(pc3)
  pc.sim <- merge(pc2,pc3,by="species")
  names(pc.sim) <- c("species","PC2_body","PC3_body")
  
  pc2 <- pc.sim$PC2_body
  pc3 <- pc.sim$PC3_body
  names(pc2) <- names(pc3) <- pc.sim$species
  
  # ancestral reconstructions PCs
  rec.pc2sim <- fastAnc(ptree2,pc.sim$PC2_body)
  rec.pc3sim <- fastAnc(ptree3,pc.sim$PC3_body)
  node.traits <- data.frame("Nodes"=names(rec.pc2sim),cbind(rec.pc2sim,rec.pc3sim))
  
  node.pairs <- data.frame(ptree2$edge,ptree2$edge.length)
  names(node.pairs) <- c("Node.Anc","Node.Des","edge.length")
  node.traitsA <- node.traits[,c(1:3)]
  names(node.traitsA) <- c("Node.Anc","pc2.Anc","pc3.Anc")
  node.pairs <- merge(node.pairs,node.traitsA,by="Node.Anc")
  node.traitsD <- node.traits[,c(1:3)]
  names(node.traitsD) <- c("Node.Des","pc2.Des","pc3.Des")
  node.pairs <- merge(node.pairs,node.traitsD,by="Node.Des")
  # Save ASR of beak simulated traits
  save(node.traits,node.pairs,pc.sim,file=paste0(outFolder,"AncestralReconstructionBody_RateHet_",sim_i,".rda"))
  cat("Body simulated: Tree",sim_i,"\n")
} # end focal tree

# 2 # Calculate distance to centroid

Nsim <- 100
ratio.BodySimCen <- numeric(Nsim)

for(sim_i in 1:Nsim){
  
  # Load results from simulated reconstructions.
  load(file=paste0(outFolder,"AncestralReconstructionBody_RateHet_",sim_i,".rda"))
  
  # Distance to centroid
  centroidSpp <- data.frame("pc2"=mean(pc.sim$PC2_body),"pc3"=mean(pc.sim$PC3_body))
  node.traits$distCen <- pointDistance(node.traits[,c("rec.pc2sim","rec.pc3sim")],centroidSpp,F) # Distance between Sps and centroid
  # Merge distance metrics with Node pairs.
  node.traitsA <- node.traits[,c("Nodes","distCen")]
  names(node.traitsA) <- c("Node.Anc","distCen.Anc")
  node.pairs <- merge(node.pairs,node.traitsA,by="Node.Anc")
  node.traitsD <- node.traits[,c("Nodes","distCen")]
  names(node.traitsD) <- c("Node.Des","distCen.Des")
  node.pairs <- merge(node.pairs,node.traitsD,by="Node.Des")
  # Distance change centroid
  node.pairs$dif.CenDist <- node.pairs$distCen.Des - node.pairs$distCen.Anc
  node.pairs$dif.CenSign <- ifelse(node.pairs$dif.CenDist>0,"Reduced (closer Arch)","Increased (dist Arch)")
  ratio.BodySimCen[sim_i] <- table(node.pairs$dif.CenSign)[2]/nrow(node.pairs) 
  cat("Tree",sim_i,"ratioCen:",ratio.BodySimCen[sim_i],"\n")
  
} # End for tree "sim_i"

# Save directionality of beak simulated traits
save(ratio.BodySimCen,file=paste0(outFolder,"BodyDirectionEvolved_RateHet_Simulated.rda"))

## Summarize results:

ratio.BeakObsCen # Observed beak
mean(ratio.BeakSimCen) # Simulated beak
quantile(ratio.BeakSimCen,c(0.025,0.975))

ratio.BodyObsCen # Observed beak
mean(ratio.BodySimCen) # Simulated beak
quantile(ratio.BodySimCen,c(0.025,0.975))

## End script 5.1 ##