## 6.2 ## Model extinction risk as a function of distance to archetypes

outFolder <-paste0(OutputFolder,'out6_ER_DR_Models/')

## fit models of extinction risk (Using BPMM) ####

data.mod <- read.table(paste0(outFolder,"ExtinctionModelsInput.csv"))
names(data.mod)

#set ndata.mod#set number of iterations, burnin and sampling frequency for mcmc sampler
Nnitt = 100000 
Nburnin = Nnitt/10
Nthin = Nnitt/1000
Nnittcorr <- Nnitt+Nburnin

## MCMCglmm models (with binomial response)
tree <- read.nexus(paste0(DataFolder,"phylo/AllBirdsHackett1_summary.tre"))

# Omit NA's in data:
names(data.mod)
data.sel <- na.omit(data.mod)

#prune the tree to only contain species in the dataset	
todrop = setdiff(tree$tip.label,data.sel$species)
length(tree$tip.label)
phylo.tree <- drop.tip(phy = tree, tip = todrop)

#set non-informative prior
prior_phylo = list(R = list(V = 1, fix=1), G = list(G1 = list(V = 1,nu = 0.002)))

# find inverse matrix of tree
Ainv  =  inverseA(phylo.tree)$Ainv

## Run MCMC models

# Extinction
phylo_mod1 = MCMCglmm(er.bin ~ distMinAR_Beak + distMinAR_Body + GenerationLengthLog + HabitatBreadth + InsularityIndex,
                      random = ~species, ginverse = list(species=Ainv), data = data.sel,
                      nitt=Nnittcorr, burnin=Nburnin, thin=Nthin, prior=prior_phylo, family = "ordinal",verbose=F)
summary(phylo_mod1)

# Diversification rate (DR)
phylo_mod2 = MCMCglmm(log.dr ~ distMinARBeak + distMinAR_Body + GenerationLengthLog + HabitatBreadth + InsularityIndex,
                      random = ~species, ginverse = list(species=Ainv), data = data.sel,
                      nitt=Nnittcorr, burnin=Nburnin, thin=Nthin, prior=prior_phylo, family = "gaussian",verbose=F)
summary(phylo_mod2)

save(phylo_mod1,phylo_mod2,file=paste0(outFolder,"PhyloMCMC_models.Rdata"))

## End script 6.2 ##