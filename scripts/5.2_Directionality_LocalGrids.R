## 5.2 ## Directionality of trait evolution (Local analysis)

outFolder <-paste0(OutputFolder,'out5_ASR_Empirical_and_Null_RateHet/')
# Create folders to save node data:
Nodes_Output <-paste0(OutputFolder,'out5_NodeTraits/')
if(!file.exists(Nodes_Output)) dir.create(Nodes_Output)

# Create grid for observed beak:
load(paste0(outFolder,"AncestralReconstructionBeakTraits_RateHet_Observed.rda"))
cellIntersectDatObs<-cellIntersect(node.pairs,nCells=30,toplot=F,printRow=F)
cellTabO<-makeCellTab(node.pairs,cellIntersectDatObs[[1]],cellIntersectDatObs[[4]])
save(cellIntersectDatObs,cellTabO,file=paste0(Nodes_Output,"BeakNodeTraits_Obs_cellIntersections.rda"))

# Create grid for observed body traits:
load(paste0(outFolder,"AncestralReconstructionBodyTraits_RateHet_Observed.rda"))
cellIntersectDatObs<-cellIntersect(node.pairs,nCells=30,toplot=FALSE)
cellTabO<-makeCellTab(node.pairs,cellIntersectDatObs[[1]],cellIntersectDatObs[[4]])
save(cellIntersectDatObs,cellTabO,file=paste0(Nodes_Output,"BodyNodeTraits_Obs_cellIntersections.rda"))

## Null models: ##

Nsim <- 100

## 1A ## Calculate Simulated cellIntersect for beak nodes: ####

# Loop 1/2 beak
for(sim_i in 1:Nsim){
  load(paste0(outFolder,"AncestralReconstructionBeak_RateHet_",sim_i,".rda"))
  cellIntersectDatNull <- cellIntersect(node.pairs,nCells=30,toplot=F,printRow=F)
  cellTabNull <- makeCellTab(node.pairs,cellIntersectDatNull[[1]],cellIntersectDatNull[[4]])
  save(cellIntersectDatNull,cellTabNull,file=paste0(Nodes_Output,"BeakNodeTraits_SimRH_",sim_i,"_cellTabNull.rda"))
  print(sim_i)
} # End loop 1/2 Beak

## 1B ## Calculate Simulated cellIntersect for beak nodes: ####

# Loop 1/2 Body
for(sim_i in 1:Nsim){
  load(paste0(outFolder,"AncestralReconstructionBody_RateHet_",sim_i,".rda"))
  cellIntersectDatNull<-cellIntersect(node.pairs,nCells=30,toplot=F,printRow=F)
  cellTabNull<-makeCellTab(node.pairs,cellIntersectDatNull[[1]],cellIntersectDatNull[[4]])
  save(cellIntersectDatNull,cellTabNull,file=paste0(Nodes_Output,"BodyNodeTraits_SimRH_",sim_i,"_cellTabNull.rda"))
  print(sim_i)
} # End Loop 1/2 Body

## 2A ## Compare directionality with null models (Beak) ####

# Loop 2/2 Beak
cellTabNullC<-data.frame()
for(sim_i in 1:Nsim){
  load(paste0(Nodes_Output,"BeakNodeTraits_SimRH_",sim_i,"_cellTabNull.rda"))
  cellTabNull$sim_i <- sim_i
  cellTabNullC <- rbind(cellTabNullC,cellTabNull)
  print(sim_i)
}  # End loop 2/2 Beak

load(paste0(Nodes_Output,"BeakNodeTraits_Obs_cellIntersections.rda"))
# use quantile regression to identify vectors of more consistent direction than expected under null model
qmod <- rq(cellTabNullC$xdegV ~ log(cellTabNullC$nb),tau = 0.05)
cellTabO$xdegVPred <- summary(qmod)$coef[1,1]+(summary(qmod)$coef[2,1]*log(cellTabO$nb))
cellTabO$RejectNull <- 0
cellTabO$RejectNull[which((cellTabO$xdegVPred-cellTabO$xdegV)>0)] <- 1
table(cellTabO$RejectNull)
save(cellIntersectDatObs,cellTabO,file=paste0(Nodes_Output,"BeakNodeTraits_Null_cellIntersections.rda"))

## 2B ## Compare directionality with null models (Body) ####

# Start Loop 2/2 Body
cellTabNullC<-data.frame()
for(sim_i in 1:Nsim){
  load(paste0(Nodes_Output,"BodyNodeTraits_SimRH_",sim_i,"_cellTabNull.rda"))
  cellTabNull$sim_i<-sim_i
  cellTabNullC<-rbind(cellTabNullC,cellTabNull)
  print(sim_i)
} # End Loop 2/2 Body

load(paste0(Nodes_Output,"BodyNodeTraits_Obs_cellIntersections.rda"))
# use quantile regression to identify vectors of more consistent direction than expected under null model 
qmod<-rq(cellTabNullC$xdegV ~ log(cellTabNullC$nb),tau = 0.05)
cellTabO$xdegVPred <- summary(qmod)$coef[1,1]+(summary(qmod)$coef[2,1]*log(cellTabO$nb))
cellTabO$RejectNull<-0;cellTabO$RejectNull[which((cellTabO$xdegVPred-cellTabO$xdegV)>0)]<-1
table(cellTabO$RejectNull)
save(cellIntersectDatObs,cellTabO,file=paste0(Nodes_Output,"BodyNodeTraits_Null_cellIntersections.rda"))

## End script 5.2 ##