## 6.1 ## Speciation and extinction rates across morphospace

# Create folder to save Extinction risk and DR models
outFolder <-paste0(OutputFolder,'out6_ER_DR_Models/')
if(!file.exists(outFolder)) dir.create(outFolder)

## Create dataset of ER, DR and confounds ####

pc.dat23beak <- read.csv(paste0(DataFolder,"data_tmp/BeakPC23.csv"))
pc.dat23body <- read.csv(paste0(DataFolder,"data_tmp/BodyPC23.csv"))
fdat <- merge(pc.dat23beak,pc.dat23body,by="species")

# Beak shape metrics
dm <- read.csv(paste0(DataFolder,"data_tmp/BeakDistanceMetrics.csv"))
# Distance to minimum archetype:
fdat$distMinAR_Beak<-dm$distMinAR[match(fdat$species,dm$species)]
# Body shape metrics:
dm <- read.csv(paste0(DataFolder,"data_tmp/BodyDistanceMetrics.csv"))
# Distance to minimum archetype:
fdat$distMinAR_Body<-dm$distMinAR[match(fdat$species,dm$species)]

# IUCN Threat data (+ confounds)
data.mod <- read.csv(paste0(DataFolder,"Data5_ExtinctionRisk_Confounds.csv"))
names(data.mod)[1] <- "species"
data.mod <- merge(fdat,data.mod,by="species")
data.mod <- data.mod[data.mod$status2022!="DD",]
data.mod$er.bin <-1
data.mod$er.bin[which(data.mod$status2022 %in% c("LC"))] <- 0

dr.data <- read.csv(paste0(DataFolder,"Data6_Diversification_DR.csv"))
dr.data$log.dr <- log(dr.data$dr)
data.mod <- merge(data.mod,dr.data,by="species")
data.mod$dr.bin <- ifelse(data.mod$log.dr>quantile(data.mod$log.dr,0.75),yes = 1,no=0)

# Log-transform bodyMass and Generation length
data.mod$BodyMassLog<-log(data.mod$BodyMass)
data.mod$GenerationLengthLog<-log(data.mod$GenerationLength)
write.table(data.mod,file = paste0(outFolder,"ExtinctionModelsInput.csv"))

## Run GAM Models with extinction risk and DR across morphospace:

# Folder to save Extinction risk and DR models (GAM models)
outFolder <-paste0(OutputFolder,'out6_ER_DR_Models/')
# Nodes ouput (for vector grids)
Nodes_Output <-paste0(OutputFolder,'out5_NodeTraits/')

# Open data extinction models
data.mod <- read.table(paste0(outFolder,"ExtinctionModelsInput.csv"))
names(data.mod)

# Add data on Invertivore gleaner specialist:
for.niche<- read.csv(file=paste0(DataFolder,"Data1_Foraging_Feeding.csv")) # Data main foraging niches
for.niche$InvGlean_specialist <- ifelse(for.niche$ForagingNiche_main=="Invertivore_glean_elevated",yes=1,no=0)
data.mod <- merge(data.mod,for.niche[,c("species","InvGlean_specialist")],by="species")

# Add data on beak task specialization:
beak.task <- read.csv(paste0(DataFolder,"Data2_BeakTaskEnrichment.csv"))
beak.task$beakTask_specialist <- apply(beak.task[,c("Crush_task","Engulf_task","Reach_task")], 1, function(row) any(row >= 60))
beak.task$beakTask_specialist <- as.numeric(beak.task$beakTask_specialist)
data.mod <- merge(data.mod,beak.task[,c("species","beakTask_specialist")],by="species")

# Add data on body task specialization
body.task <- read.csv(paste0(DataFolder,"Data3_BodyTaskEnrichment.csv"))
body.task$bodyTask_specialist <- apply(body.task[,c("Fly_task","Walk_task","Swim_task")], 1, function(row) any(row >= 60))
body.task$bodyTask_specialist <- as.numeric(body.task$bodyTask_specialist)
data.mod <- merge(data.mod,body.task[,c("species","bodyTask_specialist")],by="species")

## Prepare morphospace grid (from vector map) to predict values:

# Open beak morphospace grid
load(file=paste0(Nodes_Output,"BeakNodeTraits_Obs_cellIntersections.rda"))

# PC2 beak (x values)
xvals<-sort(unique(cellTabO$X))
xdiff<-diff(sort(unique(cellTabO$X)))[1]
xmax<-max(seq(max(xvals),ceiling(max(data.mod$PC2_beak)),xdiff))
xmin<-max(seq(abs(min(xvals)),abs(floor(min(data.mod$PC2_beak))),xdiff))*-1
PC2_beak<-seq(xmin,xmax,xdiff)

# PC3 beak (y values)
yvals<-sort(unique(cellTabO$Y))
ydiff<-diff(sort(unique(cellTabO$Y)))[1]
ymax<-max(seq(max(yvals),ceiling(max(data.mod$PC3_beak)),ydiff))
ymin<-max(seq(abs(min(yvals)),abs(floor(min(data.mod$PC3_beak))),ydiff))*-1
PC3_beak<-seq(ymin,ymax,ydiff)

# Open body morphospace grid:
load(file=paste0(Nodes_Output,"BodyNodeTraits_Obs_cellIntersections.rda"))

# PC2 body (x values)
xvals<-sort(unique(cellTabO$X))
xdiff<-diff(sort(unique(cellTabO$X)))[1]
xmax<-max(seq(max(xvals),ceiling(max(data.mod$PC2_body)),xdiff))
xmin<-max(seq(abs(min(xvals)),abs(floor(min(data.mod$PC2_body))),xdiff))*-1
PC2_body<-seq(xmin,xmax,xdiff)

# PC3 body (y values)
yvals<-sort(unique(cellTabO$Y))
ydiff<-diff(sort(unique(cellTabO$Y)))[1]
ymax<-max(seq(max(yvals),ceiling(max(data.mod$PC3_body)),ydiff))
ymin<-max(seq(abs(min(yvals)),abs(floor(min(data.mod$PC3_body))),ydiff))*-1
PC3_body<-seq(ymin,ymax,ydiff)

# Calculate number of species for beak grid:
beak_grid <- expand.grid(PC2_beak=PC2_beak,PC3_beak=PC3_beak)
beak_grid$spp <- 0
for (i in 1:nrow(data.mod)){
  # Select species
  pc2i <- data.mod$PC2_beak[i]
  pc3i <- data.mod$PC3_beak[i]
  # Select closest cell coordinates
  sel.pc2 <- which(abs(beak_grid$PC2_beak-pc2i)==min(abs(beak_grid$PC2_beak-pc2i)))
  sel.pc3 <- which(abs(beak_grid$PC3_beak-pc3i)==min(abs(beak_grid$PC3_beak-pc3i)))
  # Select row and add +1 in species count
  row.i <- sel.pc2[sel.pc2 %in% sel.pc3]
  beak_grid[row.i,]$spp <- beak_grid[row.i,]$spp+1
} # end for i

## Save body grid cell intersections:
save(cellIntersectDatObs,cellTabO,file=paste0(outFolder,"BeakNodeTraits_cellIntersections.rda"))

# Calculate number of species for body grid:
body_grid <- expand.grid(PC2_body=PC2_body,PC3_body=PC3_body)
body_grid$spp <- 0
for (i in 1:nrow(data.mod)){
  # Select species
  pc2i <- data.mod$PC2_body[i]
  pc3i <- data.mod$PC3_body[i]
  # Select closest cell coordinates
  sel.pc2 <- which(abs(body_grid$PC2_body-pc2i)==min(abs(body_grid$PC2_body-pc2i)))
  sel.pc3 <- which(abs(body_grid$PC3_body-pc3i)==min(abs(body_grid$PC3_body-pc3i)))
  # Select row and add +1 in species count
  row.i <- sel.pc2[sel.pc2 %in% sel.pc3]
  body_grid[row.i,]$spp <- body_grid[row.i,]$spp+1
} # end for i

## Run GAM models (Extinction) ####

names(data.mod)

i <- 6
btF <- gam(er.bin~te(PC2_beak,PC3_beak,k=i)+te(PC2_body,PC3_body,k=i)+
             GenerationLengthLog+InsularityIndex+HabitatBreadth,data=data.mod,family=binomial())
summary(btF)

# Predict extinction of beak PCs:
newdata<-expand.grid(PC2_beak=PC2_beak,PC3_beak=PC3_beak,
                     PC2_body=0,PC3_body=0,
                     GenerationLengthLog=0,InsularityIndex=0,HabitatBreadth=0)
bt <- gam(er.bin~te(PC2_beak,PC3_beak,k=i),data=data.mod,family=binomial())
summary(bt)
pred<-predict(bt,newdata,type="response", se.fit = T)
newdata$out<-pred[[1]]
newdata$out.se<-pred[[2]]
pred<-predict(btF,newdata,type="response", se.fit = T)
newdata$outF<-pred[[1]]
newdata$outF.se<-pred[[2]]
newdata <- cbind(newdata,"spp"=beak_grid[,3])

write.csv(newdata,file=paste0(outFolder,"GAMBeakExtinctionModelPrediction.csv"))

# Predict extinction of body PCs:
newdata<-expand.grid(PC2_beak=0,PC3_beak=0,
                     PC2_body=PC2_body,PC3_body=PC3_body,
                     GenerationLengthLog=0,InsularityIndex=0,HabitatBreadth=0)
bt <- gam(er.bin~te(PC2_body,PC3_body,k=i),data=data.mod,family=binomial())
summary(bt)
pred<-predict(bt,newdata,type="response", se.fit = T)
newdata$out<-pred[[1]]
newdata$out.se<-pred[[2]]
pred<-predict(btF,newdata,type="response", se.fit = T)
newdata$outF<-pred[[1]]
newdata$outF.se<-pred[[2]]
newdata <- cbind(newdata,"spp"=body_grid[,3])

write.csv(newdata,file=paste0(outFolder,"GAMBodyExtinctionModelPrediction.csv")) 

## GAM Models with DR ####

# New data for beak
newdata<-expand.grid(PC2_beak=PC2_beak,PC3_beak=PC3_beak)

# Gam model (DR ~ beak)
i <- 6
bt <- gam(dr.bin~te(PC2_beak,PC3_beak,k=i),data=data.mod,family=binomial())
summary(bt)
pred<-predict(bt,newdata,type="response", se.fit = T)
newdata$out<-pred[[1]]
newdata$out.se<-pred[[2]]
newdata <- cbind(newdata,"spp"=beak_grid[,3])
write.csv(newdata,file=paste0(outFolder,"GAM_BeakModelPrediction_DR.csv"))

# New data for body
newdata<-expand.grid(PC2_body=PC2_body,PC3_body=PC3_body)

# Gam model (DR ~ body)
i <- 6
bt <- gam(dr.bin~te(PC2_body,PC3_body,k=i),data=data.mod,family=binomial())
summary(bt)
pred<-predict(bt,newdata,type="response", se.fit = T)
newdata$out<-pred[[1]]
newdata$out.se<-pred[[2]]
newdata <- cbind(newdata,"spp"=body_grid[,3])

write.csv(newdata,file=paste0(outFolder,"GAM_BodyModelPrediction_DR.csv"))


## GAM Model Beak task specialization ####

# New data for beak
newdata <- expand.grid(PC2_beak=PC2_beak,PC3_beak=PC3_beak)

# Gam model (beak.MaxTask ~ beak)
i<-6
bt <- gam(beakTask_specialist~te(PC2_beak,PC3_beak,k=i),data=data.mod,family=binomial())
summary(bt)
pred<-predict(bt,newdata,type="response", se.fit = T)
newdata$out<-pred[[1]]
newdata$out.se<-pred[[2]]
newdata <- cbind(newdata,"spp"=beak_grid[,3])

write.csv(newdata,file=paste0(outFolder,"GAM_BeakModelPrediction_BeakTask.csv"))

# New data for body
newdata<-expand.grid(PC2_body=PC2_body,PC3_body=PC3_body)

# Gam model (bodyTaskSpecialist ~ body morphoscape)
i<-6
bt <- gam(bodyTask_specialist~te(PC2_body,PC3_body,k=i),data=data.mod,family=binomial())
summary(bt)
pred<-predict(bt,newdata,type="response", se.fit = T)
newdata$out<-pred[[1]]
newdata$out.se<-pred[[2]]
newdata <- cbind(newdata,"spp"=body_grid[,3])

write.csv(newdata,file=paste0(outFolder,"GAM_BodyModelPrediction_BodyTask.csv"))

## GAM Model Invert. glean specialist ####

# New data for beak
newdata<-expand.grid(PC2_beak=PC2_beak,PC3_beak=PC3_beak)

# Gam model (InvertGleanSpec ~ beak morphospace)
i<-6
bt <- gam(InvGlean_specialist~te(PC2_beak,PC3_beak,k=i),data=data.mod,family=binomial())
summary(bt)
pred<-predict(bt,newdata,type="response", se.fit = T)
newdata$out<-pred[[1]]
newdata$out.se<-pred[[2]]
newdata <- cbind(newdata,"spp"=beak_grid[,3])

write.csv(newdata,file=paste0(outFolder,"GAM_BeakModelPrediction_InvertGleanSpec.csv"))

# New data for body
newdata<-expand.grid(PC2_body=PC2_body,PC3_body=PC3_body)

# Gam model (InvertGleanSpec ~ body morphospace)
i<-6
bt <- gam(InvGlean_specialist~te(PC2_body,PC3_body,k=i),data=data.mod,family=binomial())
summary(bt)
pred<-predict(bt,newdata,type="response", se.fit = T)
newdata$out<-pred[[1]]
newdata$out.se<-pred[[2]]
newdata <- cbind(newdata,"spp"=body_grid[,3])
write.csv(newdata,file=paste0(outFolder,"GAM_BodyModelPrediction_InvertGleanSpec.csv"))

## End script 6.1 ##