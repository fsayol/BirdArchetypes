## 4.2 ## Beak and Body archetype analysis (Using three dimensions from PC123)

outFolder1 <- paste0(OutputFolder,"out4_PCBeakArchetypes_PC123/")
if(!file.exists(outFolder1)) dir.create(outFolder1)

outFolder2 <- paste0(OutputFolder,"out4_PCBodyArchetypes_PC123/")
if(!file.exists(outFolder2)) dir.create(outFolder2)

## 1A ## Run archetype analysis 3d (Beak, PC123) ####

pc.dat123 <- read.csv(paste0(DataFolder,"data_tmp/BeakPC123.csv"))
pc.dat123 <- pc.dat123[,-1]

Nrep <- 100 
Nkmin<-1
Nk <- 5
Kvals<-Nkmin:Nk
maxTrials<-10

for(i in 1:length(Kvals)){
  
  # aim to obtain Nrep estimates
  for(ii in 1:Nrep){
    
    #only run the model if the k by rep combination does not yet exist 
    outfile<-paste0(outFolder1,"ArchetypesObs_PC123Traits_RobustAlg_k_",
                    Kvals[i],"_rep_",ii,"_Beaks.rda")
    exF<-paste0(outFolder1,list.files(outFolder1))
    
    if(length(intersect(outfile,exF))==0){
      
      for(iii in 1:maxTrials){
        print(paste0(i,"_",ii,"_",iii))
        #run may fail to converge
        asAR<-NULL
        try(asAR<-robustArchetypes(data=pc.dat123, k=Kvals[i]))
        
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

## 1B ## Extract rss scores and archetypes Beak ####

rssMat<-matrix(nrow=length(Kvals),ncol=Nrep)
for(i in 1:length(Kvals)){
  
  # aim to obtain Nrep estimates
  for(ii in 1:Nrep){
    
    outfile<-paste0(outFolder1,"ArchetypesObs_PC123Traits_RobustAlg_k_",
                    Kvals[i],"_rep_",ii,"_Beaks.rda")
    exF<-paste0(outFolder1,list.files(outFolder1))
    
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

archCoords<-list()
for(i in 1:length(Kvals)){
  ii<-which.min(rssMat[i,])
  outfile<-paste0(outFolder1,"ArchetypesObs_PC123Traits_RobustAlg_k_",
                  Kvals[i],"_rep_",ii,"_Beaks.rda")
  (load(outfile))
  archCoords[[i]]<-arch
} 

save(rssMat,rssMatC,archCoords,file=paste0(OutputFolder,"out4_ArchetypesObs_RobustAlg_PC123_Beaks.rda"))


## 2A ## Run archetype analysis 3d (Body, PC123) ####

pc.dat123<-read.csv(paste0(DataFolder,"data_tmp/BodyPC123.csv"))
pc.dat123<-pc.dat123[,-1]

Nrep <- 100 
Nkmin<-1
Nk <- 5
Kvals<-Nkmin:Nk
maxTrials<-10

for(i in 1:length(Kvals)){
  
  # aim to obtain Nrep estimates
  for(ii in 1:Nrep){
    
    #only run the model if the k by rep combination does not yet exist 
    outfile<-paste0(outFolder2,"ArchetypesObs_PC123Traits_RobustAlg_k_",
                    Kvals[i],"_rep_",ii,"_Body.rda")
    exF<-paste0(outFolder2,list.files(outFolder2))
    
    if(length(intersect(outfile,exF))==0){
      
      for(iii in 1:maxTrials){
        print(paste0(i,"_",ii,"_",iii))
        #run may fail to converge
        asAR<-NULL
        try(asAR<-robustArchetypes(data=pc.dat123, k=Kvals[i]))
        
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

## 2B ## Extract rss scores and archetypes Body ####

rssMat<-matrix(nrow=length(Kvals),ncol=Nrep)
for(i in 1:length(Kvals)){
  
  # aim to obtain Nrep estimates
  for(ii in 1:Nrep){
    
    outfile<-paste0(outFolder2,"ArchetypesObs_PC123Traits_RobustAlg_k_",
                    Kvals[i],"_rep_",ii,"_Body.rda")
    exF<-paste0(outFolder2,list.files(outFolder2))
    
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

archCoords<-list()
for(i in 1:length(Kvals)){
  ii<-which.min(rssMat[i,])
  outfile<-paste0(outFolder2,"ArchetypesObs_PC123Traits_RobustAlg_k_",
                  Kvals[i],"_rep_",ii,"_Body.rda")
  (load(outfile))
  archCoords[[i]]<-arch
} 

save(rssMat,rssMatC,archCoords,file=paste0(OutputFolder,"out4_ArchetypesObs_RobustAlg_PC123_Body.rda"))

## End script 4.2 ##