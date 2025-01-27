## 2.2 ## Run robust archetypes for beak shape (irregular polygons)

outFolder <- paste0(OutputFolder,'out2_PCBeakArchetypes/')
if(!file.exists(outFolder)) dir.create(outFolder)

## Run archetype analysis (Beak shape)
pc.dat23<-read.csv(paste0(DataFolder,"data_tmp/BeakPC23.csv"))
pc.dat23<-pc.dat23[,-1]

Nrep <- 100 
Nkmin<-1
Nk <- 10
Kvals<-Nkmin:Nk
maxTrials<-10

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

save(rssMat,rssMatC,archCoords,file=paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Beaks.rda"))

## End script 2.2 ##
