## 1.0 ## Download external datasets

## 1a ## Avonet (Tobias et al. 2022): ####

url <- "https://figshare.com/ndownloader/files/34480856?private_link=b990722d72a26b5bfead"
# Specify the destination of the file
dest_file <- paste0(DataTmpFolder,"Avonet_dataset_s1.xlsx")
# Download the file
download.file(url, destfile = dest_file, mode = "wb")
# Read xlsx and select BirdTree sheet
avonet <- read_excel(paste0(DataTmpFolder,"Avonet_dataset_s1.xlsx"),sheet=4)
# Select variables and add underscore in species names:
avonet$species <- gsub(" ","_",avonet$Species3)
names(avonet)

# Add NA's for traits that have been inferred:
avonet$Traits.inferred
# Beak traits:
avonet[(str_detect(avonet$Traits.inferred, "Beak Culmen")),]$Beak.Length_Culmen <- NA
avonet[(str_detect(avonet$Traits.inferred, "Beak Width")),]$Beak.Width <- NA
avonet[(str_detect(avonet$Traits.inferred, "Beak Depth")),]$Beak.Depth <- NA
# Body traits:
avonet[(str_detect(avonet$Traits.inferred, "Tarsus Length")),]$Tarsus.Length <- NA
avonet[(str_detect(avonet$Traits.inferred, "Wing Length")),]$Wing.Length <- NA
avonet[(str_detect(avonet$Traits.inferred, "Kipp's Distance")),]$Kipps.Distance <- NA
avonet[(str_detect(avonet$Traits.inferred, "Tail Length")),]$Tail.Length <- NA

# Save avonet dataset (BirdTree version):
avonet.sel <- avonet[,c("species","Beak.Length_Culmen","Beak.Width","Beak.Depth","Tarsus.Length","Wing.Length","Kipps.Distance","Tail.Length")]
write.csv(avonet.sel,file=paste0(DataTmpFolder,"Avonet_BirdTree.csv"),row.names=F)

## 1b ## Beak scans (Cooney et al. 2017) ####

url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature21074/MediaObjects/41586_2017_BFnature21074_MOESM80_ESM.csv"
Cooney2017_data <- read.csv(url)
head(Cooney2017_data)
names(Cooney2017_data)[1] <- "species"
pc.dat12 <- Cooney2017_data[,c("species","PC1","PC2")]
# Scale and invert to facilitate interpretation (i.e., Larger +)
pc.dat12$PC1<- scale(pc.dat12$PC1)*(-1)
pc.dat12$PC2<- scale(pc.dat12$PC2)*(-1)
write.csv(pc.dat12,file=paste0(DataTmpFolder,"BeakScanPC12.csv"),row.names=F)
