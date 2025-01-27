## 4.4 ## Figure S5: Plots for 3-5 archetypes (Across PC123)

## Open archetypes Beak PC23 (To save 2D archetypes)
load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Beaks.rda"))
archs3_beak_PC23 <- data.frame(archCoords[[3]])

# Open archetypes Beak PC123 (3D beak archetypes)
load(paste0(OutputFolder,"out4_ArchetypesObs_RobustAlg_PC123_Beaks.rda"))
# 4 archetypes
archs4_beak_PC123 <- data.frame(archCoords[[4]])
# 5 archetypes
archs5_beak_PC123 <- data.frame(archCoords[[5]])

# Compute convex hulls for each dataset
archs3_beak_hull <- archs3_beak_PC23[chull(archs3_beak_PC23$PC2_beak, archs3_beak_PC23$PC3_beak), ]
archs4_beak_hull <- archs4_beak_PC123[chull(archs4_beak_PC123$PC2_beak, archs4_beak_PC123$PC3_beak), ]
archs5_beak_hull <- archs5_beak_PC123[chull(archs5_beak_PC123$PC2_beak, archs5_beak_PC123$PC3_beak), ]

# Plot 4 archetypes (3D) over 3 archetypes (in the 2D plane):
pBeakPC123_toPC23 <- ggplot(archs4_beak_PC123, aes(x=PC2_beak, y=PC3_beak)) + theme_classic() +
  geom_point(data = archs3_beak_PC23, mapping = aes(x = PC2_beak, y = PC3_beak), color="black", size =4) + 
  geom_point(data = archs4_beak_PC123, mapping = aes(x = PC2_beak, y = PC3_beak), color="blue", size =4) + 
  geom_point(data = archs5_beak_PC123, mapping = aes(x = PC2_beak, y = PC3_beak), color="red", size =4) +
  geom_polygon(data = archs3_beak_hull, aes(x = PC2_beak, y = PC3_beak), fill = NA, colour = "black") +
  geom_polygon(data = archs4_beak_hull, aes(x = PC2_beak, y = PC3_beak), fill = NA, colour = "blue") +
  geom_polygon(data = archs5_beak_hull, aes(x = PC2_beak, y = PC3_beak), fill = NA, colour = "red") +
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12))+
  labs(x = "PC2 (Beak shape)",y = "PC3 (Beak shape)")

pBeakPC123_toPC23

### Body shape (PC123): ####

# Open archetypes Body PC23 (To save 2D archetypes)
load(paste0(OutputFolder,"out2_ArchetypesObs_RobustAlg_PC23_Body.rda"))
archs3_body_PC23 <- data.frame(archCoords[[3]])
# Open archetypes Body PC123 (3D body archetypes)
load(paste0(OutputFolder,"out4_ArchetypesObs_RobustAlg_PC123_Body.rda"))
# 4 archs
archs4_body_PC123 <- data.frame(archCoords[[4]])
archs4_body_PC123 <- archs4_body_PC123[c(1,4,2,3),]
# 5 archs
archs5_body_PC123 <- data.frame(archCoords[[5]])
archs5_body_PC123 <- archs5_body_PC123[c(2,3,1,4,5),]

# Compute convex hulls for each dataset
archs3_hull <- archs3_body_PC23[chull(archs3_body_PC23$PC2_body, archs3_body_PC23$PC3_body), ]
archs4_hull <- archs4_body_PC123[chull(archs4_body_PC123$PC2_body, archs4_body_PC123$PC3_body), ]
archs5_hull <- archs5_body_PC123[chull(archs5_body_PC123$PC2_body, archs5_body_PC123$PC3_body), ]

# Plot 4 archetypes (3D) over 3 archetypes (in the 2D plane):
pBodyPC123_toPC23 <- ggplot(archs4_body_PC123, aes(x=PC2_body, y=PC3_body)) + theme_classic() +
  geom_point(data = archs3_body_PC23, mapping = aes(x = PC2_body, y = PC3_body), color="black", size =4) + 
  geom_point(data = archs4_body_PC123, mapping = aes(x = PC2_body, y = PC3_body), color="blue", size =4) + 
  geom_point(data = archs5_body_PC123, mapping = aes(x = PC2_body, y = PC3_body), color="red", size =4) + 
  geom_polygon(data = archs3_hull, aes(x = PC2_body, y = PC3_body), fill = NA, colour = "black") +
  geom_polygon(data = archs4_hull, aes(x = PC2_body, y = PC3_body), fill = NA, colour = "blue") +
  geom_polygon(data = archs5_hull, aes(x = PC2_body, y = PC3_body), fill = NA, colour = "red") +
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12))+
  labs(x = "PC2 (Body shape)",y = "PC3 (Body shape)")

pBodyPC123_toPC23

## Arrange plots:
plot_grid(pBeakPC123_toPC23, pBodyPC123_toPC23, labels=c("A","B"), nrow = 1)

# Save the plot grid as a PDF
ggsave(paste0(FigureFolder,"Figure_S5.pdf"), device = "pdf", width = 8, height = 4)

## End script 4.4 ##