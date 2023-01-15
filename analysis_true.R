###### Comments: ######



# -------------------------------------------------
# Beginning ####
# -------------------------------------------------
setwd("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\")

load("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")

library(raster); library(rgdal); library(maptools); library(dismo)

library(ggfortify); library(ggplot2); library(ggbiplot);  library(dplyr); library(corrplot)

require(HDMD); library(ade4); library(adegenet); library(adephylo); library(ecodist); library(phyloclim); library(ecospat); library(vegan)

library(phytools); library(phangorn)


# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")

# ------------------------------------------------
# Import and prepare the bioclimatic variables ####
# ------------------------------------------------

# Importing the raster files from which to compute the climatic 
bio01 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio1_a.grd")
bio02 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio2_a.grd")
bio03 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio3_a.grd")
bio04 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio4_a.grd")
bio05 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio5_a.grd")
bio06 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio6_a.grd")
bio07 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio7_a.grd")
bio08 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio8_a.grd")
bio09 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio9_a.grd")
bio10 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio10_a.grd")
bio11 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio11_a.grd")
bio12 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio12_a.grd")
bio13 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio13_a.grd")
bio14 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio14_a.grd")
bio15 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio15_a.grd")
bio16 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio16_a.grd")
bio17 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio17_a.grd")
bio18 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio18_a.grd")
bio19 <- raster("B:\\WorldClim\\WorldClim 1.4\\Mata_Atlantica\\Bioclim\\bio19_a.grd")



# Save environmental image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")



# Create the stack
bioclimStack <- stack(bio01, bio02, bio03, bio04, bio05, bio06, bio07, bio08, bio09, bio10, bio11, bio12, bio13, bio14, bio15, bio16, bio17, bio18, bio19)



# Compute the correletion between bioclimatic variables
correlation <- layerStats(bioclimStack, "pearson", na.rm = T)
correlationMatrix <- correlation$`pearson correlation coefficient`
colnames(correlationMatrix) <- c("Mean_Temp", "Diurnal_Range", "Isothermality (Bio2/Bio7)", "Temp_Seasonality (STDEV*100)", "Max_Temp_Warmest_Month", "Min_Temp_Warmest_Month","Temp_Annual_Range (Bio5-Bio6)", "Mean_Temp_Wettest_Quarter", "Mean_Temp_Driest_Quarter", "Mean_Temp_Warmest_Quarter", "Mean_Temp_Coldest_Quarter", "Annual_Precip", "Precip_Wettest_Month", "Precip_Driest_Month", "Precip_Seasonality", "Precip_Wettest_Quarter", "Precip_Driest_Quarter", "Precip_Warmest_Quarter", "Precip_Coldest_Quarter")
rownames(correlationMatrix) <- c("Mean_Temp", "Diurnal_Range", "Isothermality (Bio2/Bio7)", "Temp_Seasonality (STDEV*100)", "Max_Temp_Warmest_Month", "Min_Temp_Warmest_Month","Temp_Annual_Range (Bio5-Bio6)", "Mean_Temp_Wettest_Quarter", "Mean_Temp_Driest_Quarter", "Mean_Temp_Warmest_Quarter", "Mean_Temp_Coldest_Quarter", "Annual_Precip", "Precip_Wettest_Month", "Precip_Driest_Month", "Precip_Seasonality", "Precip_Wettest_Quarter", "Precip_Driest_Quarter", "Precip_Warmest_Quarter", "Precip_Coldest_Quarter")

correlation.row.names.alternative <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")

colnames(correlationMatrix) <- correlation.row.names.alternative
rownames(correlationMatrix) <- correlation.row.names.alternative


# Plot the correlatiogram
corrplot.mixed(correlationMatrix, lower = "square", upper = "number", tl.col="black", tl.cex = 1, win.asp = 1)



# Stacking the uncorrelated variables
uncorVar <- stack(bio01, bio02, bio03, bio04, bio07, bio12, bio13, bio14, bio15, bio18)



# ------------------------------------------------
# Import the occurrence data and conduct the extraction ####
# ------------------------------------------------
# Importing occurrence points dataframe
occurrence <- read.csv("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\ocorrencias_concatenadas.csv", header = T)
head(occurrence, 12)


# Removing duplicate records
occurrenceNOTdupli <- occurrence[!duplicated(occurrence), ]
rownames(occurrenceNOTdupli) <- NULL
head(occurrenceNOTdupli, 12)

write.csv(occurrenceNOTdupli, "occurrenceNOTdupli.csv") # cbind(groups, 

occurReady <- occurrenceNOTdupli[ , -1]
occurReady <- data.frame(occurReady[,2], occurReady[,1])
colnames(occurReady) <- c("lon", "lat")
head(occurReady, 12)



# Extract values from raster for the occurrence points
extracted <- extract(uncorVar, occurReady, method = "simple")
head(extracted, 12)

colnames(extracted) <- gsub(colnames(extracted), pattern = "_a", replacement = "")

extractedWITHcoord <- cbind(occurrenceNOTdupli, extracted)
head(extractedWITHcoord, 12)



# Save environmental image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")





# ------------------------------------------------
# Calculating the z-scores (subtracting the mean and dividing it by its standard deviation) and performing the PCA ####
# ------------------------------------------------
# Verify which rows may have NAs because of bad trimming
extracted[!complete.cases(extracted),]

# Compute the means for each variable (ROWS 88 AND 116 HAVE NAs BECAUSE OF WRONG MAP TRIMMING)
means <- apply(extracted, 2, FUN = "mean")
means

# Compute the standard deviation for each variable (ROWS 88 AND 116 HAVE NAs BECAUSE OF WRONG MAP TRIMMING)
standardDev <- apply(extracted, 2, FUN = "sd")
standardDev

# Calculate the z-scores for the variables (Remember: it results the variables as lines and the occurrences as columns)
zScores <- apply(extracted, 1, scale)
head(zScores, 10)


zScoresWITHcoord <- t(zScores)
colnames(zScoresWITHcoord) <- colnames(extracted)

zScoresWITHcoord <- cbind(zScoresWITHcoord, occurrenceNOTdupli)

head(zScoresWITHcoord, 10)




# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")



# -----------------------------------
# PCA with the z-scores to normalize the different magnitudes of the variables #### GROUPS IS HERE
# -----------------------------------
# PCA
PCA <- prcomp(zScoresWITHcoord[,c(1:10)], scale. = FALSE)
PCA
summary(PCA)
plot(PCA, main = "")
biplot(PCA)


# Plotting PCA results (by clades)
matrixPCA <- as.data.frame(PCA$x)
head(matrixPCA)

#group <- zScoresWITHcoord[complete.cases(zScoresWITHcoord),]
#write.csv(group, file = "C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\group.csv")
#groups.not.dupli <- groups[!duplicated(groups),]
#write.csv(groups.not.dupli, file = "C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\groups_not_dupli.csv")


groups <- read.csv("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\groups.csv", header = T)
groups

p <- ggplot(matrixPCA, aes(x = PC1, y = PC2, color = groups$X))+
  geom_point(size = 10, alpha=0.6) + 
  theme(legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) + 
  scale_colour_manual(values=c("goldenrod", "darkolivegreen4", "hotpink4", "lightsteelblue4", "lightgray"))
p



q <- p + stat_ellipse(aes(x = PC1, y = PC2, fill = factor(groups$X)),
                      geom = "polygon", level = 0.95, alpha = 0.10) +
  guides(color = guide_legend("Clades"), fill = guide_legend("Clades"))
q




p.biplot <- ggbiplot(PCA, obs.scale = 0.5, var.scale = 1, groups = groups$X, 
                     ellipse = T, circle = F, varname.size = 3)

p.biplot + guides(color = guide_legend("Clades"), fill = guide_legend("Clades")) +
  theme_classic() + 
  scale_colour_manual(values=c("goldenrod", "darkolivegreen4", "hotpink4", "lightsteelblue4", "lightgray"))




p.zscores.biplot <- ggbiplot(PCA, 
                             obs.scale = 0.4, 
                             var.scale = 1, 
                             var.axes = T,
                             groups = groups$X, 
                             ellipse = T, 
                             ellipse.prob = .95,
                             circle = F, 
                             labels = rownames(PCA$x), 
                             labels.size = 0, 
                             varname.size = 5) + 
  guides(color = guide_legend("Clades"), fill = guide_legend("Clades"))

p.zscores.biplot +
  scale_colour_manual(values=c("goldenrod", "darkolivegreen4", "hotpink4", "lightsteelblue4", "gray")) +
  geom_point(aes(color = groups), size = 7, alpha = 0.7) +
  theme_bw() +
  theme(line = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed")








# Plotting PCA results (by groups)
p.groups <- ggplot(matrixPCA, aes(x = PC1, y = PC2, color = groups$group))+
  geom_point(size = 10, alpha=0.6) + 
  theme(legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))
p.groups


q.groups <- p.groups + stat_ellipse(aes(x = PC1, y = PC2, fill = factor(groups$group)),
                                    geom = "polygon", level = 0.95, alpha = 0.10) +
  guides(color = guide_legend("Groups"), fill = guide_legend("Groups"))

q.groups


# Biplot
q.zscores.biplot <- ggbiplot(PCA, 
                             obs.scale = 0.28, 
                             var.scale = 1, 
                             var.axes = T,
                             groups = groups$group, 
                             ellipse = T, 
                             ellipse.prob = .95,
                             circle = F, 
                             labels = rownames(PCA$x), 
                             labels.size = 0, 
                             varname.size = 5) + 
  guides(color = guide_legend("Groups"), fill = guide_legend("Groups"))

q.zscores.biplot +
  geom_point(aes(color = groups), size = 7, alpha = 0.7) +
  theme_bw() +
  theme(line = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed")



# Save environmental image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")




# -----------------------------------
# PCA WITHOUT computing the z-scores ####
# -----------------------------------
PCAfull <- prcomp(extracted, scale = TRUE) # By setting the argument scale = TRUE, we are telling R to use a correlation matrix to calculate the eingenvector and eigenvalues. We use this method of computation when our variables differ widely in their scale from one another, haveing different magnitudes. We tend to use the covariance matrix when the variable scales are similar and the correlation matrix when variables are on different scales. https://stats.stackexchange.com/questions/53/pca-on-correlation-or-covariance
PCAfull
summary(PCAfull)
plot(PCAfull, main = "")
biplot(PCAfull)

matrixPCAfull <- as.data.frame(PCAfull$x)

pfull <- ggplot(matrixPCAfull,aes(x = PC1, y = PC2, color = groups$X))+ 
  geom_point(size = 10, alpha = 0.5) + 
  theme(legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) + 
  stat_ellipse(aes(x = PC1, y = PC2, fill = factor(groups$X)),
               geom = "polygon", level = 0.95, alpha = 0.10) +
  guides(color = guide_legend("Clades"), fill = guide_legend("Clades"))+ 
  scale_colour_manual(values=c("goldenrod", "darkolivegreen4", "hotpink4", "lightsteelblue4", "lightgray"))
pfull



# Save environmental image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")




# -----------------------------------
# Calculating the k-means in the PCA ####
# -----------------------------------
# K-means computation for 5 centers
kmeansClust <- kmeans(PCA$x, centers = 5)
kmeansClust

PCAmatrixWithKMeans <- cbind(PCA$x, cluster = factor(kmeansClust$cluster), groups = groups)

head(PCAmatrixWithKMeans)

kmean.plot <- ggplot(PCAmatrixWithKMeans, aes(x = PC1, y = PC2, shape = cluster, color = groups.X)) + 
  geom_point(size = 8, alpha = 0.5)  + 
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))
kmean.plot



# K-means computation for 4 centers
kmeansClust4 <- kmeans(PCA$x, centers = 4)
kmeansClust4

PCAmatrixWithKMeans4 <- cbind(PCA$x, cluster = factor(kmeansClust4$cluster), groups = groups)

head(PCAmatrixWithKMeans4)

kmean4.plot <- ggplot(PCAmatrixWithKMeans4, aes(x = PC1, y = PC2, shape = cluster, color = groups.X)) + 
  geom_point(size = 8, alpha = 0.5)  + 
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))
kmean4.plot



# K-means computation for 3 centers
kmeansClust3 <- kmeans(PCA$x, centers = 3)
kmeansClust3

PCAmatrixWithKMeans3 <- cbind(PCA$x, cluster = factor(kmeansClust3$cluster), groups = groups)

head(PCAmatrixWithKMeans3)

kmean3.plot <- ggplot(PCAmatrixWithKMeans3, aes(x = PC1, y = PC2, shape = cluster, color = groups.X)) + 
  geom_point(size = 8, alpha = 0.5)  + 
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))
kmean3.plot



# K-means computation for 2 centers
kmeansClust2 <- kmeans(PCA$x, centers = 2)
kmeansClust2

PCAmatrixWithKMeans2 <- cbind(PCA$x, cluster = factor(kmeansClust2$cluster), groups = groups)

head(PCAmatrixWithKMeans2)

kmean2.plot <- ggplot(PCAmatrixWithKMeans2, aes(x = PC1, y = PC2, shape = cluster, color = groups.X)) + 
  geom_point(size = 8, alpha = 0.5)  + 
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))
kmean2.plot


# Save environmental image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")





# -----------------------------------
# Calculating means for each variable for each species ####
# -----------------------------------
head(extractedWITHcoord)

actaeus <- apply(extractedWITHcoord[extractedWITHcoord$species == "actaeus",4:13], 2, mean)
albolineatus <- apply(extractedWITHcoord[extractedWITHcoord$species == "albolineatus",4:13], 2, mean)
alipioi <- apply(extractedWITHcoord[extractedWITHcoord$species == "alipioi",4:13], 2, mean)
auroguttatus <- apply(extractedWITHcoord[extractedWITHcoord$species == "auroguttatus",4:13], 2, mean)
boticario <- apply(extractedWITHcoord[extractedWITHcoord$species == "boticario",4:13], 2, mean)
brunneus <- apply(extractedWITHcoord[extractedWITHcoord$species == "brunneus",4:13], 2, mean)
bufonoides <- apply(extractedWITHcoord[extractedWITHcoord$species == "bufonoides",4:13], 2, mean)
coloratus <- apply(extractedWITHcoord[extractedWITHcoord$species == "coloratus",4:13], 2, mean)
crispus <- apply(extractedWITHcoord[extractedWITHcoord$species == "crispus",4:13], 2, mean)
curupira <- apply(extractedWITHcoord[extractedWITHcoord$species == "curupira",4:13], 2, mean)
darkside <- apply(extractedWITHcoord[extractedWITHcoord$species == "darkside",4:13], 2, mean)
didactylus <- apply(extractedWITHcoord[extractedWITHcoord$species == "didactylus",4:13], 2, mean)
ephippium <- apply(extractedWITHcoord[extractedWITHcoord$species == "ephippium",4:13], 2, mean)
ferruginus <- apply(extractedWITHcoord[extractedWITHcoord$species == "ferruginus",4:13], 2, mean)
fuscolineatus <- apply(extractedWITHcoord[extractedWITHcoord$species == "fuscolineatus",4:13], 2, mean)
garbeanus <- apply(extractedWITHcoord[extractedWITHcoord$species == "garbeanus",4:13], 2, mean)
guarani<- apply(extractedWITHcoord[extractedWITHcoord$species == "guarani",4:13], 2, mean)
hermogenesi<- apply(extractedWITHcoord[extractedWITHcoord$species == "hermogenesi",4:13], 2, mean)
izecksohni<- apply(extractedWITHcoord[extractedWITHcoord$species == "izecksohni",4:13], 2, mean)
leopardus<- apply(extractedWITHcoord[extractedWITHcoord$species == "leopardus",4:13], 2, mean)
margaritatus<- apply(extractedWITHcoord[extractedWITHcoord$species == "margaritatus",4:13], 2, mean)
mariaeterezae<- apply(extractedWITHcoord[extractedWITHcoord$species == "mariaeterezae",4:13], 2, mean)
nodoterga<- apply(extractedWITHcoord[extractedWITHcoord$species == "nodoterga",4:13], 2, mean)
olivaceus<- apply(extractedWITHcoord[extractedWITHcoord$species == "olivaceus",4:13], 2, mean)
pernix<- apply(extractedWITHcoord[extractedWITHcoord$species == "pernix",4:13], 2, mean)
pitanga<- apply(extractedWITHcoord[extractedWITHcoord$species == "pitanga",4:13], 2, mean)
pombali<- apply(extractedWITHcoord[extractedWITHcoord$species == "pombali",4:13], 2, mean)
pulex<- apply(extractedWITHcoord[extractedWITHcoord$species == "pulex",4:13], 2, mean)
quiririensis<- apply(extractedWITHcoord[extractedWITHcoord$species == "quiririensis",4:13], 2, mean)
sp1<- apply(extractedWITHcoord[extractedWITHcoord$species == "sp1",4:13], 2, mean)
sp4<- apply(extractedWITHcoord[extractedWITHcoord$species == "sp4",4:13], 2, mean)
sp5<- apply(extractedWITHcoord[extractedWITHcoord$species == "sp5",4:13], 2, mean)
sulfuratus<- apply(extractedWITHcoord[extractedWITHcoord$species == "sulfuratus",4:13], 2, mean)
toby<- apply(extractedWITHcoord[extractedWITHcoord$species == "toby",4:13], 2, mean)
tridactylus<- apply(extractedWITHcoord[extractedWITHcoord$species == "tridactylus",4:13], 2, mean)
verrucosus<- apply(extractedWITHcoord[extractedWITHcoord$species == "verrucosus",4:13], 2, mean)
vertebralis<- apply(extractedWITHcoord[extractedWITHcoord$species == "vertebralis",4:13], 2, mean)


species <- as.data.frame(extractedWITHcoord$species)


species <- species[!duplicated(species), ]
species
#write.csv(species, "species.csv")


big.matrix <- rbind(actaeus, albolineatus, alipioi, auroguttatus, boticario, brunneus, bufonoides, coloratus, crispus,  curupira, darkside, didactylus, ephippium, ferruginus, fuscolineatus, garbeanus, guarani, hermogenesi, izecksohni,leopardus,margaritatus,mariaeterezae, nodoterga, olivaceus, pernix,pitanga,   pombali,   pulex,     quiririensis, sp1, sp4, sp5, sulfuratus, toby, tridactylus, verrucosus, vertebralis)
head(big.matrix)
dim(big.matrix)

#write.csv(big.matrix, "big_matrix.csv")



# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")




# -----------------------------------
# Calculating means for the z-scores of each variable for each species ####
# -----------------------------------
head(zScoresWITHcoord)


actaeus.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "actaeus",1:10], 2, mean)
albolineatus.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "albolineatus",1:10], 2, mean)
alipioi.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "alipioi",1:10], 2, mean)
auroguttatus.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "auroguttatus",1:10], 2, mean)
boticario.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "boticario",1:10], 2, mean)
brunneus.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "brunneus",1:10], 2, mean)
bufonoides.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "bufonoides",1:10], 2, mean)
coloratus.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "coloratus",1:10], 2, mean)
crispus.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "crispus",1:10], 2, mean)
curupira.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "curupira",1:10], 2, mean)
darkside.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "darkside",1:10], 2, mean)
didactylus.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "didactylus",1:10], 2, mean)
ephippium.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "ephippium",1:10], 2, mean)
ferruginus.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "ferruginus",1:10], 2, mean)
fuscolineatus.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "fuscolineatus",1:10], 2, mean)
garbeanus.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "garbeanus",1:10], 2, mean)
guarani.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "guarani",1:10], 2, mean)
hermogenesi.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "hermogenesi",1:10], 2, mean)
izecksohni.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "izecksohni",1:10], 2, mean)
leopardus.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "leopardus",1:10], 2, mean)
margaritatus.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "margaritatus",1:10], 2, mean)
mariaeterezae.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "mariaeterezae",1:10], 2, mean)
nodoterga.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "nodoterga",1:10], 2, mean)
olivaceus.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "olivaceus",1:10], 2, mean)
pernix.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "pernix",1:10], 2, mean)
pitanga.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "pitanga",1:10], 2, mean)
pombali.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "pombali",1:10], 2, mean)
pulex.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "pulex",1:10], 2, mean)
quiririensis.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "quiririensis",1:10], 2, mean)
sp1.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "sp1",1:10], 2, mean)
sp4.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "sp4",1:10], 2, mean)
sp5.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "sp5",1:10], 2, mean)
sulfuratus.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "sulfuratus",1:10], 2, mean)
toby.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "toby",1:10], 2, mean)
tridactylus.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "tridactylus",1:10], 2, mean)
verrucosus.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "verrucosus",1:10], 2, mean)
vertebralis.zscores <- apply(zScoresWITHcoord[zScoresWITHcoord$species == "vertebralis",1:10], 2, mean)

big.matrix.zscores <- rbind(actaeus.zscores, albolineatus.zscores, alipioi.zscores, auroguttatus.zscores, boticario.zscores, brunneus.zscores, bufonoides.zscores, coloratus.zscores, crispus.zscores,  curupira.zscores, darkside.zscores, didactylus.zscores, ephippium.zscores, ferruginus.zscores, fuscolineatus.zscores, garbeanus.zscores, guarani.zscores, hermogenesi.zscores, izecksohni.zscores,leopardus.zscores, margaritatus.zscores, mariaeterezae.zscores, nodoterga.zscores, olivaceus.zscores, pernix.zscores, pitanga.zscores,   pombali.zscores,   pulex.zscores,     quiririensis.zscores, sp1.zscores, sp4.zscores, sp5.zscores, sulfuratus.zscores, toby.zscores, tridactylus.zscores, verrucosus.zscores, vertebralis.zscores)

head(big.matrix.zscores)
dim(big.matrix.zscores)

rownames(big.matrix.zscores) <- gsub(rownames(big.matrix.zscores), patter = ".zscores", replacement = "")

head(big.matrix.zscores, 10)


#write.csv(big.matrix.zscores, "big_matrix_zscores.csv")



# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")




# -----------------------------------
# PCA with the means for each species #### (DOESN't MAKE SENSE)
# -----------------------------------
PCA.big.matrix <- prcomp(big.matrix, scale = T)
plot(PCA.big.matrix, main = "")
summary(PCA.big.matrix)
biplot(PCA.big.matrix)
head(PCA.big.matrix$x)

PCA.big.matrix.table <- cbind(groups[!duplicated(groups),], PCA.big.matrix$x)
rownames(PCA.big.matrix.table) <- NULL
PCA.big.matrix.table[30,3] <- "SUL"
View(PCA.big.matrix.table)

# Plotting the results by group
b <- ggplot(PCA.big.matrix.table, aes(x = PC1, y = PC2, color = group, label = species))
b + geom_point(size = 5) + 
  stat_ellipse(aes(x = PC1, y = PC2, fill = factor(group)),
               geom = "polygon", level = 0.95, alpha = 0.10) +
  guides(color = guide_legend("Groups"), fill = guide_legend("Groups")) +
  geom_text(size=6, nudge_y = 0.5, nudge_x = 0, check_overlap = F, show.legend = F)


# Plotting the results by clade
c <- ggplot(PCA.big.matrix.table, aes(x = PC1, y = PC2, color = X, label = species))
c + geom_point(size = 5) + 
  stat_ellipse(aes(x = PC1, y = PC2, fill = factor(X)),
               geom = "polygon", level = 0.95, alpha = 0.10) +
  guides(color = guide_legend("Clades"), fill = guide_legend("Clades")) + 
  scale_colour_manual(values=c("goldenrod", "darkolivegreen4", "hotpink4", "lightsteelblue4", "gray")) +
  geom_text(size=6, nudge_y = 0.5, nudge_x = 0, check_overlap = F, show.legend = F)





# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")




# -----------------------------------
# PCA with the z-score means for each species ####
# -----------------------------------

PCA.big.matrix.zscores <- prcomp(big.matrix.zscores, scale. = F)
plot(PCA.big.matrix.zscores, main = "")
summary(PCA.big.matrix.zscores)
biplot(PCA.big.matrix.zscores)
head(PCA.big.matrix.zscores$x)

PCA.big.matrix.table.zscores <- cbind(groups[!duplicated(groups),], PCA.big.matrix.zscores$x)
rownames(PCA.big.matrix.table.zscores) <- NULL
View(PCA.big.matrix.table.zscores)
PCA.big.matrix.table.zscores[30,3] <- "SUL"


# Plotting the results by group
b.zscores <- ggplot(PCA.big.matrix.table.zscores, aes(x = PC1, y = PC2, color = group, label = species))
b.zscores + geom_point(size = 10, alpha=0.6) + 
  theme(legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) + 
  stat_ellipse(aes(x = PC1, y = PC2, fill = factor(group)),
               geom = "polygon", level = 0.95, alpha = 0.10) +
  guides(color = guide_legend("Groups"), fill = guide_legend("Groups")) +
  geom_text(size=7, nudge_y = 0.015, nudge_x = 0.0, check_overlap = F, show.legend = F)


# Plotting the results by clade
c <- ggplot(PCA.big.matrix.table.zscores, aes(x = PC1, y = PC2, color = X, label = species))
c + geom_point(size = 10, alpha=0.6) +  
  theme(legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18)) + 
  stat_ellipse(aes(x = PC1, y = PC2, fill = factor(X)),
               geom = "polygon", level = 0.95, alpha = 0.10) +
  guides(color = guide_legend("Clades"), fill = guide_legend("Clades")) + 
  scale_colour_manual(values=c("goldenrod", "darkolivegreen4", "hotpink4", "lightsteelblue4", "gray")) +
  geom_text(size=7, nudge_y = 0.015, nudge_x = 0.0, check_overlap = F, show.legend = F)



# Plotting the biplots along with the points
b.zscores.biplot <- ggbiplot(PCA.big.matrix.zscores, 
                             obs.scale = 0.6, 
                             var.scale = 1, 
                             var.axes = T,
                             groups = groups.not.dupli$X, 
                             ellipse = T, 
                             ellipse.prob = .95,
                             circle = F, 
                             labels = rownames(PCA.big.matrix.zscores$x), 
                             labels.size = 6, 
                             varname.size = 5) + 
  guides(color = guide_legend("Clades"), fill = guide_legend("Clades"))

b.zscores.biplot +
  scale_colour_manual(values=c("goldenrod", "darkolivegreen4", "hotpink4", "lightsteelblue4", "gray")) +
  geom_point(aes(color = groups.not.dupli$X), size = 7, alpha = 0.3) +
  theme_bw() +
  theme(line = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed")




# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")




# -----------------------------------
# Calculating the k-means in the PCA but for one occurrence per species #### (UNDONE)
# -----------------------------------
# K-means computation for 5 centers
kmeansClust <- kmeans(PCA$x, centers = 5)
kmeansClust

PCAmatrixWithKMeans <- cbind(PCA$x, cluster = factor(kmeansClust$cluster), groups = groups)

head(PCAmatrixWithKMeans)

kmean.plot <- ggplot(PCAmatrixWithKMeans, aes(x = PC1, y = PC2, shape = cluster, color = groups.X)) + 
  geom_point(size = 8, alpha = 0.5)  + 
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))
kmean.plot



# K-means computation for 4 centers
kmeansClust4 <- kmeans(PCA$x, centers = 4)
kmeansClust4

PCAmatrixWithKMeans4 <- cbind(PCA$x, cluster = factor(kmeansClust4$cluster), groups = groups)

head(PCAmatrixWithKMeans4)

kmean4.plot <- ggplot(PCAmatrixWithKMeans4, aes(x = PC1, y = PC2, shape = cluster, color = groups.X)) + 
  geom_point(size = 8, alpha = 0.5)  + 
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))
kmean4.plot



# K-means computation for 3 centers
kmeansClust3 <- kmeans(PCA$x, centers = 3)
kmeansClust3

PCAmatrixWithKMeans3 <- cbind(PCA$x, cluster = factor(kmeansClust3$cluster), groups = groups)

head(PCAmatrixWithKMeans3)

kmean3.plot <- ggplot(PCAmatrixWithKMeans3, aes(x = PC1, y = PC2, shape = cluster, color = groups.X)) + 
  geom_point(size = 8, alpha = 0.5)  + 
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))
kmean3.plot



# K-means computation for 2 centers
kmeansClust2 <- kmeans(PCA$x, centers = 2)
kmeansClust2

PCAmatrixWithKMeans2 <- cbind(PCA$x, cluster = factor(kmeansClust2$cluster), groups = groups)

head(PCAmatrixWithKMeans2)

kmean2.plot <- ggplot(PCAmatrixWithKMeans2, aes(x = PC1, y = PC2, shape = cluster, color = groups.X)) + 
  geom_point(size = 8, alpha = 0.5)  + 
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18))
kmean2.plot


# Save environmental image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")





# -----------------------------------
# Calculating means for each variable for each clade ####
# -----------------------------------
pernix.clade <- rbind(actaeus, albolineatus, auroguttatus, boticario, mariaeterezae, olivaceus, pernix, pombali, quiririensis, verrucosus, fuscolineatus, ferruginus)
pernix.clade.mean <- apply(pernix.clade, 2, mean)


ephippium.clade <- rbind(ephippium, garbeanus, alipioi, hermogenesi, nodoterga, pitanga, toby, vertebralis, didactylus)
ephippium.clade.mean <- apply(ephippium.clade, 2, mean)


brunneus.clade <- rbind(brunneus, curupira, izecksohni, leopardus, sp4, sp5, tridactylus)
brunneus.clade.mean <- apply(brunneus.clade, 2, mean)


sulfuratus.clade <- rbind(sulfuratus, sp1)
sulfuratus.clade.mean <- apply(sulfuratus.clade, 2, mean)


clade.matrix <- rbind(pernix.clade.mean, ephippium.clade.mean, brunneus.clade.mean, sulfuratus.clade.mean)

head(clade.matrix)




# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")

# -----------------------------------
# Calculating means for each variable for each group ####
# -----------------------------------
pernix.group <- rbind(actaeus, albolineatus, auroguttatus, boticario, brunneus, coloratus, curupira, ferruginus, fuscolineatus, izecksohni, leopardus, mariaeterezae, olivaceus, pernix, pombali, quiririensis, sp4, sp5, tridactylus, verrucosus)
pernix.group.mean <- apply(pernix.group, 2, mean)


didactylus.group <- rbind(didactylus, hermogenesi, pulex, sp1, sulfuratus)
didactylus.group.mean <- apply(didactylus.group, 2, mean)


ephippium.group <- rbind(alipioi, bufonoides, crispus, darkside, ephippium, garbeanus, guarani, margaritatus, nodoterga, pitanga, toby, vertebralis)
ephippium.group.mean <- apply(ephippium.group, 2, mean)


group.matrix <- rbind(pernix.group.mean, ephippium.group.mean, didactylus.group.mean)

group.matrix




# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")




# -----------------------------------
# Mahalanobis distances between species ####
# -----------------------------------
extractedWITHgroupsANDclades <- cbind(groups, extracted)

mahalonobis.dist.species <- pairwise.mahalanobis(extractedWITHgroupsANDclades[,4:13], grouping = extractedWITHgroupsANDclades$species)

rownames(mahalonobis.dist.species$distance) <- species
colnames(mahalonobis.dist.species$distance) <- species

corrplot(mahalonobis.dist.species$distance, is.corr = FALSE, tl.col = "black")



# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")


# -----------------------------------
# Mahalanobis distance matrix between clades ####
# -----------------------------------
mahalonobis.dist.clades <- pairwise.mahalanobis(extractedWITHgroupsANDclades[,4:13], grouping = extractedWITHgroupsANDclades$X)

colnames(mahalonobis.dist.clades$distance) <- c("Brunneus", "Ephippium", "Pernix", "Sulfuratus", "Unknown")
rownames(mahalonobis.dist.clades$distance) <- c("Brunneus", "Ephippium", "Pernix", "Sulfuratus", "Unknown")

corrplot(mahalonobis.dist.clades$distance, is.corr = FALSE, tl.col = "black")




# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")



# -----------------------------------
# Mahalanobis distance matrix between groups ####
# -----------------------------------
mahalonobis.dist.groups <- pairwise.mahalanobis(extractedWITHgroupsANDclades[,4:13], grouping = extractedWITHgroupsANDclades$group)

rownames(mahalonobis.dist.groups$distance) <- c("Didactylus", "Ephippium", "Pernix")
colnames(mahalonobis.dist.groups$distance) <- c("Didactylus", "Ephippium", "Pernix")



corrplot(mahalonobis.dist.groups$distance, is.corr = FALSE, tl.col = "black")



# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")



# -----------------------------------
# Phylogenetic distances ####
# -----------------------------------
tree <- read.tree("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\phylogenetic\\tree.nwk")
phylo.dist <- distTips(tree, tips ="all", method = "patristic", useC = TRUE)
phylo.dist

# Preparing the phylogentic matrix:
phylo.dist.matrix <- as.matrix(phylo.dist)
rownames(phylo.dist.matrix) <- gsub(rownames(as.matrix(phylo.dist)), pattern = "Brachycephalus_", replacement = "" )
colnames(phylo.dist.matrix) <- gsub(colnames(as.matrix(phylo.dist)), pattern = "Brachycephalus_", replacement = "" )
phylo.dist.matrix <- phylo.dist.matrix[-31,-31]

phylo.dist.matrix <- phylo.dist.matrix[order(rownames(phylo.dist.matrix)), order(colnames(phylo.dist.matrix))]

View(phylo.dist.matrix)


# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")


# -----------------------------------
# Mantel test: Mahalanobis per species vs phylogenetic distance ####
# -----------------------------------
# Preparing the Mahalanobis matrx:
Mahala.matrix.species <- as.matrix(mahalonobis.dist.species$distance)
Mahala.matrix.species.cleaned <- Mahala.matrix.species[-c(7,8,9,11,17,21,28,39),-c(7,8,9,11,17,21,28,39)]

View(Mahala.matrix.species.cleaned)


# Verifying if the matrices have the same dimensions
dim(phylo.dist.matrix) == dim(Mahala.matrix.species.cleaned)
rownames(phylo.dist.matrix) == rownames(Mahala.matrix.species.cleaned)


# Mantel test
mantel.MahalaVSphylogen.species <- mantel.rtest(as.dist(Mahala.matrix.species.cleaned), as.dist(phylo.dist.matrix), nrepet = 9999)
mantel.MahalaVSphylogen.species

# Plotting the Mantel histogram
plot(mantel.MahalaVSphylogen.species, main = "Simulation Histogram - Mahalanobis Distances \n vs. Phylogenetic Distances", xlab = "Simulation r-Values", ylab = "Permutation Frequency")



# Mantel correlogram
mantel.correlog.Mahala <- mantel.correlog(D.geo = as.dist(Mahala.matrix.species.cleaned), D.eco = as.dist(phylo.dist.matrix), nperm = 9999, cutoff = F)
mantel.correlog.Mahala
plot(mantel.correlog.Mahala)



# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")


# -----------------------------------
# Importing niche overlap (Schoener's D) table ####
# -----------------------------------
niche.overlap.Dtable <- read.csv("B:\\Brachycephalus\\ENMTools\\Overlap\\inphylogeny_nicheoverlap_D_output.csv", header = T)
rownames.niche.overlap.Dtable <- niche.overlap.Dtable[,1]
rownames.niche.overlap.Dtable <- gsub(pattern = ".asc", rownames.niche.overlap.Dtable, replacement = "")
rownames(niche.overlap.Dtable) <- rownames.niche.overlap.Dtable
niche.overlap.Dtable <- niche.overlap.Dtable[,-1]
colnames(niche.overlap.Dtable) <- rownames.niche.overlap.Dtable

niche.overlap.Dtable <- niche.overlap.Dtable[-c(4,9),-c(4,9)] # Removing darkside and margaritatus


corrplot.mixed(as.matrix(niche.overlap.Dtable), is.corr = F, lower = "circle", upper = "number", tl.col = "black", tl.cex = 2 , number.cex = 0.75, cl.cex = 0.8, tl.pos = "lt", mar = c(0,0,0,2))


# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")


# -----------------------------------
# Mantel test: niche overlap (Schoener's D) per species vs phylogenetic distance ####
# -----------------------------------
# Preparing the phylogenetic matrix:
phylo.dist.matrix.for.overlap <- phylo.dist.matrix[-c(2,4,5,7,10,11,14,15,16,19,21,24:30), -c(2,4,5,7,10,11,14,15,16,19,21,24:30)]


# Verifying the dimensions:
dim(phylo.dist.matrix.for.overlap) == dim(niche.overlap.Dtable)
rownames(phylo.dist.matrix.for.overlap) == rownames(niche.overlap.Dtable)


# Mantel test:
mantel.overlapDVSphylogen.species <- mantel.rtest(as.dist(niche.overlap.Dtable), as.dist(phylo.dist.matrix.for.overlap), nrepet = 9999)
mantel.overlapDVSphylogen.species
1 - mantel.overlapDVSphylogen.species$pvalue # This is the one-tailed p-value for negative Mantel correlation tests


plot(mantel.overlapVSphylogen.species, main = "Simulation Histogram - Niche Overlap (Schoener's D)", xlab = "Simulation r-Values", ylab = "Permutation Frequency")


# Mantel correlogram
mantel.correlog.Doverlap <- mantel.correlog( D.geo = as.dist(niche.overlap.Dtable), D.eco = as.dist(phylo.dist.matrix.for.overlap), nperm = 9999, cutoff = F)
mantel.correlog.Doverlap
plot(mantel.correlog.Doverlap)



# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")


# -----------------------------------
# Importing niche overlap (Warren's I) table ####
# -----------------------------------

niche.overlap.Itable <- read.csv("B:\\Brachycephalus\\ENMTools\\Overlap\\inphylogeny_nicheoverlap_I_output.csv", header = T)
rownames.niche.overlap.Itable <- niche.overlap.Itable[,1]
rownames.niche.overlap.Itable <- gsub(pattern = ".asc", rownames.niche.overlap.Itable, replacement = "")
rownames(niche.overlap.Itable) <- rownames.niche.overlap.Itable
niche.overlap.Itable <- niche.overlap.Itable[,-1]
colnames(niche.overlap.Itable) <- rownames.niche.overlap.Itable

niche.overlap.Itable <- niche.overlap.Itable[-c(4,9),-c(4,9)] # Removing darkside and margaritatus


corrplot.mixed(as.matrix(niche.overlap.Itable), is.corr = F, lower = "circle", upper = "number", tl.col = "black", tl.cex = 2 , number.cex = 0.75, cl.cex = 0.8, tl.pos = "lt", mar = c(0,0,0,2))


# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")


# -----------------------------------
# Mantel test: niche overlap (Warren's I) per species vs phylogenetic distance ####
# -----------------------------------

# Verifying the dimensions:
dim(phylo.dist.matrix.for.overlap) == dim(niche.overlap.Itable)
rownames(phylo.dist.matrix.for.overlap) == rownames(niche.overlap.Itable)


# Mantel test:
mantel.overlapIVSphylogen.species <- mantel.rtest(as.dist(niche.overlap.Itable), as.dist(phylo.dist.matrix.for.overlap), nrepet = 9999)
mantel.overlapIVSphylogen.species
1 - mantel.overlapIVSphylogen.species$pvalue # This is the one-tailed p-value for negative Mantel correlation tests


plot(mantel.overlapVSphylogen.species, main = "Simulation Histogram - Niche Overlap (Warren's I)", xlab = "Simulation r-Values", ylab = "Permutation Frequency")



# Mantel correlogram
mantel.correlog.Ioverlap <- mantel.correlog( D.geo = as.dist(niche.overlap.Itable), D.eco = as.dist(phylo.dist.matrix.for.overlap), nperm = 9999, cutoff = F)
mantel.correlog.Ioverlap
plot(mantel.correlog.Ioverlap)



# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")




# -----------------------------------
# Corrplot with D and I at the same plot ####
# -----------------------------------
########### Per species
cor_D <- niche.overlap.Dtable #cor(mydata[which(mydata$group=='A'),c(1:5)])
cor_I <- niche.overlap.Itable #cor(mydata[which(mydata$group=='B'),c(1:5)])

dim(cor_D) == dim(cor_I)

# replace the values
cor_D[lower.tri(cor_D)] <- cor_I[lower.tri(cor_I)]
cor_D <- as.matrix(cor_D)


col1 <- colorRampPalette(c("black","black", "black", "white", "blue", "green", "yellow", "orange", "red"))
corrplot(cor_D, is.corr = F, method = "square", tl.offset = 0.03, tl.col = "black", col = col1(200), diag = F)




########### Per clades
# Import D
clade.niche.overlap.D <- read.csv("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Artigo\\Molecular phylogenetics and niche conservatism of Brachycephalus\\ENMTools\\Overlap\\clades_niche_overlap_D_output.csv", header = T, row.names = 1)
clade.niche.overlap.D

# Import I
clade.niche.overlap.I <- read.csv("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Artigo\\Molecular phylogenetics and niche conservatism of Brachycephalus\\ENMTools\\Overlap\\clades_niche_overlap_I_output.csv", header = T, row.names = 1)
clade.niche.overlap.I


# Creating a new matrix with the combination of the matrix above
corr.niche.overlap.clades <- clade.niche.overlap.D
corr.niche.overlap.clades[lower.tri(corr.niche.overlap.clades)] <- clade.niche.overlap.I[lower.tri(clade.niche.overlap.I)]
corr.niche.overlap.clades <- as.matrix(corr.niche.overlap.clades)
rownames(corr.niche.overlap.clades) <- c("BRU","EPH","PER","SUL")
colnames(corr.niche.overlap.clades) <- c("BRU","EPH","PER","SUL")


# Plotting
col1 <- colorRampPalette(c("black","black", "black", "white", "blue", "green", "yellow", "orange", "red"))
corrplot(corr.niche.overlap.clades, is.corr = F, method = "square", tl.cex = 2.5, tl.pos = "d", tl.offset = 0.03, tl.col = "black", col = col1(200), diag = T)




# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")



# -----------------------------------
# PNO: predicted niche occupancy ####
# -----------------------------------

uncorVar


## Bioclimatic Variable 1
# Per species for bioclim 1 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.1 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio1_a.asc", 
                         path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.1, legend.pos = "bottomleft")

# Per clade for bioclim 1 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.1 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio1_a.asc", 
                       path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.1, legend.pos = "bottomleft")





## Bioclimatic Variable 2
# Per species for bioclim 2 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.2 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio2_a.asc", 
                         path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.2, legend.pos = "bottomleft")

# Per clade for bioclim 2 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.2 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio2_a.asc", 
                       path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.2, legend.pos = "bottomleft")








## Bioclimatic Variable 3
# Per species for bioclim 3 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.3 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio3_a.asc", 
                         path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.3, legend.pos = "bottomleft")

# Per clade for bioclim 3 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.3 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio3_a.asc", 
                       path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.3, legend.pos = "bottomleft")









## Bioclimatic Variable 4
# Per species for bioclim 4 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.4 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio4_a.asc", 
                       path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.4, legend.pos = "bottomleft")

# Per clade for bioclim 4 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.4 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio4_a.asc", 
                       path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.4, legend.pos = "bottomleft")






## Bioclimatic Variable 5
# Per species for bioclim 5 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.5 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio5_a.asc", 
                         path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.5, legend.pos = "bottomleft")

# Per clade for bioclim 5 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.5 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio5_a.asc", 
                       path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.5, legend.pos = "bottomleft")







## Bioclimatic Variable 6
# Per species for bioclim 6 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.6 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio6_a.asc", 
                         path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.6, legend.pos = "bottomleft")

# Per clade for bioclim 6 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.6 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio6_a.asc", 
                       path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.6, legend.pos = "bottomleft")








## Bioclimatic Variable 7
# Per species for bioclim 7 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.7 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio7_a.asc", 
                         path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.7, legend.pos = "bottomleft")

# Per clade for bioclim 7 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.7 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio7_a.asc", 
                       path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.7, legend.pos = "bottomleft")






## Bioclimatic Variable 8
# Per species for bioclim 8 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.8 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio8_a.asc", 
                         path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.8, legend.pos = "bottomleft")

# Per clade for bioclim 8 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.8 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio8_a.asc", 
                       path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.8, legend.pos = "bottomleft")






## Bioclimatic Variable 9
# Per species for bioclim 9 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.9 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio9_a.asc", 
                         path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.9, legend.pos = "bottomleft")

# Per clade for bioclim 9 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.9 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio9_a.asc", 
                       path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.9, legend.pos = "bottomleft")






## Bioclimatic Variable 10
# Per species for bioclim 10 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.10 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio10_a.asc", 
                          path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.10, legend.pos = "bottomleft")

# Per clade for bioclim 10 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.10 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio10_a.asc", 
                        path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.10, legend.pos = "bottomleft")






## Bioclimatic Variable 11
# Per species for bioclim 11 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.11 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio11_a.asc", 
                          path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.11, legend.pos = "bottomleft")

# Per clade for bioclim 11 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.11 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio11_a.asc", 
                        path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.11, legend.pos = "bottomleft")







## Bioclimatic Variable 12
# Per species for bioclim 12 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.12 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio12_a.asc", 
                       path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.12, legend.pos = NULL)

# Per clade for bioclim 12 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.12 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio12_a.asc", 
                     path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.12, legend.pos = "bottomleft")









## Bioclimatic Variable 13
# Per species for bioclim 13 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.13 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio13_a.asc", 
                          path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.13, legend.pos = NULL)

# Per clade for bioclim 13 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.13 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio13_a.asc", 
                        path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.13, legend.pos = "bottomleft")








## Bioclimatic Variable 14
# Per species for bioclim 14 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.14 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio14_a.asc", 
                          path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.14, legend.pos = NULL)

# Per clade for bioclim 14 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.14 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio14_a.asc", 
                        path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.14, legend.pos = "bottomleft")









## Bioclimatic Variable 15
# Per species for bioclim 15 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.15 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio15_a.asc", 
                          path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.15, legend.pos = NULL)

# Per clade for bioclim 15 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.15 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio15_a.asc", 
                        path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.15, legend.pos = "bottomleft")







## Bioclimatic Variable 16
# Per species for bioclim 16 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.16 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio16_a.asc", 
                          path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.16, legend.pos = "bottomleft")

# Per clade for bioclim 16 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.16 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio16_a.asc", 
                        path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.16, legend.pos = "bottomleft")






## Bioclimatic Variable 17
# Per species for bioclim 17 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.17 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio17_a.asc", 
                          path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.17, legend.pos = "bottomleft")

# Per clade for bioclim 17 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.17 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio17_a.asc", 
                        path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.17, legend.pos = "bottomleft")








## Bioclimatic Variable 18
# Per species for bioclim 18 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.18 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio18_a.asc", 
                       path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.18, legend.pos = NULL)

# Per clade for bioclim 18 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.18 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio18_a.asc", 
                     path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.18, legend.pos = "bottomleft")








## Bioclimatic Variable 19
# Per species for bioclim 19 (Perform this for the rest of the influent variables of the PCAs)
pno.per.species.19 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio19_a.asc", 
                          path_model = "B:/Brachycephalus/MaxEnt/Brachycephalus_Results_raw/asc")

plotPNO(pno.per.species.19, legend.pos = "bottomleft")

# Per clade for bioclim 19 (Perform this for the rest of the influent variables of the PCAs)
pno.per.clade.19 <- pno(path_bioclim = "B:/WorldClim/WorldClim 1.4/Mata_Atlantica/Bioclim/asc/bio19_a.asc", 
                        path_model = "B:/Brachycephalus/MaxEnt/Clades_Results_raw/asc")

plotPNO(pno.per.clade.19, legend.pos = "bottomleft")









# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")



# -----------------------------------
# PNO weighted mean for species ####
# -----------------------------------

# Incorporate the pno.weighted.mean function to R
pno.weighted.mean <- 						
  function(x, subset = NULL, normalize = TRUE){
    
    #x <- cbind("VAR" = as.numeric(rownames(x)), x)
    
    # subset matrix
    # ---------------------
    if (!is.null(subset))
      x <- x[, c(1, which(names(x) %in% subset))]
    
    nb <- dim(x)[2]
    wm <- vector(length = (nb - 1))
    for (i in 2:nb){
      if (normalize){
        nf <- sum(x[, i])
        wm[i - 1] <- sum(x[, 1] * (x[, i] / nf))
      }
      else
        wm[i - 1] <- sum(x[, 1] * x[, i])
    }
    names(wm) <- colnames(x)[2:nb]
    wm
  }




# Perform the mean calculation
pno.mean.bio1.species <- pno.weighted.mean(pno.per.species.1, normalize = T)
pno.mean.bio2.species <- pno.weighted.mean(pno.per.species.2, normalize = T)
pno.mean.bio3.species <- pno.weighted.mean(pno.per.species.3, normalize = T)
pno.mean.bio4.species <- pno.weighted.mean(pno.per.species.4, normalize = T)
pno.mean.bio5.species <- pno.weighted.mean(pno.per.species.5, normalize = T)
pno.mean.bio6.species <- pno.weighted.mean(pno.per.species.6, normalize = T)
pno.mean.bio7.species <- pno.weighted.mean(pno.per.species.7, normalize = T)
pno.mean.bio8.species <- pno.weighted.mean(pno.per.species.8, normalize = T)
pno.mean.bio9.species <- pno.weighted.mean(pno.per.species.9, normalize = T)
pno.mean.bio10.species <- pno.weighted.mean(pno.per.species.10, normalize = T)
pno.mean.bio11.species <- pno.weighted.mean(pno.per.species.11, normalize = T)
pno.mean.bio12.species <- pno.weighted.mean(pno.per.species.12, normalize = T)
pno.mean.bio13.species <- pno.weighted.mean(pno.per.species.13, normalize = T)
pno.mean.bio14.species <- pno.weighted.mean(pno.per.species.14, normalize = T)
pno.mean.bio15.species <- pno.weighted.mean(pno.per.species.15, normalize = T)
pno.mean.bio16.species <- pno.weighted.mean(pno.per.species.16, normalize = T)
pno.mean.bio17.species <- pno.weighted.mean(pno.per.species.17, normalize = T)
pno.mean.bio18.species <- pno.weighted.mean(pno.per.species.18, normalize = T)
pno.mean.bio19.species <- pno.weighted.mean(pno.per.species.19, normalize = T)




# Binding all the means into a single matrix
pno.means.species <- rbind(pno.mean.bio1.species, pno.mean.bio2.species, pno.mean.bio3.species, pno.mean.bio4.species, pno.mean.bio5.species, pno.mean.bio6.species, pno.mean.bio7.species, pno.mean.bio8.species ,pno.mean.bio9.species, pno.mean.bio10.species,pno.mean.bio11.species, pno.mean.bio12.species, pno.mean.bio13.species, pno.mean.bio14.species, pno.mean.bio15.species, pno.mean.bio16.species, pno.mean.bio17.species, pno.mean.bio18.species, pno.mean.bio19.species)
pno.means.species

View(pno.means.species)



t.pno.means.species <- t(pno.means.species)
t.pno.means.species[,c(1:10)] <- t.pno.means.species[,c(1:10)]/10


write.csv(round(t.pno.means.species, 2), "pno_means_species.csv")




# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")




# -----------------------------------
# PNO weighted mean for clades ####
# -----------------------------------
pno.mean.bio1.clades <- pno.weighted.mean(pno.per.clade.1, normalize = T)
pno.mean.bio2.clades <- pno.weighted.mean(pno.per.clade.2, normalize = T)
pno.mean.bio3.clades <- pno.weighted.mean(pno.per.clade.3, normalize = T)
pno.mean.bio4.clades <- pno.weighted.mean(pno.per.clade.4, normalize = T)
pno.mean.bio5.clades <- pno.weighted.mean(pno.per.clade.5, normalize = T)
pno.mean.bio6.clades <- pno.weighted.mean(pno.per.clade.6, normalize = T)
pno.mean.bio7.clades <- pno.weighted.mean(pno.per.clade.7, normalize = T)
pno.mean.bio8.clades <- pno.weighted.mean(pno.per.clade.8, normalize = T)
pno.mean.bio9.clades <- pno.weighted.mean(pno.per.clade.9, normalize = T)
pno.mean.bio10.clades <- pno.weighted.mean(pno.per.clade.10, normalize = T)
pno.mean.bio11.clades <- pno.weighted.mean(pno.per.clade.11, normalize = T)
pno.mean.bio12.clades <- pno.weighted.mean(pno.per.clade.12, normalize = T)
pno.mean.bio13.clades <- pno.weighted.mean(pno.per.clade.13, normalize = T)
pno.mean.bio14.clades <- pno.weighted.mean(pno.per.clade.14, normalize = T)
pno.mean.bio15.clades <- pno.weighted.mean(pno.per.clade.15, normalize = T)
pno.mean.bio16.clades <- pno.weighted.mean(pno.per.clade.16, normalize = T)
pno.mean.bio17.clades <- pno.weighted.mean(pno.per.clade.17, normalize = T)
pno.mean.bio18.clades <- pno.weighted.mean(pno.per.clade.18, normalize = T)
pno.mean.bio19.clades <- pno.weighted.mean(pno.per.clade.19, normalize = T)



# Binding all the means into a single matrix
pno.means.clades <- rbind(pno.mean.bio1.clades, pno.mean.bio2.clades, pno.mean.bio3.clades, pno.mean.bio4.clades, pno.mean.bio5.clades, pno.mean.bio6.clades, pno.mean.bio7.clades, pno.mean.bio8.clades ,pno.mean.bio9.clades, pno.mean.bio10.clades,pno.mean.bio11.clades, pno.mean.bio12.clades, pno.mean.bio13.clades, pno.mean.bio14.clades, pno.mean.bio15.clades, pno.mean.bio16.clades, pno.mean.bio17.clades, pno.mean.bio18.clades, pno.mean.bio19.clades)

View(pno.means.clades)

t.pno.means.clades <- t(pno.means.clades)
t.pno.means.clades[,c(1:10)] <- t.pno.means.clades[,c(1:10)]/10


write.csv(round(t.pno.means.clades, 2), "pno_means_clades.csv")




# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")



# -----------------------------------
# PCA with PNO means per species ####
# -----------------------------------
# Preparing the data
pno.means.species <- t(pno.means.species)
colnames(pno.means.species) <- gsub(colnames(pno.means.species), pattern = "pno.mean.", replacement = "")
colnames(pno.means.species) <- gsub(colnames(pno.means.species), pattern = ".species", replacement = "")


PCA.PNO.species <- prcomp(pno.means.species, scale = T)
plot(PCA.PNO.species)
summary(PCA.PNO.species)
biplot(PCA.PNO.species)


# Plotting the PCA with the biplot

PCA.PNO.biplot <- ggbiplot(PCA.PNO.species, 
                             obs.scale = 0.6, 
                             var.scale = 1, 
                             var.axes = T,
                             groups = NULL, 
                             ellipse = T, 
                             ellipse.prob = .95,
                             circle = F, 
                             labels = rownames(PCA.PNO.species$x), 
                             labels.size = 6, 
                             varname.size = 5) + 
  guides(color = guide_legend("Clades"), fill = guide_legend("Clades"))

PCA.PNO.biplot #+
  scale_colour_manual(values=c("goldenrod", "darkolivegreen4", "hotpink4", "lightsteelblue4", "gray")) +
  geom_point(aes(color = groups.not.dupli$X), size = 7, alpha = 0.3) +
  theme_bw() +
  theme(line = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed")









# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")




# -----------------------------------
# Ancestral Niche Reconstruction Bio 12 ####
# -----------------------------------

plot(tree)

colnames(pno.per.species.12)

dropped.tree <- drop.tip(phy = tree, trim.internal = T, subtree = F, interactive = T)
dropped.tree <- drop.tip(phy = dropped.tree, trim.internal = T, subtree = F, interactive = T)
# c("Brachycephalus vertebralis", "Brachycephalus toby", "Brachycephalus hermogenesi", "Brachycephalus sulfuratus", "Brachycephalus tridactylus", "Brachycephalus leopardus", "Brachycephalus sp4", "Brachycephalus izecksohni", "Brachycephalus sp5", "Brachycephalus curupira", "Brachycephalus pernix", "Brachycephalus pombali", "Brachycephalus ferruginus", "Brachycephalus auroguttatus", "Brachycephalus verrucosus", "Brachycephalus mariaeterezae", "Brachycephalus boticario", "Brachycephalus fuscolineatus", "Brachycephalus albolineatus", "Ischnocnema")
plot(dropped.tree)
dropped.tree$tip.label <- gsub(dropped.tree$tip.label, pattern = "Brachycephalus_", replacement = "")


pno.anc.bio12 <- data.frame(pno.per.species.12[,-c(5,10)])


colnames(pno.anc.bio12)

colnames(pno.anc.bio12[,c(2,10,12,4,13,3,8,9,11,6,7,5)]) == dropped.tree$tip.label
pno.anc.bio12 <- pno.anc.bio12[,c(1,2,10,12,4,13,3,8,9,11,6,7,5)]

names(pno.anc.bio12[,2:13]) == dropped.tree$tip.label
pno.anc.bio12


anc.clim.bio12 <- anc.clim(target = dropped.tree, pno = pno.anc.bio12)
plotAncClim(anc.clim.bio12)







# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")




# -----------------------------------
# Ancestral Niche Reconstruction Bio 18 ####
# -----------------------------------

plot(dropped.tree)

pno.anc.bio18 <- data.frame(pno.per.species.18[,-c(5,10)])


colnames(pno.anc.bio18)

colnames(pno.anc.bio18[,c(2,10,12,4,13,3,8,9,11,6,7,5)]) == dropped.tree$tip.label
pno.anc.bio18 <- pno.anc.bio18[,c(1,2,10,12,4,13,3,8,9,11,6,7,5)]

names(pno.anc.bio18[,2:13]) == dropped.tree$tip.label
pno.anc.bio18


anc.clim.bio18 <- anc.clim(target = dropped.tree, pno = pno.anc.bio18)
plotAncClim(anc.clim.bio18)






# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")




# -----------------------------------
# Ancestral Niche Reconstruction Bio 4 ####
# -----------------------------------

plot(dropped.tree)

pno.anc.bio4 <- data.frame(pno.per.species.4[,-c(5,10)])


colnames(pno.anc.bio4)

colnames(pno.anc.bio4[,c(2,10,12,4,13,3,8,9,11,6,7,5)]) == dropped.tree$tip.label
pno.anc.bio4 <- pno.anc.bio4[,c(1,2,10,12,4,13,3,8,9,11,6,7,5)]

names(pno.anc.bio4[,2:13]) == dropped.tree$tip.label
pno.anc.bio4


anc.clim.bio4 <- anc.clim(target = dropped.tree, pno = pno.anc.bio4)
plotAncClim(anc.clim.bio4)





# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")




# -----------------------------------
#  ####
# -----------------------------------





# Save environment image
save.image("C:\\Users\\piato\\Dropbox\\BIOLOGIA\\Brachycephalus\\Analysis_True\\R\\analysis_true.RData")

