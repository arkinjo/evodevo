
phenotypes <- read.table("pvectestphen8.dat",header=FALSE)

pca <- prcomp(phenotypes)

pcavecs <- pca$rotation
lapply(pcavecs, write, "pcavectest8.dat",append=TRUE)
