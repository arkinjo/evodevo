
phenotypes <- read.table("pvectestphen6.dat",header=FALSE)

pca <- prcomp(phenotypes)

pcavecs <- pca$rotation
lapply(pcavecs, write, "pcavectest6.dat",append=TRUE)
