
phenotypes <- read.table("m2pvectestphen.dat",header=FALSE)

pca <- prcomp(phenotypes)

pcavecs <- pca$rotation
lapply(pcavecs, write, "m2pcavectest.dat",append=TRUE)
