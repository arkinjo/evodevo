
phenotypes <- read.table("upop20211223trainphen.dat",header=FALSE)

pca <- prcomp(phenotypes)

pcavecs <- pca$rotation
print((pca$sdev[1]/pca$sdev[length(pca$sdev)])^2)

lapply(pcavecs, write, "upop20211223trainpvec.dat",append=TRUE)
