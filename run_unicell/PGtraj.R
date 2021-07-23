

for (i in c(1:200)) {
  filename <- paste("PGtraj1_",i,".dat",sep="")
  pg <- read.table(filename,header=TRUE)
  p <- pg$Phenotype
  g <- pg$Genotype
  #pdfout <- paste("PGtraj1_",i,".pdf",sep="")
  #pdf(pdfout)
  plot(g,p,xlim=c(-15,15),ylim=c(-2,2),xlab="genotype",ylab="phenotype")
  #dev.off()
}