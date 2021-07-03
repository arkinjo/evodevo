### Plot phenotypes
# Defining axes

ancphen <- read.table("phenotypes_1.dat",header=FALSE)
c.ancphen = scale(ancphen,scale=FALSE)
#c.ancphen = ancphen


#for (i in c(1:length(c.ancphen))){
#  for (j in c(1:length(c.ancphen[i]))){
#    c.ancphen[i][j] = c.ancphen[i][j] - mean(unlist(ancphen[i]))
#  }
#}

svd.ancphen <- svd(c.ancphen)
xaxis <- svd.ancphen$v[,1]

novphen <- read.table("phen5_200.dat",header=FALSE)
c.novphen = scale(novphen,scale=FALSE)
#c.novphen = novphen
#for (i in c(1:length(c.novphen))){
#  for (j in c(1:length(c.novphen[i]))){
#    c.novphen[i][j] = c.novphen[i][j] - mean(unlist(novphen[i]))
#  }
#}

svd.novphen <- svd(c.novphen)
yaxis <- svd.novphen$v[,1]

for (gen in c(1:200)){
  filename = paste("phen5_",gen,".dat",sep="")
  phen <- read.table(filename,header=FALSE)
  c.phen = scale(phen,scale=FALSE)
  #c.phen = phen
  #for (i in c(1:length(c.phen))){
  #  for (j in c(1:length(c.phen[i]))){
  #    c.phen[i][j] = c.phen[i][j] - mean(unlist(phen[i]))
  #  }
  #}
  mat.c.phen = as.matrix(c.phen)
  x = c(1:length(mat.c.phen[,1]))
  for (i in c(1:length(mat.c.phen[,1]))){
    x[i] = sum(xaxis*mat.c.phen[i,])
  }
  y = c(1:length(mat.c.phen[,1]))
  for (i in c(1:length(mat.c.phen[,1]))){
    y[i] = sum(yaxis*mat.c.phen[i,])
  }
  title = paste("Generation",gen,sep="_")
  plot(x,y,type="p",main=title,xlim=c(-1,1),ylim=c(-1,1),xlab="Ancestral principal phenotype",ylab="Novel principal phenotype")
  abline(h=0)
  abline(v=0)
}

csvs <- svd.novphen$d
plot(csvs,main="Plot of singular values",xlab="Index",ylab="Singular value")
plot(cumsum(csvs),main="Cumulative sum of singular values",xlab="Index",ylab="Sum")

#cumsum.csvs <- cumsum(csvs)
vars <- csvs*csvs
plot(vars,main="Plot of square of singular values",xlab="Index",ylab="Singular value")
plot(cumsum(vars),main="Cumulative sum of squares",xlab="Index",ylab="Sum")
