
phen <- read.table("phenotypes_1.dat",header=FALSE)

cphen = phen
for (i in c(1:length(cphen))){ #Centralize
  for (j in c(1:length(cphen[i]))){
    cphen[i][j] = cphen[i][j]-mean(unlist(phen[i]))
  }
}

svd.cphen = svd(cphen)

csvs <- svd(cphen)$d
plot(csvs,main="Plot of singular values",xlab="Index",ylab="Singular value")


v1 <- svd(cphen)$v[,1] #Axis one
v2 <- svd(cphen)$v[,2] #Axis two

#Sanity check
print(sum(v1*v1))
print(sum(v2*v2))
print(sum(v1*v2))


mat.cphen = as.matrix(cphen)

x = c(1:length(mat.cphen[,1]))
for (i in c(1:length(mat.cphen[,1]))){
  x[i] = sum(v1*mat.cphen[i,])
}

y = c(1:length(mat.cphen[,1]))
for (i in c(1:length(mat.cphen[,1]))){
  y[i] = sum(v2*mat.cphen[i,])
}

plot(x,y,type="p",main="Visualization plot test",xlim=c(-2,2),ylim=c(-2,2),xlab="First principal component",ylab="Second principle component")
abline(h=0)
abline(v=0)
