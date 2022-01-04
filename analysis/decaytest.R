
test <- read.table("upop20220101traintraj.dat",header=TRUE)

gen <- test$Generation
MSE <- test$MSE

f <- function(x,a,b,c){
  c + a*exp(-x/b)
}

a0 = MSE[1]-MSE[length(MSE)]
b0 = 10.0
c0 = MSE[length(MSE)]

expfit <- nls(MSE~f(gen,a,b,c),start=list(a=a0,b=b0,c=c0))
print(summary(expfit))

plot(gen,MSE,ylim=c(0,2),main="Exponential fit")
lines(gen,fitted(expfit),type="l",col="blue")

print("Exponential stopping time")
print(ceiling(3*coef(expfit)['b']))

g <- function(x,a,b,c,d) {
  d + a/(1+exp((x-c)/b))
}

a0 = MSE[1] - MSE[length(MSE)]
b0 = 10.0
c0 = 0.1
d0 = MSE[length(MSE)]

sigfit <- nls(MSE~g(gen,a,b,c,d),start=list(a=a0,b=b0,c=c0,d=d0))

print(summary(sigfit))
plot(gen,MSE,ylim=c(0,2),main="Sigmoidal fit")
lines(gen,fitted(sigfit),type="l",col="blue")

print("Sigmoidal stopping time")
print(ceiling(3*coef(sigfit)['b']+coef(sigfit)['c']))

