
test <- read.table("gccss6test.dat",header=TRUE)

gen <- test$Generation

plot(gen,test$MSE,type="l",ylim=c(0,0.1))
abline(h=0.04,col="red")

cuts = 200*c(0:20)
abline(v=cuts,col="red")

dp <- abs(diff(test$Obs_Plas))
de <- abs(diff(test$MSE))

#plot(abs(dp)[0:200],ylim=c(0,0.02),type="l")
plot(dp[0:200],type="l")

abline(h=0.01,col="red")
plot(de[0:200],type="l")


plot(de,ylim=c(0,0.1),type="l")
lines(dp,ylim=c(0,0.1),type="l",col="blue")
abline(h=0.01,col="red")

plot(abs(diff(test$EMA_Pl[1801:1999])),type="l")
axis(side=1,at=10*c(0:20))
abline(v=10*c(0:20),h=0.05*c(0:200),lty=3,col="grey")
axis(side=2,at=0.05*c(0:200))
abline(h=0.01,col="red")

f <- function(x,a,b,c){
  c + a*exp(-x/b)
}


gen <- seq(1:200)
MSE <- test$MSE[3401:3600]
a0 = 0.05
c0 = MSE[1]-0.05

fit <- nls(MSE~f(gen,a,b,c),start=list(a=1.5,b=6,c=0.05))


