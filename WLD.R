#In this project we study the behavior of a newely establishe Weibull-Lindley Distribution
install.packages('ggplot2')
install.packages('plyr')
install.packages('dplyr')
install.packages('ggthemes')
install.packages('MASS')
install.packages('car')
instal.packages('')
help("dgamma")
help("curve")
curve(dgamma(x,1,1),0,10)
for (i in 2:6) {curve(dgamma(x,i,1),0,10,add = TRUE,legend())
  
}
install.packages("fBasics")

#Fitting a normal distribution
X=rnorm(200,12,4)
hist(X)
density(X)
plot(density(X))
Y=(X-mean(X))/sd(X)
qqnorm(Y)
abline(0,1)
W1=rweibull(200,2.2,3)
W2=rweibull(100,0.6,3)
qqplot(W2,W1)
abline(0,1)
curve(dgamma(x,0.4,2),0,7)
for(i in 2:6){
  curve(dgamma(x,i,i-1),add = TRUE)
}
install.packages("stats4")
load(stats4)
load("stats4")
library(stats4)
mle(Y)
X=rnorm(1000,100,15)
hist(X,probability = TRUE)
XX=seq(min(X),max(X),length=100)
lines(XX,dnorm(XX,100,15))#Create pdf

Q=rgamma(1000,2,0.5)
hist(Q,probability = TRUE)
Q1=seq(min(Q),max(Q),200)
lines(Q1,dgamma(Q1,2,0.5))
W=qbinom(seq(0.1,0.9,0.1),10,1/3)
hist(W,probability = T)
curve(dgamma(x,0.98,3),0,10)
for(i in 2:5){
  curve(dgamma(x,i,3),0,10,add = TRUE)
}

curve(dbeta(x,1,11),0.1,0.9)
for(i in 2:10){
  curve(dbeta(x,i,12-i),0.1,0.9,add = TRUE,col=i)
}

curve(dnorm(x,0,5),-4,4)
for(i in 1:6){
  curve(dnorm(x,0,i),-4,4,add = TRUE)
}

curve(dnig(x,1,0.4,1),0,10)
install.packages("rmutil")
library(rmutil)
curve(dinvgauss(x,1,1),0.1,10,col="yellow",main="invgauss",
      xlab = "Time(Ages)",ylab = "Density",
      legend("topright",legend = c("shape=c(1,2,3,4,5)","scale=4")))
for(i in c(3,4,5)){
  curve(dinvgauss(x,i,1),0.1,10,col=i,add = T)
}
library(rmutil)
x2=seq(0.1,5,0.01)
y1=dinvgauss(x2,1,1)
y2=dinvgauss(x2,3,1)
y3=dinvgauss(x2,4,1)
y4=dinvgauss(x2,5,1)
xrange1=range(x2)
yrange1=range(y1,y2,y3,y4)
plot(x2,y1,xlab = "LogAges",ylab = "PDF",col="black",type = "l",lwd=2)
main="The inverse Gaussian pdf curve for different shape and scale parameters"
lines(x2,y2,col="green",lwd=2)
lines(x2,y3,col="red",lwd=2)
lines(x2,y4,col="orange",lwd=2)
legend("topright",legend=c("shape=1,scale=1","shape=3,scale=1",
                           "shape=4,scale=1","shape=5,scale=1"),
       col=c("black","green","red","orange"),lty=1)


pdf("Gamma.pdf")
curve(dgamma(x,1,4),0,3,col="yellow",
      main="Gamma density plot for different parameter values",
      xlab = "logTime(Ages)",ylab = "Density",ylim=c(0,2),lwd=2)
for(i in c(3,4,5)){
  curve(dgamma(x,i,4),0,3,col=i,lwd=c(2,2,2),add = T)
}
dev.off()

curve(dinvgauss(x,1,1),0.1,10,col="yellow",main="invgauss",
      xlab = "Time(Ages)",ylab = "Density",
      legend("topright",legend = c("shape=c(1,2,3,4,5)","scale=4")))
for(i in c(3,4,5)){
  curve(dinvgauss(x,i,1),0.1,10,col=i,add = T)
}

pdf("Gamma1.pdf")
curve(dgamma(x,4,1),0,3,col="yellow",
      main="Gamma density plot for different parameter values",
      xlab = "logTime(Ages)",ylab = "Density",ylim=c(0,2),lwd=2)
for(i in c(3,4,5)){
  curve(dgamma(x,4,i),0,3,col=i,lwd=c(2,2,2),add = T)
}
dev.off()

pdf("Invgauss5")
curve(dinvgauss(x,0.5,1),0.01,4,col="yellow",
      main="Inverse Gaussian density plot for different parameters",
      xlab = "logTime(Ages)",ylab = "Density",lwd=2,ylim=c(0,2))
for(i in c(1,2)){
  curve(dinvgauss(x,i,1),0.01,4,col=i,add = T,lwd=c(2,2,2))
}
dev.off()
P1=rnorm(200,10,2)
Q1=(P1-mean(P1)/sd(P1))
qqnorm(Q1)
abline(0,1)
P2=rgamma(200,2,1.5)
P3=rgamma(200,2,1.3)
qqplot(P2,P2)
abline(0,1)
curve(dnorm(x,8,2),0,20,
      xlab = "TTTT",ylab = "JJJJ")
for(i in 1:4){
  curve(dnorm(x,10-i,2),add = T,col=i,lwd=c(2,2,2))
}
install.packages("Newdistins")
library(AdequacyModel)
library(Newdistns)
install.packages("AdequacyModel")
install.packages("MPS")
library(MPS)
library(Newdistns)
#######Generalized distribution by Newdstis#####
x=seq(0.1,10,0.1)
#y1=dbetag(x,"gamma",a=2,b=2,shape=1,scale=1)
y1=dbetag(x,"gamma",a=2,b=3,shape=1,scale=1)
y2=dbetag(x,"gamma",a=2,b=3,shape=2,scale=2)
y3=dbetag(x,"gamma",a=2,b=3,shape=3,scale=1)
y4=dbetag(x,"gamma",a=2,b=3,shape=5,scale=2)
######Plotting the pdf of beta gamma distri#####
xrange=range(x)
yrange=range(y1,y2,y3,y4)
plot(x,y1,xlab="x",ylab="pdf",xlim = xrange,ylim = yrange,col="black",type="l")
lines(x,y2,col="red")
lines(x,y3,col="green")
lines(x,y4,col="pink")
legend("topright",legend=c("shape=1,scale=1","shape=2,scale=2",
                           "shape=3,scale=1","shape=5,scale=2"),
       col=c("black","red","green","pink"),lty = 1)




func1=function(z1,z2){
  volt=sqrt(z1^2+z2^2)
}
z1=z2=seq(1,5,length.out = 50)
z=outer(z1,z2,func1)
persp(z1,z2,z)

func2=function(c1,betta){
  h1=qweibullg(c(6/8,4/8,2/8,7/8,5/8,3/8,1/8),"lindley",beta = betta,c=c1,theta=0.5)
  K1=(h1[1]-2*h1[2]+h1[3])/(h1[1]-h1[3])
  S1=(h1[4]-h1[5]+h1[6]+h1[7])/(h1[1]-h1[3])
}
c1=betta=seq(1,5,length.out = 50)
z2=outer(c1,betta,func2)
persp(c1,betta,z2)

