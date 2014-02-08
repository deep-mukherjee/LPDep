#############################################################
#
#
# Survival: Stanford Heart Transplant data
#
##############################################################

source("LPDFun.R")
libray("survival")
data(stanford2)

S<-stanford2

S<-na.omit(S)

Y<-S[,2:3]
X<-S[,4:5]



##--(i) computing marginal scores for AGE and T5

m <- 4

LP.score <- rep(NA,ncol(X))

for(j in 1:ncol(X)){

LP.score[j] <- MLP.Dep.pval(X[,j],Y,m) 


 }




#--(ii) Spectral analysis of COH(Y,Age)

m.max<-7
TX <- Score.mat(as.matrix(X[,1]),m.max)
TY <- Score.mat(as.matrix(Y),m.max)
COH <- solve(var(TX))%*%cov(TX,TY)%*%solve(var(TY))%*%cov(TY,TX)


pdf("suvival1.pdf")

plot(svd(COH)$d,type="b",cex=1.7,ylab="Singular Values",cex.axis=1.2,cex.lab=1.24,xlab="")
dev.off()


ss<- c(expression(T[1]),expression(T[2]),expression(T[3]),expression(T[4]),expression(T[5]),expression(T[6]),expression(T[7]))

pdf("suvival2.pdf")

plot(abs(svd(COH)$v[,1]),type="b",lwd=1.6,col = "blue",,xaxt='n',lty=1,cex.axis=1.24,cex.lab=1.24,ylab="",xlab="",ylim=c(0,.87))
axis(1,at=1:7,labels=ss,cex.axis=1.4)
points(abs(svd(COH)$v[,2]), col = "red", type = 'b',lwd=1.6,lty=2)
legend("topright",c("First loading vector","Second loading vector"), lty =1:2, col = c("blue","red"),cex=1.2)
dev.off()



