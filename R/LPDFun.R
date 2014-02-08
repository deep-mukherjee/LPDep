library("orthopolynom")
library("MASS")
library("psych")  #--trace func

## ---------This function X --> FMID


Fmid <- function(X){
  U.tr <- matrix(0,nrow(X),ncol(X))
  for(j in 1:ncol(X)){
      U.tr[,j] <-  ( rank(X[,j],"keep",ties.method = c("average")) - .5)/nrow(X)  # here we have to incorporat missing observations by "keep"
  }
return(U.tr)
}


## --------This will generate the score functions
Score.fun <- function(u,m){
      U <- sort(u,index=TRUE)
	u <- U$x
	u.indx <- U$ix
       m <- min(length(unique(u ))-1, m   )
   if(m==0) {S.mat <- NULL}
   if(m > 0){
      S.mat <- as.matrix(poly(u ,df=m))  
	for(j in 1:m){
         S.mat[,j] <- S.mat[,j]/csd(S.mat[,j]) 
	  #R[j] <- cor(Y.ss,S.mat[,j])
      }
      S.mat[ u.indx ,  ] <- S.mat
    S.mat <- as.matrix(S.mat)
  }
return(S.mat)
}

#-- sd using def of n instead of (n-1)
csd <- function(x){ 
n <- length(x)
return(sqrt( ((n-1)*var(x))/n ) )}

##-------------------------------------------------------------

Score.mat <-function(X,m){
  U <- Fmid(X)
  SX <- c()
  for(j in 1:ncol(U)){
      SX <- as.matrix( cbind(SX,Score.fun(U[,j],m) ) ) 
      }
    return(SX)
   }


##------------------------------------------------------

Leg.val <- function(x,d){

     d <- min(length(unique(x))-1, d   )
     poly <-  slegendre.polynomials(d,normalized=TRUE)
     X <- matrix(NA,length(x),d)

for(j in 1:d){

    X[,j] <- predict(poly[[j+1]],x)
    }

    return(X)

}

##############################################################################
#
#
#   LP Dependece Number
#
###############################################################################


LP.Dep <-function(x,y,m.max){
Lcor <- rep(NA,3)  ##---LPINFOR,LPMax,LPCan in this order
n <- length(x)
U <- Fmid(cbind(x,y))
m1 <- min(length(unique(x))-1,m.max)
m2 <- min(length(unique(y))-1,m.max)

SX <- Score.fun(U[,1],m1)

SY <- Score.fun(U[,2],m2)

Lcor[3] <- cancor(SX,SY)$cor[1]

LP <- cov(SX, SY)

Lcor[2] <- max(abs(LP))

#--LP[abs(LP) < sqrt(2*log(n)/n)] <- 0  #--for simulation smoothing not required

Lcor[1] <- sum(LP^2)

 return(Lcor)

}

##############################################################################
#
#
#   Diagnosics for Bivariate normal or short-tailedness
#
###############################################################################

LP.BN <- function(x,y){

  ux <- ecdf(x)(x) 
  uy <- ecdf(y)(y) 
  
  R <- cor(x,y)

  T1x  <- (ux -.5)/sd(ux)
  T1y <- (uy -.5)/sd(uy)


  DD <- R - cor(x,T1x)*cor(T1x,T1y)*cor(y,T1y)

  return(DD)

}



LP.BN.Test <- function(x,y){

   B <- 500
   ss <- cov(x,y)
   Si.hat <- matrix(c(var(x),ss,ss,var(y)),byrow=TRUE,2,2)

   lp <- rep(NA,B)

   for(i in 1:B){

     T <- mvrnorm(n=length(x),mu=c(0,0),Sigma=Si.hat)

     lp[i] <- LP.BN(T[,1],T[,2])

      }


     lp.val <- LP.BN(x,y)

     pval <-  sum(lp>abs(lp.val))/length(lp)

    return(pval)

    }

##############################################################################
#
#
#   MULTIVARIATE LP Dependece Number
#
###############################################################################

MLP.Dep <-function(X,Y,m.max){

   TX <- Score.mat(as.matrix(X),m.max)
   TY <- Score.mat(as.matrix(Y),m.max)

   COH <- solve(var(TX))%*%cov(TX,TY)%*%solve(var(TY))%*%cov(TY,TX)


   return(tr(COH))


 }


MLP.Dep.pval <-function(X,Y,m.max){

   B <- 1500
   TX <- Score.mat(as.matrix(X),m.max)
   TY <- Score.mat(as.matrix(Y),m.max)

   COH <- solve(var(TX))%*%cov(TX,TY)%*%solve(var(TY))%*%cov(TY,TX)

   lpcoh<- tr(COH)

    lpcoh.null <- rep(NA,B)

   for(k in 1:B){
   
   TY <- Score.mat(as.matrix(Y[sample(1:nrow(Y)),]),m.max)
   lpcoh.null[k] <-  tr( solve(var(TX))%*%cov(TX,TY)%*%solve(var(TY))%*%cov(TY,TX)   )

               }

      p.val <- length(lpcoh.null[lpcoh.null>lpcoh])/B

   return(p.val)


 }




MLP.Dep.max <-function(X,Y,m.max){

   TX <- Score.mat(as.matrix(X),m.max)
   TY <- Score.mat(as.matrix(Y),m.max)

   COH <- solve(var(TX))%*%cov(TX,TY)%*%solve(var(TY))%*%cov(TY,TX)


   return(svd(COH)$d[1])


 }

MLP.Dep.max.pvalue <-function(X,Y,m.max){

   B <- 500
   TX <- Score.mat(as.matrix(X),m.max)
   TY <- Score.mat(as.matrix(Y),m.max)

   COH <- solve(var(TX))%*%cov(TX,TY)%*%solve(var(TY))%*%cov(TY,TX)

   lpcoh<- svd(COH)$d[1]

   lpcoh.null <- rep(NA,B)

   for(k in 1:B){
   
   TY <- Score.mat(as.matrix(Y[sample(1:nrow(Y)),]),m.max)
   lpcoh.null[k] <-  svd(solve(var(TX))%*%cov(TX,TY)%*%solve(var(TY))%*%cov(TY,TX))$d[1] 

               }

      p.val <- length(lpcoh.null[lpcoh.null>lpcoh])/B
 
   return(p.val)


 }




##############################################################################
#
#
#   LP Regression
#
###############################################################################

LP.mul.reg <-function(y,X,m.max){

  Tx <- Score.mat(as.matrix(X),m.max)
   LR <-  lm(y~Tx)
  tt <- summary(LR)$coefficients[-1,3]

 
 return(tt)

   }





LP.reg <-function(y,x){

 sx <- sort(x,index.return=TRUE)

 x <- sx$x

 y <- y[sx$ix]

 U <-Fmid(cbind(x,y))

 Tx <- Score.fun(U[,1],4)

 Ty <- Score.fun(U[,2],4)

 LR <-  lm(y~Tx)

 ZLP <- LR$coef[-1]

 pval <- summary(LR)$coefficients[-1,4]

 ZLP[pval > .01] <- 0


 yhat <- mean(y) +  Tx%*%ZLP

 par(mfrow = c(1,2), cex = 0.45)

#--generate plots

u <- ecdf(x)(x)

plot(u,y,xlab="F(X)",ylab="Y",main="LP Smooth Regression: Quantile domain",cex.lab=1.8,cex.main=1.6,cex.axis=1.5)
lines(u,yhat,col="blue",lwd=1.5)



plot(x,y,xlab="X",ylab="Y",main="LP Smooth Regression: Distribution domain",cex.lab=1.8,cex.main=1.6,cex.axis=1.5)
lines(x,yhat,col="blue",lwd=1.5)



 LR <- list()

 LR$coef <- ZLP

 LR$fit <- cbind(x,y,yhat)


return(LR)


   }

######################################################
#
#   KIC.find 
#
#######################################################


###----performs BIC regularizing 

BIC.find <- function(CR,n,m){
  CR.s <- sort(CR^2,decreasing=TRUE,index=TRUE)$x
  aa <- rep(0,length(CR.s))
  penalty <- log(n)*log(m)
  aa[1] <-CR.s[1] - penalty/n

  for(i in 2: length(CR.s)){
    aa[i] <- aa[(i-1)] + (CR.s[i] - penalty/n)
    }
CR.ind <- sort(abs(CR),decreasing=TRUE,index=TRUE)$ix
max.ind <- which(aa==max(aa))
B<- list()
B$ind <-  CR.ind[1:max.ind]
B$BIC <- aa
B$screen <- CR[CR.ind[1:max.ind]]
return(B)
}


