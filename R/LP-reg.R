#######################################
#
#   LP Regression
#
#######################################



LP.reg <-function(y,X,m){

  Tx <- Score.mat(as.matrix(X),m)
   LR <-  lm(y~Tx)
  tt <- summary(LR)$coefficients[-1,3]

 
 return(tt)

   }


#################################################
#  Turlachs Example
#################################################

Turlach<-function(n,p){

X <- matrix(runif(n*p),n,p)

y <- (X[,1] - .5)^2 + apply(X[,2:5],1,sum) + rnorm(n,sd=.05)

return(list(y=y,X=X))
}


T <- Turlach(n=500,p=50)
y <-T$y
X<-T$X

Tx <- Score.mat(as.matrix(X),2)
fit1 <- glmnet(Tx,y)
plot(fit1,label=TRUE)




