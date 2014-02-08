#######################################
#
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