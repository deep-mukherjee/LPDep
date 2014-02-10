#############################################################
#
# Construct LP-orthonormal score functions
#
# [Input X output TX]
#
##############################################################


LP.poly <- function(x,m){ #----for X random variable
    n <- length(x)
    u <- (rank(x,"keep") - .5)/n #--mid-distribution transformation
   TX <- as.matrix(poly(u ,degree= min(length(unique(u ))-1, m)  ))  

    for(j in 1:m){
         TX[,j] <- TX[,j]/(sqrt( ((n-1)*var(TX[,j]))/n ))
      }

 return(TX)
  }


Score.mat <-function(X,m){  #----for X random vector
  X<- as.matrix(X)
  TX <- c()

  for(j in 1:ncol(X)){
      TX <- as.matrix( cbind(TX,LP.poly(X[,j],m) ) ) 
      }

    return(TX)
   }

