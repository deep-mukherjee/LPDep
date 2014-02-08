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
   
   TS <- as.matrix(poly(u ,degree= min(length(unique(u ))-1, m)  ))  

    for(j in 1:m){
         TS[,j] <- TS[,j]/(sqrt( ((n-1)*var(TS[,j]))/n ))
      }
 return(TS)

  }


Score.mat <-function(X,m){  #----for X random vector
  X<- as.matrix(X)
  SX <- c()
  for(j in 1:ncol(U)){
      SX <- as.matrix( cbind(SX,LP.poly(X[,j],m) ) ) 
      }
    return(SX)
   }

