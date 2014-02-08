#############################################################
#
# Construct LP-orthonormal score functions
#
# [Input X output TX]
#
##############################################################


LP.poly <- function(x,m){

    n <- length(x)
  
    u <- (rank(x,"keep") - .5)/n #--mid-distribution transformation
   
   TS <- as.matrix(poly(u ,degree= min(length(unique(u ))-1, m)  ))  

    for(j in 1:m){
         TS[,j] <- TS[,j]/(sqrt( ((n-1)*var(TS[,j]))/n ))
      }

 return(TS)


  }