########################################################
#
# Computes LP-comoment of (X,Y), where (X,Y)
# 
#
########################################################


 LP-comoment<-function(x,y,m){

  TX <- LP.poly(x,m) #--LP-Transformation X --> TX
  TY <- LP.poly(y,m) #--LP-Transformation Y --> TY

 return(cov(TX,TY))  #--Cross-cov

   }





  LP.smooth <-function(x,y,m,method=c("AIC","BIC","KIC")){  #--output smooth-LP-comoment matrix

   TX <- LP.poly(x,m) 
   TY <- LP.poly(y,m)
   n <- length(x)
   
   CR <- c(cov(TX,TY))
   CR.s <- sort(CR^2,decreasing=TRUE,index=TRUE)$x
   aa <- rep(0,length(CR.s))

   if(method=="AIC") penalty <- 2
   if(method=="BIC") penalty <- log(n)
   if(method=="KIC") penalty <- log(n)*log(max(nrow(LP),ncol(LP)))


   aa[1] <-CR.s[1] - penalty/n
  for(i in 2: length(CR.s)){
    aa[i] <- aa[(i-1)] + (CR.s[i] - penalty/n)
    }
 # plot(aa,type="b",ylab=method,cex.axis=1.2,cex.lab=1.2)
  LP[LP^2<CR.s[which(aa==max(aa))]] <- 0
  return(LP)
   }
  