########################################################
#
#
#   LP Dependece Number Algorithm
#
########################################################


LP.Dep <-function(x,y,m){

	Lcor <- rep(NA,3)  ##---LPINFOR,LPMax,LPCan in this order

	SX <- LP.poly(x,m)

	SY <- LP.poly(y,m)

	Lcor[3] <- cancor(SX,SY)$cor[1]

	LP <- cov(SX, SY)

	Lcor[2] <- max(abs(LP))


	Lcor[1] <- sum(LP^2)

 return(Lcor)

}



MLP.Dep <-function(X,Y,m){  #--For X,Y random vectors

   TX <- Score.mat(as.matrix(X),m)
   TY <- Score.mat(as.matrix(Y),m)

   COH <- solve(var(TX))%*%cov(TX,TY)%*%solve(var(TY))%*%cov(TY,TX)


   return(tr(COH))


 }


