#####################################################################
#
#
# Zelterman Sparse Table data analysis
#
#
#####################################################################

library("psych") # convert table to matrix
source("LPDFun.R")


M<- Z
rownames(M) <-seq(1:nrow(Z))
colnames(M) <-seq(1:ncol(Z))

 n<- margin.table(Z)

M <- table2matrix(M) 

 U <-Fmid(M)

 Tx <- Score.fun(U[,1],4)

 Ty <- Score.fun(U[,2],4)

 LP <- cov(Tx,Ty)