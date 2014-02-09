#########################################################
#  LP-CCA
#  Example: breast cancer data: gene expression + DNA
#
##########################################################

library("PMA")

data(breastdata)
attach(breastdata)
Y<-dna <- t(dna)
X<- rna <- t(rna)
m <- 1
TX <- Score.mat(as.matrix(X),m)

perm.out <- CCA.permute(x=TX,z=Y[,chrom==1],typex="standard", typez="ordered",nperms=5,penaltyxs=seq(.02,.2,length=20))

print(perm.out)

ppx <- .029  #--such that we select ~20 non-zero u's to compare with the result of Witten et al. (2009)


out <- CCA(x=TX,z=Y[,chrom==1], typex="standard", typez="ordered",penaltyx=ppx,
v=perm.out$v.init, penaltyz=perm.out$bestpenaltyz, xnames=substr(genedesc,1,20),
znames=paste("Pos", sep="", nuc[chrom==1]))

#---details of the top 20 genes which are associated with DNA at chrom=1


ind <- sort(abs(out$u),decreasing=TRUE,index.return=TRUE)$ix[1:20]
breastdata$genenames[ind]
breastdata$genechr[ind]
breastdata$genedesc[ind]


#--run the same module with m=2, to see whether we have any hidden nonlinear multivariate effect.

