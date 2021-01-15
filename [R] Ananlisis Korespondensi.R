data=matrix(c(33,25,71,30),2,2)
data

data=matrix(c(82,11,6,18,16,26),2,3)
data

data=matrix(c(156,124,77,14,20,11,2,5,7),3,3)
data

data=matrix(c(178,570,138,138,648,252,108,442,252),3,3)
rownames(data)=c("Less than HS","HS/junior college","Bachelor/Graduate")
colnames(data)=c("Fundamentalist","Moderate","Liberal")
data

data=matrix(c(714,730,498,221,33,425,68,17,320,813,1072,142,284,276,325,188),4,4)
data

data.awal=read.csv("putussekolahdki2018.csv",header=T)
data.awal
data=data.awal[,-1]
rownames(data)=data.awal[,1]
data

chisq.test(data)
library(ca)
hasil=ca(data)
plot(hasil)
hasil

korespondensi=function(data){
  i=nrow(data)
  j=ncol(data)
  k=min(i,j)-1
  
  P=data/sum(data)
  r=as.matrix(rowSums(P))
  R=diag(rowSums(P))
  c=as.matrix(colSums(P))
  C=diag(colSums(P))
  S=solve(sqrt(R))%*%as.matrix((P-r%*%t(c)))%*%solve(sqrt(C))
  
  if (k==1){
    if (i <= j){
      SSt=S%*%t(S)
      lamda=round(eigen(SSt)$value,5)
      A=lamda[which(lamda!=0)]
      D=diag(sqrt(A)[1:k])
      U=eigen(SSt)$vector[,1:k]
      V=t(solve(D)%*%t(U)%*%S)
    }
    
    else if (i > j){
      StS=t(S)%*%S
      lamda=round(eigen(StS)$value,5)
      A=lamda[which(lamda!=0)]
      D=sqrt(A)[1]
      V=eigen(StS)$vector[,1]
      U=S%*%V%*%solve(D)
    }
  }
  
  else if(k>1){
    if (i <= j){
      SSt=S%*%t(S)
      lamda=round(eigen(SSt)$value,5)
      A=lamda[which(lamda!=0)]
      D=diag(sqrt(A)[1:k])
      U=eigen(SSt)$vector[,1:k]
      V=t(solve(D)%*%t(U)%*%S)
    }
    
    else if (i > j){
      StS=t(S)%*%S
      lamda=round(eigen(StS)$value,5)
      A=lamda[which(lamda!=0)]
      D=diag(sqrt(A)[1:k])
      V=eigen(StS)$vector[,1:k]
      U=S%*%V%*%solve(D)
    }
  }
  
  Y=solve(sqrt(R))%*%U%*%D
  Z=solve(sqrt(C))%*%V%*%D
  satu=matrix(c(1),1,k)
  tau=solve(satu%*%as.matrix(A))%*%t(as.matrix(A))
  
  
  hasil=list("Estimasi koordinat utama dari baris (Y)"=Y, 
             "Estimasi koordinat utama dari kolom (Z)"=Z,
             "Persentase Keragaman"=tau)
  
  print(hasil)
}

korespondensi(data)

install.packages("rgl")
library(ca)
library(rgl)
hasil=ca(data)
plot(hasil)
plot3d.ca(hasil)


library("FactoMineR")
library("factoextra")
res.ca=CA(data,ncp=2)
fviz_ca_biplot(res.ca)
data(housetasks)
res.ca <- CA(data, graph = FALSE)
fviz_ca_biplot(res.ca, repel = TRUE)

dist(data)
