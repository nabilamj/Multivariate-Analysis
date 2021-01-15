data=read.csv("diskriminan1.csv",header=T)
data

Dis2=function(data){
  
  Y=data[,1]
  k1=data[which(Y==1),-1]
  k2=data[which(Y==2),-1]
  
  n=nrow(data)
  p=ncol(data)-1
  n1=nrow(k1)
  n2=nrow(k2)
  mu1=colMeans(k1)
  mu2=colMeans(k2)
  mus=mu1-mu2
  S1=cov(k1)
  S2=cov(k2)
  Spl=((n1-1)*S1+(n2-1)*S2)/(n1+n2-2)
  Splinv=solve(Spl)
  
  a=Splinv%*%mus
  
  k1b=matrix(0,n1,p)
  k2b=matrix(0,n2,p)
  for (j in (1:p)){
    k1b[,j] = a[j]*k1[,j]
    k2b[,j] = a[j]*k2[,j]
  }
  z1=matrix(0,n1,1)
  for (i in (1:n1)){
    z1[i,]=sum(k1b[i,])
  }
  z2=matrix(0,n2,1)
  for (i in (1:n2)){
    z2[i,]=sum(k2b[i,])
  }
  
  muz1=mean(z1)
  muz2=mean(z2)
  m=(n1*muz2+n2*muz1)/(n1+n2)
  
  kat1=matrix(0,n1,1)
  for (i in (1:n1)){
    kat1[i]=if(z1[i] >= m) 1 else 2
  }
  kat2=matrix(0,n2,1)
  for (i in (1:n2)){
    kat2[i]=if(z2[i] >= m) 1 else 2
  }
  
  b=length(kat1[which(kat1==1)])
  c=length(kat1[which(kat1==2)])
  d=length(kat2[which(kat2==1)])
  e=length(kat2[which(kat2==2)])
  
  Kes=matrix(c(b,c,d,e),2,2,byrow = T)
  dimnames(Kes)=list(Estimasi=c("K1","K2"),Observasi=c("K1","K2"))
  
  hitratio=(b+e)/n
  
  hasil=list("Koefisien Kombinasi Linier Fisher"=a,"Matriks Estimasi dan Observasi"=Kes,
             "Hit Ratio"=hitratio)
  print(hasil)
}

Dis2(data)
