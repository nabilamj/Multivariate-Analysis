data=read.csv("DATA PCA.csv",header=T,sep=",")

FactorAnalysis=function(data, standardize=FALSE){
  
  if(standardize == TRUE){
    data = scale(data)
  }
  
  S = cov(data)
  
  E = eigen(S)$values
  A = diag(E)
  V = eigen(S)$vector
  
  trace_S=sum(diag(S))
  prop=E/trace_S
  q = length(prop)
  propcum = c()
  for (i in 1:q){
    propcum[i]=sum(prop[1:i])
  }
  
  L=V%*%sqrt(A)
  
  miu=colMeans(data)
  n=nrow(data)
  p=ncol(data)
  Xc=matrix(0,n,p)
  for (i in (1:n)){
    for (j in (1:p)){
      Xc[i,j]=data[i,j]-miu[j]
    }
  }
  
  F=Xc%*%solve(S)%*%L
  
  plot(eigen_val, main="Scree Plot", type="o")
  
  hasil = list("Matriks Varkov/Korelasi"=S, "Eigen Value"=E, 
               "Eigen Vector"=V, "Proporsi Faktor"=prop, 
               "Proporsi Kumulatif Faktor"=propcum, 
               "Loading Factor"=L, "Factor Score"=F)
  print(hasil)
}

FactorAnalysis(data,standardize = T)
  
  