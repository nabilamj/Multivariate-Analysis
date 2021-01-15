data=read.csv("DATA PCA.csv",header=T)

x1=as.matrix(data$X1)
x2=as.matrix(data$X2)
x3=as.matrix(data$X3)
x4=as.matrix(data$X4)
x5=as.matrix(data$X5)
x6=as.matrix(data$X6)
x7=as.matrix(data$X7)

library(mvnTest)
mvn(data,mvnTest="royston")
R.test(data)
bartlett.test(data)

PCA=function(data, standardize=FALSE){
  if(standardize == TRUE){
    data = scale(data)
  }
  
  S = cov(data)
  
  eigen_val = eigen(S)$values
  
  eigen_vec = eigen(S)$vector
  
  n = length(eigen_val)
  prop = c()
  for (i in (1:n)){
    prop[i] = (eigen_val[i])/sum(eigen_val)
  }
  
  q = length(prop)
  propcum = c()
  for (i in 1:q){
    propcum[i]=sum(prop[1:i])
  }
  
  p = nrow(eigen_vec)
  corr = matrix(0,p,p)
  for (i in (1:p)){
    for (j in (1:p)){
      corr[i,j] = (eigen_vec[i,j]*sqrt(eigen_val[i]))/S[j,j]
    }
  }
  
  hasil = list("Matriks Varkov/Korelasi"=S,"Eigen Value"=eigen_val,"Eigen Vector"=eigen_vec,"Proporsi"=prop,"Proporsi Kumulatif"=propcum,"Matriks Korelasi Y dan X"=corr)
  print(hasil)
}

PCA(data, standardize = F)
