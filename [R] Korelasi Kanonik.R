canocor = function(x,y, standardize = FALSE){
  p=ncol(y)
  q=ncol(x)
  k = min(p,q)
  data=as.matrix(data.frame(y,x))
  
  if(standardize == TRUE){
    data = scale(data)
  }
  S = cov(data)
  s11 = S[1:p, 1:p]
  s12 = S[1:p, -(1:p)]
  s21 = S[-(1:p), 1:p]
  s22 = S[-(1:p),-(1:p)]
  msrs11=eigen(s11)$vectors%*%diag(1/sqrt(eigen(s11)$value))%*%t(eigen(s11)$vectors)
  msrs22=eigen(s22)$vectors%*%diag(1/sqrt(eigen(s22)$value))%*%t(eigen(s22)$vectors)
  #msr artinya min square root matriks
  rhomax=msrs11%*%s12%*%solve(s22)%*%s21%*%msrs11
  d=eigen(rhomax)$vectors
  c=msrs22%*%s21%*%msrs11%*%d
  
  # corr=diag((t(c)%*%msrs22%*%s21%*%msrs11%*%d)/(sqrt(t(c)%*%c)%*%sqrt(t(d)%*%d)))
  rhosq=solve(s11)%*%s12%*%solve(s22)%*%s21
  ev_rhosq=eigen(rhosq)$values
  for (i in 1:k){
    corr[i]=sqrt(ev_rhosq[i])
  }
  
  hasil_u = matrix(0,p,k)
  for(i in (1:k)){
    hasil_u[,i] = msrs11%*%d[,i]
  }
  
  #hasil_v = matrix(0,k,k)
  #for (i in (1:k)){
  #  b1 = (solve(s22))%*%s21%*%hasil_u[,i]
  #  b2 = t(b1)%*%s22%*%b1
  #  b11 = 1/sqrt(b2[1])*b1
  #  hasil_v[,i] = c(b11)
  #}
  
  b1=matrix(0,q,k)
  for (i in (1:k)){
    b1[,i] = (solve(s22))%*%s21%*%hasil_u[,i]
  }
  
  b2 = matrix(0,k,1)
  for (i in (1:k)){
    b = b1[,i]
    b2[i,] = t(b)%*%s22%*%b
  }
  
  b11 = matrix(0,q,k)
  for (j in (1:k)){
    b11[,j] = (1/sqrt(b2[j,]))*b1[,j]
  }
  
  return(list("Korelasi Kanonik"=corr, "Koefisien x"=b11, "Koefisien y"=hasil_u))
}
