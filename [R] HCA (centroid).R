HCA=function(data, method=c("single","complete","average"),k){
  if(method=="single"){
    single=function(data){
      n=nrow(data)
      jarak=dist(data,method="euclidean")
      proses=function(jarak){
        D=as.matrix(jarak)
        min1=which(D==min(jarak),arr.ind=TRUE)
        b1=min1[1,1]
        k1=min1[1,2]
        m=nrow(D)
        d1=matrix(nrow=m,ncol=1)
        for (i in 1:m){
          d1[i]=min(cbind(D[b1,i],D[k1,i]))
        }
        d1=d1[-c(b1,k1),]
        D2=D[-c(b1,k1),-c(b1,k1)]
        mat1=cbind(d1,D2)
        
        D=as.matrix(rbind(cbind(0,t(d1)),mat1))
        jarak=as.dist(D)
        
      }
      
      tabel=list(NA)
      for (i in 1:(n-2)){
        tabel[[i]]=matrix(nrow=(n-i),ncol=(n-i))
      }
      
      tabel[[1]]=as.matrix(proses(jarak))
      for (i in 1:(n-3)){
        tabel[[i+1]]=as.matrix(proses(as.dist(tabel[[i]])))
      }
      print(tabel[[n-2]])
    }
    single(data)
    jarak = dist(data,method="euclidean")
    dendrogram = hclust(jarak, method = "single")
    plot(dendrogram, cex = 0.5, hang = -1)
    rect.hclust(dendrogram, k = k, border = 2:k+1)
  }
  
  else if(method=="complete") {
    complete=function(data){
      n=nrow(data)
      jarak=dist(data,method="euclidean")
      proses=function(jarak){
        D=as.matrix(jarak)
        min1=which(D==min(jarak),arr.ind=TRUE)
        b1=min1[1,1]
        k1=min1[1,2]
        m=nrow(D)
        d1=matrix(nrow=m,ncol=1)
        for (i in 1:m){
          d1[i]=max(cbind(D[b1,i],D[k1,i]))
        }
        d1=d1[-c(b1,k1),]
        D2=D[-c(b1,k1),-c(b1,k1)]
        mat1=cbind(d1,D2)
        
        D=as.matrix(rbind(cbind(0,t(d1)),mat1))
        jarak=as.dist(D)
        
      }
      
      tabel=list(NA)
      for (i in 1:(n-2)){
        tabel[[i]]=matrix(nrow=(n-i),ncol=(n-i))
      }
      
      tabel[[1]]=as.matrix(proses(jarak))
      for (i in 1:(n-3)){
        tabel[[i+1]]=as.matrix(proses(as.dist(tabel[[i]])))
      }
      print(tabel[[n-2]])
    }
    complete(data)
    jarak = dist(data,method="euclidean")
    dendrogram = hclust(jarak, method = "complete")
    plot(dendrogram, cex = 0.5, hang = -1)
    rect.hclust(dendrogram, k = k, border = 2:k+1)
  }
  
  else if(method=="average"){
    average=function(data){
      n=nrow(data)
      jarak=dist(data,method="euclidean")
      proses=function(jarak){
        D=as.matrix(jarak)
        min1=which(D==min(jarak),arr.ind=TRUE)
        b1=min1[1,1]
        k1=min1[1,2]
        m=nrow(D)
        d1=matrix(nrow=m,ncol=1)
        for (i in 1:m){
          d1[i]=mean(cbind(D[b1,i],D[k1,i]))
        }
        d1=d1[-c(b1,k1),]
        D2=D[-c(b1,k1),-c(b1,k1)]
        mat1=cbind(d1,D2)
        
        D=as.matrix(rbind(cbind(0,t(d1)),mat1))
        jarak=as.dist(D)
        
      }
      
      tabel=list(NA)
      for (i in 1:(n-2)){
        tabel[[i]]=matrix(nrow=(n-i),ncol=(n-i))
      }
      
      tabel[[1]]=as.matrix(proses(jarak))
      for (i in 1:(n-3)){
        tabel[[i+1]]=as.matrix(proses(as.dist(tabel[[i]])))
      }
      print(tabel[[n-2]])
    }
    average(data)
    jarak = dist(data,method="euclidean")
    dendrogram = hclust(jarak, method = "average")
    plot(dendrogram, cex = 0.5, hang = -1)
    rect.hclust(dendrogram, k = k, border = 2:k+1)
  }
  
  else if(method=="centroid"){
    centroid=function(data){
      n=nrow(data)
      jarak=dist(data,method="euclidean")
      proses=function(data,jarak){
        data=as.matrix(data)
        D=as.matrix(jarak)
        min1=which(D==min(jarak),arr.ind=TRUE)
        b1=min1[1,1]
        k1=min1[1,2]
        
        d1=(data[b1,]+data[k1,])/2
        data1=data[-c(b1,k1),]
        
        data=rbind(d1,data1)
        jarak=dist(data,method="euclidean")
      }
      
      tabel=list(NA)
      for (i in 1:(n-2)){
        tabel[[i]]=matrix(nrow=(n-i),ncol=(n-i))
      }
      
      tabel[[1]]=as.matrix(proses(data,jarak))
      for (i in 1:(n-3)){
        tabel[[i+1]]=as.matrix(proses(as.dist(tabel[[i]]),dist(tabel[[i]])))
      }
      print(tabel[[n-2]])
    }
    centroid(data)
    jarak = dist(data,method="euclidean")
    dendrogram = hclust(jarak, method = "centroid")
    plot(dendrogram, cex = 0.5, hang = -1)
    rect.hclust(dendrogram, k = k, border = 2:k+1)
  }
}
HCA(data,method="single",k=5)
HCA(data,method="complete",k=5)
HCA(data,method="average",k=5)
HCA(data,method="centroid",k=5)
