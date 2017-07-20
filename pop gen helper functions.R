# Helper functions for pop gen stats --------------------------------------

getp<-function(data, pops){
  pqloci<-lapply(data, getpqmatrix, pops)
  return(pqloci)
}

getpqmatrix<-function(data, pops){
  popdat<-split(data, pops)  
  genomat<-lapply(popdat, count.geno)
  p<-sapply(genomat, function(x) (x[[1]]*2+x[[2]])/((x[[1]]+x[[2]]+x[[3]])*2))
  q<-sapply(genomat, function(x) (x[[3]]*2+x[[2]])/((x[[1]]+x[[2]]+x[[3]])*2))
  pq<-rbind(p, q)
  return(pq)
}

count.geno<-function(locus){
  nhom1<-sum(locus==0)
  nhom2<-sum(locus==2)
  nhet<-sum(locus==1)
  results<-cbind(nhom1, nhet, nhom2)
  colnames(results)<-c("hom1","het", "hom2")
  return(results)
}

calc.Hobs<-function(input){
  HoTable<-lapply(input, count.geno)
  n<-nrow(input)
  HoSums<-sapply(HoTable, function(x) x[[2]]/n)
  return(HoSums)
}

calc.Hobs.pop<-function(input, pops){
  datasplit<-split(input, pops)
  popHo<-lapply(datasplit, calc.Hobs)
  return(popHo)
}

is.polymorphic<-function(locus){
  nhom1<-sum(locus==0)
  nhom2<-sum(locus==2)
  nhet<-sum(locus==1)
  if ((nhom1==0 & nhet==0) | (nhom2==0 & nhet==0)){
    result = "mono"
  } else {
    result = "poly"
  }
  return(result)
}

count.polymorphic<-function(input){
  record.poly<-sapply(input, is.polymorphic)
  poly<-sum(record.poly=="poly")
  return(poly)
}

calc.Hexp<-function(input){
  HeTable<-lapply(input, count.geno)
  n<-nrow(input)
  p2<-sapply(HeTable, function(x) ((x[[1]]*2+x[[2]])/(n*2))^2)
  q2<-sapply(HeTable, function(x) ((x[[3]]*2+x[[2]])/(n*2))^2)
  Hexp<-1-p2-q2
  return(Hexp)
}

calc.Hexp.pop<-function(input, pops){
  datasplit<-split(input, pops)
  popHexp<-lapply(datasplit, calc.Hexp)
  return(popHexp)
}

### get p and q ###
getpq<-function(genomat, npop){
  plist<-vector()
  qlist<-vector()
  for (i in 1:length(genomat)){
    h1<-genomat[[i]][[1]]
    h2<-genomat[[i]][[2]]
    h3<-genomat[[i]][[3]]
    p<-(h1*2+h2)/((h1+h2+h3)*2)
    q<-(h3*2+h1)/((h1+h2+h3)*2) 
    plist[i]<-(p/npop)^2
    qlist[i]<-(q/npop)^2
  }
  return(list("p" = plist, "q" = qlist))
}
