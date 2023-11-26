

#require foreach and doparallel for parallel processing
if (!require("foreach")) install.packages("foreach")
library(foreach)
if (!require("doParallel")) install.packages("doParallel")
library(doParallel)

numCores <- detectCores()
#registering clusters, can set a smaller number using numCores-1 

registerDoParallel(numCores)

#require randtoolbox for random number generations
if (!require("randtoolbox")) install.packages("randtoolbox")
library(randtoolbox)
#require Rfast for faster computation
if (!require("Rfast")) install.packages("Rfast")
library(Rfast)
if (!require("gnorm")) install.packages("gnorm")
library(gnorm)
if (!require("stats")) install.packages("stats")
library(stats)
if (!require("np")) install.packages("np")
library(np)
if (!require("nloptr")) install.packages("nloptr")
library(nloptr)

factdivide<-function(n1,n2){
  decin1<-n1-floor(n1)
  if(decin1==0){decin1=1}
  decin2<-n2-floor(n2)
  if(decin2==0){decin2=1}
  n1seq<-seq(decin1,n1, by=1)
  n2seq<-seq(decin2,n2,by=1)
  all<-list(n1seq,n2seq)
  maxlen <- max(lengths(all))
  all2 <- as.data.frame(lapply(all, function(lst) c(lst, rep(1, maxlen - length(lst)))))
  division<-all2[,1]/all2[,2]
  answer<-exp(sum(log(division)))*(gamma(decin1)/gamma(decin2))
  return(answer)
}
unbiasedsd<-function (x){
  n<-length(x)
  if(n==1){
    return(1000000)
  }
  sd1<-sd(x)
  if(n==2){
    c1<-0.7978845608
  }else if (n==3){
    c1<-0.8862269255 
  }else if (n==4){
    c1<-0.9213177319 
  }else{
    c1<-sqrt(2/(n-1))*factdivide(n1=((n/2)-1),n2=(((n-1)/2)-1))
  }
  listall<-sd1*c1
  (listall)
}
correctfactor<-function (n){
  if(n==2){
    c1<-0.7978845608
  }else if (n==3){
    c1<-0.8862269255 
  }else if (n==4){
    c1<-0.9213177319 
  }else{
    c1<-sqrt(2/(n-1))*factdivide(n1=((n/2)-1),n2=(((n-1)/2)-1))
  }
  listall<-c1
  (listall)
}
#load the deterministic simulation functions of 9 common unimodal distributions
dsexp<-function (uni,location=1,scale=1) {
  sample1<-qexp(uni,rate=scale)
  sample1
}
dsRayleigh<-function (uni,location=1,scale=1) {
  sample1 <- scale * sqrt(-2 * log((uni)))
  sample1[scale <= 0] <- NaN
  rev(sample1)
}
dsnorm<-function (uni,location=0,scale=1) {
  sample1<-qnorm(uni,mean =location,sd=scale)
  sample1
}
dsLaplace<-function (uni,location=0,scale=1) {
  sample1<-location - sign(uni - 0.5) * scale * (log(2) + ifelse(uni < 0.5, log(uni), log1p(-uni)))
  sample1
}
dslogis<-function (uni,location=0,scale=1) {
  sample1<-qlogis(uni,location=location,scale=scale)
  sample1
}
dsPareto<-function (uni,shape,scale=1) {
  sample1 <- scale*((uni))^(-1/shape)
  sample1[scale <= 0] <- NaN
  sample1[shape <= 0] <- NaN
  rev(sample1)
}
dslnorm<-function (uni,location=0,scale) {
  sample1 <- qlnorm(uni,meanlog=location,sdlog = scale)
  sample1
}
dsgamma<-function (uni,shape,scale = 1) {
  sample1<-qgamma(uni,shape=shape,scale=scale)
  sample1
}
dsWeibull<-function (uni,shape, scale = 1){
  sample1<-qweibull(uni,shape=shape, scale = scale)
  sample1
}
dsgnorm<-function (uni,location,shape, scale = 1){
  sample1<-qgnorm(p=uni, mu = location, alpha = scale, beta = shape)
  sample1
}

dsbeta<-function (uni,shape1,shape2) {
  sample1<-qbeta(uni,shape1=shape1,shape2=shape2)
  sample1
}

#moments for checking the accuracy of bootstrap 
moments<-function (x){
  n<-length(x)
  m1<-mean(x)
  var1<-(sum((x - m1)^2)/n)
  tm1<-(sum((x - m1)^3)/n)
  fm1<-(sum((x - m1)^4)/(n))
  listall<-c(mean=m1,variance=var1,tm=tm1,fm=fm1)
  (listall)
}
unbiasedmoments<-function (x){
  n<-length(x)
  m1<-mean(x)
  var1<-sd(x)^2
  var2<-(sum((x - m1)^2)/n)
  tm1<-(sum((x - m1)^3)/n)*(n^2/((n-1)*(n-2)))
  fm1<-(sum((x - m1)^4)/n)
  ufm1<--3*var2^2*(2*n-3)*n/((n-1)*(n-2)*(n-3))+(n^2-2*n+3)*fm1*n/((n-1)*(n-2)*(n-3))
  listall<-c(mean=m1,variance=var1,tm=tm1,fm=ufm1)
  (listall)
}

standardizedmoments<-function (x){
  n<-length(x)
  m1<-mean(x)
  var1<-sd(x)^2
  var2<-(sum((x - m1)^2)/n)
  tm1<-(sum((x - m1)^3)/n)*(n^2/((n-1)*(n-2)))
  fm1<-(sum((x - m1)^4)/n)
  ufm1<--3*var2^2*(2*n-3)*n/((n-1)*(n-2)*(n-3))+(n^2-2*n+3)*fm1*n/((n-1)*(n-2)*(n-3))
  listall<-c(mean=m1,variance=var1,skewness=tm1/((var1)^(3/2)),kurtosis=ufm1/((var1)^(2)))
  (listall)
}

biasedmoments_expected<-function (n,targetm,targetvar,targettm,targetfm){
  m1<-targetm
  var1=targetvar/(n/(n-1))
  
  tm1<-targettm/(n^2/((n-1)*(n-2)))
  
  fm1<-(targetfm+3*(targetvar/(n/(n-1)))^2*(2*n-3)*n/((n-1)*(n-2)*(n-3)))/((n^2-2*n+3)*n/((n-1)*(n-2)*(n-3)))
  
  listall<-c(mean=m1,variance=var1,tm=tm1,fm=fm1)
  (listall)
}
weighted_SE<-function(x,weights){
  i <- !is.na(x)
  weights <- weights[i]
  x <- x[i]
  n_eff1 <- (sum(weights))^2/(sum(weights^2))
  wSE1<-sqrt((n_eff1/(n_eff1-1) * (sum(weights*(x-weighted.mean(x,weights))^2)/sum(weights)))/n_eff1)
  return(wSE1)
}
weighted_SD<-function(x,weights){
  i <- !is.na(x)
  weights <- weights[i]
  x <- x[i]
  SD1 =sqrt(sum(weights*(x-weighted.mean(x,weights))^2)/sum(weights))
  return(SD1)
}

random_sequences_process<-function(function1=NULL,expect1=NULL,expect2=NULL,samplesize=NULL,seed1=NULL){

  dataframe_Gaussian<-c()

  #first calculate the expected bias based on unbiased moments
  Results1<-c(0,expect1)
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  dataframe_exp<-c()
  
  Results1<-c(0,expect2)
  
  dataframe_exp<-rbind(dataframe_exp,Results1)
  
  length2<-samplesize
  
  set.seed(seed1)
  ran1<-runif(length2)
  uniran1<-c(ran1)
  
  x<-dsnorm(uni=uniran1,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(1,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  x<-dsexp(uni=uniran1,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(1,function1(x=sortedx))
  
  dataframe_exp<-rbind(dataframe_exp,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran2<-c(ran1)
  
  x<-dsnorm(uni=uniran2,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(2,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  x<-dsexp(uni=uniran2,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(2,function1(x=sortedx))
  
  dataframe_exp<-rbind(dataframe_exp,Results1)
  
  length2<-samplesize
  ran1<-runif(length2)
  uniran3<-c(ran1)
  
  x<-dsnorm(uni=uniran3,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(3,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  x<-dsexp(uni=uniran3,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(3,function1(x=sortedx))
  
  dataframe_exp<-rbind(dataframe_exp,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran4<-c(ran1)
  
  x<-dsnorm(uni=uniran4,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(4,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  x<-dsexp(uni=uniran4,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(4,function1(x=sortedx))
  
  dataframe_exp<-rbind(dataframe_exp,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran5<-c(ran1)
  
  x<-dsnorm(uni=uniran5,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(5,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  x<-dsexp(uni=uniran5,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(5,function1(x=sortedx))
  
  dataframe_exp<-rbind(dataframe_exp,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran6<-c(ran1)
  
  x<-dsnorm(uni=uniran6,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(6,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  x<-dsexp(uni=uniran6,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(6,function1(x=sortedx))
  
  dataframe_exp<-rbind(dataframe_exp,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran7<-c(ran1)
  
  x<-dsnorm(uni=uniran7,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(7,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  x<-dsexp(uni=uniran7,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(7,function1(x=sortedx))
  
  dataframe_exp<-rbind(dataframe_exp,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran8<-c(ran1)
  
  x<-dsnorm(uni=uniran8,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(8,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  x<-dsexp(uni=uniran8,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(8,function1(x=sortedx))
  
  dataframe_exp<-rbind(dataframe_exp,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran9<-c(ran1)
  
  x<-dsnorm(uni=uniran9,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(9,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  x<-dsexp(uni=uniran9,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(9,function1(x=sortedx))
  
  dataframe_exp<-rbind(dataframe_exp,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran10<-c(ran1)
  
  x<-dsnorm(uni=uniran10,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(10,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  x<-dsexp(uni=uniran10,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(10,function1(x=sortedx))
  
  dataframe_exp<-rbind(dataframe_exp,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran11<-c(ran1)
  
  x<-dsnorm(uni=uniran11,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(11,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  x<-dsexp(uni=uniran11,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(11,function1(x=sortedx))
  
  dataframe_exp<-rbind(dataframe_exp,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran12<-c(ran1)
  
  x<-dsnorm(uni=uniran12,location = 0,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(12,function1(x=sortedx))
  
  dataframe_Gaussian<-rbind(dataframe_Gaussian,Results1)
  
  x<-dsexp(uni=uniran12,scale =1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(12,function1(x=sortedx))
  
  dataframe_exp<-rbind(dataframe_exp,Results1)
  
  dataframe1<-rbind(dataframe_Gaussian,dataframe_exp)
  
  return(dataframe1)
}
distribution_sequences_process<-function(function1=NULL,disfunction=NULL,location1=NULL,scale1=NULL,samplesize=NULL,seed1=NULL){
  
  dataframe_distribution1<-c()
  
  length2<-samplesize
  
  set.seed(seed1)
  ran1<-runif(length2)
  uniran1<-c(ran1)
  
  x<-disfunction(uni=uniran1,location = location1,scale =scale1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(1,function1(x=sortedx))
  
  dataframe_distribution1<-rbind(dataframe_distribution1,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran1<-c(ran1)
  
  x<-disfunction(uni=uniran1,location = location1,scale =scale1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(2,function1(x=sortedx))
  
  dataframe_distribution1<-rbind(dataframe_distribution1,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran1<-c(ran1)
  
  x<-disfunction(uni=uniran1,location = location1,scale =scale1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(3,function1(x=sortedx))
  
  dataframe_distribution1<-rbind(dataframe_distribution1,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran1<-c(ran1)
  
  x<-disfunction(uni=uniran1,location = location1,scale =scale1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(4,function1(x=sortedx))
  
  dataframe_distribution1<-rbind(dataframe_distribution1,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran1<-c(ran1)
  
  x<-disfunction(uni=uniran1,location = location1,scale =scale1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(5,function1(x=sortedx))
  
  dataframe_distribution1<-rbind(dataframe_distribution1,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran1<-c(ran1)
  
  x<-disfunction(uni=uniran1,location = location1,scale =scale1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(6,function1(x=sortedx))
  
  dataframe_distribution1<-rbind(dataframe_distribution1,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran1<-c(ran1)
  
  x<-disfunction(uni=uniran1,location = location1,scale =scale1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(7,function1(x=sortedx))
  
  dataframe_distribution1<-rbind(dataframe_distribution1,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran1<-c(ran1)
  
  x<-disfunction(uni=uniran1,location = location1,scale =scale1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(8,function1(x=sortedx))
  
  dataframe_distribution1<-rbind(dataframe_distribution1,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran1<-c(ran1)
  
  x<-disfunction(uni=uniran1,location = location1,scale =scale1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(9,function1(x=sortedx))
  
  dataframe_distribution1<-rbind(dataframe_distribution1,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran1<-c(ran1)
  
  x<-disfunction(uni=uniran1,location = location1,scale =scale1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(10,function1(x=sortedx))
  
  dataframe_distribution1<-rbind(dataframe_distribution1,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran1<-c(ran1)
  
  x<-disfunction(uni=uniran1,location = location1,scale =scale1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(11,function1(x=sortedx))
  
  dataframe_distribution1<-rbind(dataframe_distribution1,Results1)
  
  length2<-samplesize
  
  ran1<-runif(length2)
  uniran1<-c(ran1)
  
  x<-disfunction(uni=uniran1,location = location1,scale =scale1)
  
  sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  Results1<-c(12,function1(x=sortedx))
  
  dataframe_distribution1<-rbind(dataframe_distribution1,Results1)
  
  return(dataframe_distribution1)
}

eval_f <- function(x){
  return ((x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
             x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
           +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
           -dataframe1[1,2])^2+(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                  x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                -dataframe1[1,3])^2+(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                       x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                     +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                     -dataframe1[1,4])^2+(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                            x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                          +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                          -dataframe1[1,5])^2+(x[1]*dataframe_Process[5,1]+x[2]*dataframe_Process[5,2]+x[3]*dataframe_Process[5,3]+x[4]*dataframe_Process[5,4]+
                                                                                                 x[5]*dataframe_Process[5,5]+x[6]*dataframe_Process[5,6]+x[7]*dataframe_Process[5,7]+x[8]*dataframe_Process[5,8]
                                                                                               +x[9]*dataframe_Process[5,9]+x[10]*dataframe_Process[5,10]+x[11]*dataframe_Process[5,11]+x[12]*dataframe_Process[5,12]
                                                                                               -dataframe1[14,2])^2+(x[1]*dataframe_Process[6,1]+x[2]*dataframe_Process[6,2]+x[3]*dataframe_Process[6,3]+x[4]*dataframe_Process[6,4]+
                                                                                                                       x[5]*dataframe_Process[6,5]+x[6]*dataframe_Process[6,6]+x[7]*dataframe_Process[6,7]+x[8]*dataframe_Process[6,8]
                                                                                                                     +x[9]*dataframe_Process[6,9]+x[10]*dataframe_Process[6,10]+x[11]*dataframe_Process[6,11]+x[12]*dataframe_Process[6,12]
                                                                                                                     -dataframe1[14,3])^2+(x[1]*dataframe_Process[7,1]+x[2]*dataframe_Process[7,2]+x[3]*dataframe_Process[7,3]+x[4]*dataframe_Process[7,4]+
                                                                                                                                             x[5]*dataframe_Process[7,5]+x[6]*dataframe_Process[7,6]+x[7]*dataframe_Process[7,7]+x[8]*dataframe_Process[7,8]
                                                                                                                                           +x[9]*dataframe_Process[7,9]+x[10]*dataframe_Process[7,10]+x[11]*dataframe_Process[7,11]+x[12]*dataframe_Process[7,12]
                                                                                                                                           -dataframe1[14,4])^2+(x[1]*dataframe_Process[8,1]+x[2]*dataframe_Process[8,2]+x[3]*dataframe_Process[8,3]+x[4]*dataframe_Process[8,4]+
                                                                                                                                                                   x[5]*dataframe_Process[8,5]+x[6]*dataframe_Process[8,6]+x[7]*dataframe_Process[8,7]+x[8]*dataframe_Process[8,8]
                                                                                                                                                                 +x[9]*dataframe_Process[8,9]+x[10]*dataframe_Process[8,10]+x[11]*dataframe_Process[8,11]+x[12]*dataframe_Process[8,12]
                                                                                                                                                                 -dataframe1[14,5])^2)
}
eval_grad_f <-function(x){
  return(c(2*dataframe_Process[1,1]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,1]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,1]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,1]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5])+2*dataframe_Process[5,1]*(x[1]*dataframe_Process[5,1]+x[2]*dataframe_Process[5,2]+x[3]*dataframe_Process[5,3]+x[4]*dataframe_Process[5,4]+
                                                                                                                                                                                                                       x[5]*dataframe_Process[5,5]+x[6]*dataframe_Process[5,6]+x[7]*dataframe_Process[5,7]+x[8]*dataframe_Process[5,8]
                                                                                                                                                                                                                     +x[9]*dataframe_Process[5,9]+x[10]*dataframe_Process[5,10]+x[11]*dataframe_Process[5,11]+x[12]*dataframe_Process[5,12]
                                                                                                                                                                                                                     -dataframe1[14,2])+2*dataframe_Process[6,1]*(x[1]*dataframe_Process[6,1]+x[2]*dataframe_Process[6,2]+x[3]*dataframe_Process[6,3]+x[4]*dataframe_Process[6,4]+
                                                                                                                                                                                                                                                                    x[5]*dataframe_Process[6,5]+x[6]*dataframe_Process[6,6]+x[7]*dataframe_Process[6,7]+x[8]*dataframe_Process[6,8]
                                                                                                                                                                                                                                                                  +x[9]*dataframe_Process[6,9]+x[10]*dataframe_Process[6,10]+x[11]*dataframe_Process[6,11]+x[12]*dataframe_Process[6,12]
                                                                                                                                                                                                                                                                  -dataframe1[14,3])+2*dataframe_Process[7,1]*(x[1]*dataframe_Process[7,1]+x[2]*dataframe_Process[7,2]+x[3]*dataframe_Process[7,3]+x[4]*dataframe_Process[7,4]+
                                                                                                                                                                                                                                                                                                                 x[5]*dataframe_Process[7,5]+x[6]*dataframe_Process[7,6]+x[7]*dataframe_Process[7,7]+x[8]*dataframe_Process[7,8]
                                                                                                                                                                                                                                                                                                               +x[9]*dataframe_Process[7,9]+x[10]*dataframe_Process[7,10]+x[11]*dataframe_Process[7,11]+x[12]*dataframe_Process[7,12]
                                                                                                                                                                                                                                                                                                               -dataframe1[14,4])+2*dataframe_Process[8,1]*(x[1]*dataframe_Process[8,1]+x[2]*dataframe_Process[8,2]+x[3]*dataframe_Process[8,3]+x[4]*dataframe_Process[8,4]+
                                                                                                                                                                                                                                                                                                                                                              x[5]*dataframe_Process[8,5]+x[6]*dataframe_Process[8,6]+x[7]*dataframe_Process[8,7]+x[8]*dataframe_Process[8,8]
                                                                                                                                                                                                                                                                                                                                                            +x[9]*dataframe_Process[8,9]+x[10]*dataframe_Process[8,10]+x[11]*dataframe_Process[8,11]+x[12]*dataframe_Process[8,12]
                                                                                                                                                                                                                                                                                                                                                            -dataframe1[14,5]),
           2*dataframe_Process[1,2]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,2]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,2]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,2]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5])+2*dataframe_Process[5,2]*(x[1]*dataframe_Process[5,1]+x[2]*dataframe_Process[5,2]+x[3]*dataframe_Process[5,3]+x[4]*dataframe_Process[5,4]+
                                                                                                                                                                                                                       x[5]*dataframe_Process[5,5]+x[6]*dataframe_Process[5,6]+x[7]*dataframe_Process[5,7]+x[8]*dataframe_Process[5,8]
                                                                                                                                                                                                                     +x[9]*dataframe_Process[5,9]+x[10]*dataframe_Process[5,10]+x[11]*dataframe_Process[5,11]+x[12]*dataframe_Process[5,12]
                                                                                                                                                                                                                     -dataframe1[14,2])+2*dataframe_Process[6,2]*(x[1]*dataframe_Process[6,1]+x[2]*dataframe_Process[6,2]+x[3]*dataframe_Process[6,3]+x[4]*dataframe_Process[6,4]+
                                                                                                                                                                                                                                                                    x[5]*dataframe_Process[6,5]+x[6]*dataframe_Process[6,6]+x[7]*dataframe_Process[6,7]+x[8]*dataframe_Process[6,8]
                                                                                                                                                                                                                                                                  +x[9]*dataframe_Process[6,9]+x[10]*dataframe_Process[6,10]+x[11]*dataframe_Process[6,11]+x[12]*dataframe_Process[6,12]
                                                                                                                                                                                                                                                                  -dataframe1[14,3])+2*dataframe_Process[7,2]*(x[1]*dataframe_Process[7,1]+x[2]*dataframe_Process[7,2]+x[3]*dataframe_Process[7,3]+x[4]*dataframe_Process[7,4]+
                                                                                                                                                                                                                                                                                                                 x[5]*dataframe_Process[7,5]+x[6]*dataframe_Process[7,6]+x[7]*dataframe_Process[7,7]+x[8]*dataframe_Process[7,8]
                                                                                                                                                                                                                                                                                                               +x[9]*dataframe_Process[7,9]+x[10]*dataframe_Process[7,10]+x[11]*dataframe_Process[7,11]+x[12]*dataframe_Process[7,12]
                                                                                                                                                                                                                                                                                                               -dataframe1[14,4])+2*dataframe_Process[8,2]*(x[1]*dataframe_Process[8,1]+x[2]*dataframe_Process[8,2]+x[3]*dataframe_Process[8,3]+x[4]*dataframe_Process[8,4]+
                                                                                                                                                                                                                                                                                                                                                              x[5]*dataframe_Process[8,5]+x[6]*dataframe_Process[8,6]+x[7]*dataframe_Process[8,7]+x[8]*dataframe_Process[8,8]
                                                                                                                                                                                                                                                                                                                                                            +x[9]*dataframe_Process[8,9]+x[10]*dataframe_Process[8,10]+x[11]*dataframe_Process[8,11]+x[12]*dataframe_Process[8,12]
                                                                                                                                                                                                                                                                                                                                                            -dataframe1[14,5]),
           2*dataframe_Process[1,3]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,3]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,3]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,3]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5])+2*dataframe_Process[5,3]*(x[1]*dataframe_Process[5,1]+x[2]*dataframe_Process[5,2]+x[3]*dataframe_Process[5,3]+x[4]*dataframe_Process[5,4]+
                                                                                                                                                                                                                       x[5]*dataframe_Process[5,5]+x[6]*dataframe_Process[5,6]+x[7]*dataframe_Process[5,7]+x[8]*dataframe_Process[5,8]
                                                                                                                                                                                                                     +x[9]*dataframe_Process[5,9]+x[10]*dataframe_Process[5,10]+x[11]*dataframe_Process[5,11]+x[12]*dataframe_Process[5,12]
                                                                                                                                                                                                                     -dataframe1[14,2])+2*dataframe_Process[6,3]*(x[1]*dataframe_Process[6,1]+x[2]*dataframe_Process[6,2]+x[3]*dataframe_Process[6,3]+x[4]*dataframe_Process[6,4]+
                                                                                                                                                                                                                                                                    x[5]*dataframe_Process[6,5]+x[6]*dataframe_Process[6,6]+x[7]*dataframe_Process[6,7]+x[8]*dataframe_Process[6,8]
                                                                                                                                                                                                                                                                  +x[9]*dataframe_Process[6,9]+x[10]*dataframe_Process[6,10]+x[11]*dataframe_Process[6,11]+x[12]*dataframe_Process[6,12]
                                                                                                                                                                                                                                                                  -dataframe1[14,3])+2*dataframe_Process[7,3]*(x[1]*dataframe_Process[7,1]+x[2]*dataframe_Process[7,2]+x[3]*dataframe_Process[7,3]+x[4]*dataframe_Process[7,4]+
                                                                                                                                                                                                                                                                                                                 x[5]*dataframe_Process[7,5]+x[6]*dataframe_Process[7,6]+x[7]*dataframe_Process[7,7]+x[8]*dataframe_Process[7,8]
                                                                                                                                                                                                                                                                                                               +x[9]*dataframe_Process[7,9]+x[10]*dataframe_Process[7,10]+x[11]*dataframe_Process[7,11]+x[12]*dataframe_Process[7,12]
                                                                                                                                                                                                                                                                                                               -dataframe1[14,4])+2*dataframe_Process[8,3]*(x[1]*dataframe_Process[8,1]+x[2]*dataframe_Process[8,2]+x[3]*dataframe_Process[8,3]+x[4]*dataframe_Process[8,4]+
                                                                                                                                                                                                                                                                                                                                                              x[5]*dataframe_Process[8,5]+x[6]*dataframe_Process[8,6]+x[7]*dataframe_Process[8,7]+x[8]*dataframe_Process[8,8]
                                                                                                                                                                                                                                                                                                                                                            +x[9]*dataframe_Process[8,9]+x[10]*dataframe_Process[8,10]+x[11]*dataframe_Process[8,11]+x[12]*dataframe_Process[8,12]
                                                                                                                                                                                                                                                                                                                                                            -dataframe1[14,5]),
           2*dataframe_Process[1,4]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,4]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,4]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,4]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5])+2*dataframe_Process[5,4]*(x[1]*dataframe_Process[5,1]+x[2]*dataframe_Process[5,2]+x[3]*dataframe_Process[5,3]+x[4]*dataframe_Process[5,4]+
                                                                                                                                                                                                                       x[5]*dataframe_Process[5,5]+x[6]*dataframe_Process[5,6]+x[7]*dataframe_Process[5,7]+x[8]*dataframe_Process[5,8]
                                                                                                                                                                                                                     +x[9]*dataframe_Process[5,9]+x[10]*dataframe_Process[5,10]+x[11]*dataframe_Process[5,11]+x[12]*dataframe_Process[5,12]
                                                                                                                                                                                                                     -dataframe1[14,2])+2*dataframe_Process[6,4]*(x[1]*dataframe_Process[6,1]+x[2]*dataframe_Process[6,2]+x[3]*dataframe_Process[6,3]+x[4]*dataframe_Process[6,4]+
                                                                                                                                                                                                                                                                    x[5]*dataframe_Process[6,5]+x[6]*dataframe_Process[6,6]+x[7]*dataframe_Process[6,7]+x[8]*dataframe_Process[6,8]
                                                                                                                                                                                                                                                                  +x[9]*dataframe_Process[6,9]+x[10]*dataframe_Process[6,10]+x[11]*dataframe_Process[6,11]+x[12]*dataframe_Process[6,12]
                                                                                                                                                                                                                                                                  -dataframe1[14,3])+2*dataframe_Process[7,4]*(x[1]*dataframe_Process[7,1]+x[2]*dataframe_Process[7,2]+x[3]*dataframe_Process[7,3]+x[4]*dataframe_Process[7,4]+
                                                                                                                                                                                                                                                                                                                 x[5]*dataframe_Process[7,5]+x[6]*dataframe_Process[7,6]+x[7]*dataframe_Process[7,7]+x[8]*dataframe_Process[7,8]
                                                                                                                                                                                                                                                                                                               +x[9]*dataframe_Process[7,9]+x[10]*dataframe_Process[7,10]+x[11]*dataframe_Process[7,11]+x[12]*dataframe_Process[7,12]
                                                                                                                                                                                                                                                                                                               -dataframe1[14,4])+2*dataframe_Process[8,4]*(x[1]*dataframe_Process[8,1]+x[2]*dataframe_Process[8,2]+x[3]*dataframe_Process[8,3]+x[4]*dataframe_Process[8,4]+
                                                                                                                                                                                                                                                                                                                                                              x[5]*dataframe_Process[8,5]+x[6]*dataframe_Process[8,6]+x[7]*dataframe_Process[8,7]+x[8]*dataframe_Process[8,8]
                                                                                                                                                                                                                                                                                                                                                            +x[9]*dataframe_Process[8,9]+x[10]*dataframe_Process[8,10]+x[11]*dataframe_Process[8,11]+x[12]*dataframe_Process[8,12]
                                                                                                                                                                                                                                                                                                                                                            -dataframe1[14,5]),
           2*dataframe_Process[1,5]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,5]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,5]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,5]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5])+2*dataframe_Process[5,5]*(x[1]*dataframe_Process[5,1]+x[2]*dataframe_Process[5,2]+x[3]*dataframe_Process[5,3]+x[4]*dataframe_Process[5,4]+
                                                                                                                                                                                                                       x[5]*dataframe_Process[5,5]+x[6]*dataframe_Process[5,6]+x[7]*dataframe_Process[5,7]+x[8]*dataframe_Process[5,8]
                                                                                                                                                                                                                     +x[9]*dataframe_Process[5,9]+x[10]*dataframe_Process[5,10]+x[11]*dataframe_Process[5,11]+x[12]*dataframe_Process[5,12]
                                                                                                                                                                                                                     -dataframe1[14,2])+2*dataframe_Process[6,5]*(x[1]*dataframe_Process[6,1]+x[2]*dataframe_Process[6,2]+x[3]*dataframe_Process[6,3]+x[4]*dataframe_Process[6,4]+
                                                                                                                                                                                                                                                                    x[5]*dataframe_Process[6,5]+x[6]*dataframe_Process[6,6]+x[7]*dataframe_Process[6,7]+x[8]*dataframe_Process[6,8]
                                                                                                                                                                                                                                                                  +x[9]*dataframe_Process[6,9]+x[10]*dataframe_Process[6,10]+x[11]*dataframe_Process[6,11]+x[12]*dataframe_Process[6,12]
                                                                                                                                                                                                                                                                  -dataframe1[14,3])+2*dataframe_Process[7,5]*(x[1]*dataframe_Process[7,1]+x[2]*dataframe_Process[7,2]+x[3]*dataframe_Process[7,3]+x[4]*dataframe_Process[7,4]+
                                                                                                                                                                                                                                                                                                                 x[5]*dataframe_Process[7,5]+x[6]*dataframe_Process[7,6]+x[7]*dataframe_Process[7,7]+x[8]*dataframe_Process[7,8]
                                                                                                                                                                                                                                                                                                               +x[9]*dataframe_Process[7,9]+x[10]*dataframe_Process[7,10]+x[11]*dataframe_Process[7,11]+x[12]*dataframe_Process[7,12]
                                                                                                                                                                                                                                                                                                               -dataframe1[14,4])+2*dataframe_Process[8,5]*(x[1]*dataframe_Process[8,1]+x[2]*dataframe_Process[8,2]+x[3]*dataframe_Process[8,3]+x[4]*dataframe_Process[8,4]+
                                                                                                                                                                                                                                                                                                                                                              x[5]*dataframe_Process[8,5]+x[6]*dataframe_Process[8,6]+x[7]*dataframe_Process[8,7]+x[8]*dataframe_Process[8,8]
                                                                                                                                                                                                                                                                                                                                                            +x[9]*dataframe_Process[8,9]+x[10]*dataframe_Process[8,10]+x[11]*dataframe_Process[8,11]+x[12]*dataframe_Process[8,12]
                                                                                                                                                                                                                                                                                                                                                            -dataframe1[14,5]),
           2*dataframe_Process[1,6]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,6]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,6]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,6]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5])+2*dataframe_Process[5,6]*(x[1]*dataframe_Process[5,1]+x[2]*dataframe_Process[5,2]+x[3]*dataframe_Process[5,3]+x[4]*dataframe_Process[5,4]+
                                                                                                                                                                                                                       x[5]*dataframe_Process[5,5]+x[6]*dataframe_Process[5,6]+x[7]*dataframe_Process[5,7]+x[8]*dataframe_Process[5,8]
                                                                                                                                                                                                                     +x[9]*dataframe_Process[5,9]+x[10]*dataframe_Process[5,10]+x[11]*dataframe_Process[5,11]+x[12]*dataframe_Process[5,12]
                                                                                                                                                                                                                     -dataframe1[14,2])+2*dataframe_Process[6,6]*(x[1]*dataframe_Process[6,1]+x[2]*dataframe_Process[6,2]+x[3]*dataframe_Process[6,3]+x[4]*dataframe_Process[6,4]+
                                                                                                                                                                                                                                                                    x[5]*dataframe_Process[6,5]+x[6]*dataframe_Process[6,6]+x[7]*dataframe_Process[6,7]+x[8]*dataframe_Process[6,8]
                                                                                                                                                                                                                                                                  +x[9]*dataframe_Process[6,9]+x[10]*dataframe_Process[6,10]+x[11]*dataframe_Process[6,11]+x[12]*dataframe_Process[6,12]
                                                                                                                                                                                                                                                                  -dataframe1[14,3])+2*dataframe_Process[7,6]*(x[1]*dataframe_Process[7,1]+x[2]*dataframe_Process[7,2]+x[3]*dataframe_Process[7,3]+x[4]*dataframe_Process[7,4]+
                                                                                                                                                                                                                                                                                                                 x[5]*dataframe_Process[7,5]+x[6]*dataframe_Process[7,6]+x[7]*dataframe_Process[7,7]+x[8]*dataframe_Process[7,8]
                                                                                                                                                                                                                                                                                                               +x[9]*dataframe_Process[7,9]+x[10]*dataframe_Process[7,10]+x[11]*dataframe_Process[7,11]+x[12]*dataframe_Process[7,12]
                                                                                                                                                                                                                                                                                                               -dataframe1[14,4])+2*dataframe_Process[8,6]*(x[1]*dataframe_Process[8,1]+x[2]*dataframe_Process[8,2]+x[3]*dataframe_Process[8,3]+x[4]*dataframe_Process[8,4]+
                                                                                                                                                                                                                                                                                                                                                              x[5]*dataframe_Process[8,5]+x[6]*dataframe_Process[8,6]+x[7]*dataframe_Process[8,7]+x[8]*dataframe_Process[8,8]
                                                                                                                                                                                                                                                                                                                                                            +x[9]*dataframe_Process[8,9]+x[10]*dataframe_Process[8,10]+x[11]*dataframe_Process[8,11]+x[12]*dataframe_Process[8,12]
                                                                                                                                                                                                                                                                                                                                                            -dataframe1[14,5]),
           2*dataframe_Process[1,7]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,7]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,7]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,7]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5])+2*dataframe_Process[5,7]*(x[1]*dataframe_Process[5,1]+x[2]*dataframe_Process[5,2]+x[3]*dataframe_Process[5,3]+x[4]*dataframe_Process[5,4]+
                                                                                                                                                                                                                       x[5]*dataframe_Process[5,5]+x[6]*dataframe_Process[5,6]+x[7]*dataframe_Process[5,7]+x[8]*dataframe_Process[5,8]
                                                                                                                                                                                                                     +x[9]*dataframe_Process[5,9]+x[10]*dataframe_Process[5,10]+x[11]*dataframe_Process[5,11]+x[12]*dataframe_Process[5,12]
                                                                                                                                                                                                                     -dataframe1[14,2])+2*dataframe_Process[6,7]*(x[1]*dataframe_Process[6,1]+x[2]*dataframe_Process[6,2]+x[3]*dataframe_Process[6,3]+x[4]*dataframe_Process[6,4]+
                                                                                                                                                                                                                                                                    x[5]*dataframe_Process[6,5]+x[6]*dataframe_Process[6,6]+x[7]*dataframe_Process[6,7]+x[8]*dataframe_Process[6,8]
                                                                                                                                                                                                                                                                  +x[9]*dataframe_Process[6,9]+x[10]*dataframe_Process[6,10]+x[11]*dataframe_Process[6,11]+x[12]*dataframe_Process[6,12]
                                                                                                                                                                                                                                                                  -dataframe1[14,3])+2*dataframe_Process[7,7]*(x[1]*dataframe_Process[7,1]+x[2]*dataframe_Process[7,2]+x[3]*dataframe_Process[7,3]+x[4]*dataframe_Process[7,4]+
                                                                                                                                                                                                                                                                                                                 x[5]*dataframe_Process[7,5]+x[6]*dataframe_Process[7,6]+x[7]*dataframe_Process[7,7]+x[8]*dataframe_Process[7,8]
                                                                                                                                                                                                                                                                                                               +x[9]*dataframe_Process[7,9]+x[10]*dataframe_Process[7,10]+x[11]*dataframe_Process[7,11]+x[12]*dataframe_Process[7,12]
                                                                                                                                                                                                                                                                                                               -dataframe1[14,4])+2*dataframe_Process[8,7]*(x[1]*dataframe_Process[8,1]+x[2]*dataframe_Process[8,2]+x[3]*dataframe_Process[8,3]+x[4]*dataframe_Process[8,4]+
                                                                                                                                                                                                                                                                                                                                                              x[5]*dataframe_Process[8,5]+x[6]*dataframe_Process[8,6]+x[7]*dataframe_Process[8,7]+x[8]*dataframe_Process[8,8]
                                                                                                                                                                                                                                                                                                                                                            +x[9]*dataframe_Process[8,9]+x[10]*dataframe_Process[8,10]+x[11]*dataframe_Process[8,11]+x[12]*dataframe_Process[8,12]
                                                                                                                                                                                                                                                                                                                                                            -dataframe1[14,5]),
           2*dataframe_Process[1,8]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,8]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,8]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,8]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5])+2*dataframe_Process[5,8]*(x[1]*dataframe_Process[5,1]+x[2]*dataframe_Process[5,2]+x[3]*dataframe_Process[5,3]+x[4]*dataframe_Process[5,4]+
                                                                                                                                                                                                                       x[5]*dataframe_Process[5,5]+x[6]*dataframe_Process[5,6]+x[7]*dataframe_Process[5,7]+x[8]*dataframe_Process[5,8]
                                                                                                                                                                                                                     +x[9]*dataframe_Process[5,9]+x[10]*dataframe_Process[5,10]+x[11]*dataframe_Process[5,11]+x[12]*dataframe_Process[5,12]
                                                                                                                                                                                                                     -dataframe1[14,2])+2*dataframe_Process[6,8]*(x[1]*dataframe_Process[6,1]+x[2]*dataframe_Process[6,2]+x[3]*dataframe_Process[6,3]+x[4]*dataframe_Process[6,4]+
                                                                                                                                                                                                                                                                    x[5]*dataframe_Process[6,5]+x[6]*dataframe_Process[6,6]+x[7]*dataframe_Process[6,7]+x[8]*dataframe_Process[6,8]
                                                                                                                                                                                                                                                                  +x[9]*dataframe_Process[6,9]+x[10]*dataframe_Process[6,10]+x[11]*dataframe_Process[6,11]+x[12]*dataframe_Process[6,12]
                                                                                                                                                                                                                                                                  -dataframe1[14,3])+2*dataframe_Process[7,8]*(x[1]*dataframe_Process[7,1]+x[2]*dataframe_Process[7,2]+x[3]*dataframe_Process[7,3]+x[4]*dataframe_Process[7,4]+
                                                                                                                                                                                                                                                                                                                 x[5]*dataframe_Process[7,5]+x[6]*dataframe_Process[7,6]+x[7]*dataframe_Process[7,7]+x[8]*dataframe_Process[7,8]
                                                                                                                                                                                                                                                                                                               +x[9]*dataframe_Process[7,9]+x[10]*dataframe_Process[7,10]+x[11]*dataframe_Process[7,11]+x[12]*dataframe_Process[7,12]
                                                                                                                                                                                                                                                                                                               -dataframe1[14,4])+2*dataframe_Process[8,8]*(x[1]*dataframe_Process[8,1]+x[2]*dataframe_Process[8,2]+x[3]*dataframe_Process[8,3]+x[4]*dataframe_Process[8,4]+
                                                                                                                                                                                                                                                                                                                                                              x[5]*dataframe_Process[8,5]+x[6]*dataframe_Process[8,6]+x[7]*dataframe_Process[8,7]+x[8]*dataframe_Process[8,8]
                                                                                                                                                                                                                                                                                                                                                            +x[9]*dataframe_Process[8,9]+x[10]*dataframe_Process[8,10]+x[11]*dataframe_Process[8,11]+x[12]*dataframe_Process[8,12]
                                                                                                                                                                                                                                                                                                                                                            -dataframe1[14,5]),
           2*dataframe_Process[1,9]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                       x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                     +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                     -dataframe1[1,2])+2*dataframe_Process[2,9]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                   x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                 +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                 -dataframe1[1,3])+2*dataframe_Process[3,9]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                               x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                             +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                             -dataframe1[1,4])+2*dataframe_Process[4,9]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                           x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                         +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                         -dataframe1[1,5])+2*dataframe_Process[5,9]*(x[1]*dataframe_Process[5,1]+x[2]*dataframe_Process[5,2]+x[3]*dataframe_Process[5,3]+x[4]*dataframe_Process[5,4]+
                                                                                                                                                                                                                       x[5]*dataframe_Process[5,5]+x[6]*dataframe_Process[5,6]+x[7]*dataframe_Process[5,7]+x[8]*dataframe_Process[5,8]
                                                                                                                                                                                                                     +x[9]*dataframe_Process[5,9]+x[10]*dataframe_Process[5,10]+x[11]*dataframe_Process[5,11]+x[12]*dataframe_Process[5,12]
                                                                                                                                                                                                                     -dataframe1[14,2])+2*dataframe_Process[6,9]*(x[1]*dataframe_Process[6,1]+x[2]*dataframe_Process[6,2]+x[3]*dataframe_Process[6,3]+x[4]*dataframe_Process[6,4]+
                                                                                                                                                                                                                                                                    x[5]*dataframe_Process[6,5]+x[6]*dataframe_Process[6,6]+x[7]*dataframe_Process[6,7]+x[8]*dataframe_Process[6,8]
                                                                                                                                                                                                                                                                  +x[9]*dataframe_Process[6,9]+x[10]*dataframe_Process[6,10]+x[11]*dataframe_Process[6,11]+x[12]*dataframe_Process[6,12]
                                                                                                                                                                                                                                                                  -dataframe1[14,3])+2*dataframe_Process[7,9]*(x[1]*dataframe_Process[7,1]+x[2]*dataframe_Process[7,2]+x[3]*dataframe_Process[7,3]+x[4]*dataframe_Process[7,4]+
                                                                                                                                                                                                                                                                                                                 x[5]*dataframe_Process[7,5]+x[6]*dataframe_Process[7,6]+x[7]*dataframe_Process[7,7]+x[8]*dataframe_Process[7,8]
                                                                                                                                                                                                                                                                                                               +x[9]*dataframe_Process[7,9]+x[10]*dataframe_Process[7,10]+x[11]*dataframe_Process[7,11]+x[12]*dataframe_Process[7,12]
                                                                                                                                                                                                                                                                                                               -dataframe1[14,4])+2*dataframe_Process[8,9]*(x[1]*dataframe_Process[8,1]+x[2]*dataframe_Process[8,2]+x[3]*dataframe_Process[8,3]+x[4]*dataframe_Process[8,4]+
                                                                                                                                                                                                                                                                                                                                                              x[5]*dataframe_Process[8,5]+x[6]*dataframe_Process[8,6]+x[7]*dataframe_Process[8,7]+x[8]*dataframe_Process[8,8]
                                                                                                                                                                                                                                                                                                                                                            +x[9]*dataframe_Process[8,9]+x[10]*dataframe_Process[8,10]+x[11]*dataframe_Process[8,11]+x[12]*dataframe_Process[8,12]
                                                                                                                                                                                                                                                                                                                                                            -dataframe1[14,5]),
           2*dataframe_Process[1,10]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                        x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                      +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                      -dataframe1[1,2])+2*dataframe_Process[2,10]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                     x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                   +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                   -dataframe1[1,3])+2*dataframe_Process[3,10]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                                  x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                                +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                                -dataframe1[1,4])+2*dataframe_Process[4,10]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                               x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                             +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                             -dataframe1[1,5])+2*dataframe_Process[5,10]*(x[1]*dataframe_Process[5,1]+x[2]*dataframe_Process[5,2]+x[3]*dataframe_Process[5,3]+x[4]*dataframe_Process[5,4]+
                                                                                                                                                                                                                            x[5]*dataframe_Process[5,5]+x[6]*dataframe_Process[5,6]+x[7]*dataframe_Process[5,7]+x[8]*dataframe_Process[5,8]
                                                                                                                                                                                                                          +x[9]*dataframe_Process[5,9]+x[10]*dataframe_Process[5,10]+x[11]*dataframe_Process[5,11]+x[12]*dataframe_Process[5,12]
                                                                                                                                                                                                                          -dataframe1[14,2])+2*dataframe_Process[6,10]*(x[1]*dataframe_Process[6,1]+x[2]*dataframe_Process[6,2]+x[3]*dataframe_Process[6,3]+x[4]*dataframe_Process[6,4]+
                                                                                                                                                                                                                                                                          x[5]*dataframe_Process[6,5]+x[6]*dataframe_Process[6,6]+x[7]*dataframe_Process[6,7]+x[8]*dataframe_Process[6,8]
                                                                                                                                                                                                                                                                        +x[9]*dataframe_Process[6,9]+x[10]*dataframe_Process[6,10]+x[11]*dataframe_Process[6,11]+x[12]*dataframe_Process[6,12]
                                                                                                                                                                                                                                                                        -dataframe1[14,3])+2*dataframe_Process[7,10]*(x[1]*dataframe_Process[7,1]+x[2]*dataframe_Process[7,2]+x[3]*dataframe_Process[7,3]+x[4]*dataframe_Process[7,4]+
                                                                                                                                                                                                                                                                                                                        x[5]*dataframe_Process[7,5]+x[6]*dataframe_Process[7,6]+x[7]*dataframe_Process[7,7]+x[8]*dataframe_Process[7,8]
                                                                                                                                                                                                                                                                                                                      +x[9]*dataframe_Process[7,9]+x[10]*dataframe_Process[7,10]+x[11]*dataframe_Process[7,11]+x[12]*dataframe_Process[7,12]
                                                                                                                                                                                                                                                                                                                      -dataframe1[14,4])+2*dataframe_Process[8,10]*(x[1]*dataframe_Process[8,1]+x[2]*dataframe_Process[8,2]+x[3]*dataframe_Process[8,3]+x[4]*dataframe_Process[8,4]+
                                                                                                                                                                                                                                                                                                                                                                      x[5]*dataframe_Process[8,5]+x[6]*dataframe_Process[8,6]+x[7]*dataframe_Process[8,7]+x[8]*dataframe_Process[8,8]
                                                                                                                                                                                                                                                                                                                                                                    +x[9]*dataframe_Process[8,9]+x[10]*dataframe_Process[8,10]+x[11]*dataframe_Process[8,11]+x[12]*dataframe_Process[8,12]
                                                                                                                                                                                                                                                                                                                                                                    -dataframe1[14,5]),
           2*dataframe_Process[1,11]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                        x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                      +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                      -dataframe1[1,2])+2*dataframe_Process[2,11]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                     x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                   +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                   -dataframe1[1,3])+2*dataframe_Process[3,11]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                                  x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                                +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                                -dataframe1[1,4])+2*dataframe_Process[4,11]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                               x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                             +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                             -dataframe1[1,5])+2*dataframe_Process[5,11]*(x[1]*dataframe_Process[5,1]+x[2]*dataframe_Process[5,2]+x[3]*dataframe_Process[5,3]+x[4]*dataframe_Process[5,4]+
                                                                                                                                                                                                                            x[5]*dataframe_Process[5,5]+x[6]*dataframe_Process[5,6]+x[7]*dataframe_Process[5,7]+x[8]*dataframe_Process[5,8]
                                                                                                                                                                                                                          +x[9]*dataframe_Process[5,9]+x[10]*dataframe_Process[5,10]+x[11]*dataframe_Process[5,11]+x[12]*dataframe_Process[5,12]
                                                                                                                                                                                                                          -dataframe1[14,2])+2*dataframe_Process[6,11]*(x[1]*dataframe_Process[6,1]+x[2]*dataframe_Process[6,2]+x[3]*dataframe_Process[6,3]+x[4]*dataframe_Process[6,4]+
                                                                                                                                                                                                                                                                          x[5]*dataframe_Process[6,5]+x[6]*dataframe_Process[6,6]+x[7]*dataframe_Process[6,7]+x[8]*dataframe_Process[6,8]
                                                                                                                                                                                                                                                                        +x[9]*dataframe_Process[6,9]+x[10]*dataframe_Process[6,10]+x[11]*dataframe_Process[6,11]+x[12]*dataframe_Process[6,12]
                                                                                                                                                                                                                                                                        -dataframe1[14,3])+2*dataframe_Process[7,11]*(x[1]*dataframe_Process[7,1]+x[2]*dataframe_Process[7,2]+x[3]*dataframe_Process[7,3]+x[4]*dataframe_Process[7,4]+
                                                                                                                                                                                                                                                                                                                        x[5]*dataframe_Process[7,5]+x[6]*dataframe_Process[7,6]+x[7]*dataframe_Process[7,7]+x[8]*dataframe_Process[7,8]
                                                                                                                                                                                                                                                                                                                      +x[9]*dataframe_Process[7,9]+x[10]*dataframe_Process[7,10]+x[11]*dataframe_Process[7,11]+x[12]*dataframe_Process[7,12]
                                                                                                                                                                                                                                                                                                                      -dataframe1[14,4])+2*dataframe_Process[8,11]*(x[1]*dataframe_Process[8,1]+x[2]*dataframe_Process[8,2]+x[3]*dataframe_Process[8,3]+x[4]*dataframe_Process[8,4]+
                                                                                                                                                                                                                                                                                                                                                                      x[5]*dataframe_Process[8,5]+x[6]*dataframe_Process[8,6]+x[7]*dataframe_Process[8,7]+x[8]*dataframe_Process[8,8]
                                                                                                                                                                                                                                                                                                                                                                    +x[9]*dataframe_Process[8,9]+x[10]*dataframe_Process[8,10]+x[11]*dataframe_Process[8,11]+x[12]*dataframe_Process[8,12]
                                                                                                                                                                                                                                                                                                                                                                    -dataframe1[14,5]),
           2*dataframe_Process[1,12]*(x[1]*dataframe_Process[1,1]+x[2]*dataframe_Process[1,2]+x[3]*dataframe_Process[1,3]+x[4]*dataframe_Process[1,4]+
                                        x[5]*dataframe_Process[1,5]+x[6]*dataframe_Process[1,6]+x[7]*dataframe_Process[1,7]+x[8]*dataframe_Process[1,8]
                                      +x[9]*dataframe_Process[1,9]+x[10]*dataframe_Process[1,10]+x[11]*dataframe_Process[1,11]+x[12]*dataframe_Process[1,12]
                                      -dataframe1[1,2])+2*dataframe_Process[2,12]*(x[1]*dataframe_Process[2,1]+x[2]*dataframe_Process[2,2]+x[3]*dataframe_Process[2,3]+x[4]*dataframe_Process[2,4]+
                                                                                     x[5]*dataframe_Process[2,5]+x[6]*dataframe_Process[2,6]+x[7]*dataframe_Process[2,7]+x[8]*dataframe_Process[2,8]
                                                                                   +x[9]*dataframe_Process[2,9]+x[10]*dataframe_Process[2,10]+x[11]*dataframe_Process[2,11]+x[12]*dataframe_Process[2,12]
                                                                                   -dataframe1[1,3])+2*dataframe_Process[3,12]*(x[1]*dataframe_Process[3,1]+x[2]*dataframe_Process[3,2]+x[3]*dataframe_Process[3,3]+x[4]*dataframe_Process[3,4]+
                                                                                                                                  x[5]*dataframe_Process[3,5]+x[6]*dataframe_Process[3,6]+x[7]*dataframe_Process[3,7]+x[8]*dataframe_Process[3,8]
                                                                                                                                +x[9]*dataframe_Process[3,9]+x[10]*dataframe_Process[3,10]+x[11]*dataframe_Process[3,11]+x[12]*dataframe_Process[3,12]
                                                                                                                                -dataframe1[1,4])+2*dataframe_Process[4,12]*(x[1]*dataframe_Process[4,1]+x[2]*dataframe_Process[4,2]+x[3]*dataframe_Process[4,3]+x[4]*dataframe_Process[4,4]+
                                                                                                                                                                               x[5]*dataframe_Process[4,5]+x[6]*dataframe_Process[4,6]+x[7]*dataframe_Process[4,7]+x[8]*dataframe_Process[4,8]
                                                                                                                                                                             +x[9]*dataframe_Process[4,9]+x[10]*dataframe_Process[4,10]+x[11]*dataframe_Process[4,11]+x[12]*dataframe_Process[4,12]
                                                                                                                                                                             -dataframe1[1,5])+2*dataframe_Process[5,12]*(x[1]*dataframe_Process[5,1]+x[2]*dataframe_Process[5,2]+x[3]*dataframe_Process[5,3]+x[4]*dataframe_Process[5,4]+
                                                                                                                                                                                                                            x[5]*dataframe_Process[5,5]+x[6]*dataframe_Process[5,6]+x[7]*dataframe_Process[5,7]+x[8]*dataframe_Process[5,8]
                                                                                                                                                                                                                          +x[9]*dataframe_Process[5,9]+x[10]*dataframe_Process[5,10]+x[11]*dataframe_Process[5,11]+x[12]*dataframe_Process[5,12]
                                                                                                                                                                                                                          -dataframe1[14,2])+2*dataframe_Process[6,12]*(x[1]*dataframe_Process[6,1]+x[2]*dataframe_Process[6,2]+x[3]*dataframe_Process[6,3]+x[4]*dataframe_Process[6,4]+
                                                                                                                                                                                                                                                                          x[5]*dataframe_Process[6,5]+x[6]*dataframe_Process[6,6]+x[7]*dataframe_Process[6,7]+x[8]*dataframe_Process[6,8]
                                                                                                                                                                                                                                                                        +x[9]*dataframe_Process[6,9]+x[10]*dataframe_Process[6,10]+x[11]*dataframe_Process[6,11]+x[12]*dataframe_Process[6,12]
                                                                                                                                                                                                                                                                        -dataframe1[14,3])+2*dataframe_Process[7,12]*(x[1]*dataframe_Process[7,1]+x[2]*dataframe_Process[7,2]+x[3]*dataframe_Process[7,3]+x[4]*dataframe_Process[7,4]+
                                                                                                                                                                                                                                                                                                                        x[5]*dataframe_Process[7,5]+x[6]*dataframe_Process[7,6]+x[7]*dataframe_Process[7,7]+x[8]*dataframe_Process[7,8]
                                                                                                                                                                                                                                                                                                                      +x[9]*dataframe_Process[7,9]+x[10]*dataframe_Process[7,10]+x[11]*dataframe_Process[7,11]+x[12]*dataframe_Process[7,12]
                                                                                                                                                                                                                                                                                                                      -dataframe1[14,4])+2*dataframe_Process[8,12]*(x[1]*dataframe_Process[8,1]+x[2]*dataframe_Process[8,2]+x[3]*dataframe_Process[8,3]+x[4]*dataframe_Process[8,4]+
                                                                                                                                                                                                                                                                                                                                                                      x[5]*dataframe_Process[8,5]+x[6]*dataframe_Process[8,6]+x[7]*dataframe_Process[8,7]+x[8]*dataframe_Process[8,8]
                                                                                                                                                                                                                                                                                                                                                                    +x[9]*dataframe_Process[8,9]+x[10]*dataframe_Process[8,10]+x[11]*dataframe_Process[8,11]+x[12]*dataframe_Process[8,12]
                                                                                                                                                                                                                                                                                                                                                                    -dataframe1[14,5])) )
}
eval_g_eq<-function(x){
  constr<-c((1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10]))^2)
  grad<-c(-2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10])),
          -2*(1-(x[12]+x[11]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9]+x[10]))
  )
  return(list("constraints"=constr,"jacobian"=grad))
}
#Pool adjacent violators (PAV) Algorithm
MonotonicSmoothing<-function(x,y,reverse1=FALSE,auto=FALSE){
  n<-length(y)
  if (auto){
    if((max(y)-y[1])>=(y[1]-min(y))){
      reverse1<-FALSE
    }else{
      reverse1<-TRUE
    }
  }
  
  if(reverse1){
    y <- rev(y)
  }
  w <- rep(1,n)
  rep1 <- rep(1,n)
  repeat {
    Indicator <- TRUE
    i <- 1
    while(i < n) {
      if(y[i] > y[i+1]) {
        Indicator <- FALSE
        w_w1 <- w[i] + w[i+1]
        sum1 <- (w[i]*y[i] + w[i+1]*y[i+1])/w_w1
        y[i+1] <- sum1
        w[i+1] <- w_w1
        y <- y[-i]
        w <- w[-i]
        rep1[i+1] <- rep1[i] + rep1[i+1]
        rep1 <- rep1[-i]
        n <- n-1
      }
      i <- i+1
    }
    if(Indicator) {
      break
    }
  }
  y  <- rep(y,rep1)
  if(reverse1){
    y <- rev(y)
  }
  return(data.frame(x=x,y=y))
}

simulatedbatch_FSD_calibration<-foreach(batchnumber =c(c(seq(from=16, to=400, by=8),8*round(2^seq(from=6, to=10, by=0.18)))), .combine = 'rbind') %dopar% {
  library(Rfast)
  if (!require("foreach")) install.packages("foreach")
  library(foreach)
  if (!require("doParallel")) install.packages("doParallel")
  library(doParallel)
  #registering clusters, can set a smaller number using numCores-1 
  
  #require randtoolbox for random number generations
  if (!require("randtoolbox")) install.packages("randtoolbox")
  library(randtoolbox)
  #require Rfast for faster computation
  if (!require("Rfast")) install.packages("Rfast")
  library(Rfast)
  if (!require("gnorm")) install.packages("gnorm")
  library(gnorm)
  
  n <- batchnumber
  
  samplesize=batchnumber
  
  batchsize<-ceiling(100/sqrt(samplesize))
  
  seed1=1
  
  criteria<-1e-10
  
  uniran<-c()
  repeat{
    dataframe1<-random_sequences_process(function1 = moments,expect1=biasedmoments_expected(n=samplesize,targetm=0,targetvar=1,targettm=0,targetfm=3),
                                         expect2=biasedmoments_expected(n=samplesize,targetm=1,targetvar=1,targettm=2,targetfm=9),
                                         samplesize=samplesize,seed1=seed1)
    
    dataframe_Process<-t(cbind(dataframe1[2:13,2:5],dataframe1[15:26,2:5]))
    
    x0<-rep(1,12)
    lb<-rep(0,12)
    ub<-rep(1,12)
    opts <- list("algorithm"="NLOPT_LD_SLSQP",
                 "xtol_rel"=1.0e-10,"maxeval" = 1000)
    
    library(nloptr)
    results1 <- nloptr( x0=x0,
                        eval_f=eval_f,
                        lb=lb,
                        ub=ub,
                        eval_grad_f=eval_grad_f,
                        eval_g_eq=eval_g_eq,
                        opts=opts)
    
    objective=results1$objective
    
    weights=results1$solution
    
    if(objective<criteria){
      uniran<-rbind(uniran,c(samplesize=samplesize,seed1=seed1,objective=results1$objective,weights=results1$solution))
    }
    if (!is.null(nrow(uniran))){
    if (nrow(uniran)==batchsize){
      break
    }}
    seed1=seed1+1
  }
  uniran
}
simulatedbatch_FSD_calibration<-as.data.frame(simulatedbatch_FSD_calibration)
saveRDS(simulatedbatch_FSD_calibration, paste("Calibrated_Monte_Carlo_Gaussian_exp.Rds", sep = ",")) 
# 
# #filter the error larger than 1e^-4
# simulatedbatch_FSD_calibration_1e4<-simulatedbatch_FSD_calibration[simulatedbatch_FSD_calibration$objective<1e-4,]
# 
# sample_group <- function(x, n) {
#   if (length(x) <= n){
#     return(x)}else{
#   x[x %in% sample(x, length(x)/5)]}
# }
# 
# simulatedbatch_FSD_calibration_1e4<-simulatedbatch_FSD_calibration_1e4[unlist(lapply(split(1:nrow(simulatedbatch_FSD_calibration_1e4), simulatedbatch_FSD_calibration_1e4$samplesize), sample_group, n = 8)), ]


#test with the standard deviation of Gaussian distribution
simulatedbatch_Gaussian_calibrated_Monte_sd<-foreach(batchnumber =c(c(seq(from=16, to=400, by=8),8*round(2^seq(from=6, to=10, by=0.18)))), .combine = 'rbind') %dopar% {
  library(Rfast)
  if (!require("foreach")) install.packages("foreach")
  library(foreach)
  if (!require("doParallel")) install.packages("doParallel")
  library(doParallel)
  #registering clusters, can set a smaller number using numCores-1 
  
  #require randtoolbox for random number generations
  if (!require("randtoolbox")) install.packages("randtoolbox")
  library(randtoolbox)
  #require Rfast for faster computation
  if (!require("Rfast")) install.packages("Rfast")
  library(Rfast)
  if (!require("gnorm")) install.packages("gnorm")
  library(gnorm)
  
  n <- batchnumber
  
  samplesize=batchnumber
  
  weights1<-simulatedbatch_FSD_calibration[simulatedbatch_FSD_calibration$samplesize==samplesize,]
  
  SEbataches1<-c()
  dataframe2<-c()
  for(batch1 in (1:nrow(weights1))){
    seed2<-as.numeric(weights1[batch1,2])
    weight2<-as.numeric(weights1[batch1,4:15])
    dataframe2<-distribution_sequences_process(function1=sd,disfunction=dsnorm,location1=0,scale1=1,samplesize=samplesize,seed1=seed2)
    weightedmean<-c(n,weighted.mean(dataframe2[1:12,2],weight2))
    SEbataches1<-rbind(SEbataches1,weightedmean)
  }
  
  SEbataches2<-c()
  for (batch1 in 1:nrow(weights1)){
    uniran1<-runif(n)
    x<-c(dsnorm(uni=uniran1,location=0,scale=1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    x<-c()
    sdx<-c(n,sd(x=sortedx),correctfactor(n))
    sortedx<-c()
    SEbataches2<-rbind(SEbataches2,sdx)
  }
  cbind(SEbataches1,SEbataches2)
}

#A standard dataframe
Finite_Gaussian_sd <- data.frame(Size = simulatedbatch_Gaussian_calibrated_Monte_sd[,1],
                                 Calibrated_Monte=simulatedbatch_Gaussian_calibrated_Monte_sd[,2],
                                 
                                 Monte=simulatedbatch_Gaussian_calibrated_Monte_sd[,4],
                                 Exact=simulatedbatch_Gaussian_calibrated_Monte_sd[,5]
)

rownames(Finite_Gaussian_sd)<-c()

Finite_Gaussian_sd_Calibrated_Monte<-tapply(Finite_Gaussian_sd$Calibrated_Monte, Finite_Gaussian_sd$Size, mean)

Finite_Gaussian_sd_Monte<-tapply(Finite_Gaussian_sd$Monte, Finite_Gaussian_sd$Size, mean)

Finite_Gaussian_sd <- Finite_Gaussian_sd[!duplicated(Finite_Gaussian_sd$Size), ]

Finite_Gaussian_sd$Calibrated_Monte<-Finite_Gaussian_sd_Calibrated_Monte

Finite_Gaussian_sd$Monte<-Finite_Gaussian_sd_Monte  

Calibrated_Monte_MonotonSmooth1<-MonotonicSmoothing(x=Finite_Gaussian_sd$Size,y=Finite_Gaussian_sd$Calibrated_Monte,auto=TRUE)

Finite_Gaussian_sd$Calibrated_Monte_mono<-Calibrated_Monte_MonotonSmooth1$y

Monte_MonotonSmooth1<-MonotonicSmoothing(x=Finite_Gaussian_sd$Size,y=Finite_Gaussian_sd$Monte,auto=TRUE)

Finite_Gaussian_sd$Monte_mono<-Monte_MonotonSmooth1$y

Allresults_Monte101 <- Finite_Gaussian_sd[Finite_Gaussian_sd$Size %in% c(c(seq(from=16, to=400, by=8))),]#,8*round(2^seq(from=6, to=9, by=0.18))),]

theme_set(theme_bw(base_size=16))
gfg_plot <- ggplot(Finite_Gaussian_sd, aes(Size)) +
  geom_line(aes(y = Monte), color = "red", size=1) +
  geom_line(aes(y = Calibrated_Monte), color = "green", size=1)+geom_line(aes(y = Monte_mono), color = "blue", size=1)+
  geom_line(aes(y = Calibrated_Monte_mono), color = "black", size=1)+geom_line(aes(y = Exact), color = "purple", size=1)


gfg_plot
gfg_plot <- ggplot(Allresults_Monte101, aes(Size)) +
  geom_line(aes(y = Monte), color = "red", size=1) +
  geom_line(aes(y = Calibrated_Monte), color = "green", size=1)+geom_line(aes(y = Monte_mono), color = "blue", size=1)+
  geom_line(aes(y = Calibrated_Monte_mono), color = "black", size=1)+geom_line(aes(y = Exact), color = "purple", size=1)

gfg_plot

Finite_Gaussian_sd$Diff_Calibrated_Monte<-abs(Finite_Gaussian_sd$Exact-Finite_Gaussian_sd$Calibrated_Monte)
Finite_Gaussian_sd$Diff_Monte<-abs(Finite_Gaussian_sd$Exact-Finite_Gaussian_sd$Monte)

Finite_Gaussian_sd$Diff_Calibrated_Monte_mono<-abs(Finite_Gaussian_sd$Exact-Finite_Gaussian_sd$Calibrated_Monte_mono)
Finite_Gaussian_sd$Diff_Monte_mono<-abs(Finite_Gaussian_sd$Exact-Finite_Gaussian_sd$Monte_mono)

Finite_Gaussian_sd$Improvement1<-Finite_Gaussian_sd$Diff_Calibrated_Monte_mono/Finite_Gaussian_sd$Diff_Monte

Finite_Gaussian_sd$Improvement2<-Finite_Gaussian_sd$Diff_Calibrated_Monte_mono/Finite_Gaussian_sd$Diff_Calibrated_Monte

mean(Finite_Gaussian_sd$Improvement1)

median(Finite_Gaussian_sd$Improvement2)

colMeans(Finite_Gaussian_sd)
#test with the median of exponential distribution
simulatedbatch_exp_calibrated_Monte_median<-foreach(batchnumber =c(c(seq(from=16, to=400, by=8),8*round(2^seq(from=6, to=10, by=0.18)))), .combine = 'rbind') %dopar% {
  library(Rfast)
  if (!require("foreach")) install.packages("foreach")
  library(foreach)
  if (!require("doParallel")) install.packages("doParallel")
  library(doParallel)
  #registering clusters, can set a smaller number using numCores-1 
  
  #require randtoolbox for random number generations
  if (!require("randtoolbox")) install.packages("randtoolbox")
  library(randtoolbox)
  #require Rfast for faster computation
  if (!require("Rfast")) install.packages("Rfast")
  library(Rfast)
  if (!require("gnorm")) install.packages("gnorm")
  library(gnorm)
  
  n <- batchnumber
  
  samplesize=batchnumber
  
  weights1<-simulatedbatch_FSD_calibration[simulatedbatch_FSD_calibration$samplesize==samplesize,]
  
  SEbataches1<-c()
  dataframe2<-c()
  for(batch1 in (1:nrow(weights1))){
    seed2<-as.numeric(weights1[batch1,2])
    weight2<-as.numeric(weights1[batch1,4:15])
    dataframe2<-distribution_sequences_process(function1=median,disfunction=dsexp,location1=0,scale1=1,samplesize=samplesize,seed1=seed2)
    weightedmean<-c(n,weighted.mean(dataframe2[1:12,2],weight2))
    SEbataches1<-rbind(SEbataches1,weightedmean)
  }
  
  SEbataches2<-c()
  for (batch1 in 1:nrow(weights1)){
    uniran1<-runif(n)
    x<-c(dsexp(uni=uniran1, scale=1))
    sortedx<-Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    x<-c()
    medianx<-c(n,median(x=sortedx))
    sortedx<-c()
    SEbataches2<-rbind(SEbataches2,medianx)
  }
  cbind(SEbataches1,SEbataches2)
}

medianexp<- read.csv(paste("medianexp.csv", sep = ","))
medianexp<-data.frame(Size = (1:10000),
                      Exact =medianexp)
medianexp <- medianexp[medianexp$Size %in% c(c(seq(from=16, to=400, by=8)),8*round(2^seq(from=6, to=10, by=0.18))),]


#A standard dataframe
Finite_exp_median <- data.frame(Size = simulatedbatch_exp_calibrated_Monte_median[,1],
                                 Calibrated_Monte=simulatedbatch_exp_calibrated_Monte_median[,2],
                                 
                                 Monte=simulatedbatch_exp_calibrated_Monte_median[,4]
)

rownames(Finite_exp_median)<-c()

Finite_exp_median_Calibrated_Monte<-tapply(Finite_exp_median$Calibrated_Monte, Finite_exp_median$Size, mean)

Finite_exp_median_Monte<-tapply(Finite_exp_median$Monte, Finite_exp_median$Size, mean)

Finite_exp_median <- Finite_exp_median[!duplicated(Finite_exp_median$Size), ]

Finite_exp_median$Calibrated_Monte<-Finite_exp_median_Calibrated_Monte

Finite_exp_median$Monte<-Finite_exp_median_Monte  
Finite_exp_median$Exact<-medianexp$medianexp


Calibrated_Monte_MonotonSmooth1<-MonotonicSmoothing(x=Finite_exp_median$Size,y=Finite_exp_median$Calibrated_Monte,auto=TRUE)

Finite_exp_median$Calibrated_Monte_mono<-Calibrated_Monte_MonotonSmooth1$y

Monte_MonotonSmooth1<-MonotonicSmoothing(x=Finite_exp_median$Size,y=Finite_exp_median$Monte,auto=TRUE)

Finite_exp_median$Monte_mono<-Monte_MonotonSmooth1$y


Allresults_Monte101 <- Finite_exp_median[Finite_exp_median$Size %in% c(c(seq(from=16, to=400, by=8))),]#,8*round(2^seq(from=6, to=9, by=0.18))),]

theme_set(theme_bw(base_size=16))
gfg_plot <- ggplot(Finite_exp_median, aes(Size)) +
  geom_line(aes(y = Monte), color = "red", size=1) +
  geom_line(aes(y = Calibrated_Monte), color = "green", size=1)+geom_line(aes(y = Monte_mono), color = "blue", size=1)+
  geom_line(aes(y = Calibrated_Monte_mono), color = "black", size=1)+geom_line(aes(y = Exact), color = "purple", size=1)


gfg_plot
gfg_plot <- ggplot(Allresults_Monte101, aes(Size)) +
  geom_line(aes(y = Monte), color = "red", size=1) +
  geom_line(aes(y = Calibrated_Monte), color = "green", size=1)+geom_line(aes(y = Monte_mono), color = "blue", size=1)+
  geom_line(aes(y = Calibrated_Monte_mono), color = "black", size=1)+geom_line(aes(y = Exact), color = "purple", size=1)

gfg_plot

Finite_exp_median$Diff_Calibrated_Monte<-abs(Finite_exp_median$Exact-Finite_exp_median$Calibrated_Monte)
Finite_exp_median$Diff_Monte<-abs(Finite_exp_median$Exact-Finite_exp_median$Monte)

Finite_exp_median$Diff_Calibrated_Monte_mono<-abs(Finite_exp_median$Exact-Finite_exp_median$Calibrated_Monte_mono)
Finite_exp_median$Diff_Monte_mono<-abs(Finite_exp_median$Exact-Finite_exp_median$Monte_mono)

Finite_exp_median$Improvement1<-Finite_exp_median$Diff_Calibrated_Monte_mono/Finite_exp_median$Diff_Monte

Finite_exp_median$Improvement2<-Finite_exp_median$Diff_Calibrated_Monte_mono/Finite_exp_median$Diff_Calibrated_Monte

mean(Finite_exp_median$Improvement1)

median(Finite_exp_median$Improvement2)

colMeans(Finite_exp_median)

