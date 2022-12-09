#' @docType data
#' @title A illustration dataset
#' @name twotypeVehicle
#' @description A dataset used to illustrate the performance of SMAP and logistic method.
#' @examples
#' \dontrun{
#' data(twotypeVehicle)
#' attach(twotypeVehicle)
#' }
NULL

#' @title A SMAP function to estimate the weights for each sub model using R
#' @description  A SMAP function to estimate the weights for each sub model using R
#' @param data which data you want to use to build a SMAP model
#' @param p the number of variables
#' @param p_dis the number of discrete variables
#' @param n_D the sample size
#' @return the weights for each sub-model by SMAP and improved SMAP
#' @examples 
#' \dontrun{
#'      data<-data(twotypeVehicle)
#'      w<-SMAP_weights{data,10,3}
#' }
#' @import mvtnorm
#' @import splines2
#' @import nnet
#' @import quadprog
#' @import splitTools
#' @import stats
#' @import stats4
#' @import boot
#' @import DAAG
#' @import bootstrap
#' @export
SMAP_weights<-function(data,p,p_dis,n_D){
  #编写函数Design.matrix
  Design.matrix <- function(U, x, kn, degree,p,n)
  {
    q=kn+degree+1
    B=bSpline(U,df=q,intercept=TRUE,degree=degree)          
    Num<-p*(q)
    xx<-matrix(0,n,Num)
    Bb<-array(0,dim=c(p,Num,n)) 
    for (j in 1:n)  
    {
      for (i in 1:p)
      {Bb[i,((i-1)*(q)+1):(i*(q)),j]<-B[j,] }
      xx[j,]<-x[j,]%*%Bb[,,j]
    }
    design = list(Bb = Bb, xx = xx,B=B)
    return(design)
  } 
  
  ###数据输入
  realdata<-data
  J=2#二分类预测问题
  HitRite<-0
  HitRite1<-0
  HAT.W.VEC_B=matrix(0,1,p-p_dis)
  HAT.W.VEC_B_new=matrix(0,1,p-p_dis)
  ##############改进后的权重分析###############
  ############model averaging varying coefficient model#(SMAP)
  ####estimate model weights################
  ####B-spline################
  #数据集划分(训练集：测试集：验证集=8:1:1)
  train.rows<-sample(n_D,n_D*0.6)
  valid.rows<-sample(setdiff(seq(1:n_D),train.rows),n_D*0.2)
  test.rows<-setdiff(seq(1:n_D),union(train.rows,valid.rows))
  train<-realdata[train.rows,]
  valid<-realdata[valid.rows,]
  test<-realdata[test.rows,]
  
  n_tr=dim(train)[1];n_val=dim(valid)[1];n_test=dim(test)[1]
  x_tr=train[,1:p]
  y_tr_sta=train[,(p+1)]
  y_tr_sta[which(y_tr_sta==0)]=2
  y_tr=matrix(0,n_tr,2)
  y_tr[which(y_tr_sta==1),1]=1
  y_tr[which(y_tr_sta==2),2]=1
  #y_sta_tr=y_sta_D[out]
  x_val=valid[,1:p]
  y_val_sta=valid[,(p+1)]
  y_val_sta[which(y_val_sta==0)]=2
  y_val=matrix(0,n_val,2)
  y_val[which(y_val_sta==1),1]=1
  y_val[which(y_val_sta==2),2]=1
  x_test=test[,1:p]
  y_test_sta=test[,(p+1)]
  y_test_sta[which(y_test_sta==0)]=2
  y_test=matrix(0,n_test,2)
  y_test[which(y_test_sta==1),1]=1
  y_test[which(y_test_sta==2),2]=1
  
  n=n_D-n_test-n_val
  kn=floor((n)^(1/5));degree=1;q=kn+degree+1
  h=1*(((J-1)*p)/n)^(0.2)
  
  
  #############bspline改进前###########################
  hat.P=array(0,dim=c(n_tr,2,p-p_dis))
  hat.P.test=array(0,dim=c(n_test,2,p-p_dis))
  for (j in 1:(p-p_dis)){
    x_sub=Design.matrix(x_tr[,j+p_dis], as.matrix(cbind(1,x_tr[,-(j+p_dis)])),kn,degree,p=p,n=n_tr)$xx
    x_test_sub=Design.matrix(x_test[,j+p_dis],as.matrix(cbind(1,x_test[,-(j+p_dis)])),kn, degree,p=p,n=n_test)$xx
    lm0=multinom(y_tr~x_sub+0)
    hat.P[,,j]=lm0$fitted.values
    
    hat.mu1=x_test_sub%*%(coefficients(lm0)[1,])
    hat.p1=1/(1+exp(hat.mu1))
    hat.p2=exp(hat.mu1)/(1+exp(hat.mu1))
    hat.P.test[,,j]=cbind(hat.p1,hat.p2)
  }
  
  hat.P.A=y.vec=NULL 
  for(c0 in 1:n_tr) {
    hat.P.A=rbind(hat.P.A,hat.P[c0,,])
    y.vec=c(y.vec,y_tr[c0,])
  }
  lm.opt=solve.QP(Dmat=t(hat.P.A)%*%hat.P.A, 
                  dvec=t(y.vec)%*%hat.P.A, Amat=cbind(rep(1,p-p_dis),diag(p-p_dis)),
                  bvec=c(1,rep(0,p-p_dis)), meq=1) 
  w.opt=lm.opt$solution
  HAT.W.VEC_B[1,]=w.opt
  
  Hat.w.P=matrix(0,n_test,2)
  for(i in 1:n_test){
    for(j in 1:J){
      Hat.w.P[i,j]=sum(hat.P.test[i,j,]*w.opt) 
    }
  }
  Hat.y=apply(Hat.w.P,1,which.max)
  HitRate=mean(Hat.y==y_test_sta)
  
  
  #############bspline改进后###########################
  hat.P.val=array(0,dim=c(n_val,2,p-p_dis))
  hat.P.test=array(0,dim=c(n_test,2,p-p_dis))
  for (j in 1:(p-p_dis)){
    x_sub=Design.matrix(x_tr[,j+p_dis], as.matrix(cbind(1,x_tr[,-(j+p_dis)])),kn,degree,p=p,n=n_tr)$xx
    x_val_sub=Design.matrix(x_val[,j+p_dis],as.matrix(cbind(1,x_val[,-(j+p_dis)])),kn,degree,p=p,n=n_val)$xx
    x_test_sub=Design.matrix(x_test[,j+p_dis],as.matrix(cbind(1,x_test[,-(j+p_dis)])),kn, degree,p=p,n=n_test)$xx
    lm0=multinom(y_tr~x_sub+0)
    #hat.P[,,j]=lm0$fitted.values
    
    hat.mu1_val=x_val_sub%*%(coefficients(lm0)[1,])
    hat.p1_val=1/(1+exp(hat.mu1_val))
    hat.p2_val=exp(hat.mu1_val)/(1+exp(hat.mu1_val))
    hat.P.val[,,j]=cbind(hat.p1_val,hat.p2_val)
    
    hat.mu1=x_test_sub%*%(coefficients(lm0)[1,])
    hat.p1=1/(1+exp(hat.mu1))
    hat.p2=exp(hat.mu1)/(1+exp(hat.mu1))
    hat.P.test[,,j]=cbind(hat.p1,hat.p2)
  }
  
  hat.P.A=y.vec=NULL 
  for(c0 in 1:n_val) {
    hat.P.A=rbind(hat.P.A,hat.P.val[c0,,])
    y.vec=c(y.vec,y_val[c0,])
  }
  
  lm.opt=solve.QP(Dmat=t(hat.P.A)%*%hat.P.A, 
                  dvec=t(y.vec)%*%hat.P.A, Amat=cbind(rep(1,(p-p_dis)),diag(p-p_dis)),
                  bvec=c(1,rep(0,p-p_dis)), meq=1) 
  w.opt=lm.opt$solution
  HAT.W.VEC_B_new[1,]=w.opt
  
  Hat.w.P=matrix(0,n_test,2)
  for(i in 1:n_test){
    for(j in 1:J){
      Hat.w.P[i,j]=sum(hat.P.test[i,j,]*w.opt) 
    }
  }
  Hat.y=apply(Hat.w.P,1,which.max)
  HitRate1=mean(Hat.y==y_test_sta)
  
  results<-list(rbind(HAT.W.VEC_B,HAT.W.VEC_B_new),rbind(HitRate,HitRate1))
  names(results)<-c("Weights","Accuracy")
  return(results)
}

