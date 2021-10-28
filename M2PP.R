rm(list = ls())

library(e1071)
library(caret)
library(randomForest)
library(Rcpp)
library(entropy)

source("calAUPR.R")
cppSN <- lapply("fastKgipMat.cpp", sourceCpp, verbose = FALSE)

LD <- read.table("data/target-disease-association matrix.txt")
LD=as.matrix(LD)

LM <- read.table("data/target-drug-association matrix.txt")
LM=as.matrix(LM)

DM <- read.table("data/disease-drug-association matrix.txt")
DM=as.matrix(DM)

DD1 <- read.table("data/disease MESH similarity.txt")
DD1=as.matrix(DD1)

DD2 <- read.table("data/disease DO similarity.txt")
DD2=as.matrix(DD2)

LL1 <-read.table("data/target secquence similarity.txt")
LL1=as.matrix(LL1)

LL2 <- read.table("data/target functional similarity.txt")
LL2=as.matrix(LL2)

set.seed(123)
kfold=5
times=5
nL=nrow(LD)
nD=ncol(LD)
DD3=0.5*DD1+0.5*DD2
LL3=0.5*LL1+0.5*LL2
exdrug=which(colSums(LM)==0)
LM=LM[,-exdrug]
DM=DM[,-exdrug]
nM=ncol(LM)
MM_D <- exp(-fastKgipMat(t(DM), 1))
MM_L <- exp(-fastKgipMat(t(LM), 1))
AUC=matrix(0,nrow = times,ncol = kfold)
AUPR=matrix(0,nrow = times,ncol = kfold)
Label_All=as.vector(LD)
positive_number=sum(LD)
positive_set_index=which(Label_All==1)
index0=which(Label_All==0)
for(time in 1:times){
  negative_set_index=sample(x=index0,size=length(positive_set_index))
  folds_positive=createFolds(y=positive_set_index,k=kfold)
  folds_negative=createFolds(y=negative_set_index,k=kfold)
  Fold=list()
  for(i in 1:kfold)
  {
    Fold[[i]]=union(positive_set_index[folds_positive[[i]]],negative_set_index[folds_negative[[i]]])
  }
  for(i in 1:kfold){
    test_positive=positive_set_index[folds_positive[[i]]]
    LD_temp=as.vector(LD)
    LD_temp[test_positive]=0
    LD_temp=matrix(LD_temp,nrow = nL,ncol = nD)
    DD3=0.5*DD1+0.5*DD2
    for(a in 1:nL){
      for(b in 1:nL){
        index_gi=which(LD[a,]!=0)
        index_gj=which(LD[b,]!=0)
        if((length(index_gi)!=0)&(length(index_gj)!=0)){
          sum1=0
          for(k in 1:length(index_gi)){
            sum1=sum1+max(DD3[index_gi[k],index_gj])
          }
          sum2=0
          for(k in 1:length(index_gj)){
            sum2=sum2+max(DD3[index_gj[k],index_gi])
          }
          LL3[a,b]=(sum1+sum2)/(length(index_gi)+length(index_gj))
        }else{
          LL3[a,b]=0
        }
      }
    }
    write.table(LL3,file=".txt",row.names = FALSE, col.names = FALSE)
    LL3=read.table(".txt")
    unlink(".txt")
    LL3=as.matrix(LL3)
    LL4 <- exp(-fastKgipMat(LD_temp, 1))
    LL=matrix(0,nrow = nL,ncol = nL)
    mea1=c()
    aa1=0
    mea2=c()
    aa2=0
    for(a in 1:nL){
      if(entropy(LL3[a,],unit = "log2")<=entropy(LL4[a,],unit = "log2")){
        aa1=aa1+1
        mea1[aa1]=a
      }else{
        aa2=aa2+1
        mea2[aa2]=a
      }
    }
    LL[mea1,mea1]=LL3[mea1,mea1]
    LL[mea2,mea2]=LL4[mea2,mea2]
    for(a in 1:aa1){
      if(entropy(LL3[mea1[a],mea2],unit = "log2")<=entropy(LL4[mea1[a],mea2],unit = "log2"))
        LL[mea1[a],mea2]= LL3[mea1[a],mea2]
      else
        LL[mea1[a],mea2]= LL4[mea1[a],mea2]
    }
    LL[mea2,mea1]=t(LL[mea1,mea2])
    for(a in 1:nD){
      for(b in 1:nD){
        index_gi=which(LD[,a]!=0)
        index_gj=which(LD[,b]!=0)
        if((length(index_gi)!=0)&(length(index_gj)!=0)){
          sum1=0
          for(k in 1:length(index_gi)){
            sum1=sum1+max(LL1[index_gi[k],index_gj])
          }
          sum2=0
          for(k in 1:length(index_gj)){
            sum2=sum2+max(LL1[index_gj[k],index_gi])
          }
          DD3[a,b]=(sum1+sum2)/(length(index_gi)+length(index_gj))
        }
        else{
          DD3[a,b]=0
        }
      }
    }
    write.table(DD3,file=".txt",row.names = FALSE, col.names = FALSE)
    DD3=read.table(".txt")
    unlink(".txt")
    DD3=as.matrix(DD3)
    DD4 <- exp(-fastKgipMat(t(LD_temp), 1))
    DD=matrix(0,nrow = nD,ncol = nD)
    mea1=c()
    aa1=0
    mea2=c()
    aa2=0
    for(a in 1:nD){
      if(entropy(DD3[a,],unit = "log2")<=entropy(DD4[a,],unit = "log2")){
        aa1=aa1+1
        mea1[aa1]=a
      }else{
        aa2=aa2+1
        mea2[aa2]=a
      }
    }
    DD[mea1,mea1]=DD3[mea1,mea1]
    DD[mea2,mea2]=DD4[mea2,mea2]
    for(a in 1:aa1){
      if(entropy(DD3[mea1[a],mea2],unit = "log2")<=entropy(DD4[mea1[a],mea2],unit = "log2"))
        DD[mea1[a],mea2]= DD3[mea1[a],mea2]
      else
        DD[mea1[a],mea2]= DD4[mea1[a],mea2]
    }
    DD[mea2,mea1]=t(DD[mea1,mea2])
    Feature=matrix(0,nrow = nL*nD,ncol = 9)
    LLD_0=matrix(0,nrow = nL,ncol = nD)
    LDD_0=matrix(0,nrow = nL,ncol = nD)
    for(a in 1:nL){
      for(b in 1:nD){
        index1=which((LL[a,]!=0)&(LD_temp[,b]==0))
        LLD_0[a,b]=mean(LL[a,index1])
        index2=which((LD_temp[a,]==0)&(DD[,b]!=0))
        LDD_0[a,b]=mean(DD[index2,b])
      }
    }
    Feature[,1]=as.matrix(as.vector(LLD_0))
    Feature[,2]=as.matrix(as.vector(LDD_0))
    diag(DD)=0
    diag(LL)=0
    maxscore=matrix(0,nrow = nL,ncol = nD)
    for(a in 1:nL){
      for(b in 1:nD){
        maxg=max(LL[a,])
        indexmaxg=which(LL[a,]==maxg)
        maxd=max(DD[,b])
        indexmaxd=which(DD[,b]==maxd)
        maxscore[a,b]=sum(maxg*LD_temp[indexmaxg,b])+sum(maxd*LD_temp[a,indexmaxd])+sum(maxg*maxd*LD_temp[indexmaxg,indexmaxd])
      }
    }
    Feature[,3]=as.matrix(as.vector(maxscore))
    diag(DD)=1
    diag(LL)=1
    LM_score=LM%*%MM_L
    for(a in 1:nL){
      for(b in 1:nM){
        if(LM_score[a,b]!=0)
          LM_score[a,b]=LM_score[a,b]/length(which((LM[a,]==1)&(MM_L[,b]!=0)))
      }
    }
    DM_score=DM%*%MM_D
    for(a in 1:nD){
      for(b in 1:nM){
        if(DM_score[a,b]!=0)
          DM_score[a,b]=DM_score[a,b]/length(which((DM[a,]==1)&(MM_D[,b]!=0)))
      }
    }
    infer=matrix(0,nrow = nL,ncol = nD)
    for(a in 1:nL){
      for(b in 1:nD){
        indexdrug=which((LM[a,]==1)&(DM[b,]==1))
        if(length(indexdrug)!=0)
          infer[a,b]=max(DM_score[b,indexdrug]/LM_score[a,indexdrug])
      }
    }
    Feature[,4]=as.matrix(as.vector(infer))
    LL_int=as.vector(LL)
    LL_int[which(LL_int!=0)]=1
    LL_int=matrix(LL_int,nrow = nL,ncol = nL)
    DD_int=as.vector(DD)
    DD_int[which(DD_int!=0)]=1
    DD_int=matrix(DD_int,nrow = nD,ncol = nD)
    Fea1=(LL%*%LD_temp)/(LL_int%*%LD_temp)
    Fea1[is.na(Fea1)]=0
    Fea2=(LD_temp%*%DD)/(LD_temp%*%DD_int)
    Fea2[is.na(Fea2)]=0
    Fea3=(LL%*%LL%*%LD_temp)/(LL_int%*%LL_int%*%LD_temp)
    Fea3[is.na(Fea3)]=0
    Fea4=(LL%*%LD_temp%*%DD)/(LL_int%*%LD_temp%*%DD_int)
    Fea4[is.na(Fea4)]=0
    Fea6=(LD_temp%*%DD%*%DD)/(LD_temp%*%DD_int%*%DD_int)
    Fea6[is.na(Fea6)]=0
    Feature[,5]=as.matrix(as.vector(Fea1))
    Feature[,6]=as.matrix(as.vector(Fea2))
    Feature[,7]=as.matrix(as.vector(Fea3))
    Feature[,8]=as.matrix(as.vector(Fea4))
    Feature[,9]=as.matrix(as.vector(Fea6))
    center <- sweep(Feature, 2, apply(Feature, 2, min),'-')
    R <- apply(Feature, 2, max) - apply(Feature,2,min)
    Feature<- sweep(center, 2, R, "/")
    train=unlist(Fold[-i])
    Label_train=as.matrix(as.vector(LD_temp))[train,]
    model=randomForest(x=Feature[train,],y=Label_train,ntree = 500,importance = FALSE)
    test=setdiff(c(1:(nL*nD)),unlist(Fold[-i]))
    prob=predict(model,Feature[test,])
    Label_test=as.matrix(as.vector(LD))[test,]
    result=calAUPR(Label_test,prob)
    print(result)
    AUC[time,i]=result[1]
    AUPR[time,i]=result[2]
  }
  print(paste0(time," time 5-fold CV: auroc mean = ",mean(AUC[time,]),"; aupr mean =",mean(AUPR[time,])))
}
print(paste0("Five times 5-fold CV: mean AUC = ",mean(AUC),"; mean AUPR = ",mean(AUPR)))