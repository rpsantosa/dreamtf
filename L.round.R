
# ##########################################################################
#1-denase
#2-makefasta 
#3-meme scores
#4-make r objects
#5-analysis
#############################################################################
#!/usr/bin/env r

args = commandArgs(trailingOnly=TRUE)
tf<-args[1]
if (length(args)==0) {
  stop("usage: Rscript  /home/ricardo/hd/projects/dream_tf_competition/rscript/L.round.R   'ATF1'
", call.=FALSE)
} 

library(xgboost)
library(ranger)
require(biovizBase)
require(ggbio)
require(R.utils)
require(gdata)
require(plyr)
library(GenomicRanges)
library(ShortRead)  #clean
library(rtracklayer)
library(gdata)
library(Biostrings)
library(data.table)
library(caret)
  pwriteup<-'/home/ricardo/hd/projects/dream_tf_competition/data/writeup/'
  ptf2<-'/home/ricardo/hd/projects/dream_tf_competition/data/TF2/'
  tfs<-read.xls(paste0(pwriteup,"tfs.xls"))
  nk=100
  source(paste0(pwriteup,'functions.R')) 
  xgbtrain<-function(dd,dxg){
    t_trainx<-setdiff(t_train,'SK-N-SH')
    train.xg<-f_combine(dd[t_trainx],1:90,1:length(t_trainx))
    test.xg<-f_combine(dd[t_trainx],91:100,1:length(t_trainx))
    train.xgd<- xgb.DMatrix(as.matrix(train.xg[,-ncol(train.xg),with=F]),label=as.matrix(as.numeric(train.xg[,bind])))
    test.xgd<- xgb.DMatrix(as.matrix(test.xg[,-ncol(test.xg),with=F]),label=as.matrix(as.numeric(test.xg[,bind])))
    watchlist <- list(train = train.xgd, eval = test.xgd)
    eta<-.479  #0.01
    max.depth<-31  #40
    colsample.bytree<-.9779
    maxd<- 5.2
    param <- list(objective = "binary:logistic",eval_metric= 'auc',
                  eta=eta,
                  subsample=0.632,
                  'maximize'=T,
                  colsample_bytree=colsample.bytree,
                  max_depth=max.depth,
                  max_delta_step=maxd
    )
    bst <-xgb.train(param=param, data =train.xgd, nrounds=15,watchlist,early.stop.round=3)
    bstp<-predict(bst, dxg)
    # id<-getinfo(dxg,'nrow')
    save(bst,file=paste0('bst.RData'))
    return(bstp)
  }
  xgbtrainf<-function(dd,dxg){
    # total<-ldply(dd[t_train],rbind)
    # total<-total[,-1];#total<-total[,c(1:7)]
    t_trainx<-setdiff(t_train,'SK-N-SH')
    train.xg<-f_combine(dd[t_trainx],1:.9*nk,1:length(t_trainx))
    test.xg<-f_combine(dd[t_trainx],.9*nk:nk,1:length(t_trainx))
    train.xgd<- xgb.DMatrix(as.matrix(train.xg[,-ncol(train.xg),with=F]),label=as.matrix(as.numeric(train.xg[,bind])))
    test.xgd<- xgb.DMatrix(as.matrix(test.xg[,-ncol(test.xg),with=F]),label=as.matrix(as.numeric(test.xg[,bind])))
    #out<-sparse.model.matrix(bind ~.-1, data = total)
    watchlist <- list( train = train.xgd, eval = test.xgd)
    #out<- xgb.DMatrix(as.matrix(total),label=as.matrix(total$bind))
    #outr<-total$bind
    eta<-0.01 #0.01
    max.depth<-14 #40
    colsample.bytree<-.6867#.6
    maxd<- 20 #3.5518
    param <- list(objective = "binary:logistic",eval_metric= 'auc',
                  eta=eta,
                  subsample=0.632,
                  'maximize'=T,
                  colsample_bytree=colsample.bytree,
                  max_depth=max.depth,
                  max_delta_step=maxd
    )
    bstf <-xgb.train(param=param, data =train.xgd, nrounds=20,watchlist,early.stop.round=3)
    bstp<-predict(bstf, dxg)
    #id<-getinfo(dxg,'nrow')
    save(bstf,file=paste0('bstf.RData'))
    return(bstp)
  }
  rftrain<-function(dd,dat.testsub){
    t_trainx<-setdiff(t_train,'SK-N-SH')
    dat.train<-f_combine(dd[t_trainx],1,1:length(t_trainx))
    # ww<-ranger(dependent.variable.name = "bind",
    #            data = dat.train,  importance='impurity',min.node.size=ns,
    #            write.forest=T,num.trees = ntree,
    #            mtry=mtry)
    # splitselectw<-ff(importance(outx))
    out<-function(i,mtry=2,ns=1,ntree=200){
      t_trainx<-setdiff(t_train,'SK-N-SH')
      dat.train<-setDF(f_combine(dd[t_trainx],i,1:length(t_trainx)))
      outx<-ranger(dependent.variable.name = "bind",
                   data = dat.train,  importance='impurity',min.node.size=ns,
                   write.forest=T,num.trees = ntree,
                   mtry=mtry)
      print(i)
      return(outx)
    }
    ra<-lapply(1:nk,out)
    pred<-sapply(1:nk,function(i){out<-predict(ra[[i]],dat.testsub)$predictions})
    predxm<-apply(pred,1,mean)
    save(ra,file='ra.RData')
    return(predxm)
  }
  
  DIR_TF<-paste0(ptf2,tf,'/')
  con_submission<-paste0(pwriteup,'ladder_regions.blacklistfiltered.bed')
  con_final<-paste0(pwriteup,'test_regions.blacklistfiltered.bed.gz')
  j<-which(tfs[,1]==tf)
  setwd(DIR_TF)
  t_train<-strsplit(as.character(tfs[j,2]),',')[[1]];t_train<-sub('[[:space:]]','',t_train)
  t_test<-strsplit(as.character(tfs[j,3]),',')[[1]];t_test<-sub('[[:space:]]','',t_test)
  res<-load_train_and_createfolds(tf,t_train,kk=nk,t_test,firstn=F)
  dd<-res[[1]];dd<-lapply(dd,function(x){x[,sapply(x,class)=='numeric']})
  
  dfolds<-res[[2]]
  if(length(t_test)!=0){
    for(r in 1:length(t_test)){
      dat.testsub<-dd[length(t_train)+r][[1]]
      dat.testf<-dd[length(c(t_train,t_test))+r][[1]]
      dtestsub <- xgb.DMatrix(as.matrix(dat.testsub))
      dtestf <- xgb.DMatrix(as.matrix(dat.testf))
      dd1<-lapply(dd[t_train],
                  function(x){subset(x,select=c('yavg','ymax','y75','bind'))})
      
      bstp<-xgbtrain(dd,dtestsub)
      rfp<-rftrain(dd,dat.testsub)
      bstpf<-xgbtrainf(dd1,dtestf)
      #rfpf<-rftrainf(dd1,dat.testf)
      fcs<-function(bstp,rfp,bstpf){
        load(file=paste0('chip_dnase','.',tf,'.',t_test[r],'.RData'))
        sub<-import(con<-(gzfile(con_submission)),  format = "BED") 
        mcols(sub)<-data.frame(value=1:length(sub))
        xfd<-data.frame(value=chip_dnase$value ,pred2=(bstp+rfp))
        xff<-data.frame(value=1:length(bstpf),pred=(bstpf))
        xfd<-setDT(xfd);xff<-setDT(xff)
        xfa<-merge(xff,xfd,by='value',all.x=T,all.y=F)
        setDF(xfa)
        xfa[is.na(xfa)]<-0;xfa$submit<-apply(xfa[,c(2,3)],1,sum)
        xf3<-data.frame(mold(sub),pred2=ff(xfa$submit))
        xf4<-xf3[,c(1,2,3,8)]
        xf4$start<-xf4$start-1;#xf4$pred<-format(xf4$pred,scientific = T,digits = 2)
        xf5<-xf4;
        xf5$start<-as.integer(xf5$start)
        xf5$end<-as.integer(xf5$end)
        filename=paste0('L.',tf,'.',t_test[r],'.tab')
        write.table(xf5,file=filename,quote=F,col.names=F,row.names=F,sep='\t')
        if(file.exists(paste0(filename,'.gz'))){file.rename(paste0(filename,'.gz'),paste0(filename,'.gz.old'))}
        gzip(filename, destname=sprintf("%s.gz", filename),skip=F)
        return(NULL)
      }
      fcs(bstp,rfp,bstpf)
      filename<-paste0('L.',tf,'.',t_test[r],'.tab.gz')
      #loadsyn(filename)
    }
  }else{
    print('There is no tissue to test')
  }
