
###############################################################################
#  copy all data in the base directory. 
#
      # training_data.ChIPseq.tar
      # training_data.DNASE_wo_bams.tar
      # training_data.RNAseq.tar
      # training_data.annotations.tar
#
#  Extract all files
#  In annotations directory, do:
#  Use bedtools to extract fasta from bed coordinates:
#  bedtools getfasta -fi hg19.genome.fa -bed test_regions.blacklistfiltered.bed -fo test_regions.blacklistfiltered.fa
#  bedtools getfasta -fi hg19.genome.fa -bed ladder_regions.blacklistfiltered.bed -fo ladder_regions.blacklistfiltered.fa
# * background file for ama suite 
#  in the writeup directory:
#  fasta-get-markov ../annotations/hg19.genome.fa hg19markov.bkg
#
#   
#
###############################################################################

# execute execbash.sh on the tf directory
tf<-'EGR1';
source('aux_funcs.R')
##################################### run main process##################
run_execbash_sh_on_tf_folder_after_this(tf)
# execute execbash.sh on the tf directory
preprocess_writeup(tf)
####################################run machine learning process
dd<-load_features(tf)
# take the 90/100 folds of each data/tissue
dat.train<-f_combinew(dd[train],1:90,kk=100)  

# to verify if number of lines i dat.train is ok --------------------------
#dat.train<-f_combinew(dd[train],1)
#sum(sapply(dd[train],function(x)dim(x)[[1]])/20)
# end ---------------------------------------------------------------------

# prepare to train --------------------------------------------------------
s<-function(x){ifelse(x=='B',1,0)}
trainRf<-setDF(f_combinew(dd[train],1,kk=100));
trainRF<-trainRf[,.SD,.SDcols=-c('score','index_nona',train[1])]
ladderSub<-dd[leaderboard]  
testXg<-f_combinew(dd[train],91:100,kk=100)  
trainXg<-f_combinew(dd[train],1:90,kk=100)  
trainXg<- xgb.DMatrix(as.matrix(trainXg[,.SD,.SDcols=-c('score','index_nona',train[1])])
                        ,label=as.matrix(s(trainXg[,.SD,.SDcols=train[1]][[1]])))
testXg<- xgb.DMatrix(as.matrix(testXg[,.SD,.SDcols=-c('score','index_nona',train[1])])
                      ,label=as.matrix(s(testXg[,.SD,.SDcols=train[1]][[1]])))                       

# end ---------------------------------------------------------------------




a<-s(trainXg[,.SD,.SDcols=train[1]][[1]])

xgbtrain<-function(trainXg,testXg){
  watchlist <- list(train = trainXg, eval = testXg)
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
  bst <-xgb.train(param=param, data =trainXg, nrounds=15,watchlist,early.stop.round=3)
  bstp<-predict(bst, dxg)
  # id<-getinfo(dxg,'nrow')
  save(bst,file=paste0('bst.RData'))
  return(bstp)
}
xgbtrainf<-function(trainRf,ladderSub){
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
 out<-function(i,mtry=2,ns=1,ntree=200){
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
