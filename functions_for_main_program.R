#######################set the tf and base directory##############
tf<-'FOXA1'
base<-'~/hd/projects/dream_tf_competition/data'
writeup <- file.path(base,'writeup')
subDir <- file.path(base,'writeup','results')
if (!file_test("-d",writeup)){
  dir.create(file.path(writeup))
}
if (!file_test("-d",subDir)){
  dir.create(file.path(subDir))
}
if (!file_test("-d",file.path(subDir,tf))){
  dir.create(file.path(subDir,tf))
}
setwd(writeup)
###################################################################

ff<-function(x){(x-min(x))/(max(x)-min(x))}
load_lables_tsv<-function(tf){
  con_chipseq_label_tf<-file.path(base,'ChIPseq/labels',paste0(tf,'.train.labels.tsv.gz'))
  x<-data.table::fread(paste0('gzip -dc ',con_chipseq_label_tf))
  chip<-makeGRangesFromDataFrame(x,keep.extra.columns = T,starts.in.df.are.0based=T)
  names(mcols(chip))<-sub('\\.','-',names(mcols(chip)))   # set MCF.7 to MCF-7 for example
  #mcols(chip)<-data.frame(mcols(chip),value=1:length(chip)) # to indices purposes
  return(chip)
}
remove_na_save_index<-function(){
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(GenomicRanges)
  library(data.table)
  library(DNAshapeR)
  train<-file.path( annotationDir, 'train_regions.blacklistfiltered.bed.gz')
  test<-file.path( annotationDir,'test_regions.blacklistfiltered.bed.gz')
  ladder<-file.path( annotationDir,'ladder_regions.blacklistfiltered.bed.gz')
  fremovena<-function(x){
    ## remove sequences where there are any "N"
    chip<-makeGRangesFromDataFrame(x,keep.extra.columns = F,starts.in.df.are.0based=T)
    chip$index_nona<-1:nrow(x)
    xb<-getSeq(Hsapiens,chip)
    db<-data.table(letterFrequency(xb,  c('A','C','T','G','N')))
    #db<-data.table(oligonucleotideFrequency(xb, width=1));
    #if(ncol(x)>3){db<-cbind(db,x[,-c(1,2,3),with=F])}
    #h<-db[,A!=0 & T!=0 & C!=0 & G!=0]
    #if(ncol(x)>3){db<-cbind(db,x[,-c(1,2,3),with=F])}
    h<-db[,A!=0 & T!=0 & C!=0 & G!=0 & N ==0 ]
    return(chip[h])
  }
  
  xtest<-data.table::fread(paste0('gzip -dc ',test));names(xtest)<-c('chr','start','end')
  x<-fremovena(xtest);index_nona<-x$index_nona
  rm(xtest);gc()
  getFasta(x,	Hsapiens,	width	= 200,	filename	= file.path(annotationDir,'test_nona.fa'))
  save(index_nona,file=file.path(annotationDir,'test_nona.RData'))
  
  xladder<-data.table::fread(paste0('gzip -dc ',ladder));names(xladder)<-c('chr','start','end')
  x<-fremovena(xladder);index_nona<-x$index_nona
  rm(xladder);gc()
  getFasta(x,	Hsapiens,	width	= 200,	filename	= file.path(annotationDir,'ladder_nona.fa'))
  save(index_nona,file=file.path(annotationDir,'ladder_nona.RData'))
  
  xtrain<-data.table::fread(paste0('gzip -dc ',train));names(xtrain)<-c('chr','start','end')
  x<-fremovena(xtrain);index_nona<-x$index_nona
  rm(xtrain);gc()
  getFasta(x,	Hsapiens,	width	= 200,	filename	= file.path(annotationDir,'train_nona.fa'))
  save(index_nona,file=file.path(annotationDir,'train_nona.RData'))
  rm(x);gc()
}
scoresh<-function(tf,filesh,filescore,filename){
  text<-paste0(memeAma,'/ama  --o-format gff ' ,
               '  ',writeup,'/pwm_plus/',tf,'.txt',
               '  ',file.path(annotationDir,filename),
               '  ',writeup,'/hg19markov.bkg > ',
               '  ',tfDir,'/',filescore,'.',tf,'.txt')
  fileConn<-file.path(tfDir,paste0(filesh,'.',tf,'.sh'))
  writeLines(text, fileConn)
}
read_train_leaderboard_test<-function(tf){
  library(gdata)
  library(data.table)
  base<-'~/hd/projects/dream_tf_competition/data'
  writeup <- file.path(base,'writeup')
  tfs<-read.xls(file.path(writeup,"tfs.xls"))
  setkey(setDT(tfs),F.Name)
  test<-strsplit(as.character(gsub('\xa0','',tfs[tf,Final.Submission.Cell.Types])),',')[[1]]
  leaderboard<-strsplit(as.character(tfs[tf,Leaderboard.Cell.Types]),',')[[1]]
  train<-strsplit(as.character(tfs[tf,Training.Cell.Types]),',')[[1]]
  test<-sub('[[:space:]]','',test);#test<-sub('-','.',test)
  leaderboard<-sub('[[:space:]]','',leaderboard);#leaderboard<-sub('-','.',leaderboard)
  train<-sub('[[:space:]]','',train);  #train<-sub('-','.',train)
  return(list(train,leaderboard,test))
}
run_execbash_sh_on_tf_folder_after_this<-function(tf){
  aux<-read_train_leaderboard_test(tf)
  train<-aux[[1]];leaderboard<-aux[[2]];test<-aux[[3]]
  if(!identical(test, character(0))){
    for(e in test){
      scoresh(tf=tf,'score_test','test','test_nona.fa')
    }
  }
  if(!identical(leaderboard, character(0))){
    for(e in leaderboard){
      scoresh(tf=tf,'score_ladder','ladder','ladder_nona.fa')
    }
  }
  if(!identical(train, character(0))){
    for(e in train){
      scoresh(tf=tf,'score_train','train','train_nona.fa')
    }
  }
  file.copy(file.path(writeup,'execbash.sh'),file.path(tfDir),overwrite = T)
}
# execute execbash.sh on the tf directory

faux<-function(x){  
  ll<-table(x$V10)[[1]]
  np<-dim(x)[1]/ll;y<-x[1:ll,]
  if(np>1){
    for(i in 1:(np-1)){
      y<-cbind(y,x[(i*ll+1):(ll+i*ll),])
    }
    y<-y[,c(-seq(3,3*np,3),-seq(1,1+3*np,3))]
  }else{
    y<-x$V6
  }
  return(y)
}
feature_tf_test<-function(tf,test){
  con_test<-file.path(base,'annotations/test_regions.blacklistfiltered.bed.gz')
  x<-data.table::fread(paste0('gzip -dc ',con_test));names(x)<-c('chr','start','stop')
  grange_test<-makeGRangesFromDataFrame(x,keep.extra.columns = T,starts.in.df.are.0based=T) 
  load(file.path(annotationDir,'test_nona.RData'))  # load index with no 'N's - index_nona
  grange_test<-grange_test[index_nona];grange_test$index_nona<-index_nona
  score<-fread(paste0(tfDir,'/test.',tf,'.txt'), data.table=T,
               colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
  rm(x);gc()
  for( e in test){
    #file to the final submission: train and held out
    con_dnase_peak<-paste0(base,'/essential_training_data/DNASE/peaks/conservative/',
                           'DNASE.',e,'.conservative.narrowPeak.gz')
    con_dnase_fc<-paste0(base,'/essential_training_data/DNASE/fold_coverage_wiggles/',
                         'DNASE.',e,'.fc.signal.bigwig')
    extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")
    #load files:
    aux <- file.path(con_dnase_fc);bwf <- BigWigFile(aux)
    dnasefc <- import(bwf)
    dnasepeak <- import(gzfile(con_dnase_peak), format = "BED", extraCols = extraCols_narrowPeak)
    
    #overlaps:
    over<-findOverlaps(grange_test,dnasepeak)
    x<-data.table(data.frame(qh=over@queryHits,sh=over@subjectHits,
                             mcols(dnasepeak[over@subjectHits])[-c(1,4,6)]))  # remove c('name','peak')
    #indextest=grange_test[unique(over@queryHits)]$index_nona) 
    
    ##### test variables#########
    # qplot(pValue, signalValue, colour = bind, shape = bind, 
    #       data = x)
    # p <- ggplot(x[,.(pValue,signalValue,bind)], aes(pValue,signalValue))
    #p + geom_boxplot(aes(colour =bind))
    ######test     extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")
    ########
    x2<-x[, lapply(.SD,max), by=qh,.SDcols = -c('sh')]            # ==> toooo fast!
    x2m<-x[, lapply(.SD,mean), by=qh,.SDcols = -c('sh')] ;names(x2m)<-paste0(names(x2m),'m')
    #chip_dnase<-chip[unique(over@queryHits)] 
    
    grange_test_peak<-grange_test[unique(over@queryHits)]
    x3<-cbind(x2[,.(score,signalValue,qValue)],
              x2m[,.(scorem,signalValuem,qValuem)])
    mcols(grange_test_peak)<-cbind(mcols(grange_test_peak),as.data.frame(x3))
    rm(x3,x2m,x2,x);gc()
    
    over<-findOverlaps(grange_test_peak,dnasefc)
    x<-data.table(qh=over@queryHits,sh=over@subjectHits,maxfc=dnasefc[over@subjectHits]$score)
    x1<-x[, lapply(.SD,max), by=qh,.SDcols = -c('sh')]  ;
    x1m<-x[, lapply(.SD,mean), by=qh,.SDcols = -c('sh')]
    grange_test_peak_fc<-grange_test_peak;
    rm(grange_test_peak);gc()
    grange_test_peak_fc$maxfc<-x1[,maxfc];  grange_test_peak_fc$maxfcm<-x1m[,maxfc]
    featurefcpeak<-setDT(data.frame(mcols(grange_test_peak_fc)))
    rm(x,x1,x1m);gc()
    
    #merge with motif scores
    
    #setDT(score)
    motifsc<-score[,V6] #faux(score)
    dfa<-data.table(motifsc=motifsc,index_nona=index_nona)
    feature<-merge(featurefcpeak,dfa,by='index_nona')
    save(feature,
         file=file.path(tfDir,paste0('feature_test_',e,'.RData')))
  }
}  
feature_tf_train<-function(tf,train){
  grange_train_labels<-load_lables_tsv(tf)
  load(file.path(annotationDir,'train_nona.RData'))  # load index with no 'N's - index_nona
  grange_train_labels<-grange_train_labels[index_nona];grange_train_labels$index_nona<-index_nona
  score<-fread(file.path(tfDir,paste0('train.',tf,'.txt')), data.table=T,
               colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
  #setDT(score)
  for( e in train){
    ex<-sub('-','\\.',e) # in  labelfiles, e is written with '.'
    grange_train_labels_loop<-grange_train_labels[,c(e,'index_nona')]
    #remove ambiguous
    a<-setDT(as.data.frame(mcols(grange_train_labels_loop)))
    al<-a[,grep('A',get(ex),invert = T)] # lines dont have 'A'
    index_nona_noA<-a[al,][,index_nona]
    
    # aux<-data.table(bind=mcols(grange_train_labels_loop)[[1]]);idnoA<-aux[,bind]!='A'
    # grange_train_labels_loop<-grange_train_labels_loop[idnoA]
    
    con_dnase_peak<-paste0(base,'/essential_training_data/DNASE/peaks/conservative/',
                           'DNASE.',e,'.conservative.narrowPeak.gz')
    con_dnase_fc<-paste0(base,'/essential_training_data/DNASE/fold_coverage_wiggles/',
                         'DNASE.',e,'.fc.signal.bigwig')
    extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")
    #load files:
    aux <- file.path(con_dnase_fc);bwf <- BigWigFile(aux)
    dnasefc <- import(bwf)
    dnasepeak <- import(gzfile(con_dnase_peak), format = "BED", extraCols = extraCols_narrowPeak)
    
    #overlaps:
    over<-findOverlaps(grange_train_labels_loop,dnasepeak)
    x<-data.table(data.frame(qh=over@queryHits,sh=over@subjectHits,
                             mcols(dnasepeak[over@subjectHits])[-c(1,4,6)]))  # remove c('name','pValue','peak')
    #indextrain=grange_train_labels_loop[unique(over@queryHits)]$index_nona) 
    
    ##### train variables#########
    # qplot(pValue, signalValue, colour = bind, shape = bind, 
    #       data = x)
    # p <- ggplot(x[,.(pValue,signalValue,bind)], aes(pValue,signalValue))
    #p + geom_boxplot(aes(colour =bind))
    ######train     extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")
    ########
    x2<-x[, lapply(.SD,max), by=qh,.SDcols = -c('sh')]            # ==> toooo fast!
    x2m<-x[, lapply(.SD,mean), by=qh,.SDcols = -c('sh')] ;names(x2m)<-paste0(names(x2m),'m')
    #chip_dnase<-chip[unique(over@queryHits)] 
    
    grange_train_labels_peak<-grange_train_labels_loop[unique(over@queryHits)]
    x3<-cbind(x2[,.(score,signalValue,qValue)],
              x2m[,.(scorem,signalValuem,qValuem)])
    mcols(grange_train_labels_peak)<-cbind(mcols(grange_train_labels_peak),as.data.frame(x3))
    rm(x3,x2m,x2,x);gc()
    
    over<-findOverlaps(grange_train_labels_peak,dnasefc)
    x<-data.table(qh=over@queryHits,sh=over@subjectHits,maxfc=dnasefc[over@subjectHits]$score)
    x1<-x[, lapply(.SD,max), by=qh,.SDcols = -c('sh')]  ;
    x1m<-x[, lapply(.SD,mean), by=qh,.SDcols = -c('sh')]
    grange_train_labels_peak_fc<-grange_train_labels_peak;
    rm(grange_train_labels_peak);gc()
    grange_train_labels_peak_fc$maxfc<-x1[,maxfc];  grange_train_labels_peak_fc$maxfcm<-x1m[,maxfc]
    featurefcpeak<-setDT(data.frame(mcols(grange_train_labels_peak_fc)))
    rm(dnasefc,dnasepeak,x,x1,
       x1m,grange_train_labels_peak_fc);gc()
    
    #merge with motif scores
    
    motifsc<-score[,V6]#faux(score)
    dfa<-data.table(motifsc=motifsc,index_nona=index_nona)
    feature<-merge(featurefcpeak,dfa,by='index_nona')
    feature<-feature[index_nona %in% index_nona_noA,]
    save(feature,
         file=file.path(tfDir,paste0('feature_train_',e,'.RData')))
  }
}
feature_tf_ladder<-function(tf,leaderboard){
  con_ladder<-file.path(base,'annotations/ladder_regions.blacklistfiltered.bed.gz')
  x<-data.table::fread(paste0('gzip -dc ',con_ladder));names(x)<-c('chr','start','stop')
  grange_ladder<-makeGRangesFromDataFrame(x,keep.extra.columns = T,starts.in.df.are.0based=T) 
  load(file.path(annotationDir,'ladder_nona.RData'))  # load index with no 'N's - index_nona
  grange_ladder<-grange_ladder[index_nona];grange_ladder$index_nona<-index_nona
  score<-fread(paste0(tfDir,'/ladder.',tf,'.txt'), data.table=T,
               colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
  rm(x);gc()
  for( e in leaderboard){
    con_dnase_peak<-paste0(base,'/essential_training_data/DNASE/peaks/conservative/',
                           'DNASE.',e,'.conservative.narrowPeak.gz')
    con_dnase_fc<-paste0(base,'/essential_training_data/DNASE/fold_coverage_wiggles/',
                         'DNASE.',e,'.fc.signal.bigwig')
    extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")
    #load files:
    aux <- file.path(con_dnase_fc);bwf <- BigWigFile(aux)
    dnasefc <- import(bwf,format = 'bigWig')
    dnasepeak <- import(gzfile(con_dnase_peak), format = "BED", extraCols = extraCols_narrowPeak)
    
    #overlaps:
    over<-findOverlaps(grange_ladder,dnasepeak)
    x<-data.table(data.frame(qh=over@queryHits,sh=over@subjectHits,
                             mcols(dnasepeak[over@subjectHits])[-c(1,4,6)]))  # remove c('name','peak')
    #indexladder=grange_ladder[unique(over@queryHits)]$index_nona) 
    
    ##### ladder variables#########
    # qplot(pValue, signalValue, colour = bind, shape = bind, 
    #       data = x)
    # p <- ggplot(x[,.(pValue,signalValue,bind)], aes(pValue,signalValue))
    #p + geom_boxplot(aes(colour =bind))
    ######test     extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")
    ########
    x2<-x[, lapply(.SD,max), by=qh,.SDcols = -c('sh')]            # ==> toooo fast!
    x2m<-x[, lapply(.SD,mean), by=qh,.SDcols = -c('sh')] ;names(x2m)<-paste0(names(x2m),'m')
    #chip_dnase<-chip[unique(over@queryHits)] 
    
    grange_ladder_peak<-grange_ladder[unique(over@queryHits)]
    x3<-cbind(x2[,.(score,signalValue,qValue)],
              x2m[,.(scorem,signalValuem,qValuem)])
    mcols(grange_ladder_peak)<-cbind(mcols(grange_ladder_peak),as.data.frame(x3))
    rm(x3,x2m,x2,x);gc()
    
    over<-findOverlaps(grange_ladder_peak,dnasefc)
    x<-data.table(qh=over@queryHits,sh=over@subjectHits,maxfc=dnasefc[over@subjectHits]$score)
    x1<-x[, lapply(.SD,max), by=qh,.SDcols = -c('sh')]  ;
    x1m<-x[, lapply(.SD,mean), by=qh,.SDcols = -c('sh')]
    grange_ladder_peak_fc<-grange_ladder_peak;rm(grange_ladder_peak);gc()
    grange_ladder_peak_fc$maxfc<-x1[,maxfc];  grange_ladder_peak_fc$maxfcm<-x1m[,maxfc]
    featurefcpeak<-setDT(data.frame(mcols(grange_ladder_peak_fc)))
    
    #merge with motif scores
    
    #setDT(score)
    motifsc<-score[,V6]#faux(score)
    dfa<-data.table(motifsc=motifsc,index_nona=index_nona)
    feature<-merge(featurefcpeak,dfa,by='index_nona')
    save(feature,
         file=file.path(tfDir,paste0('feature_ladder_',e,'.RData')))
  }
}  
load_features<-function(tf){
  subDir <- file.path(base,'writeup','results')
  tfDir<-file.path(subDir,tf)
  aux<-read_train_leaderboard_test(tf)
  train<-aux[[1]];leaderboard<-aux[[2]];test<-aux[[3]]
  if(!identical(train, character(0))){
    faux<-function(x){
      load(file.path(tfDir,paste0('feature_train_',x,'.RData')))
      return(feature)}
    auxtrain<-lapply(train,faux); names(auxtrain)<-train
    dd<-auxtrain
  }
  if(!identical(leaderboard, character(0))){
    faux<-function(x){
      load(file.path(tfDir,paste0('feature_ladder_',x,'.RData')))
      return(feature)} 
    auxladder<-lapply(leaderboard,faux); names(auxladder)<-leaderboard
    dd<-c(dd,auxladder)
  }
  if(!identical(test, character(0))){
    faux<-function(x){
      load(file.path(tfDir,paste0('feature_test_',x,'.RData')))
      return(feature)} 
    auxtest<-lapply(test,faux);names(auxtest)<-test
    dd<-c(dd,auxtest)
  }
  #dd<-c(auxtrain,auxladder,auxtest)
  return(dd)
}
preprocess_writeup<-function(tf){
  #remove_na_save_index()
  if(!identical(test, character(0))){
    feature_tf_test(tf,test)
  }
  if(!identical(leaderboard, character(0))){
    feature_tf_ladder(tf,leaderboard)
  }
  if(!identical(train, character(0))){
    feature_tf_train(tf,train)
  }
}
load_train_and_createfolds<-function(tf,t_train,kk=20,t_test,firstn=F){
  set.seed(2423)
  library(ranger)
  e<-c(t_train,t_test);
  dd<-vector('list',length(e));folds<-vector('list',length(t_train))
  for(i in 1:length(e)){
    load(file=paste0('dflinm_allscores.',tf,'.',e[i],'.RData'))
    dflinm<-dflinm[,setdiff(colnames(dflinm),'value')]
    id<-sapply(dflinm[1,],class);ids<-id=='numeric'
    #dflinm[,c('yavg','ymax','y75','maxfc','maxfcm')]<-log10(dflinm[,c('yavg','ymax','y75','maxfc','maxfcm')])
    dflinm[,ids]<-log10(dflinm[,ids])
    xxx<-dflinm[,ids]
    aux<-xxx== -Inf
    xxx[aux]<- -10
    dflinm[,ids]<-xxx 
    dflinm[,ids]<-scale(dflinm[,ids])
    if(i  <= length(t_train)){ 
      dflinm$bind<-ifelse(dflinm$bind=='B',1,0)
    }
    dd[[i]]<-as.data.frame(dflinm)
    if(i  > length(t_train)){ 
      names(dd)[i]<-paste0(e[i],'_dna')
    }else{
      names(dd)[i]<-e[i]
    }
  }
  if(length(t_test)>0){
    for(i in 1:length(t_test)){
      load(file=paste0('dflinm_full.',tf,'.',t_test[i],'.RData'))
      dflinm[,c('yavg','ymax','y75')]<-log10(dflinm[,c('yavg','ymax','y75')])
      xxx<-dflinm[,c('yavg','ymax','y75')]
      aux<-xxx== -Inf
      xxx[aux]<- -10
      dflinm[,c('yavg','ymax','y75')]<-xxx   
      dflinm[,c('yavg','ymax','y75')]<-scale(dflinm[,c('yavg','ymax','y75')])
      dflinm<-dflinm[,setdiff(colnames(dflinm),'svalue')]
      dd<-c(dd,list(dflinm))
      names(dd)[length(dd)]<-paste0(t_test[i],'_4')
    }
  }
  for(i in 1:length(t_train)){
    dat0<-dd[[i]]#[,c(imp,'bind')]
    dfold<-createFolds(dat0$bind, k = kk, list = TRUE, returnTrain = FALSE)
    folds[[i]]<-dfold
  }
  if(length(t_tf)!=0){
    for( et in t_tf ){
      load(file=paste0('dflinm_final_dnase',tf,'.',et,'.RData'))
      dflinm<-dflinm[,setdiff(colnames(dflinm),'value')]
      id<-sapply(dflinm[1,],class);ids<-id=='numeric'
      #dflinm[,c('yavg','ymax','y75','maxfc','maxfcm')]<-log10(dflinm[,c('yavg','ymax','y75','maxfc','maxfcm')])
      dflinm[,ids]<-log10(dflinm[,ids])
      xxx<-dflinm[,ids]
      aux<-xxx== -Inf
      xxx[aux]<- -10
      dflinm[,ids]<-xxx 
      dflinm[,ids]<-scale(dflinm[,ids])
      dd<-c(dd,list(dflinm))
      names(dd)[length(dd)]<-paste0(et,'_dna')
      
      load(file=paste0('dflinm_final.',tf,'.',et,'.RData'))
      dflinm<-dflinm[,setdiff(colnames(dflinm),'svalue')]
      dflinm[,c('yavg','ymax','y75')]<-log10(dflinm[,c('yavg','ymax','y75')])
      xxx<-dflinm[,c('yavg','ymax','y75')]
      aux<-xxx== -Inf
      xxx[aux]<- -10
      dflinm[,c('yavg','ymax','y75')]<-xxx   
      dflinm[,c('yavg','ymax','y75')]<-scale(dflinm[,c('yavg','ymax','y75')])
      dd<-c(dd,list(dflinm))
      names(dd)[length(dd)]<-paste0(et,'_4')
    }
  }
  # if(firstn){
  #   dat.train<-dd[[1]][folds[[1]][[1]],];
  #   rranger<-ranger(dependent.variable.name = "bind", data = dat.train,
  #           importance='impurity',write.forest=T,num.trees = 100,mtry=ncol(dat.train)-1)
  #   imp<-ranger::importance((rranger));
  #   imp<-names(imp[order(imp,decreasing = T)])[1:firstn]
  # } else {imp<-setdiff(colnames(dd[[1]]),'bind')}
  # dd[1:(length(t_train))]<-lapply(1:(length(t_train)),function(i,x){x[[i]]<-x[[i]][,c(imp,'bind')]},x=dd)
  # if(length(t_test)>0){
  #   dd[(length(t_train)+1):length(e)]<-lapply((length(t_train)+1):length(e),
  #                 function(i,x){x[[i]]<-x[[i]][,imp]},x=dd)
  # };
  print(names(dd))
  return(list(dd=dd,folds=folds))
}
createF<-function(x,dd,kk){
  y<-x
  x<-sub('-','\\.',x)
  a<-createFolds(dd[y][[1]][,get(x)], k=kk , list = TRUE, returnTrain = FALSE)
  return(a)
}
f_combinew<-function(dd,f20,kk=20){
  dfolds<-lapply(train,createF,dd=dd,kk=kk);names(dfolds)<-train
  #dd: list of data
  #f20: interval. which folds from k=20 to get? 1:k folds
  out<-lapply(seq_along(dd),function(i){
    ids<-unlist(dfolds[[i]][f20])
    dd[[i]][ids,]})
  outx<-rbindlist(out)
}
s<-function(x){ifelse(x=='B',1,0)}
xgbtrain<-function(e){
  testXg<-f_combinew(dd[train],91:100,kk=100);setnames(testXg,train[1],'bind')  
  trainXg<-f_combinew(dd[train],1:90,kk=100);setnames(trainXg,train[1],'bind')  
  trainXg<- xgb.DMatrix(as.matrix(trainXg[,.SD,.SDcols=-c('score','index_nona','bind')])
                        ,label=as.matrix(s(trainXg[,bind])))
  testXg<- xgb.DMatrix(as.matrix(testXg[,.SD,.SDcols=-c('score','index_nona','bind')])
                       ,label=as.matrix(s(testXg[,bind])))    
  eSub<-xgb.DMatrix(as.matrix(dd[e][[1]][,.SD,.SDcols=-c('score','index_nona')]))  

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
  bst <-xgb.train(param=param, data =trainXg, nrounds=15,watchlist,early_stop_round=3)
  xgscore<-predict(bst, eSub)
  # id<-getinfo(dxg,'nrow')
  save(xgscore,file=file.path(tfDir,paste0('xgscore_',e,'.RData')))
  return(xgscore)
}
rftrain<-function(e){
  nk=100
  # trainRf<-f_combinew(dd[train],1,kk=100);setnames(trainRf,train[1],'bind')
  # trainRf<-trainRf[,.SD,.SDcols=-c('score','index_nona')]
  # trainRf[,bind:=s(bind)]
  eSub<-dd[e][[1]][,.SD,.SDcols=-c('score','index_nona')] 
  out<-function(i,mtry=2,ns=1,ntree=200){
    set.seed(12)
    trainRf<-f_combinew(dd[train],i,kk=100);setnames(trainRf,train[1],'bind')
    trainRf<-trainRf[,.SD,.SDcols=-c('score','index_nona')]
    trainRf[,bind:=s(bind)]
    outx<-ranger(dependent.variable.name = "bind",
                 data = trainRf,  importance='impurity',min.node.size=ns,
                 write.forest=T,num.trees = ntree,
                 mtry=mtry)
    print(i)
    return(outx)
  }
  ra<-lapply(1:nk,out)
  pred<-sapply(1:nk,function(i){out<-predict(ra[[i]], eSub)$predictions})
  rfscore<-apply(pred,1,mean)
  save(rfscore,file=file.path(tfDir,paste0('rfscore_',e,'.RData')))
  save(ra,file=file.path(tfDir,paste0('ra_',e,'.RData')))
  return(rfscore)
}
fcs<-function(xgscore,rfscore,e){
  if(any(e==leaderboard)){
    ladderfile<-file.path( annotationDir,'ladder_regions.blacklistfiltered.bed.gz')
    score<-fread(file.path(tfDir,paste0('ladder.',tf,'.txt')), data.table=T,
                 colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
    load(file.path(annotationDir,'ladder_nona.RData'))  # load index with no 'N's - index_nona
    xladder<-fread(paste0('gzip -dc ',ladderfile));names(xladder)<-c('chr','start','end')
    sub<-makeGRangesFromDataFrame(xladder,keep.extra.columns = F,starts.in.df.are.0based=T)
    filename=paste0('L.',tf,'.',e,'.tab')
  }
  if(any(e==test)){
    testfile<-file.path( annotationDir,'test_regions.blacklistfiltered.bed.gz')
    score<-fread(file.path(tfDir,paste0('test.',tf,'.txt')), data.table=T,
                 colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
    load(file.path(annotationDir,'test_nona.RData'))  # load index with no 'N's - index_nona
    xtest<-fread(paste0('gzip -dc ',testfile));names(xtest)<-c('chr','start','end')
    sub<-makeGRangesFromDataFrame(xtest,keep.extra.columns = F,starts.in.df.are.0based=T)
    filename=paste0('F.',tf,'.',e,'.tab')
  }
  dtsub<-data.frame(mold(sub),index=1:length(sub),final=0);setDT(dtsub)
  index_ladder<-dd[[e]][,index_nona]
  #dtsub[index %in% index_nona,final:=ff(log(score[,V6]))/10]
  dtsub[index %in% index_ladder,final:=ff(xgscore+rfscore)]
  dtsub<-dtsub[,-c(4,5,6,7),with=F]
  dtsub[,start:=start-1]

  filex<-file.path(tfDir,filename)
  #write.table(dtsub,file=filex,quote=F,col.names=F,row.names=F,sep='\t')
  fwrite(dtsub,file=filex,quote=F,col.names=F,row.names=F,sep='\t')
  if(file.exists(paste0(filex,'.gz'))){file.rename(paste0(filex,'.gz'),paste0(filex,'.gz.old'))}
  gzip(filex, destname=file.path(paste0(filex,'.gz')),skip=F)
  return(NULL)
}

################################### library and paths set###############
#download here, in the base directory these dat: annotation, DNASE, RNAseq, CHIPseq
#set directory of meme suit (ama)
tfDir<-file.path(subDir,tf)
annotationDir<-file.path(base,'annotations')
gencodev19<-file.path(base,'annotations/gencode.v19.annotation.gtf.gz')
memeAma<-file.path('~/hd/meme/bin')
library(xgboost)
library(ranger)
library(biovizBase)
library(ggbio)
library(R.utils)
library(gdata)
library(plyr)
library(GenomicRanges)
library(ShortRead)  #clean
library(rtracklayer)
library(gdata)
library(Biostrings)
library(data.table)
library(caret)
#source(file.path(writeup,"functions.R"))
library(gdata)
tfs<-read.xls(file.path(writeup,"tfs.xls"))
setkey(setDT(tfs),F.Name)
test<-strsplit(as.character(gsub('\xa0','',tfs[tf,Final.Submission.Cell.Types])),',')[[1]]
leaderboard<-strsplit(as.character(tfs[tf,Leaderboard.Cell.Types]),',')[[1]]
train<-strsplit(as.character(tfs[tf,Training.Cell.Types]),',')[[1]]

test<-sub('[[:space:]]','',test)
leaderboard<-sub('[[:space:]]','',leaderboard)
train<-sub('[[:space:]]','',train)

print(train)
print(leaderboard)
print(test)


