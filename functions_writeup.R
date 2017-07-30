
ff<-function(x){(x-min(x))/(max(x)-min(x))}
evalerror <- function(preds, dtrain) {
  library(PRROC)
  labels <- getinfo(dtrain, "label")
  pr<- pr.curve( scores.class0=preds,  weights.class0=labels)
  return(list(metric = "aupr", value = pr$auc.davis.goadrich))
}
froc<-function(gs,pr){
  library(PRROC)
  prec<-ifelse(pr=='B',1,0);gs<-ifelse(gs=='B',1,0)
  roc <- roc.curve( scores.class0=prec,  weights.class0=gs)
  pr<- pr.curve( scores.class0=prec,  weights.class0=gs)
  return(data.frame(roc=roc$auc,pr=pr$auc.davis.goadrich))
}
froc_regression<-function(gs,pr){
  library(PRROC)
  prec<-pr;gs=gs
  roc <- roc.curve( scores.class0=prec,  weights.class0=gs)
  pr<- pr.curve( scores.class0=prec,  weights.class0=gs)
  return(data.frame(roc=roc$auc,pr=pr$auc.davis.goadrich))
}
f_combine<-function(dd,f,t_train){
  out<-lapply(seq_along(t_train),function(i){
    ids<-unlist(dfolds[[i]][f])
    dd[[i]][ids,]})
  outx<-rbindlist(out)
}
preprocess<-function(tf){
  #download here, in the base directory these dat: annotation, DNASE, RNAseq, CHIPseq
  #set directory of meme suit (ama)
  base<-'/home/ricardo/hd/projects/dream_tf_competition/data'
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
  con_chipseq_label_tf<-file.path(base,'ChIPseq/labels',paste0(tf,'.train.labels.tsv.gz'))
  tfDir<-file.path(subDir,tf)
  #setwd(DIR_TF)
  chiplabelnona<-file.path(base,'annotations/chiplabel_nona.RData')
  chipladdernona<-file.path(base,'annotations/chipladder_nona.RData')
  chiptestnona<-file.path(base,'annotations/chiptest_nona.RData')
  
  gencodev19<-file.path(base,'annotations/gencode.v19.annotation.gtf.gz')
  pladder<-file.path(base,'annotations/ladder_regions.blacklistfiltered.bed.gz')
  ptest<-file.path(base,'annotations/test_regions.blacklistfiltered.bed.gz')
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
  tfs<-read.xls(file.path(writeup,"tfs.xls"))
  source(file.path(writeup,"functions.R"))
  j<-which(tfs[,1]==tf)
  t_train<-strsplit(as.character(tfs[j,2]),',')[[1]];t_train<-sub('[[:space:]]','',t_train)
  t_test<-strsplit(as.character(tfs[j,3]),',')[[1]];t_test<-sub('[[:space:]]','',t_test)
  t_tf<-strsplit(as.character(gsub('\xa0','',tfs[j,4])),',')[[1]]
  
  print(t_train)
  print(t_test)
  print(t_tf)
  make_dnase_over_train_and_test(tf,t_train,t_test)
  x<-data.table::fread(paste0('gzip -dc ',con_chipseq_label_tf))
  chip<-makeGRangesFromDataFrame(x,keep.extra.columns = T,starts.in.df.are.0based=T)
  rm(x);gc()
  nslots<-names(mcols(chip))
  nslots<-gsub('.','-',nslots,fixed=T)
  load(chiplabelnona)
  chip<-chip[value];chip$value<-value
  
  gtf<-import(gzfile(gencodev19)) #2619444
  gtf<-gtf[gtf$gene_type=='protein_coding' & gtf$transcript_status!='PUTATIVE' & 
             gtf$transcript_status!='PUTATIVE']
  pt<-promoters(gtf[gtf$type=='gene'])
  #motifs_scores(tf)
  #mscores(); ##exec bash *.sh in the tf directory
  #for i in *.sh; do bash $i;done
  mscoreswithin(tf)
  mscoresacross(tf);file.copy(from='execbash.sh',to=tfDir)
  for(e in t_train){
    mcols(chip)<-data.frame(mcols(chip))[,c('value',e)]
    extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")
    con_dnase_peak<-paste0(base,'/essential_training_data/DNASE/peaks/conservative/',
                           'DNASE.',e,'.conservative.narrowPeak.gz')
    con_dnase_fc<-paste0(base,'/essential_training_data/DNASE/fold_coverage_wiggles/',
                         'DNASE.',e,'.fc.signal.bigwig')
    aux <- file.path(con_dnase_fc);bwf <- BigWigFile(aux)
    dnasefc <- import(bwf)
    dnasepeak <- import(c1<-gzfile(con_dnase_peak), format = "BED", extraCols = extraCols_narrowPeak)
    dflinm<-fannot(chip,dnasepeak,dnasefc)
    ####################### motifscores#############
    xlab1<-fread(paste0(tfDir,'/meme_label_avg1.',e,'.txt'), data.table=FALSE,
                colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
    xlab2<-fread(paste0(tfDir,'/meme_label_avg2.',e,'.txt'), data.table=FALSE,
                colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
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
    # nn<-c('xlab1','xlab2','xlad1','xlad2','xtest1','xtest2')
    # ni<-sapply(nn,exists);nn[ni]
    y1<-faux(xlab1)
    y2<-faux(xlab2)
    load(chiplabelnona)
    mot<-data.table(y1=y1,y2=y2,value=value)
    motm<-merge(dflinm,mot,by='value')
    setDF(motm);id<-motm[,e]!='A';dflinm<-motm[id,];names(dflinm)[2]<-'bind'
    save(dflinm,file=paste0(tfDir,'/dflinm_label.',tf,'.',e,'.RData'))
    dflinm4<-data.frame(y1=y1,y2=y2,bind=data.frame(mcols(chip))[,e])
    save(dflinm4,file=paste0(tfDir,'/dflinm_label4.',tf,'.',e,'.RData'))
    rm(dflinm,mot,motm,xlab1,xlab2,dnasefc,dnasepeak);gc()
  }
  for(e in t_test){
    if(length(t_test)!=0){
      load(chiptestnona)
      x<-data.table::fread(paste0('gzip -dc ',pladder));names(x)<-c('chr','start','stop')
      chip<-makeGRangesFromDataFrame(x,keep.extra.columns = T,starts.in.df.are.0based=T)
      chip<-chip[value];chip$value<-value
      extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")
      con_dnase_peak<-paste0(base,'/essential_training_data/DNASE/peaks/conservative/',
                             'DNASE.',e,'.conservative.narrowPeak.gz')
      con_dnase_fc<-paste0(base,'/essential_training_data/DNASE/fold_coverage_wiggles/',
                           'DNASE.',e,'.fc.signal.bigwig')
      aux <- file.path(con_dnase_fc);bwf <- BigWigFile(aux)
      dnasefc <- import(bwf)
      dnasepeak <- import(c1<-gzfile(con_dnase_peak), format = "BED", extraCols = extraCols_narrowPeak)
      dflinm<-fannot(chip,dnasepeak,dnasefc)
      xlad1<-fread(paste0(tfDir,'/meme_ladder_avg1.',e,'.txt'), data.table=FALSE,
                     colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
      xlad2<-fread(paste0(tfDir,'/meme_ladder_avg2.',e,'.txt'), data.table=FALSE,
                     colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
      y1<-faux(xlad1)
      y2<-faux(xlad2)
      mot<-data.table(y1=y1,y2=y2,value=value)
      dflinm<-merge(dflinm,mot,by='value')
      save(dflinm,file=paste0(tfDir,'/dflinm_ladder.',tf,'.',e,'.RData'))
      dflinm4<-data.frame(y1=y1,y2=y2)
      save(dflinm4,file=paste0(tfDir,'/dflinm_label4.',tf,'.',e,'.RData'))
      rm(dflinm,mot,dnasefc,dnasepeak,y1,y2,xlad1,xlad2);gc()
    }
  }
  for(e in t_tf){
    if(length(t_tf)!=0){
      load(chipladdernona)
      x<-data.table::fread(paste0('gzip -dc ',pladder));names(x)<-c('chr','start','stop')
      chip<-makeGRangesFromDataFrame(x,keep.extra.columns = T,starts.in.df.are.0based=T)
      chip<-chip[value];chip$value<-value
      extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")
      con_dnase_peak<-paste0(base,'/essential_training_data/DNASE/peaks/conservative/',
                             'DNASE.',e,'.conservative.narrowPeak.gz')
      con_dnase_fc<-paste0(base,'/essential_training_data/DNASE/fold_coverage_wiggles/',
                           'DNASE.',e,'.fc.signal.bigwig')
      aux <- file.path(con_dnase_fc);bwf <- BigWigFile(aux)
      dnasefc <- import(bwf)
      dnasepeak <- import(c1<-gzfile(con_dnase_peak), format = "BED", extraCols = extraCols_narrowPeak)
      dflinm<-fannot(chip,dnasepeak,dnasefc)
      xtest1<-fread(paste0(tfDir,'/meme_test_avg1.',e,'.txt'), data.table=FALSE,
                      colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
      xtest2<-fread(paste0(tfDir,'/meme_test_avg2.',e,'.txt'), data.table=FALSE,
                      colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
      y1<-faux(xtest1)
      y2<-faux(xtest2)
      mot<-data.table(y1=y1,y2=y2,value=value)
      dflinm<-merge(dflinm,mot,by='value')
      save(dflinm,file=paste0(tfDir,'/dflinm_test.',tf,'.',e,'.RData'))
      rm(dflinm,mot,motm,xlab1,xlab2,dnasefc,dnasepeak,y1,y2,xtest1,xtest2);gc()
    }
  }
  return(NULL)
} 
fannot<-function(chip_dnase,tissue){
  ###################annotations inclusion########################################
  gtf<-import(gzfile(paste0(base,'/annotations/gencode.v19.annotation.gtf.gz'))) #2619444
  gtf<-gtf[gtf$gene_type=='protein_coding' & gtf$transcript_status!='PUTATIVE' & 
             gtf$transcript_status!='PUTATIVE' & gtf$type=='gene']
  #pt<-promoters(gtf[gtf$type=='gene'])
  ####################### genes: dist, after and before#######
  fo<-follow(chip_dnase, gtf);   # return index in gtf that follow in chip_dnase
  pr<-precede(chip_dnase, gtf);
  dfp<-data.table(fo_gtf=fo,pr_gtf=pr,id_dnase= 1:length(fo))
  ################test#######
  table(is.na(fo)) # FALSE   TRUE  870793    848 
  table(is.na(pr))  # FALSE   TRUE 870669    972 
  
  setkey(dfp,fo_gtf);dfp[.(NA), fo_gtf := 0L]; 
  setkey(dfp,pr_gtf); dfp[.(NA), pr_gtf := 0L]
  #faux<-function(x){ifelse(x==0,0,gtf[x]$gene_name)}
  faux<-function(x){
    ifelse(x==0,0,gtf[x]$gene_name)
  }
  dfp[,fon:=faux(fo_gtf)]
  dfp[,prn:=faux(pr_gtf)]
  dnearest<-distanceToNearest(chip_dnase,gtf)
  ffonpr<-function(x){
    folowDistance<-distanceToNearest(chip_dnase[i],gtf[dfp[i,fo_gtf]])
    preceeDistance<-distanceToNearest(chip_dnase,gtf)
  }
  ##################the fon and per is problematic and meaningless.
  dfp[,dnearest:=x@elementMetadata$distance];#dfp<-dfp[order(id_dnase)]
  #load(paste0(base,'/RNAseq/',tissue,'.RData'))
  rna<-fread(paste0(base,'/RNAseq/','gene_expression.',tissue,'.biorep1.tsv'))
  rna[,"transcript_id(s)":=NULL]
  table_conv<-as.data.table(mcols(gtf)[c('gene_id','gene_name')])
  table_conv[,gene_id:=unique(gene_id)] #resolve for multiples gene_ids
  rna<-merge(rna, table_conv,by='gene_id')
  #xx<-as.data.table(mcols(genes_gtf))
  rna_dfp<-merge(dfp,rna,by.x='fon',by.y='gene_name')
  xx1<-merge(dfp,xx[,c(2:7),with=F],by.x='fon',by.y='gene_name',all.x=T,all.y=F,sort=F, allow.cartesian=T)
  xx1<-xx1[!duplicated(id)]
  xx2<-merge(xx1,xx[,c(2:7),with=F],by.x='prn',by.y='gene_name',all.x=T,all.y=F,sort=F, allow.cartesian=T)
  xx2<-xx2[!duplicated(id)]
  xx2[,paste0(names(xx[,c(3:7),with=F]),'.y') :=
        lapply(.SD,function(x){ifelse(is.na(x),0,x)}),.SDcols=paste0(names(xx[,c(3:7),with=F]),'.y')]
  xx2[,c('fon','prn','pr','fo'):=list(NULL,NULL,NULL,NULL)]
  a<-distanceToNearest(chip_dnase, gtf )
  ##########
  fprep<-function(gtf){
    over<-findOverlaps(gtf,chip_dnase)
    x<-data.table(value=mcols(chip_dnase[over@subjectHits])$value,type=gtf[over@queryHits]$type)
    xc <- dcast(x, value ~type,fun.aggregate = length)
    return(xc)
  }
  z<-fprep(gtf);if('Selenocysteine' %in% names(z)){z[,Selenocysteine:= NULL]}
  z1<-fprep(pt);names(z1)<-c('value','promoter2k')
  db<-data.table(as.data.frame(mcols(chip_dnase)))
  xdb<-merge(db,z,by='value',all.x=T,all.y=F)
  xdb[, setdiff(names(z),'value') :=lapply(.SD, function(x){ ifelse(is.na(x),0,x)}),.SDcols=setdiff(names(z),'value') ]
  xdb<-merge(xdb,z1,by='value',all.x=T,all.y=F)
  xdb[, c('promoter2k') :=sapply(promoter2k, function(x){ ifelse(is.na(x),0,x)})]
  return(data.frame(xdb,dist=a@elementMetadata$distance,xx2,stringsAsFactors = F))
} ## finish later 24/08

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
run_execbash_sh_on_tf_folder_after_this<-function(tf){
  tfs<-read.xls(file.path(writeup,"tfs.xls"))
  setkey(setDT(tfs),F.Name)
  test<-strsplit(as.character(gsub('\xa0','',tfs[tf,Final.Submission.Cell.Types])),',')[[1]]
  leaderboard<-strsplit(as.character(tfs[tf,Leaderboard.Cell.Types]),',')[[1]]
  train<-strsplit(as.character(tfs[tf,Training.Cell.Types]),',')[[1]]
  
  test<-sub('[[:space:]]','',test);#test<-sub('-','.',test)
  leaderboard<-sub('[[:space:]]','',leaderboard);#leaderboard<-sub('-','.',leaderboard)
  train<-sub('[[:space:]]','',train);  #train<-sub('-','.',train)

  
  # print(train)
  # print(leaderboard)
  # print(test)
  # 
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
    score<-fread(paste0(tfDir,'/test.',tf,'.txt'), data.table=T,
                 colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
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
  for( e in train){
    grange_train_labels_loop<-grange_train_labels[,c(e,'index_nona')]
    #remove ambiguous
    e<-data.table(bind=mcols(grange_train_labels_loop)[[1]]);idnoA<-e[,bind]!='A'
    grange_train_labels_loop<-grange_train_labels_loop[idnoA]
    
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
    score<-fread(file.path(tfDir,paste0('train.',tf,'.txt')), data.table=T,
                 colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
    #setDT(score)
    motifsc<-score[,V6]#faux(score)
    dfa<-data.table(motifsc=motifsc,index_nona=index_nona)
    feature<-merge(featurefcpeak,dfa,by='index_nona')
    
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
    score<-fread(paste0(tfDir,'/ladder.',tf,'.txt'), data.table=T,
                 colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
    #setDT(score)
    motifsc<-score[,V6]#faux(score)
    dfa<-data.table(motifsc=motifsc,index_nona=index_nona)
    feature<-merge(featurefcpeak,dfa,by='index_nona')
    save(feature,
         file=file.path(tfDir,paste0('feature_ladder_',e,'.RData')))
  }
}  


# execute execbash.sh on the tf directory
tf<-'MAX'
################################### library and paths set###############
#download here, in the base directory these dat: annotation, DNASE, RNAseq, CHIPseq
#set directory of meme suit (ama)
base<-'/home/ricardo/hd/projects/dream_tf_competition/data'
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

################################### library and paths set###############
run_execbash_sh_on_tf_folder_after_this(tf)
# execute execbash.sh on the tf directory

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
  
  