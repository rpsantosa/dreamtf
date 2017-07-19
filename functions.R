
prepare_test_final<-function(){
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(GenomicRanges)
  library(data.table)
  library(DNAshapeR)
  annot<-'/home/ricardo/hd/projects/dream_tf_competition/data/annotations/'
  p<-'/home/ricardo/hd/projects/dream_tf_competition/data/ChIPseq/labels/ARID3A.train.labels.tsv.gz'
  xlabel<-fread(paste0('gzip -dc ',p))
  con_final<-'/home/ricardo/hd/projects/dream_tf_competition/data/annotations/test_regions.blacklistfiltered.bed.gz'
  xtest<-data.table::fread(paste0('gzip -dc ',con_final));names(xtest)<-c('chr','start','end')
  con_submission<-'/home/ricardo/hd/projects/dream_tf_competition/data/annotations/t/ladder_regions.blacklistfiltered.bed.gz'
  xladder<-data.table::fread(paste0('gzip -dc ',con_submission));names(xladder)<-c('chr','start','end')
  fremovena<-function(x){
    chip<-makeGRangesFromDataFrame(x,keep.extra.columns = F,starts.in.df.are.0based=T)
    chip$value=1:nrow(x)
    xb<-getSeq(Hsapiens,chip)
    db<-data.table(oligonucleotideFrequency(xb, width=1));
    if(ncol(x)>3){db<-cbind(db,x[,-c(1,2,3),with=F])}
    h<-db[,A!=0 & T!=0 & C!=0 & G!=0]
    return(chip[h])
  }
  
  chipgt<-fremovena(xtest);value<-chipgt$value
  getFasta(chipgt,	Hsapiens,	width	= 200,	filename	= paste0(annot,'chiptest_nona.fa'))
  save(value,file='/home/ricardo/hd/projects/dream_tf_competition/data/annotations/chiptest_nona.RData')
  rm(xtest,chipgt);gc()
  
  chipgtl<-fremovena(xladder);value<-chipgtl$value
  getFasta(chipgtl,	Hsapiens,	width	= 200,	filename	= paste0(annot,'chipladder_nona.fa'))
  save(value,file='/home/ricardo/hd/projects/dream_tf_competition/data/annotations/chipladder_nona.RData')
  rm(xladder,chipgtl);gc()
  
  chiplabel<-fremovena(xlabel);value<-chiplabel$value
  getFasta(chiplabel,	Hsapiens,	width	= 200,	filename	= paste0(annot,'chiplabel_nona.fa'))
  save(value,file='/home/ricardo/hd/projects/dream_tf_competition/data/annotations/chiplabel_nona.RData')
  rm(xlabel,chiplabel);gc()
  
}
make_dnase_over_train_and_test<-function(tf,t_train,t_test){
  for(e in t_train){
    print(e)
  dnase_train(tf,e)
  }
  if(length(t_test)!=0){
    for(e in t_test){
      print(e)
      dnase_test(tf,e)
    }
  }
}
make_motifs_scores_train_and_test<-function(tf,t_train,t_test){
  for(e in t_train){
    make_scores(tf,e,trainn=T)
  }
  if(length(t_test)!=0){
    for(e in t_test){
      make_scores(tf,e,trainn=F)
      make_scores_full(tf,e)
    }
  }
}
make_scores<-function(tf,e,trainn){ 
  xavg<-fread(paste0('meme_avg.',tf,'.',e,'.txt'), data.table=FALSE,
                   colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
  xmax<-fread(paste0('meme_plus.',tf,'.',e,'.txt'), data.table=FALSE,
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
  yavg<-faux(xavg)
  ymax<-faux(xmax)
  y75<-(ymax*yavg)
 
  load(paste0('chip_dnase.',tf,'.',e,'.RData'))
   aux<- data.frame(mcols(chip_dnase))
  dflinm<-data.frame(yavg,ymax,y75,aux)
  id_bind<-which(dflinm[1,]=='U' | dflinm[1,]=='B' | dflinm[1,]=='A')
  id_e<-which(colnames(dflinm)==gsub('-','.',e))
  if(trainn==T){
    ids<-setdiff(1:ncol(dflinm),id_bind);idsv<-setdiff(ids,max(id_bind)+1)
    dflinm<-dflinm[,c(max(id_bind)+1,idsv,id_e)];names(dflinm)[ncol(dflinm)]<-'bind'
    # value          V6     V6.1 score signalValue   qValue    maxfc bind
    # 1  2708 0.800076000 1.278280  1000    16.74641 63.55915  3.41085    U
    # 2  2709 0.589210000 1.195440  1000    16.74641 63.55915 11.16279    U
  } else {
    ids<-setdiff(colnames(dflinm),'value')
    dflinm<-dflinm[,c('value',ids)]
  }
  save(dflinm,file=paste0('dflinm_allscores.',tf,'.',e,'.RData'))
  return(NULL)
}
make_scores_full<-function(tf,e){ #e=tissue_test 
   con_submission<-'/home/ricardo/hd/projects/dream_tf_competition/data/annotations/t/ladder_regions.blacklistfiltered.bed.gz'
   ####dist genes
   x<-data.table::fread(paste0('gzip -dc ',con_submission));names(x)<-c('chr','start','end')
   chip<-makeGRangesFromDataFrame(x,keep.extra.columns = T,starts.in.df.are.0based=T)
   gtf<-import(gzfile('/home/ricardo/hd/projects/dream_tf_competition/data/annotations/gencode.v19.annotation.gtf.gz')) #2619444
   gtf<-gtf[gtf$gene_type=='protein_coding' & gtf$transcript_status!='PUTATIVE' & 
              gtf$transcript_status!='PUTATIVE']
   pt<-promoters(gtf[gtf$type=='gene'])
   fo<-follow(chip, gtf[gtf$type=='gene'] );
   pr<-precede(chip, gtf[gtf$type=='gene']);
   dfp<-data.table(fo=fo,pr=pr);dfp[,id:=1:nrow(dfp)]
   setkey(dfp,fo);dfp[.(NA), fo := 0L]; 
   setkey(dfp,pr); dfp[.(NA), pr := 0L]
   faux<-function(x){ifelse(x==0,0,gtf[gtf$type=='gene'][x]$gene_name)}
   dfp[,fon:=faux(fo)]
   dfp[,prn:=faux(pr)];dfp<-dfp[order(id)]
   load(paste0('/home/ricardo/hd/projects/dream_tf_competition/data/RNAseq/',e,'.RData'))
   xx<-as.data.table(mcols(genes_gtf))
   xx1<-merge(dfp,xx[,c(2:7),with=F],by.x='fon',by.y='gene_name',all.x=T,all.y=F,sort=F, allow.cartesian=T)
   xx1<-xx1[!duplicated(id)]
   xx2<-merge(xx1,xx[,c(2:7),with=F],by.x='prn',by.y='gene_name',all.x=T,all.y=F,sort=F, allow.cartesian=T)
   #xx2<-xx2[!duplicated(id)][,c('TPM.x','TPM.y'),with=F];names(xx2)<-c('follow','precede')
   xx2<-xx2[!duplicated(id)]
   xx2[,paste0(names(xx[,c(3:7),with=F]),'.y') :=
         lapply(.SD,function(x){ifelse(is.na(x),0,x)}),.SDcols=paste0(names(xx[,c(3:7),with=F]),'.y')]
   xx2[,c('fon','prn','pr','fo'):=list(NULL,NULL,NULL,NULL)]
   xx2[,paste0(names(xx[,c(3:7),with=F]),'.x') :=
         lapply(.SD,function(x){ifelse(is.na(x),0,x)}),.SDcols=paste0(names(xx[,c(3:7),with=F]),'.x')]
   
   a<-distanceToNearest(chip, gtf[gtf$type=='gene'] )
   xx2[,id:=NULL]
   mcols(chip)<-data.frame(value=1:length(chip),xx2)
   ##########
   fprep<-function(gtf){
     over<-findOverlaps(gtf,chip)
     x<-data.table(value=mcols(chip[over@subjectHits])$value,type=gtf[over@queryHits]$type)
     xc <- dcast(x, value ~type,fun.aggregate = length)
     return(xc)
   }
   z<-fprep(gtf);z[,Selenocysteine:= NULL]
   z1<-fprep(pt);names(z1)<-c('value','promoter2k')
   db<-data.table(as.data.frame(mcols(chip)))
   xdb<-merge(db,z,by='value',all.x=T,all.y=F)
   xdb[, setdiff(names(z),'value') :=lapply(.SD, function(x){ ifelse(is.na(x),0,x)}),.SDcols=setdiff(names(z),'value') ]
   xdb<-merge(xdb,z1,by='value',all.x=T,all.y=F)
   xdb[, c('promoter2k') :=sapply(promoter2k, function(x){ ifelse(is.na(x),0,x)})]
   #mcols(chip)<-data.frame(xdb,dist=a@elementMetadata$distance,xx2,stringsAsFactors = F)
   rm(aux,db,x1m,x2,z,gtf,pt,xx2,z,z1,dfp);gc()
  
   ####
    xavg<-fread(paste0('meme_full_avg.',e,'.txt'), data.table=FALSE,
                     colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
    xmax<-fread(paste0('meme_full_plus.',e,'.txt'), data.table=FALSE,
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
    yavg<-faux(xavg)
    ymax<-faux(xmax)
    y75<-(ymax*yavg)
    dflinm<-data.frame(yavg,ymax,y75,dist=a@elementMetadata$distance,xdb)
    save(dflinm,file=paste0('dflinm_full.',tf,'.',e,'.RData'))
    return(NULL)
}
make_scores_final<-function(tf,e){
  #con_submission<-'/home/ricardo/hd/projects/dream_tf_competition/data/annotations/t/ladder_regions.blacklistfiltered.bed.gz'
  con_final<-'/home/ricardo/hd/projects/dream_tf_competition/data/annotations/test_regions.blacklistfiltered.bed.gz'
  ####dist genes & annotations
  x<-data.table::fread(paste0('gzip -dc ',con_final));names(x)<-c('chr','start','end')
  chip<-makeGRangesFromDataFrame(x,keep.extra.columns = T,starts.in.df.are.0based=T)
  gtf<-import(gzfile('/home/ricardo/hd/projects/dream_tf_competition/data/annotations/gencode.v19.annotation.gtf.gz')) #2619444
  gtf<-gtf[gtf$gene_type=='protein_coding' & gtf$transcript_status!='PUTATIVE' & 
             gtf$transcript_status!='PUTATIVE']
  fo<-follow(chip, gtf[gtf$type=='gene'] );
  pr<-precede(chip, gtf[gtf$type=='gene']);
  dfp<-data.table(fo=fo,pr=pr);dfp[,id:=1:nrow(dfp)]
  setkey(dfp,fo);dfp[.(NA), fo := 0L]; 
  setkey(dfp,pr); dfp[.(NA), pr := 0L]
  faux<-function(x){ifelse(x==0,0,gtf[gtf$type=='gene'][x]$gene_name)}
  dfp[,fon:=faux(fo)]
  dfp[,prn:=faux(pr)];dfp<-dfp[order(id)]
  load(paste0('/home/ricardo/hd/projects/dream_tf_competition/data/RNAseq/',e,'.RData'))
  xx<-as.data.table(mcols(genes_gtf))
  xx1<-merge(dfp,xx[,c(2:7),with=F],by.x='fon',by.y='gene_name',all.x=T,all.y=F,sort=F, allow.cartesian=T)
  xx1<-xx1[!duplicated(id)]
  xx2<-merge(xx1,xx[,c(2:7),with=F],by.x='prn',by.y='gene_name',all.x=T,all.y=F,sort=F, allow.cartesian=T)
  #xx2<-xx2[!duplicated(id)][,c('TPM.x','TPM.y'),with=F];names(xx2)<-c('follow','precede')
  xx2<-xx2[!duplicated(id)]
  xx2[,paste0(names(xx[,c(3:7),with=F]),'.y') :=
        lapply(.SD,function(x){ifelse(is.na(x),0,x)}),.SDcols=paste0(names(xx[,c(3:7),with=F]),'.y')]
  xx2[,c('fon','prn','pr','fo'):=list(NULL,NULL,NULL,NULL)]
  xx2[,paste0(names(xx[,c(3:7),with=F]),'.x') :=
        lapply(.SD,function(x){ifelse(is.na(x),0,x)}),.SDcols=paste0(names(xx[,c(3:7),with=F]),'.x')]
  
  a<-distanceToNearest(chip, gtf[gtf$type=='gene'] )
  mcols(chip)<-data.frame(value=1:length(chip),xx2)
  ##########
  fprep<-function(gtf){
    over<-findOverlaps(gtf,chip)
    x<-data.table(value=mcols(chip[over@subjectHits])$value,type=gtf[over@queryHits]$type)
    xc <- dcast(x, value ~type,fun.aggregate = length)
    return(xc)
  }
  z<-fprep(gtf);z[,Selenocysteine:= NULL]
  z1<-fprep(pt);names(z1)<-c('value','promoter2k')
  db<-data.table(as.data.frame(mcols(chip)))
  xdb<-merge(db,z,by='value',all.x=T,all.y=F)
  xdb[, setdiff(names(z),'value') :=lapply(.SD, function(x){ ifelse(is.na(x),0,x)}),.SDcols=setdiff(names(z),'value') ]
  xdb<-merge(xdb,z1,by='value',all.x=T,all.y=F)
  xdb[, c('promoter2k') :=sapply(promoter2k, function(x){ ifelse(is.na(x),0,x)})]
  xdb[,'id' := NULL]
  #mcols(chip)<-data.frame(xdb,dist=a@elementMetadata$distance,xx2,stringsAsFactors = F)
  rm(aux,db,x1m,x2,z,gtf,pt,xx2,z,z1,dfp);gc()
  ####
  xavg<-fread(paste0('meme_finalsub_avg.',tf,'.',e,'.txt'), data.table=FALSE,
                   colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
  xmax<-fread(paste0('meme_finalsub_plus.',tf,'.',e,'.txt'), data.table=FALSE,
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
  yavg<-faux(xavg)
  ymax<-faux(xmax)
  y75<-(ymax*yavg)
  dflinm<-data.frame(svalue=1:length(yavg),yavg,ymax,y75,dist=a@elementMetadata$distance,xdb)
  save(dflinm,file=paste0('dflinm_final.',tf,'.',e,'.RData'))
  return(NULL)
}
make_scores_final_dnase<-function(tf,e){
  xavg<-fread(paste0('meme_finalsub_dnase_avg.',tf,'.',e,'.txt'), data.table=FALSE,
                   colClasses=c("character",'NULL','NULL','NULL','NULL','numeric','NULL','NULL','NULL','character'))
  xmax<-fread(paste0('meme_finalsub_dnase_plus.',tf,'.',e,'.txt'), data.table=FALSE,
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
  yavg<-faux(xavg)
  ymax<-faux(xmax)
  y75<-(ymax*yavg)
  load(paste0('final_dnase.',tf,'.',e,'.RData'))
  aux<- data.frame(mcols(final_dnase))
  dflinm<-data.frame(yavg,ymax,y75,aux)
  ids<-setdiff(colnames(dflinm),'value')
  dflinm<-dflinm[,c('value',ids)]
  save(dflinm,file=paste0('dflinm_final_dnase',tf,'.',e,'.RData'))
  return(NULL)
}

dnase_train<-function(tf,tissue_train){
  # DIR_TF<-paste0(con_results,tf,'/')
  # #con_chipseq_label_tf<-paste0('/home/ricardo/hd/projects/dream_tf_competition/data/ChIPseq/labels/',tf,'.train.labels.tsv.gz.1.gz')
  # con_chipseq_label_tf<-paste0('/home/ricardo/hd/projects/dream_tf_competition/data/ChIPseq/labels/',tf,'.train.labels.tsv.gz')
  # #nslots<-strsplit(readLines(con<-gzfile(con_chipseq_labelno1),n=1)[[1]],'\t')[[1]][-c(1,2,3)]
  x<-data.table::fread(paste0('gzip -dc ',con_chipseq_label_tf))
  chip<-makeGRangesFromDataFrame(x,keep.extra.columns = T,starts.in.df.are.0based=T)
  rm(x);gc()
  nslots<-names(mcols(chip))
  nslots<-gsub('.','-',nslots,fixed=T)
  #extraCols_bind<-rep('character',length(nslots));names(extraCols_bind)<-nslots
  extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")
  tissue_train<-t_train
  con_dnase_peak_train<-paste0(base,'/essential_training_data/DNASE/peaks/conservative/',
                               'DNASE.',tissue_train,'.conservative.narrowPeak.gz')
  con_dnase_fc_train<-paste0(base,'/essential_training_data/DNASE/fold_coverage_wiggles/',
                             'DNASE.',tissue_train,'.fc.signal.bigwig')
  aux <- file.path(con_dnase_fc_train);bwf <- BigWigFile(aux)
  dnasefc <- import(bwf)
  dnasepeak <- import(c1<-gzfile(con_dnase_peak_train), format = "BED", extraCols = extraCols_narrowPeak)
  #chip<-import(con<-gzfile(con_chipseq_label_tf),  format = "BED",extraCols = extraCols_bind) 
  mcols(chip)<-data.frame(mcols(chip),value=1:length(chip)) # to indices purposes
  #on.exit(close(aux))
  
  over<-findOverlaps(chip,dnasepeak)
  chip_dnase<-chip[unique(over@queryHits)] 
  x<-data.frame(qh=over@queryHits,sh=over@subjectHits,mcols(dnasepeak[over@subjectHits]));x<-x[,-c(3,8)]
  x<-as.data.table(x)
  x2<-x[, lapply(.SD,max), by=qh]            # ==> toooo fast!
  x2m<-x[, lapply(.SD,mean), by=qh];names(x2m)<-paste0(names(x2m),'m')
  x3<-cbind(subset(x2,select=c(4,5,6)),subset(x2m,select=c(4,5,6)))
  mcols(chip_dnase)<-cbind(mcols(chip_dnase),as.data.frame(x3))
  rm(chip,dnasepeak,x3,x2m,x);gc()
  
  over<-findOverlaps(chip_dnase,dnasefc)
  x<-data.frame(qh=over@queryHits,sh=over@subjectHits,maxfc=dnasefc[over@subjectHits]$score)
  x<-as.data.table(x)
  x1<-x[, lapply(.SD,max), by=qh]  ;x1m<-x[, lapply(.SD,mean), by=qh]
  chip_dnase$maxfc<-x1[,maxfc];  chip_dnase$maxfcm<-x1m[,maxfc]
  id<-which(nslots == tissue_train)
  #filter only B and Us to speed process
  da<-mcols(chip_dnase)[id][[1]]=='A';chip_dnase<-chip_dnase[!da]
  rm(dnasefc,x,x1,over);gc()
  
  ###################annotations inclusion########################################
  gtf<-import(gzfile(file.path(base,'/annotations/gencode.v19.annotation.gtf.gz'))) #2619444
  gtf<-gtf[gtf$gene_type=='protein_coding' & gtf$transcript_status!='PUTATIVE' & 
             gtf$transcript_status!='PUTATIVE']
  pt<-promoters(gtf[gtf$type=='gene'])
  ####dist genes
  fo<-follow(chip_dnase, gtf[gtf$type=='gene'] );
  pr<-precede(chip_dnase, gtf[gtf$type=='gene']);
  dfp<-data.table(fo=fo,pr=pr);dfp[,id:=1:nrow(dfp)]
  setkey(dfp,fo);dfp[.(NA), fo := 0L]; 
  setkey(dfp,pr); dfp[.(NA), pr := 0L]
  faux<-function(x){ifelse(x==0,0,gtf[gtf$type=='gene'][x]$gene_name)}
  dfp[,fon:=faux(fo)]
  dfp[,prn:=faux(pr)];dfp<-dfp[order(id)]
  load(paste0(base,'/RNAseq/',tissue_train,'.RData'))
  xx<-as.data.table(mcols(genes_gtf))
  xx1<-merge(dfp,xx[,c(2:7),with=F],by.x='fon',by.y='gene_name',all.x=T,all.y=F,sort=F, allow.cartesian=T)
  xx1<-xx1[!duplicated(id)]
  xx2<-merge(xx1,xx[,c(2:7),with=F],by.x='prn',by.y='gene_name',all.x=T,all.y=F,sort=F, allow.cartesian=T)
  xx2<-xx2[!duplicated(id)]
  xx2[,paste0(names(xx[,c(3:7),with=F]),'.y') :=
        lapply(.SD,function(x){ifelse(is.na(x),0,x)}),.SDcols=paste0(names(xx[,c(3:7),with=F]),'.y')]
    xx2[,c('fon','prn','pr','fo'):=list(NULL,NULL,NULL,NULL)]
  a<-distanceToNearest(chip_dnase, gtf[gtf$type=='gene'] )
  ##########
  fprep<-function(gtf){
    over<-findOverlaps(gtf,chip_dnase)
    x<-data.table(value=mcols(chip_dnase[over@subjectHits])$value,type=gtf[over@queryHits]$type)
    xc <- setDT(dcast(x, value ~type,fun.aggregate = length))
    return(xc)
  }
  z<-fprep(gtf);if('Selenocysteine' %in% names(z)){z[,Selenocysteine:= NULL]}
  z1<-fprep(pt);names(z1)<-c('value','promoter2k')
  db<-data.table(as.data.frame(mcols(chip_dnase)))
  xdb<-merge(db,z,by='value',all.x=T,all.y=F)
  xdb[, setdiff(names(z),'value') :=lapply(.SD, function(x){ ifelse(is.na(x),0,x)}),.SDcols=setdiff(names(z),'value') ]
  xdb<-merge(xdb,z1,by='value',all.x=T,all.y=F)
  xdb[, c('promoter2k') :=sapply(promoter2k, function(x){ ifelse(is.na(x),0,x)})]
  mcols(chip_dnase)<-data.frame(xdb,dist=a@elementMetadata$distance,xx2,stringsAsFactors = F)
  rm(aux,db,x1m,x2,xdb,z,gtf,pt,xx2,z1,a,dfp);gc()
  #####################end anotations inclusion
  save(chip_dnase,file=paste0(tfDir,'/chip_dnase','.',tf,'.',tissue_train,'.RData'))
  #export.bed(chip_dnase,con =file(paste0(DIR_TF,'chip_dnase','.',tf,'.',tissue_train,'.bed')))
  require(BSgenome.Hsapiens.UCSC.hg19);require(DNAshapeR)
  getFasta(chip_dnase,	Hsapiens,	width	= 200,	filename	= paste0(tfDir,'/chip_dnase','.',tf,'.',tissue_train,'.fa'))
  return(NULL)
}
dnase_test<-function(tf,tissue_test){
  #DIR_TF<-paste0(con_results,tf,'/')
  con_submission<-paste0(base,'/annotations/t/ladder_regions.blacklistfiltered.bed')
  con_dnase_peak_test<-paste0(base,'/essential_training_data/DNASE/peaks/conservative/',
                              'DNASE.',tissue_test,'.conservative.narrowPeak.gz')
  con_dnase_fc_test<-paste0(base,'/essential_training_data/DNASE/fold_coverage_wiggles/',
                            'DNASE.',tissue_test,'.fc.signal.bigwig')
  chip<-import(con<-gzfile(con_submission),  format = "BED") 

  mcols(chip)<-data.frame(mcols(chip),value=1:length(chip)) # to indices purposes
  extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")
  dnasepeak <- import(gzfile(con_dnase_peak_test),
                      format = "BED", extraCols = extraCols_narrowPeak)
  aux <- file.path(con_dnase_fc_test);bwf <- BigWigFile(aux)
  dnasefc <- import(bwf)

  over<-findOverlaps(chip,dnasepeak)
  chip_dnase<-chip[unique(over@queryHits)] 
  x<-data.frame(qh=over@queryHits,sh=over@subjectHits,mcols(dnasepeak[over@subjectHits]));x<-x[,-c(3,8)]
  setDT(x)
  x2<-x[, lapply(.SD,max), by=qh]            # ==> toooo fast!
  x2m<-x[, lapply(.SD,mean), by=qh];names(x2m)<-paste0(names(x2m),'m')
  x3<-cbind(subset(x2,select=c(4,5,6)),subset(x2m,select=c(4,5,6)))
  mcols(chip_dnase)<-cbind(mcols(chip_dnase),as.data.frame(x3))
  rm(chip,dnasepeak,x3,x2m,x);gc()
  
  over<-findOverlaps(chip_dnase,dnasefc)
  x<-data.frame(qh=over@queryHits,sh=over@subjectHits,maxfc=dnasefc[over@subjectHits]$score)
  setDT(x)
  x1<-x[, lapply(.SD,max), by=qh]  ;x1m<-x[, lapply(.SD,mean), by=qh]
  chip_dnase$maxfc<-x1[,maxfc];  chip_dnase$maxfcm<-x1m[,maxfc]
  rm(dnasefc,x,x1,over);gc()
  
  ###################annotations inclusion########################################
  gtf<-import(gzfile(file.path(base,'/annotations/gencode.v19.annotation.gtf.gz'))) #2619444
  gtf<-gtf[gtf$gene_type=='protein_coding' & gtf$transcript_status!='PUTATIVE' & 
             gtf$transcript_status!='PUTATIVE']
  pt<-promoters(gtf[gtf$type=='gene'])
  ####dist genes
  fo<-follow(chip_dnase, gtf[gtf$type=='gene'] );
  pr<-precede(chip_dnase, gtf[gtf$type=='gene']);
  dfp<-data.table(fo=fo,pr=pr);dfp[,id:=1:nrow(dfp)]
  setkey(dfp,fo);dfp[.(NA), fo := 0L]; 
  setkey(dfp,pr); dfp[.(NA), pr := 0L]
  faux<-function(x){ifelse(x==0,0,gtf[gtf$type=='gene'][x]$gene_name)}
  dfp[,fon:=faux(fo)]
  dfp[,prn:=faux(pr)];dfp<-dfp[order(id)]
  load(paste0(base,'/RNAseq/',tissue_test,'.RData'))
  xx<-as.data.table(mcols(genes_gtf))
  xx1<-merge(dfp,xx[,c(2:7),with=F],by.x='fon',by.y='gene_name',all.x=T,all.y=F,sort=F, allow.cartesian=T)
  xx1<-xx1[!duplicated(id)]
  xx2<-merge(xx1,xx[,c(2:7),with=F],by.x='prn',by.y='gene_name',all.x=T,all.y=F,sort=F, allow.cartesian=T)
  #xx2<-xx2[!duplicated(id)][,c('TPM.x','TPM.y'),with=F];names(xx2)<-c('follow','precede')
  xx2<-xx2[!duplicated(id)]
  xx2[,paste0(names(xx[,c(3:7),with=F]),'.y') :=
        lapply(.SD,function(x){ifelse(is.na(x),0,x)}),.SDcols=paste0(names(xx[,c(3:7),with=F]),'.y')]
    xx2[,c('fon','prn','pr','fo'):=list(NULL,NULL,NULL,NULL)]
  a<-distanceToNearest(chip_dnase, gtf[gtf$type=='gene'] )
  ##########
  fprep<-function(gtf){
    over<-findOverlaps(gtf,chip_dnase)
    x<-data.table(value=mcols(chip_dnase[over@subjectHits])$value,type=gtf[over@queryHits]$type)
    xc <- setDT(dcast(x, value ~type,fun.aggregate = length))
    return(xc)
  }
  z<-fprep(gtf);if('Selenocysteine' %in% names(z)){z[,Selenocysteine:= NULL]}
  z1<-fprep(pt);names(z1)<-c('value','promoter2k')
  db<-data.table(as.data.frame(mcols(chip_dnase)))
  xdb<-merge(db,z,by='value',all.x=T,all.y=F)
  xdb[, setdiff(names(z),'value') :=lapply(.SD, function(x){ ifelse(is.na(x),0,x)}),.SDcols=setdiff(names(z),'value') ]
  xdb<-merge(xdb,z1,by='value',all.x=T,all.y=F)
  xdb[, c('promoter2k') :=sapply(promoter2k, function(x){ ifelse(is.na(x),0,x)})]
  mcols(chip_dnase)<-data.frame(xdb,dist=a@elementMetadata$distance,xx2,stringsAsFactors = F)
  rm(aux,db,x1m,x2,xdb,z,gtf,pt,xx2,z1,a,dfp);gc()
  #####################end anotations inclusion

  save(chip_dnase,file=paste0(tfDir,'/chip_dnase','.',tf,'.',tissue_test,'.RData'))
  #export.bed(chip_dnase,con =file(paste0(DIR_TF,'chip_dnase','.',tf,'.',tissue_test,'.bed')))
  require(BSgenome.Hsapiens.UCSC.hg19);require(DNAshapeR)
  getFasta(chip_dnase,	Hsapiens,	width	= 200,	filename	= paste0(tfDir,'/chip_dnase','.',tf,'.',tissue_test,'.fa'))
  return(NULL)
}
dnase_final<-function(tf,e){ #e=tissue_subm
  con_final<-'/home/ricardo/hd/projects/dream_tf_competition/data/annotations/test_regions.blacklistfiltered.bed.gz'
  DIR_TF<-paste0('/home/ricardo/hd/projects/dream_tf_competition/data/TF2/',tf,'/')
  con_dnase_peak_test<-paste0('/home/ricardo/hd/projects/dream_tf_competition/data/essential_training_data/DNASE/peaks/conservative/',
                              'DNASE.',e,'.conservative.narrowPeak.gz')
  con_dnase_fc_test<-paste0('/home/ricardo/hd/projects/dream_tf_competition/data/essential_training_data/DNASE/fold_coverage_wiggles/',
                            'DNASE.',e,'.fc.signal.bigwig')
  final<-import(con<-gzfile(con_final) , format = "BED")

  mcols(final)<-data.frame(mcols(final),value=1:length(final)) # to indices purposes
  extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")
  dnasepeak <- import(gzfile(con_dnase_peak_test),
                      format = "BED", extraCols = extraCols_narrowPeak)
  aux <- file.path(con_dnase_fc_test);bwf <- BigWigFile(aux)
  dnasefc <- import(bwf)

  over<-findOverlaps(final,dnasepeak)
  final_dnase<-final[unique(over@queryHits)] 
  x<-data.table(qh=over@queryHits,sh=over@subjectHits,as.data.frame(mcols(dnasepeak[over@subjectHits])));x<-x[,-c(3),with=F]
  x2<-x[, lapply(.SD,max), by=qh]            # ==> toooo fast!
  x2m<-x[, lapply(.SD,mean), by=qh];names(x2m)<-paste0(names(x2m),'m')
  x3<-data.table(x2[,c(4,5,6),with=F],x2m[,c(4,5,6),with=F])
  mcols(final_dnase)<-cbind(as.data.frame(mcols(final_dnase)),as.data.frame(x3))
  rm(final,dnasepeak,x,x2,x3,x2m);gc()

  over<-findOverlaps(final_dnase,dnasefc)
  final_dnase<-final_dnase[unique(over@queryHits)] 
  x<-data.table(qh=over@queryHits,sh=over@subjectHits,maxfc=dnasefc[over@subjectHits]$score)
  x1<-x[, lapply(.SD,max), by=qh]  ;x1m<-x[, lapply(.SD,mean), by=qh]
  final_dnase$maxfc<-x1[,maxfc];  final_dnase$maxfcm<-x1m[,maxfc]
  rm(dnasefc,x,x1,over);gc()
  
  
  ###################annotations inclusion########################################
  gtf<-import(gzfile('/home/ricardo/hd/projects/dream_tf_competition/data/annotations/gencode.v19.annotation.gtf.gz')) #2619444
  gtf<-gtf[gtf$gene_type=='protein_coding' & gtf$transcript_status!='PUTATIVE' & 
             gtf$transcript_status!='PUTATIVE']
  pt<-promoters(gtf[gtf$type=='gene'])
  ####dist genes
  fo<-follow(final_dnase, gtf[gtf$type=='gene'] );
  pr<-precede(final_dnase, gtf[gtf$type=='gene']);
  dfp<-data.table(fo=fo,pr=pr);dfp[,id:=1:nrow(dfp)]
  setkey(dfp,fo);dfp[.(NA), fo := 0L]; 
  setkey(dfp,pr); dfp[.(NA), pr := 0L]
  faux<-function(x){ifelse(x==0,0,gtf[gtf$type=='gene'][x]$gene_name)}
  dfp[,fon:=faux(fo)]
  dfp[,prn:=faux(pr)];dfp<-dfp[order(id)]
  load(paste0('/home/ricardo/hd/projects/dream_tf_competition/data/RNAseq/',e,'.RData'))
  xx<-as.data.table(mcols(genes_gtf))
  xx1<-merge(dfp,xx[,c(2:7),with=F],by.x='fon',by.y='gene_name',all.x=T,all.y=F,sort=F, allow.cartesian=T)
  xx1<-xx1[!duplicated(id)];
  xx2<-merge(xx1,xx[,c(2:7),with=F],by.x='prn',by.y='gene_name',all.x=T,all.y=F,sort=F, allow.cartesian=T)
  #xx2<-xx2[!duplicated(id)][,c('TPM.x','TPM.y'),with=F];names(xx2)<-c('follow','precede')
  xx2<-xx2[!duplicated(id)];
  xx2[,paste0(names(xx[,c(3:7),with=F]),'.y') :=
        lapply(.SD,function(x){ifelse(is.na(x),0,x)}),.SDcols=paste0(names(xx[,c(3:7),with=F]),'.y')]

  xx2[,c('fon','prn','pr','fo'):=list(NULL,NULL,NULL,NULL)]
  a<-distanceToNearest(final_dnase, gtf[gtf$type=='gene'] )
  ##########
  fprep<-function(gtf){
    over<-findOverlaps(gtf,final_dnase)
    x<-data.frame(value=mcols(final_dnase[over@subjectHits])$value,type=gtf[over@queryHits]$type,
                  stringsAsFactors = FALSE)#,eid=gtf[over@queryHits]$exon_id)
    x<-setDT(x)
    xc <- dcast(x, value ~type,fun.aggregate = length)
    return(xc)
  }
  z<-fprep(gtf);z[,Selenocysteine:= NULL]
  z1<-fprep(pt);names(z1)<-c('value','promoter2k')
  db<-data.table(as.data.frame(mcols(final_dnase)))
  xdb<-merge(db,z,by='value',all.x=T,all.y=F)
  xdb[, setdiff(names(z),'value') :=lapply(.SD, function(x){ ifelse(is.na(x),0,x)}),.SDcols=setdiff(names(z),'value') ]
  xdb<-merge(xdb,z1,by='value',all.x=T,all.y=F)
  xdb[, c('promoter2k') :=sapply(promoter2k, function(x){ ifelse(is.na(x),0,x)})]
  mcols(final_dnase)<-data.frame(xdb,dist=a@elementMetadata$distance,xx2,stringsAsFactors = F)
  rm(aux,db,x1m,x2,xdb,z,gtf,pt,xx2,z,z1,a,dfp);gc()
  
  #####################end anotations inclusion
 
  save(final_dnase,file=paste0(DIR_TF,'final_dnase','.',tf,'.',e,'.RData'))
  #export.bed(final_dnase,con =file(paste0(DIR_TF,'final_dnase','.',tf,'.',e,'.bed')))
  require(BSgenome.Hsapiens.UCSC.hg19);require(DNAshapeR)
  getFasta(final_dnase,	Hsapiens,	width	= 200,	filename	= paste0(DIR_TF,'final_dnase','.',tf,'.',e,'.fa'))
  return(NULL)
}
load_pwm_yaml<-function(){
  x<-yaml.load_file('/home/ricardo/hd/projects/dream_tf_competition/baseline/DREAM_invivo_tf_binding_prediction_challenge_baseline-master/models.yaml')
  pwm<-vector('list',length(x));nn<-vector(length=length(x))
  for(i in 1:length(x)){
    for(kk in 1:length(x[[i]])){
      if(class(x[[i]][[kk]])=='list'){
        lr<-length(x[[i]][[kk]])
        pwm[[i]]<-matrix(unlist(x[[i]][[kk]]),nrow=lr,byrow=T);nn[i]<-x[[i]]$tf_name
      }
    }
  }
  names(pwm)<-nn
  return(pwm)
}
load_pwmfile<-function(){ 
  path=paste0('/home/ricardo/hd/projects/dream_tf_submission/pwm_encode/')
  ff<-dir(path=path,pattern='motif');
  pwms_ff<-vector('list',length=length(ff))
  for (i in 1:length(ff)){
    a<-read.table(paste0(path,ff[i]));#a<-t(a);rownames(a)<-c('A','C','G','T');  
    pwms_ff[[i]]<-as.matrix(a)
  }
  return(pwms_ff)
}
rnaseqdata<-function(tissue){
  HepG20<-fread(paste0('/home/ricardo/hd/projects/dream_tf_competition/data/RNAseq/gene_expression.',tissue,'.biorep1.tsv'),header=T,data.table=F) #57820     7
  hepg20<-HepG20[,-2]; #57820
  HepG21<-fread(paste0('/home/ricardo/hd/projects/dream_tf_competition/data/RNAseq/gene_expression.',tissue,'.biorep2.tsv'),header=T,data.table=F) #57820     7
  hepg21<-HepG21[,-2]; #57820
  a<-(hepg20[,-1]+hepg21[,-1])/2;a<-data.frame(gene_id=hepg20$gene_id,a)
  gtf<-import(gzfile('/home/ricardo/hd/projects/dream_tf_competition/data/annotations/gencode.v19.annotation.gtf.gz')) #2619444 
  gtf<-gtf[gtf$gene_type=='protein_coding' & gtf$transcript_status!='PUTATIVE' & 
             gtf$transcript_status!='PUTATIVE']
  genes_gtf<-gtf[gtf$type=='gene']
  
  aux<-as.data.frame(mcols(genes_gtf));auxt<-aux[,c('gene_id','gene_name')]
  aux1<-merge(auxt,a,by='gene_id',sort=F)
  
  mcols(genes_gtf)<-aux1
  save(genes_gtf,file=paste0('/home/ricardo/hd/projects/dream_tf_competition/data/RNAseq/',tissue,'.RData'))
  rm(a,gtf,genes_gtf);gc()
  return(NULL)
}
rna<-function(path_trna){
  #path_trna<-'/home/ricardo/hd/projects/dream_tf_competition/data/RNAseq/'
  ff<-list.files(path_trna,pattern = 'biorep1.tsv')
  ff1<-sapply(ff,function(x){strsplit(x,'[.]')[[1]][2]})
  sapply(ff1,rnaseqdata)
}
make_yaml_motif<-function(tf){  # to use scpd2meme
  pwm_yaml<-load_pwm_yaml()
  aa<-(t(pwm_yaml[[tf]]));rownames(aa)<-c('A','C','G','T')
  p<-toPWM(aa)
  t<-data.frame(p$pfm)
  names(t)<-rep(paste0('>',tf),ncol(t))
  write.table(t,file=paste0('/home/ricardo/hd/projects/dream_tf_submission/pwm_baseline/',tf,'yaml','.txt'),quote=F)
}
load_train_and_createfolds<-function(tf,t_train,kk=20,t_test,firstn=F){
  DIR_TF<-paste0('/home/ricardo/hd/projects/dream_tf_competition/data/TF2/',tf,'/')
  setwd(DIR_TF)
  t_train<-setdiff(t_train,'SK-N-SH')
  k<-which(tfs[,1]==tf)
  t_tf<-strsplit(as.character(gsub('\xa0','',tfs[k,4])),',')[[1]]
  
  ##the last dd is the test set.
  set.seed(2423)
  require(ranger)
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

ff<-function(x){(x-min(x))/(max(x)-min(x))}
evalerror <- function(preds, dtrain) {
  require(PRROC)
  labels <- getinfo(dtrain, "label")
  pr<- pr.curve( scores.class0=preds,  weights.class0=labels)
  return(list(metric = "aupr", value = pr$auc.davis.goadrich))
}
froc<-function(gs,pr){
  require(PRROC)
  prec<-ifelse(pr=='B',1,0);gs<-ifelse(gs=='B',1,0)
  roc <- roc.curve( scores.class0=prec,  weights.class0=gs)
  pr<- pr.curve( scores.class0=prec,  weights.class0=gs)
  return(data.frame(roc=roc$auc,pr=pr$auc.davis.goadrich))
}
froc_regression<-function(gs,pr){
  require(PRROC)
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
# loadsyn <- function (filex) {
#   out <- tryCatch( 
#     {
#       eval_leaderboard<-7071644
#       eval_final<-7212779
#       eval<-ifelse(strsplit(filex ,'.',fixed=T)[[1]][1]=='F',eval_final,eval_leaderboard)
#       #id<-strsplit(filex)
#       filename<-filex
#       library(synapseClient)   # Synapse Client, 1.13-4
#       synapseLogin(username = 'maximus', password = 'xxxx')
#       project<-Project(name="ENCODE _TF_submissions")
#       project<-synStore(project)
#       pid<-propertyValue(project, "id")
#       evaluation<-synGetEvaluation(eval) # 7212779
#       #file<-File(path=paste0('/home/ricardo/hd/projects/dream_tf_competition/data/TF1/',tf,'/',filename), parentId=pid)
#       file<-File(path=filex, parentId=pid)
#       entity<-synStore(file )
#       submit(evaluation, entity,'maximus')
#       print(paste0(eval,' ',filex))
#     }, 
#     error = function(cond){
#       message('não pode carregar o arquivo')
#       return(NA)}
#   )
#   return(out)
# }
### this function, preprocess, substitute make_scores*
paths<-function(tf,env = parent.env()){
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
  gencodev19<-file.path(base,'annotations/gencode.v19.annotation.gtf.gz')
  pladder<-file.path(base,'annotations/ladder_regions.blacklistfiltered.bed.gz')
  ptest<-file.path('/home/ricardo/hd/projects/dream_tf_competition/data/annotations/test_regions.blacklistfiltered.bed.gz')
  memeAma<-file.path('~/hd/meme/bin')
  
}
###
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
fannot<-function(chip,dnasepeak,dnasefc){
  over<-findOverlaps(chip,dnasepeak)
  chip_dnase<-chip[unique(over@queryHits)] 
  x<-data.table(qh=over@queryHits,sh=over@subjectHits,as.data.frame(mcols(dnasepeak[over@subjectHits])));x<-x[,-c(3),,with=F]
  x2<-x[, lapply(.SD,max), by=qh]            # ==> toooo fast!
  x2m<-x[, lapply(.SD,mean), by=qh];names(x2m)<-paste0(names(x2m),'m')
  x3<-cbind(subset(x2,select=c(4,5,6)),subset(x2m,select=c(4,5,6)))
  mcols(chip_dnase)<-cbind(mcols(chip_dnase),as.data.frame(x3))
  rm(chip,dnasepeak,x3,x2m,x);gc()
  
  over<-findOverlaps(chip_dnase,dnasefc)
  x<-data.table(qh=over@queryHits,sh=over@subjectHits,maxfc=dnasefc[over@subjectHits]$score)
  x1<-x[, lapply(.SD,max), by=qh]  ;x1m<-x[, lapply(.SD,mean), by=qh]
  chip_dnase$maxfc<-x1[,maxfc];  chip_dnase$maxfcm<-x1m[,maxfc]
  if(is.element(e,t_tf)){
    id<-which(nslots == e)
    #filter only B and Us to speed process
    da<-mcols(chip_dnase)[id][[1]]=='A';chip_dnase<-chip_dnase[!da]
  }
  rm(dnasefc,x,x1,over);gc()
  
  ###################annotations inclusion########################################
  
  ####dist genes
  fo<-follow(chip_dnase, gtf[gtf$type=='gene'] );
  pr<-precede(chip_dnase, gtf[gtf$type=='gene']);
  dfp<-data.table(fo=fo,pr=pr);dfp[,id:=1:nrow(dfp)]
  setkey(dfp,fo);dfp[.(NA), fo := 0L]; 
  setkey(dfp,pr); dfp[.(NA), pr := 0L]
  faux<-function(x){ifelse(x==0,0,gtf[gtf$type=='gene'][x]$gene_name)}
  dfp[,fon:=faux(fo)]
  dfp[,prn:=faux(pr)];dfp<-dfp[order(id)]
  load(paste0(base,'/RNAseq/',e,'.RData'))
  xx<-as.data.table(mcols(genes_gtf))
  xx1<-merge(dfp,xx[,c(2:7),with=F],by.x='fon',by.y='gene_name',all.x=T,all.y=F,sort=F, allow.cartesian=T)
  xx1<-xx1[!duplicated(id)]
  xx2<-merge(xx1,xx[,c(2:7),with=F],by.x='prn',by.y='gene_name',all.x=T,all.y=F,sort=F, allow.cartesian=T)
  xx2<-xx2[!duplicated(id)]
  xx2[,paste0(names(xx[,c(3:7),with=F]),'.y') :=
        lapply(.SD,function(x){ifelse(is.na(x),0,x)}),.SDcols=paste0(names(xx[,c(3:7),with=F]),'.y')]
  xx2[,c('fon','prn','pr','fo'):=list(NULL,NULL,NULL,NULL)]
  a<-distanceToNearest(chip_dnase, gtf[gtf$type=='gene'] )
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
}
#this function to make scores give origin to other 2: 
# across-cell type and within-cell type
################################make_scoressh.R #in the writeup folder
#within-cell type
mscoreswithin<-function(tf){
  for(e in t_test){
    if(length(t_test)!=0){
  pwr<-file.path(paste0(writeup,'/'))
  annot<-file.path(base,'annotations/')
###ladder
  text<-paste0(memeAma,"/ama --o-format gff  --motif ",tf,
               " ",pwr,"pwm_baseline/all.meme ",
               " ",annot,"chipladder_nona.fa ",
               " ",pwr,"hg19markov2.bkg",
               " > meme_ladder_avg1.",e,".txt" )
  fileConn<-file(paste0(tfDir,"/scores_ladder_avg1.",e,".sh"))
  writeLines(text, fileConn)
  close(fileConn)
  text<-paste0(memeAma,"/ama --o-format gff ",  
               "  ",pwr,"pwm_plus/",tf,".txt ",
               "  ",annot,"chipladder_nona.fa ",
               "  ",pwr,"hg19markov2.bkg ",
               "  > meme_ladder_avg2.",e,".txt" )
  fileConn<-file(paste0(tfDir,"/scores_ladder_avg2.",e,".sh"))
  writeLines(text, fileConn)
  close(fileConn)
    }
  }
####label  
 for(e in t_train){
  text<-paste0(memeAma,"/ama --o-format gff  --motif ",tf,
               " ",pwr,"pwm_baseline/all.meme ",
               " ",annot,"chiplabel_nona.fa ",
               " ",pwr,"hg19markov2.bkg ",
               " >  meme_label_avg1.",e,".txt" )
  fileConn<-file(paste0(tfDir,"/scores_label_avg1.",e,".sh"))
  writeLines(text, fileConn)
  close(fileConn)
  text<-paste0(memeAma,"/ama --o-format gff ",  
               " ",pwr,"pwm_plus/",tf,".txt ",
               " ",annot,"chiplabel_nona.fa ",
               " ",pwr,"hg19markov2.bkg ",
               " > meme_label_avg2.",e,".txt " )
  fileConn<-file(paste0(tfDir,"/scores_label_avg2.",e,".sh"))
  writeLines(text, fileConn)
  close(fileConn)
  }
###test
  for(e in t_tf){
  if(length(t_tf)!=0){
    text<-paste0(memeAma,"/ama --o-format gff  --motif ",tf,
                 " ",pwr,"pwm_baseline/all.meme ",
                 " ",annot,"chiptest_nona.fa ",
                 " ",pwr,"hg19markov2.bkg ",
                 " > meme_test_avg1.",e,".txt " )
    fileConn<-file(paste0(tfDir,"/scores_test_avg1.",e,".sh"))
    writeLines(text, fileConn)
    close(fileConn)
    text<-paste0(memeAma,"/ama --o-format gff ",  
                 " ",pwr,"pwm_plus/",tf,".txt ",
                 " ",annot,"chiptest_nona.fa ",
                 " ",pwr,"hg19markov2.bkg ",
                 " > meme_test_avg2.",e,".tx t" )
    fileConn<-file(paste0(tfDir,"/scores_test_avg2.",e,".sh"))
    writeLines(text, fileConn)
    }
  }
}
#across-cell type
mscoresacross<-function(tf){
 # for(tf in tfs[,1]){
    j<-which(tfs[,1]==tf)
    t_train<-strsplit(as.character(tfs[j,2]),',')[[1]];t_train<-sub('[[:space:]]','',t_train)
    t_test<-strsplit(as.character(tfs[j,3]),',')[[1]];t_test<-sub('[[:space:]]','',t_test)  
    #DIR_TF<-paste0('../TF2/',tf,'/')
    if(length(t_test)!=0){
      for(e in t_test){
        text<-paste0(memeAma,"/ama --o-format gff  --motif ",tf,'  ',
                     writeup,"/pwm_baseline/all.meme ",
                     writeup,"/ladder_regions.blacklistfiltered.bed.gz.fasta ",
                     writeup,"/hg19markov.bkg",
                     " > meme_full_avg.",e,".txt" )
        fileConn<-file(paste0(tfDir,"/scores_full_avg.",e,".sh"))
        writeLines(text, fileConn)
        close(fileConn)
        # text<-paste0("~/hd/meme/bin/ama --o-format gff  --scoring max-odds  --motif ",tf,
        #              " ./pwm_baseline/all.meme ",
        #              " ladder_regions.blacklistfiltered.bed.gz.fasta ",
        #              " hg19markov.bkg",
        #              " > meme_full_max.",e,".txt" )
        # fileConn<-file(paste0(DIR_TF,"scores_full_max.",e,".sh"))
        # writeLines(text, fileConn)
        # close(fileConn)
      }
    }
    for(e in c(t_train,t_test)){
      filename<-paste0(tfDir,'/chip_dnase','.',tf,'.',e,'.fa')
      text<-paste0(memeAma,'/ama --o-format gff --motif ',tf,'  ',
                   writeup,"/pwm_baseline/all.meme ",
                   filename,' ',
                   writeup,'/hg19markov.bkg ',
                   ' > meme_avg.',tf,'.',e,'.txt')
      fileConn<-file(paste0(tfDir,"/scores_avg.",e,".sh"))
      writeLines(text, fileConn)
      close(fileConn)
      #   text<-paste0('~/hd/meme/bin/ama --o-format gff --scoring max-odds --motif ',tf,
      #                " ./pwm_baseline/all.meme ",
      #                DIR_TF,filename,' ',
      #              ' hg19markov.bkg ',
      #              ' > meme_max.',tf,'.',e,'.txt')
      # fileConn<-file(paste0(DIR_TF,"scores_max.",e,".sh"))
      # writeLines(text, fileConn)
      # close(fileConn)
    }
 # }
  
#  for(i in 1:nrow(tfs)){  ## for the full 
    t_tf<-strsplit(as.character(gsub('\xa0','',tfs[j,4])),',')[[1]]
    if(length(t_tf)!=0){
      for( e in t_tf ){
        #tf<-tfs[i,1]
        #DIR_TF<-paste0('../TF2/',tf,'/')
        filename<-'test_regions.blacklistfiltered.bed.gz.fasta'
        # text<-paste0('~/hd/meme/bin/ama --o-format gff --scoring max-odds  --motif ',tf,
        #              " ./pwm_baseline/all.meme ",
        #              filename,' ',
        #              ' hg19markov.bkg > ',
        #              DIR_TF,'meme_finalsub_max.',tf,'.',e,'.txt')
        # fileConn<-file(paste0(DIR_TF,"scores_finalsub_max.",tf,'.',e,".sh"))
        # writeLines(text, fileConn)
        # close(fileConn) 
        text<-paste0(memeAma,'/ama --o-format gff  --motif ',tf,'  ',
                     writeup,"/pwm_baseline/all.meme ",
                     filename,' ',
                     writeup,'/hg19markov.bkg > ',
                     tfDir,'/meme_finalsub_avg.',tf,'.',e,'.txt')
        fileConn<-file(paste0(tfDir,"/scores_finalsub_avg.",tf,".",e,".sh"))
        writeLines(text, fileConn)
        close(fileConn) 
      }
    }
 # }
  
  #for(i in 1:nrow(tfs)){  #for the dnase
    t_tf<-strsplit(as.character(gsub('\xa0','',tfs[j,4])),',')[[1]]
    if(length(t_tf)!=0){
      for( e in t_tf ){
       # tf<-tfs[i,1]
        #DIR_TF<-paste0('../TF2/',tf,'/')
        filename<-paste0(tfDir,'/final_dnase','.',tf,'.',e,'.fa')
        # text<-paste0('~/hd/meme/bin/ama --o-format gff --scoring max-odds  --motif ',tf,
        #              " ./pwm_baseline/all.meme ",
        #              filename,' ',
        #              ' hg19markov.bkg > ',
        #              DIR_TF,'meme_finalsub_dnase_max.',tf,'.',e,'.txt')
        # fileConn<-file(paste0(DIR_TF,"scores_finalsub_dnase_max.",tf,'.',e,".sh"))
        # writeLines(text, fileConn)
        # close(fileConn) 
        text<-paste0(memeAma,'/ama --o-format gff  --motif ',tf,'  ',
                     writeup,"/pwm_baseline/all.meme ",
                       filename,' ',
                     writeup,'/hg19markov.bkg > ',
                     tfDir,'/meme_finalsub_dnase_avg.',tf,'.',e,'.txt')
        fileConn<-file(paste0(tfDir,"/scores_finalsub_dnase_avg.",tf,".",e,".sh"))
        writeLines(text, fileConn)
        close(fileConn) 
      }
    }
  #}
}
preprocess_writeup<-function(tf){
  tf<-'E2F6'
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
  tfs<-read.xls(file.path(writeup,"tfs.xls"))
  source(file.path(writeup,"functions.R"))
  j<-which(tfs[,1]==tf)
  t_train<-strsplit(as.character(tfs[j,2]),',')[[1]];t_train<-sub('[[:space:]]','',t_train)
  t_test<-strsplit(as.character(tfs[j,3]),',')[[1]];t_test<-sub('[[:space:]]','',t_test)
  t_tf<-strsplit(as.character(gsub('\xa0','',tfs[j,4])),',')[[1]]
  
  print(t_train)
  print(t_test)
  print(t_tf)
  ################################### library and paths set###############
   #to train
  load_lables_tsv<-function(){
    x<-data.table::fread(paste0('gzip -dc ',con_chipseq_label_tf))
    chip<-makeGRangesFromDataFrame(x,keep.extra.columns = T,starts.in.df.are.0based=T) 
    mcols(chip)<-data.frame(mcols(chip),value=1:length(chip)) # to indices purposes
    return(chip)
  }
  extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",qValue = "numeric", peak = "integer")
  tissue_train<-t_train
  clabels<-load_lables_tsv()
  for(tissue_train in t_train){
    con_dnase_peak_train<-paste0(base,'/essential_training_data/DNASE/peaks/conservative/',
                                 'DNASE.',tissue_train,'.conservative.narrowPeak.gz')
    con_dnase_fc_train<-paste0(base,'/essential_training_data/DNASE/fold_coverage_wiggles/',
                               'DNASE.',tissue_train,'.fc.signal.bigwig')
    aux <- file.path(con_dnase_fc_train);bwf <- BigWigFile(aux)
    dnasefc <- import(bwf)
    dnasepeak <- import(gzfile(con_dnase_peak_train), format = "BED", extraCols = extraCols_narrowPeak)
    over<-findOverlaps(chip,dnasepeak)
    chip_dnase<-chip[unique(over@queryHits)] 
    x<-setDT(data.frame(qh=over@queryHits,sh=over@subjectHits,
              mcols(dnasepeak[over@subjectHits])[-1],# remove 'name'
             bind=factor(unlist(mcols(chip[over@queryHits])[tissue_train])))
             )
    x<-x[bind!='A']    #remove ambiguos
    ##### test variables#########
    # qplot(pValue, signalValue, colour = bind, shape = bind, 
    #       data = x)
    # p <- ggplot(x[,.(pValue,signalValue,bind)], aes(pValue,signalValue))
    #p + geom_boxplot(aes(colour =bind))
    ######test variables########
    x2<-x[, lapply(.SD,max), by=qh]            # ==> toooo fast!
    x2m<-x[, lapply(.SD,mean), by=qh];names(x2m)<-paste0(names(x2m),'m')
    x3<-cbind(subset(x2,select=c(4,5,6)),subset(x2m,select=c(4,5,6)))
    mcols(chip_dnase)<-cbind(mcols(chip_dnase),as.data.frame(x3))
    rm(chip,dnasepeak,x3,x2m,x);gc()
    
    over<-findOverlaps(chip_dnase,dnasefc)
    x<-data.table(qh=over@queryHits,sh=over@subjectHits,maxfc=dnasefc[over@subjectHits]$score)
    x1<-x[, lapply(.SD,max), by=qh]  ;x1m<-x[, lapply(.SD,mean), by=qh]
    chip_dnase$maxfc<-x1[,maxfc];  chip_dnase$maxfcm<-x1m[,maxfc]
    id<-which(nslots == tissue_train)
    #filter only B and Us to speed process
    da<-mcols(chip_dnase)[id][[1]]=='A';chip_dnase<-chip_dnase[!da]
    chipf<- fannot(chip_dnase,dnasepeak,dnasefc)
                     
    # add_dnase_peak<-function(chip,tissue_train){
    # add_dnase_foldcoverage<-function(chip,tissue_train){
  }
  add_annotation<-function(){
  pwm_scores<-function(){
    
  clabels<-load_lables_tsv()
    
      
  merge_all<-function(){
  #leaderboard round
  load_ladder_leaderboard()
  pwm_scores()
  add_dnase_peak()
  add_dnase_foldcoverage()
  add_annotation()
  merge_all()
  #final submission round
  load_test()
  pwm_scores()
  add_dnase_peak()
  add_dnase_foldcoverage()
  add_annotation()
  merge_all()
  
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
