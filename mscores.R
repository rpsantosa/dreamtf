
library(gdata)
tfs<-read.xls("/home/ricardo/hd/projects/dream_tf_competition/data/tfs.xls") 
motifs_scores<-function(tf){
  j<-which(tfs[,1]==tf)
  t_train<-strsplit(as.character(tfs[j,2]),',')[[1]];t_train<-sub('[[:space:]]','',t_train)
  t_test<-strsplit(as.character(tfs[j,3]),',')[[1]];t_test<-sub('[[:space:]]','',t_test)
  t_tf<-strsplit(as.character(gsub('\xa0','',tfs[j,4])),',')[[1]]
#for(tf in tfs[,1]){
  j<-which(tfs[,1]==tf)
  t_train<-strsplit(as.character(tfs[j,2]),',')[[1]];t_train<-sub('[[:space:]]','',t_train)
  t_test<-strsplit(as.character(tfs[j,3]),',')[[1]];t_test<-sub('[[:space:]]','',t_test)  
  DIR_TF<-paste0('../TF2/',tf,'/')
  if(length(t_test)!=0){
    for(e in t_test){
      text<-paste0("~/hd/meme/bin/ama --o-format gff ", 
                   " ./pwm_plus/",tf,".txt",
                   " ladder_regions.blacklistfiltered.bed.gz.fasta ",
                   " hg19markov.bkg",
                   " > meme_full_plus.",e,".txt" )
      fileConn<-file(paste0(DIR_TF,"scores_full_plus.",e,".sh"))
      writeLines(text, fileConn)
      close(fileConn)
      
    }
  }
  for(e in c(t_train,t_test)){
    filename<-paste0('chip_dnase','.',tf,'.',e,'.fa')
    text<-paste0('~/hd/meme/bin/ama --o-format gff ',
                 " ./pwm_plus/",tf,".txt ",
                 DIR_TF,filename,' ',
                 ' hg19markov.bkg ',
                 ' > meme_plus.',tf,'.',e,'.txt')
    fileConn<-file(paste0(DIR_TF,"scores_plus.",e,".sh"))
    writeLines(text, fileConn)
    close(fileConn)
  }
#}

#for(i in 1:nrow(tfs)){  ## for the full 
  t_tf<-strsplit(as.character(gsub('\xa0','',tfs[i,4])),',')[[1]]
  if(length(t_tf)!=0){
    for( e in t_tf ){
      tf<-tfs[i,1]
      DIR_TF<-paste0('../TF2/',tf,'/')
      filename<-'test_regions.blacklistfiltered.bed.gz.fasta'
      text<-paste0('~/hd/meme/bin/ama --o-format gff ' ,
                   "./pwm_plus/",tf,".txt ",
                   filename,' ',
                   ' hg19markov.bkg > ',
                   DIR_TF,'meme_finalsub_plus.',tf,'.',e,'.txt')
      fileConn<-file(paste0(DIR_TF,"scores_finalsub_plus.",tf,'.',e,".sh"))
      writeLines(text, fileConn)
      close(fileConn) 
      
    }
  }
#}

#for(i in 1:nrow(tfs)){  
  t_tf<-strsplit(as.character(gsub('\xa0','',tfs[i,4])),',')[[1]]
  if(length(t_tf)!=0){
    for( e in t_tf ){
      tf<-tfs[i,1]
      DIR_TF<-paste0('../TF2/',tf,'/')
      filename<-paste0(DIR_TF,'final_dnase','.',tf,'.',e,'.fa')
      text<-paste0('~/hd/meme/bin/ama --o-format gff ',
                   " ./pwm_plus/",tf,".txt ",
                   filename,' ',
                   ' hg19markov.bkg > ',
                   DIR_TF,'meme_finalsub_dnase_plus.',tf,'.',e,'.txt')
      fileConn<-file(paste0(DIR_TF,"scores_finalsub_dnase_plus.",tf,'.',e,".sh"))
      writeLines(text, fileConn)
      close(fileConn) 
    }
  }
#}
}