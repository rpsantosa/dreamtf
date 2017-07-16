
# ##########################################################################
# make_dnase
# source make_scoressh.R
# bash execbash.sh 
# make_scores
#-analysis
#############################################################################
# make the fasta file ladder_regions.blacklistfiltered.bed.gz.fasta  #2GB
# make the fasta file test_regions.blacklistfiltered.bed.gz.fasta     #13.7GB
# save in the writeup directory

# download the tar files in https://www.synapse.org/#!Synapse:syn6131484/wiki/402043
# in a directory named data. Extract all there.
# change the paths in functions.R accordingly, the functions dnase_train,dnase_test,dnase_final

# set the paths appropriately before run each loop for.
# here i create the TF2 folder in the same level than writeup folder 
# all the functions are in the writeup folder.

# 1- set the meme path in make_scoressh.R and run in the writeup directory
# After the first loop, the meme scores have to be calculated
# based on fasta files of each TF. To do that, the file make_scoressh.R
# creates a .sh linux script in each folder (each TF) to do that. 
# It calculates the average and maximum motifs scores values


# 2- run execbash.sh (modify if necessary. It entries in each folder and exec
# all .sh files). It calcules the meme motifs scores for each TF for each tissue.

# 3-the second and third loops prepare the data for the analysis,
# the trainning, leaderboard and final submission sets

# After the data are ready, the script F.round.R gives the final submittion file
# in the TF folder, typing: 
  
  
# Rscript  F.round.R   'ATF2'




library(GenomicRanges)
library(ShortRead)  #clean
library(rtracklayer)
library(gdata)
library(Biostrings)
library(data.table)
library(caret)
library(reshape2)
tf<-'ARID3A'
paths(tf)
setwd('/home/ricardo/hd/projects/dream_tf_competition/data/writeup')
path_to_tfsname<-"tfs.xls"
con_results<-'/home/ricardo/hd/projects/dream_tf_competition/data/writeup/results/'
tfs<-read.xls(path_to_tfsname) 
#take 12 hours
for(tf in (tfs[1,1])){
  rm(list=setdiff(ls(),c('tf','tfs')));gc()
  source('/home/ricardo/hd/projects/dream_tf_competition/data/writeup/functions.R')
  options( warn = -1 )
  # con_tf2<-'/home/ricardo/hd/projects/dream_tf_competition/data/TF2/'
  #con_tf2<-'/home/ricardo/hd/projects/dream_tf_competition/data/writeup/results/'
  DIR_TF<-paste0(con_results,tf,'/')
  dir.create(DIR_TF)
  setwd(DIR_TF)
  j<-which(tfs[,1]==tf)
  t_train<-strsplit(as.character(tfs[j,2]),',')[[1]];t_train<-sub('[[:space:]]','',t_train)
  t_test<-strsplit(as.character(tfs[j,3]),',')[[1]];t_test<-sub('[[:space:]]','',t_test)
  print(t_train)
  print(t_test)
  make_dnase_over_train_and_test(tf,t_train,t_test)
  t_tf<-strsplit(as.character(gsub('\xa0','',tfs[j,4])),',')[[1]]
  if(length(t_tf)!=0){
    for( e in t_tf ){
      rm(list=setdiff(ls(),c('tf','tfs','j','e','t_tf')));gc()
      source('../../writeup/functions.R')
      tf<-tfs[j,1]
      dnase_final(tf,e)
    }
  }
}
#this take 8 hours (the meme scores using ama from meme suit)

## now run make_scoressh.R to create the .sh scripts that will use the pwm motifs to 
## to create scores. After, write bash execbash.sh from the writeup directory to
## activate the scripts written by make_scoressh.R

#this is fast ( 20 min) 
#Goes printing on output the TF and tissues used
for(tf in tfs[,1]){
  rm(list=setdiff(ls(),c('tf','tfs')));gc()
  source('/home/ricardo/hd/projects/dream_tf_competition/data/writeup/functions.R')
  DIR_TF<-paste0('/home/ricardo/hd/projects/dream_tf_competition/data/TF2/',tf,'/')
  setwd(DIR_TF)
  j<-which(tfs[,1]==tf)
  t_train<-strsplit(as.character(tfs[j,2]),',')[[1]];t_train<-sub('[[:space:]]','',t_train)
  t_test<-strsplit(as.character(tfs[j,3]),',')[[1]];t_test<-sub('[[:space:]]','',t_test)
  print(tf)
  print(t_train)
  print(t_test)
  make_motifs_scores_train_and_test(tf,t_train,t_test)
}
## here the data for the final submission.
## Goes printing on output the TF and tissues used
## Take some time because the files are big (~3GB)
## ~2 hours
for(i in 1:nrow(tfs)){
  t_tf<-strsplit(as.character(gsub('\xa0','',tfs[i,4])),',')[[1]]
  if(length(t_tf)!=0){
    for( e in t_tf ){
      rm(list=setdiff(ls(),c('tf','tfs','i','e','t_tf')));gc()
      source('../../writeup/functions.R')
      tf<-tfs[i,1]
      DIR_TF<-paste0('/home/ricardo/hd/projects/dream_tf_competition/data/TF2/',tf,'/')
      setwd(DIR_TF)
      make_scores_final(tf,e)
      make_scores_final_dnase(tf,e)
      print(as.character(tf))
      print(as.character(e))
    }
  }
}

