
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
source('functions_for_main_program.R')
##################################### run main process##################
run_execbash_sh_on_tf_folder_after_this(tf)
# execute execbash.sh on the tf directory!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
preprocess_writeup(tf)
####################################run machine learning process
dd<-load_features(tf)
# to verify if number of lines i dat.train is ok --------------------------
# dat.train<-f_combinew(dd[train],1)
# sum(sapply(dd[train],function(x)dim(x)[[1]])/20)
# end ---------------------------------------------------------------------

#create file to submite leaderboard and test if exist --------------------------------------------------------
for(e in c(leaderboard,test)){
  xgscore<-xgbtrain(e)
  rfscore<-rftrain(e)
  fcs(xgscore,rfscore,e)
}
#end  
# 
# load(file=file.path(tfDir,paste0('rfscore_',leaderboard[i],'.RData')))
# load(file=file.path(tfDir,paste0('xgscore_',leaderboard[i],'.RData')))


# end ---------------------------------------------------------------------


