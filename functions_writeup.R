
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
dfolds<-lapply(train,createF,dd=dd,kk=20);names(dfolds)<-train
dat.train<-f_combinew(dd[trainx],1,1:1e5,dfolds,trainx)

