ENCODE-DREAM in vivo Transcription Factor Binding Site Prediction Challenge
================

Instructions
------------

Install meme suite to use **ama** to calculate tf scores along the DNA.

Define a **base** directory and, on top of that, a **writeup** directory.

Copy all the github files in the **base/writeup/** directory.

Copy these data in the **base** directory:

``` bash
training_data.ChIPseq.tar
training_data.DNASE_wo_bams.tar
training_data.RNAseq.tar
training_data.annotations.tar
```

#### Extract all files

Under the ***writeup*** dirctory, at the beginning of ***'functions\_for\_main\_program.R'*** set the variables **tf**, for what transcription factor will be assessed and the **base** directory.

Then source the file 'functions\_for\_main\_program.R'

``` r
source('functions_for_main_program.R')
```

As it is stated in the **main.R** execute the function **run\_execbash\_sh\_on\_tf\_folder\_after\_this(tf)** and, using the terminal(I used ubuntu xenial), go to the folder **'base/writeup/results/tf'** where the **tf** was defined in the beginning.

Execute execbash.sh on the tf directory, to get the file with the TF score using the meme suit

In annotations directory, do: Use bedtools to extract fasta from bed coordinates:

``` bash
 bedtools getfasta -fi hg19.genome.fa -bed test_regions.blacklistfiltered.bed -fo test_regions.blacklistfiltered.fa
```

``` bash
 bedtools getfasta -fi hg19.genome.fa -bed ladder_regions.blacklistfiltered.bed -fo ladder_regions.blacklistfiltered.fa
```

#### Get background file for ama suite in the writeup directory:

``` bash

 fasta-get-markov ../annotations/hg19.genome.fa hg19markov.bkg
```

Running Machine Learning and Genereting File to Submit
------------------------------------------------------

At this point, everything shoud work fine, generating the file to submit in the respective **tf** directory:

``` r
preprocess_writeup(tf)
dd<-load_features(tf)
#machine learnig section
for(i in seq_along(leaderboard)){
  xgscore<-xgbtrain(i)
  rfscore<-rftrain(i)
  fcs(xgscore,rfscore,i)
}
```

Name: maximus

Ricardo Paixao dos Santos

Prof. Katlin Brauer Massirer, Ph.D.

Lab of RNA and microRNA Regulation in Disease,

Center for Molecular Biology and Genetic Engineering

University of Campinas, UNICAMP

Av Candido Rondon, 400 Campinas,13083-875, Brazil
