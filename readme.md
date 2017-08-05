ENCODE-DREAM in vivo Transcription Factor Binding Site Prediction Challenge
================

Instructions
------------

**copy all data in the base directory.**

training\_data.ChIPseq.tar training\_data.DNASE\_wo\_bams.tar training\_data.RNAseq.tar training\_data.annotations.tar Extract all files In annotations directory, do: Use bedtools to extract fasta from bed coordinates: bedtools getfasta -fi hg19.genome.fa -bed test\_regions.blacklistfiltered.bed -fo test\_regions.blacklistfiltered.fa bedtools getfasta -fi hg19.genome.fa -bed ladder\_regions.blacklistfiltered.bed -fo ladder\_regions.blacklistfiltered.fa \* background file for ama suite in the writeup directory: fasta-get-markov ../annotations/hg19.genome.fa hg19markov.bkg

Including Code
--------------

You can include R code in the document as follows:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

Including Plots
---------------

You can also embed plots, for example:

![](readme_files/figure-markdown_github-ascii_identifiers/pressure-1.png)

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
