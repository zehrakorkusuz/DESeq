# DESeq

This is the repository of the replication of .... code to discover DESeq library in R.

# Overview

- Reading in table of counts
- Adding annotation
- Filtering lowly expressed genes
- Quality control
- Normalisation for composition bias
- Differential expression analysis
- Testing relative to a threshold
- Visualisation
- Gene set testing


# Datasets
Data files are available from: https://figshare.com/s/1d788fd384d33e913a2a. You should download the files listed below and place them into a folder called data in your working directory.

Data files:

GSE60450_Lactation-GenewiseCounts.txt
SampleInfo.txt
SampleInfo_Corrected.txt
Data files were originally obtained from:
ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60450/suppl/GSE60450_Lactation-GenewiseCounts.txt.gz
http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.c2.all.v7.1.entrez.rds
http://bioinf.wehi.edu.au/MSigDB/v7.1/Mm.h.all.v7.1.entrez.rds

This material has been inspired by the following resources:
http://www.statsci.org/smyth/pubs/QLedgeRPreprint.pdf (Lun, Chen, and Smyth 2016)
http://monashbioinformaticsplatform.github.io/RNAseq-DE-analysis-with-R/RNAseq_DE_analysis_with_R.html


# Resources

https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html#Overview
