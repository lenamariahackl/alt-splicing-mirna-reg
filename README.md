# Exon-specific differences in microRNA regulation due to alternative splicing
A Implementation of a pipeline to 
* predict miRNA - mRNA binding sites from sequence using [TarPmiR](http://hulab.ucf.edu/research/projects/miRNA/TarPmiR/)
* predict miRNA - mRNA interaction from [TCGA](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga) expression using Elastic net regression
* combine the two and analyse the influence of alternative splicing on miRNA regulation.

## Setup

For Python 3 but should also work for Python 2 :) To run this project, install it locally using:
```
git clone git@github.com:lenamariahackl/alt-splicing-mirna-reg.git
cd alt-splicing-mirna-reg
pip3 install -r requirements.txt
```
## Overview
![Detailed overview of the implemented steps. While the processing on the left side had to be executed per TCGA cancer type to be analyzed, the right side was executed
only once.](/figures/pipeline.PNG)

The implementation is divided into three separate jupyter notebooks following the three phases of prediction of miRNA binding sites from sequence data using TarPmiR (Phase 1), from expression data using elastic net regression (Phase 2) and the final analysis (Phase 3).

### Phase 1: Prediction of binding sites from sequence using TarPmiR
Human miRNA sequence information is downloaded from [miRBase](http://www.mirbase.org/) (release 22.1) and mRNA transcript sequences, as well as the exon and chromosome position annotation from [Ensembl](https://www.ensembl.org/) release 100 (GRCh38.p13 assembly). TarPmiR predicts binding sites on these sequences with position relative to the transcript, which we translate to be relative to the chromosome using exon start and end positions provided by Ensembl. We also add exon id, chromosome and strand. Then, we drop all binding sites spanning more than one exon and filter for binding sites with binding probability > 0.8. Then we reduce the binding site information to (miRNA,exon)-level binding information.

### Phase 2: Prediction of binding sites from expression using Elastic net regression
We used the [Xena repository](https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) to download the transcript expression (unit log2(tpm + 0.001)) and miRNA mature strand expression (unit log2(norm_value + 1)) data from TCGA Pan-Cancer (PANCAN). Until now cancer types Kidney Chromophobe (KICH) and Lower Grade Glioma data (LGG) were investigated, but it is easily possible to instead analyse any of the [33 TCGA cancer types](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga/studied-cancers) with available data. The exon expression is being filtered for disease and sample type. Then we calculate exon expression from transcript expression and calulate an alternative splicing filter to only keep genes, that show alternative splicing. Finally, we drop all genes with low variance between samples. For miRNA expression we filter for disease and sample type and apply the variance filter. 
To find the correlation between exon and miRNA expression, we use elastic net regression and try to predict miRNA expression from exon expression. Afterwards we filter for good models with low RMSE and, per miRNA, for genes where at least one of the exons has a negative coefficient.
![Example input data format for one miRNA elastic net regression model.](/figures/elastic_net.PNG)

### Phase 3: Finding common binding sites and analysis
We combine filtered (exon, miRNA) binding pairings from both TarPmiR predictions and elastic net regression. With the unbiased elastic net regression we control for the false positives in the TarPmiR predictions. We only kept (miRNA, gene) pairs, where the exon with the most negative coefficient was also predicted in TarPmiR as binding to this miRNA. We again train regression models, this time per (miRNA, gene) pair. We filter for good models and miRNA binding genes. 

## Execution
### Phase 1: Prediction of binding sites from sequence using TarPmiR
Parallel execution of TarPmiR for the 249.740 mRNA sequences and 2.656 human miRNA sequences on a cluster of 80 cores with 1.5 TB RAM takes about 10 days. The output files were read into python and combined into pandas DataFrame with 983.499.270 rows of around 130 GB in memory. The data is too big to be uploaded to git, but can be requested on lenamariahackl@gmail.com.

### Phase 2: Prediction of binding sites from expression using Elastic net regression
In the first notebook cell both path and disease, that we want to analyse the expression of, can be changed. Then the pipeline can be executed step-by-step. The training of the miRNA-level elastic net regression models requires around 6 hours per 100 models depending on the amount of input. Figures are produced and automatically saved in path/disease/'plots'.

### Phase 3: Finding common binding sites and analysis
In the first notebook cell both path and disease, that we want to analyse the expression of, can be changed. The user can decide as well whether miRNA level regression or (miRNA, gene) level regression should be analyzed. Then the pipeline can be executed step-by-step. To plot the percentage of coding vs noncoding, coding information is needed. This information needs to be obtained from [BioMart](https://www.ensembl.org/biomart/). Figures are produced and automatically saved in path/disease/'plots'.

## Project description
Alternative splicing enables the expression of several different transcripts from a single gene. Alternative splicing can produce a messenger RNA transcript unaffected by microRNA regulation by excluding an exon with a microRNA target site. In diseases such as cancer, epigenetic dysregulation including miRNA dysregulation plays a significant role. Epigenetic therapy focuses on treating a disease by targeting the regulatory mechanisms. Better understanding of the epigenetic processes and their changes in cancer is necessary to aid in the development of better treatments. 
While there is previous work focusing on the diversity of 3' UTRs affecting miRNA efficiency and investigating the genome-wide impact of alternative splicing on microRNA binding sites on transcript-level, we study the inﬂuence of alternative splicing on miRNA regulation on exon-level. We predict miRNA binding sites on mRNA from their sequence using TarPmiR. Furthermore, we seek evidence for miRNA binding correlations in exon and miRNA expression from cancerous tissue using Elastic net regression and then filter out false positive TarPmiR binding site predictions by comparison. We test our hypothesis, that alternative splicing has a signiﬁcant impact on miRNA regulation, on the examples of The Cancer Genome Atlas Kidney Chromophobe Carcinoma and Brain Lower Grade Glioma data sets and are able to show an inﬂuence given the filtering steps we applied. Finally, we also show that coding regions encompass more miRNA binding sites than previously thought. Overall we conclude that miRNAs can't sufficiently regulate transcripts, that splice out the binding exon.
