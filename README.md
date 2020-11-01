# alt-splicing-mirna-reg
Exon-specific differences in microRNA regulation due to alternative splicing

## Setup
Python 3.8.5
pyarrow 
cython 0.29.21

Also install liftOver (I used kent source version 405 from 13-Oct-2020) (using wget http://hgdownload.soe.ucsc.edu/admin/exe/). The LiftOver program requires a UCSC-generated over.chain file as input (download using wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz).