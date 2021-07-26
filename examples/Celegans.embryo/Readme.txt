
Readme

1. Obtain the data by running an R script:
library("devtools")
library("Seurat")
library("SeuratData")

devtools::install_github('satijalab/seurat-data')

InstallData("celegans.embryo")
# http://seurat.nygenome.org/src/contrib/celegans.embryo.SeuratData_0.1.0.tar.gz

write.table(GetAssayData(celegans.embryo), "Celegans.embryo.exp.byGene.txt", sep="\t", quote = FALSE)

2. Next, running Bash and Perl script:
# contents in Celegans.embryo.exp.byGene.txt
# AAACCTGCAAGACGTG.300.1.1	AAACCTGGTGTGAATA.300.1.1	...
# WBGene00010957	5	10	...

# vi Celegans.embryo.exp.byGene.txt, by adding 'CellID'
# CellID	AAACCTGCAAGACGTG.300.1.1	AAACCTGGTGTGAATA.300.1.1	...
# WBGene00010957	5	10	...

perl transform_matrix.pl -matrix Celegans.embryo.exp.byGene.txt -output Celegans.embryo.exp.byCell.txt

less -S Celegans.embryo.exp.byGene.txt|awk '{print $1}'|sort -u | perl -ne 'if(!($_=~/CellID/)){print $_}'> gene.id

perl ../generate_matrix.pl -reformat_matrix_byCellID Celegans.embryo.exp.byCell.txt -features  gene.id -output_prefix Celegans.embryo.data

gzip Celegans.embryo.data_features.tsv
gzip Celegans.embryo.data_matrix.mtx
gzip Celegans.embryo.data_barcodes.tsv
gzip Celegans.embryo.exp.byGene.txt
gzip Celegans.embryo.exp.byCell.txt

mkdir step01
cd step01
ln -s ../Celegans.embryo.data_features.tsv.gz features.tsv.gz
ln -s ../Celegans.embryo.data_matrix.mtx.gz matrix.mtx.gz
ln -s ../Celegans.embryo.data_barcodes.tsv.gz barcodes.tsv.gz

3. Using the three files for following analysis:
Celegans.embryo.data_features.tsv.gz
Celegans.embryo.data_matrix.mtx.gz
Celegans.embryo.data_barcodes.tsv.gz

Performing analysis using R, which were shown in two detailed files (02.Cel_step01.R and 02.Cel_step02.R).
