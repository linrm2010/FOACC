
Readme

1. Download data from NCBI
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR825/SRR8257100/SRR8257100.sra
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR825/SRR8257101/SRR8257101.sra

# Download SRA Toolkit (https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software) and use ‘fastq-dump’ to split sra files: 
fastq-dump --split-files SRR8257100.sra
fastq-dump --split-files SRR8257101.sra

2. Running Cell Ranger for analysis
perl change_10x_chromium_lib_sequenced_data_format.pl -read1fq  SRR8257100_1.fastq.gz  -read2fq  SRR8257100_2.fastq.gz  -output_prefix  YSZT_1_L101
perl change_10x_chromium_lib_sequenced_data_format.pl -read1fq  SRR8257101_1.fastq.gz  -read2fq  SRR8257101_2.fastq.gz  -output_prefix  YSZT_1_L102

# running Cell Ranger v2.1.1 on May 13, 2019
count --id=RWR1 --fastqs=/data/Arabidopsis_root_Ryu/WT_rep1/10x_format --transcriptome=/data/install/cell_ranger/cellranger-2.1.1/reference/Arabidopsis_thaliana_TAIR10_update/TAIR10 --localcores=20 --localmem=60 --jobmode=local --expect-cells=6000 --uiport=3600 --disable-ui

# obtain three files:
matrix.mtx.gz
genes.tsv.gz
barcodes.tsv.gz

# cd step01
# in step01 direction
ln -s ../matrix.mtx.gz matrix.mtx.gz
ln -s ../genes.tsv.gz features.tsv.gz
ln -s ../barcodes.tsv.gz barcodes.tsv.gz

The matrix.mtx.gz is split into four files of 'matrix.mtx.gz.001', 'matrix.mtx.gz.002', 'matrix.mtx.gz.003' and 'matrix.mtx.gz.004' due to its bigger than 25 Mb for upload on Github.

3. Two R scripts ('02.AthRoot_step01.R' and '02.AthRoot_step02.R') will be used for analysis.
