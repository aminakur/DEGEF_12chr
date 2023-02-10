#!/bin/bash

### Bash script for running RNAseq pipeline and DEGEF
### Run as batch job: sbatch --gres=lscratch:800 --mem=256g --cpus-per-task=28 --time=48:00:00 wrapper.sh

# load in modules we need to run the pipeline on biowulf
# sratoolkit contains fasterq-dump for SRA downloads
# fastqc and cutadapt are read quality control tools
# hicexplorer contains the tools required for processing HiC data
# subread contains featureCounts function for RNAseq analysis

# module load sratoolkit samtools python R # general tools
# module load fastqc cutadapt	# QC tools

# make sure you have hisat2 & subread tools installed & loaded
# module load hisat2 subread # RNAseq specific tools ()


### Set algorithm runtime parameters ###

# logistics define working variables
workingDir=/lscratch/$SLURM_JOBID/
tag=Yang24hHiC_redo1/
outDir=/data/gaogf/Houston/$tag/

# grab factor data format file of RNAseq bam files to use
factorData=~/Houston/InputFiles/Yang4hRNA.txt

# pass in reference genome
refGenome=/fdb/hisat/hg38/genome

# pass in featureCounts gene transcripts annotation file
gtfFile=/fdb/rsem/Homo_sapiens.GRCh38.86.gtf # hg38

# pass in number of cores available for wrapper to use
nCores=$SLURM_CPUS_PER_TASK



# build output directory structure
mkdir -p $outDir/RnaSeq/QC/

echo " --->>>                                    <<<--- "
echo " --->>> Running RNAseq Processing Pipeline <<<--- "
echo " --->>>                                    <<<--- "

cd $outDir/RnaSeq/


echo "Counting read pairs aligned to each transcript..."
# gtfFile=/fdb/rsem/Homo_sapiens.GRCh38.86.gtf # hg38
gtfFile=/fdb/rsem/Mus_musculus.GRCm38.82.gtf # mm10
echo $gtfFile
while read s; do
	outName=${s}_pairCounts.txt
	trimmedOut=${s}_pairCountsTrimmed.txt
	featureCounts -T $nCores -a $gtfFile -F GTF -t exon -g gene_id -t 'exon' -g 'gene_id' 
		-s 0 -Q 0 -p --countReadPairs --minOverlap 1 --fracOverlap 0 --fracOverlapFeature 0
		-C -o $outName ${s}.bam
	cut -f 1,7 $outName > $trimmedOut
done <$RnaSeqAccList

mkdir -p QC/featureCounts
cp *.txt.summary QC/featureCounts


echo "Concatting all columns of readpair counts together into 1 file..."
python ~/Houston/join_columns.py _pairCountsTrimmed.txt pairCounts_concat.txt


echo "Performing differential gene expression analysis using limma-vroom..."
# create annotation data file (file with ENTREZID, gene symbol, & gene name/description)
Rscript  ~/Houston/annotate_IDs.R -i pairCounts_concat.txt -o annodata.txt -h

mkdir -p limmaTreat_output/
Rscript ~/Houston/limma_voom.R -R limmaTreat_report.html -o limmaTreat_output -m pairCounts_concat.txt -f $factorData -a annodata.txt -D 'activated-resting' -c 0.5 -s 2 -P d,c,b,x,m,h,s -l 0.58 -p 0.01 -d 'BH' -G 10 -T -n 'TMM' -b
cp limma_voom/limma-voom_activated-resting.tsv ../ # or whatever its name is


# Run DEGEF
python DEGEF.py [runtime variables]