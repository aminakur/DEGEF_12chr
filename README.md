# IN THIS FORK
1. DEGEF.py was modified for n=12 genomes;
2. Fixed the warning triggered by df_es_c being a view not a copy of df_es;
3. Input file requirements included.
   
# DEGEF
Differentially Expressed Gene Expression Finder (DEGEF) discovers genomic regions of enriched differential gene expression from differential gene expression results.

# Installation
We recommend running DEGEF using Python-3.7 or greater. Required python packages include argparse, numpy, pandas, random, scipy, statsmodels, and matplotlib.

To install DEGEF, run the following code:

```
git clone --recursive https://github.com/gaog94/DEGEF
```

Then add the pulled directory to PYTHONPATH.

```
export PYTHONPATH="/path/to/DEGEF/directory/:$PYTHONPATH"
```


DEGEF is able to run using differential gene expression results in  table form as input. The differential expression table should include the following columns: xyz, xyz, xyz.

To generate these results, run a differential gene expression finder like DESeq2 or limma-voom.

featureCounts is available as part of the Rsubread/Subread package, which may be obtained at https://bioconductor.org/packages/release/bioc/html/Rsubread.html.

limma-voom may be obtained at https://bioconductor.org/packages/release/bioc/html/limma.html.

# Description
DEGEF (rhymes with "VEGF") stands for Differentially Expressed Gene Expression Finder and is an algorithm that searches for genomic areas of enriched differential gene expression. Using differential gene expression results from common differential expression analysis pipelines such as limma-voom or DESeq2, DEGEF scans the genome searching for regions with greater differential expression than expected by chance.

Details of implementation are included in our manuscript "Co-localization of Clusters of TCR-Regulated Differentially Expressed Genes with Areas of TAD Rearrangement" (in preparation).


# Usage

### Input files generation - DEGEF_input_files.ipynb

### Create and activate virtual environment

```
cd DEGEF/
python3 -m venv venv
pip install argparse numpy pandas random scipy statsmodels matplotlib
source /DEGEF/venv/bin/activate
```

### Running DEGEF from the Command Line - recommended
After generating a differential gene expression results table using differential expression analysis software (e.g., limma-voom, DESeq2, etc.), run the following:

```
python DEGEF.py -i filename -o outDir -l refLoci -c refSizes -m metric -d dir -w ws -s ss -p pval -f lfc -x max -n nBoot -q fdr
```
The following input parameters are used:
* ``-i`` ``--inputDEGFile`` ``filename`` is the input tsv file outputted by a differential expression analysis pipeline.
* ``-o`` ``--outputDirectory`` ``outDir`` is the output directory where DEGEF will write its output files.
* ``-l`` ``--refTranscriptLoci`` ``refLoci`` is the path to the reference file containing transcripts' genomic loci. (Default: "../ReferenceFiles/hg38_ensemble_loci_modified.tsv")
* ``-c`` ``--refChromSizes`` ``refSizes`` is the path to the reference file containing chromosome sizes. (Default: "../ReferenceFiles/chrom_hg38.sizes")
* ``-m`` ``--scoringMetric`` ``metric`` is the scoring metric ("count" vs "significance" vs "foldchange") to use for calculating raw scores. (Default: "foldchange")
* ``-d`` ``--scoringDirection`` ``dir`` is the directionality of expression ("up" vs "dn" vs "mx") that DEGEF will test for enrichment in. (Default: "up")
* ``-w`` ``--windowSize`` ``ws`` is the size of the genomic window used to compute enrichment scores. (Default: 1e6)
* ``-s`` ``--stepSize`` ``ss`` is the number of basepairs the genomic window shifts by when computing enrichment scores. (Default: 2e5)
* ``-p`` ``--limmaPValueThreshold`` ``pval`` is the adjusted p-value threshold used to classify genes as up/down/not regulated when computing raw scores. (Default: 0.05)
* ``-f`` ``--limmaLogFoldChangeThreshold`` ``lfc`` is the log-fold-change threshold used to classify genes as up/down/not regulated when computing raw scores. (Default: 0.58)
* ``-x`` ``--maxGeneInWindow`` ``max`` is the maximum number of genes within a window before DEGEF will switch from bootstrapping to using the Central Limit Theorem (CLT) to compute the null distribution of enrichment scores. (Default: 30)
* ``-n`` ``--nBootstraps`` ``nBoot`` is an integer representing the number of bootstraps to conduct when building null distributions for the enrichment score. (Default: 6400)
* ``-q`` ``--fdrPeakThreshold`` ``fdr`` represents the False Discovery Rate (FDR) threshold to use to define DEGEF enrichment peaks. (Default: 0.05)


<!--
### Running Differential Gene Expression Analysis and DEGEF from the Command Line
```
./DEGEF_wrapper.sh -i filename -o tag -l refLoci 
```

# Outputs
DEGEF generates a number of output files in the user-specified output directory:
* something.tsv
* somethingelse.txt
* parameters.txt
* Per Chromosome Plots
-->

# Attributions
KIRCLE is developed and maintained by Galen Gao -- Laboratory of Molecular Immunology -- National Heart Lung and Blood Institute (NHLBI), National Institutes of Health (NIH), Bethesda, MD.
