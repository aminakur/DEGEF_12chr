#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 12:47:47 2022

@author: galen

Version 0.2 uses the entire genome, rather than just the current chromosome, in
order to characterize the null distribution of raw scores. This decreases the
effect of chromosomal differences in raw score distribution on significance 
calculations. (i.e. In DEGEF v0.1, if chr21 had a different proportion of 
DEGs & raw scores compared to the entire genome, this would have skewed p-value
estimates and thus potentially misidentified significant peaks).
"""

import sys
import os
import tempfile
import shutil
import argparse

import numpy as np
import pandas as pd

from random import choices, seed
from scipy.stats import norm
from scipy.signal import find_peaks
from statsmodels.stats.multitest import multipletests

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

def load_data(fnDEG, fnRef, entrez=False):
    """Given filenames of DEG file (limma output), all genes in reference, and
    size of chromosome sizes, returns a dataframes of differentially expressed
    genes, reference."""
    # read in differentially expressed genes
    df_deg = pd.read_table(fnDEG, index_col=0)
    if entrez: # only focus on genes w/an entrezID (20232 -> 15556)
        df_deg = df_deg[~df_deg["ENTREZID"].isnull()]
    
    # read in gene locus mapping
    df_ref = pd.read_table(fnRef, index_col=0)
    df_ref = df_ref.join(df_deg[["SYMBOL", "logFC", "P.Value", "adj.P.Val"]], \
                         how="inner")

    return df_ref

def compute_raw_score(df_de, direction="up", metric="count", p=0.05, logFC=0.58):
    """Given a matrix of differentially expressed genes df_de, computes a set of
    raw scores (R-scores) for each gene. Outputs pandas series of R-scores.
    
    Metric may be one of:
        count: each gene"s R-score is the indicator function such that R(g)=1
               if g is a DEG according to p & logFC and R(g)=0 otherwise.
               final score.
        significance: each gene"s R-score is -log10(p-value) if the gene is 
                      upregulated and 0 otherwise.
        reads: each gene"s R-score is the difference in number of reads between
               the 2 groups.
    """
    
    # make sure df_de is sorted the way we want: 1st by chr and then by TSS
    df_de = df_de.sort_values(["Chr", "TSS"])
    
    # degs will be used if metric == "count"; dirGs will be used if 
    # metric == "significance" or "foldchange"
    if direction == "up":
        degs = df_de[(df_de["adj.P.Val"]<=p) & (df_de["logFC"]>=logFC)].index
        dirGs = df_de[df_de["logFC"] > 0].index
    elif direction == "dn":
        degs = df_de[(df_de["adj.P.Val"]<=p) & (df_de["logFC"]<=-logFC)].index
        dirGs = df_de[df_de["logFC"] < 0].index
    elif direction == "mx":
        degs = df_de[(df_de["adj.P.Val"]<=p) & (np.abs(df_de["logFC"])>=logFC)].index
        dirGs = df_de[np.abs(df_de["logFC"]) > 0].index

    newIdx = [r["Chr"] + ":" + str(int(r["TSS"])) for i,r in df_de.iterrows()]

    # depending on the metric being used, compute raw score from differential
    # expression results
    if metric == "count":
        # scoring function is 1 if gene is up/dn; 0 otherwise
        outputs = [1 if x in list(degs) else 0 for x in df_de.index]
        df_out = pd.DataFrame(outputs, index=newIdx, columns=["RS"])
        
    elif metric == "significance":
        # raw score is -log10(p) if logFC matches up/dn; 0 otherwise
        df_zeros = df_de[~df_de.index.isin(dirGs)]
        df_positives = df_de[df_de.index.isin(dirGs)]
        df_zeros["RS"] = 0
        df_positives["RS"] = -np.log10(df_positives["adj.P.Val"])
        
        # re-sort, smart-index, & strip-nonRS-info from output dataframe
        df_out = pd.concat([df_zeros, df_positives]).sort_values(["Chr","TSS"])
        df_out.index = newIdx
        df_out = df_out[["RS"]]
        
    elif metric == "foldchange":
        # raw score is log2(foldchange) if logFC matches up/dn; 0 otherwise
        df_zeros = df_de[~df_de.index.isin(dirGs)]
        df_positives = df_de[df_de.index.isin(dirGs)]
        df_zeros["RS"] = 0
        df_positives["RS"] = np.abs(df_positives["logFC"])
        
        # re-sort, smart-index, & strip-nonRS-info from output dataframe
        df_out = pd.concat([df_zeros, df_positives]).sort_values(["Chr","TSS"])
        df_out.index = newIdx
        df_out = df_out[["RS"]]
        
    else:
        raise ValueError("Metric must be 'count', 'significance', or 'foldchange'")
    
    return df_out

def compute_null_distributions(scores, metric="count", maxWindow=10, nboot=6400):
    """Given a set of scores, compute null distributions. Output as dictionary
    of arrays where key-value pairs are [# genes in window]: [array of scores].
    Default 6400 bootstraps in each null distribution.
    
    Use nullDistributions[0] to store mu and sigma of empirical data (for use
    w/central limit theorem for n_genes larger than maxWindow).
    Alternatively, if metric = "counts" then return a single-element dictionary
    where nullDistributions[0] = probability that any gene is up/downregulated"""
    seed(42)
    
    def bootstrap(scores, n, nboot=1000):
        """Bootstrap sum of n scores from 'scores'"""
        bootstraps = []
        for i in range(nboot):
            ES = np.sum(choices(np.array(scores), k=n))
            bootstraps.append(ES)
        return bootstraps

    if metric == "count":
        # consider each gene"s DE status as a Bernouli trial with each window
        # being a series of such trials. Thus we return an empirical estimate
        # of p, the probability with which any given gene is up/downregulated.
        nullDistributions = {0: np.average(scores)}
    else:
        # use bootstrapping/CLT-combo to compute null distributions
        # nullDistributions[0] is reserved for storing mu & sigma for CLT
        nullDistributions = {0: np.array([np.mean(scores), np.std(scores)])}
        windowPops = np.array(range(maxWindow))+1
    
        for p in windowPops:
            nullDist = bootstrap(scores, p, nboot=nboot)
            nullDistributions[p] = np.array(nullDist)

    return nullDistributions

def binomialApproximation_significance(score, n, p):
    """Given a set of scores, approximate the null distribution as a gaussian
    distribution over the binomial distribution. """
    mu = n * p
    sigma = np.sqrt(n * p * (1-p))
    dist = norm(mu, sigma)
    return dist.sf(score)

def significance(score, n, nulls):
    """Given a null distribution nullDistr (numpy array of score values) and a
    score, returns the fraction of scores in the null distribution that are
    greater than or equal to the given score."""
    if n < len(nulls):
        nullDistr = nulls[n]
        p = ((nullDistr >= score).sum()+0.1) / len(nullDistr)
    else:
        # use central limit theorem to estimate null distribution:
        # mu is empirical average; sigma is empirical std
        mu = n * nulls[0][0]
        sigma = np.sqrt(n) * nulls[0][1]
        normDist = norm(mu, sigma)
        p = normDist.sf(score)

    return p

def compute_enrichment_score(df_rs, chrSize, nullDists, windowSize=1e6, stepSize=5e4, metric="count"):
    """Given a dataframe of enrichment scores, a windowsize and a stepsize,
    computes a series of enrichment scores for all windows in this chromosome."""
    # use hyperparameters to define arrays of start/end points
    starts = np.array(range(0, int(chrSize - windowSize), int(stepSize)))
    ends = starts + windowSize
    
    # compute scores, nGenesInWindow, and significances
    scores, genesInWindow, pvals = [], [], []
    for start, end in zip(starts, ends):
        df_window = df_rs[(df_rs.index > start) & (df_rs.index < end)]
        
        # compute outputs
        n = len(df_window)
        score = np.sum(df_window["RS"])
        if metric == "count":
            pval = binomialApproximation_significance(score, n, nullDists[0])
        elif metric in ["significance", "foldchange"]:
            pval = significance(score, n, nullDists)

        # handle empty windows specially; store all values in lists
        if n > 0:
            scores.append(score)
            pvals.append(pval)
        else:
            scores.append(0)
            pvals.append(1)
        
        genesInWindow.append(n)
        
    # put output together
    df_out = pd.DataFrame([starts, ends, scores, genesInWindow, pvals], \
                          index=["start", "end", "ES", "n", "p"]).T
    return df_out

def peak_region_finder(df_es, fdrThresh=0.05, windowSize=1e6):
    """Given dataframe of enrichment scores & significances, run a peak-finder
    that identifies all peaks in the enrichment score. A peak is defined as a
    contiguous series of windows that surpass fdrThresh."""
    
    peakStarts = np.array(df_es[df_es["FDR"] <= fdrThresh]["start"])
    
    if len(peakStarts) > 0:
        # pick out all start spots that are within one window-length of each other
        # as these peaks will contain the same genes in them in the output file
        consecutives = np.split(peakStarts, \
                                np.where(np.diff(peakStarts) >= windowSize)[0]+1)
    
        # store peaks as list of tuples
        bases, summitFDRs, summitPs = [], [], []
        for peak in consecutives:
            lpeak, rpeak = peak[0], peak[-1] + windowSize
    
            df_x = df_es[(df_es["start"]>=lpeak) & (df_es["end"]<=rpeak)]
            summitFDR, summitP = np.min(df_x["FDR"]), np.min(df_x["p"])
            bases.append([lpeak, rpeak])
            summitFDRs.append(summitFDR)
            summitPs.append(summitP)
        
        # handle outputting max (min) significance/FDR levels of peak regions
        regionNames = [str(int(z[0])) + "-" + str(int(z[1])) for z in bases]
        peakSigs = pd.DataFrame([summitPs, summitFDRs], columns=regionNames, \
                                index=["peakP", "peakFDR"]).T

    else:
        bases = []
        peakSigs = pd.DataFrame([], columns=["peakP", "peakFDR"])
    
    return bases, peakSigs

def get_genes_in_zones(df_ref, peakZones, c, direction="up", p=0.05, logFC=0.58):
    """Given dataframe of DE analysis output and list of peak zones, returns a
    dictionary of all genes within each peak zone. Differentially expressed
    genes (defined by direction, p, and logFC) are marked with an asterisk."""
    peakGs = {}
    for z in peakZones:
        # define df_z as all dataframe entries with TSS within current peakZone
        df_z = df_ref[df_ref["Chr"] == c]
        df_z = df_z[(df_z["TSS"] > z[0]) & (df_z["TSS"] < z[1])]
        
        # df_dz is our subset of df_z with all up/down genes (per p & logFC)
        if direction == "up":
            df_dz = df_z[(df_z["adj.P.Val"]<=p) & (df_z["logFC"]>=logFC)]
        elif direction == "dn":
            df_dz = df_z[(df_z["adj.P.Val"]<=p) & (df_z["logFC"]<=-logFC)]
        elif direction == "mx":
            df_dz = df_z[(df_z["adj.P.Val"]<=p) & (np.abs(df_z["logFC"])<=logFC)]

        # pull names of up/dn genes (if SYMBOL unavailable, use ENSG name)
        dgenes = list(df_dz[~df_dz["SYMBOL"].isnull()]["SYMBOL"])
        densgs = list(df_dz[df_dz["SYMBOL"].isnull()].index)

        # pull names of non-up/dn genes (if SYMBOL unavailable, use ENSG name)        
        ogenes = list(df_z[(~df_z["SYMBOL"].isnull()) & \
                           (~df_z["SYMBOL"].isin(dgenes))]["SYMBOL"])
        oensgs = list(df_z[(df_z["SYMBOL"].isnull()) & \
                           (~df_z.index.isin(densgs))].index)
        
        # append result for this region/peakZone to output dictionary
        region = c + ":" + str(int(z[0])) + "-" + str(int(z[1]))
        peakGs[region] = [x+"*" for x in dgenes] + [x+"*" for x in densgs] + \
                            ogenes + oensgs

    return peakGs
    
def write_peakGenes_to_file(peakGenes, peakSigs, fname):
    """Given dictionary of genes within each peak and the significance of each
    peak in peakSigs, write these results to a file named fname."""
    # pre-compute maximum number of genes in a single peak
    peakGeneNums = [len(x) for x in peakGenes.values()]
    if len(peakGeneNums) > 0:
        maxGeneLen = np.max(peakGeneNums)
    else:
        maxGeneLen = 0

    with open(fname, "w") as f:
        # write names of regions in first row of file
        for region in peakGenes.keys():
            f.write(region + "\t")
        f.write("\n")
        
        # write significances associated with each peak in second row of file
        for region in peakGenes.keys():
            fdr = peakSigs.loc[region, "peakFDR"]
            f.write(str(fdr) + "\t")
        f.write("\n")
        
        # write genes in the remaining lines of the file
        for i in range(maxGeneLen):
            for gs in peakGenes.values():
                if i < len(gs):
                    f.write(gs[i] + "\t")
                else:
                    f.write("\t")
            f.write("\n")

    print("Output written to file.")
    
def write_params_to_file(args, outDirName):
    argnames = ["DEG filename", "Reference filename", "ChrSize filename", \
                "direction", "metric", "p", "logFC", "maxWindow", "windowSize",\
                "stepSize", "fdrThresh"]
    df_p = pd.DataFrame(args, index=argnames)
    df_p.to_csv(outDirName+"/runtime_arguments.txt", sep="\t", header=None)    

def visualize_results(df_ref, df_rs, df_es, df_pk, c, chrSize, direction="up", \
                      p=0.05, logFC=0.58, fdrThresh=0.05, outname="plot.png"):
    """Given the output dataframes from running DEGEF, visualize the results
    on a single chromosome c."""
    # pick out parts of input files relevant to the chromosome c    
    df_in = df_ref[df_ref["Chr"] == c]
    df_rs_c = df_rs[df_rs.index.str.split(":").str[0] == c]
    df_es_c = df_es[df_es.index == c].copy()
    df_pk_c = df_pk[df_pk.index.str.split(":").str[0] == c]
    
    # plotting parameters
    fig, axs = plt.subplots(nrows=5, ncols=1, sharex=True, figsize=(10,10))
    scale = 1e6
    
    # set up color/label scheme and define the DEGs
    if direction  == "up":
        label = "Upregulated Gene"
        lineColor = "#ca0020"
        peakColor = "#f4a582"
        df_de = df_in[(df_in["adj.P.Val"] <= p) & (df_in["logFC"] >= logFC)]
    elif direction == "dn":
        label = "Downregulated Gene"
        lineColor = "#0571b0"
        peakColor = "#92c5de"
        df_de = df_in[(df_in["adj.P.Val"] <= p) & (df_in["logFC"] <= -logFC)]
    elif direction == "mx":
        label = "Differentially Expressed Gene"
        lineColor = "#7b3294"
        peakColor = "#c2a5cf"
        df_de = df_in[(df_in["adj.P.Val"] <= p) & (np.abs(df_in["logFC"]) >= logFC)]

    
    # plot DEG TSS locations (and background genes that are not DE)
    df_bg = df_in[~df_in.index.isin(df_de.index)]
    axs[0].plot(df_de["TSS"]/scale, np.ones(len(df_de)), linestyle="", \
                marker="o", alpha=0.8, color=lineColor, label=label)
    axs[0].plot(df_bg["TSS"]/scale, np.zeros(len(df_bg)), linestyle="", \
                marker="o", alpha=0.8, color="#bababa", label="Non-"+label)    
    axs[0].legend(loc="center right")
    axs[0].set_ylim((-0.33, 1.33))
    axs[0].set_ylabel("DEGs")
    
    # Aesthetics: Clean up y-axis information in DEG plot (top/1st plot)
    axs[0].set_yticks([])
    axs[0].spines.left.set_visible(False)

    # plot raw scores
    df_rs_c.index = df_rs_c.index.str.split(":").str[-1].astype(int)
    axs[1].plot(df_rs_c.index/scale, df_rs_c["RS"], color=lineColor, alpha=1)
    axs[1].set_ylabel("Raw Score", fontsize=14)

    # plot enrichment scores
    df_es_c.loc[:,"mid"] = (df_es_c["start"] + df_es_c["end"]) / 2
    axs[2].plot(df_es_c["mid"]/scale, df_es_c["ES"], color=lineColor)
    axs[2].set_ylabel("Enrichment Score", fontsize=14)
    
    # plot uncorrected enrichment significances (p)
    axs[3].plot(df_es_c["mid"]/scale, -np.log10(df_es_c["p"]), color=lineColor)
    axs[3].set_ylabel(r"$-\log_{10}(p)$", fontsize=14)
    
    # plot corrected enrichment FDR
    logFDRThresh = -np.log10(fdrThresh)
    axs[4].plot(df_es_c["mid"]/scale, -np.log10(df_es_c["FDR"]), color=lineColor)
    axs[4].plot([0, chrSize/scale], [logFDRThresh, logFDRThresh], "k--")
    axs[4].set_ylabel(r"$-\log_{10}(FDR)$", fontsize=14)
    
    # plot peaks that surpass FDR threshold
    ymin, ymax = axs[4].get_ylim()
    for i in df_pk_c.index:
        l = int(i.split(":")[-1].split("-")[0])/scale
        r = int(i.split("-")[-1])/scale
        pts = np.array([[l, ymin], [l, ymax*1.2], [r, ymax*1.2], [r, ymin]])
        poly = Polygon(pts, closed=False, fc=peakColor, alpha=0.75)
        axs[4].add_patch(poly)
    axs[4].set_ylim((ymin, ymax*1.2))

    # Additional plot label functions (set xlabel & title)
    axs[-1].set_xlabel("Genomic Position (Mb)", fontsize=15)
    axs[0].set_title("Chromosome " + c, fontsize=18)
    
    # Aesthetics: despine right & top; move left & bottom out by 10 pts
    for ax in axs:
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        ax.spines.left.set_position(("outward", 10))
        ax.spines.bottom.set_position(("outward", 10))
    
    # save figure and then clear
    plt.savefig(outname, bbox_inches="tight")
    plt.close()


def run_enrichment_finder(fnDEG, fnRef, fnSize, direction="up", metric="count",\
                          p=0.05, logFC=0.58, maxWindow=10, nboot=6400,\
                          windowSize=1e6, stepSize=5e5, fdrThresh=0.05,\
                          writeToDisk=False, outDirName="xyz"):
    """Main function. Run enrichment finder across the whole genome using the
    specified input parameters.
        - fnDEG: output file from limma-voom with TREAT
        - fnRef: bedfile of genes under consideration in current genome build
        - fnSize: file of chromosome sizes
        - direction ("up"/"down"): look for clusters of up/downregulated genes
        - metric ("count"/"significance"/"foldchange"): raw-score metric
        - p: p-value threshold used to compute raw-score
        - logFC: log-fold-change threshold used to compute raw-score
        - maxWindow: max number of genes within window before switching from 
                    using bootstrapping to CLT to compute null distribution
        - nboot: number of bootstraps in each bootstrapped null distribution
        - windowSize: size (in bps) of sliding genomic window for computing ES
        - stepSize: sliding genomic window moves stepSize bps when computing ES
        - fdrThresh: FDR threshold below which enrichment peaks are significant
        - writeToDisk (bool): whether or not to write to outputs to file
        - outDirName: name of directory in which to save DEGEF outputs
    """
    
    # load reference data
    print(" --->>> Loading Reference Data <<<--- ")
    df_ref = load_data(fnDEG, fnRef)
    
    df_sizes = pd.read_table(fnSize, index_col=0, header=None)
    df_sizes.columns = ["length"]
    df_sizes.index = df_sizes.index.str[3:]

    
    # compute raw scores
    print(" --->>> Computing Raw Scores For Each Gene <<<--- ")
    df_rs = compute_raw_score(df_ref, direction=direction, metric=metric, \
                              p=p, logFC=logFC)


    # build null distributions of raw scores
    # considering all genome-wide R-scores (not just those on the local chr)
    print(" --->>> Building Null Distributions of Raw Scores <<<--- ")
    nullDists = compute_null_distributions(df_rs["RS"], metric=metric, \
                                           maxWindow=maxWindow, nboot=nboot)


    # compute enrichment scores, enrichment significances
    print(" --->>> Computing Enrichment Scores & Significances <<<--- ")
    esTot = []
    for c in mainChrs:
        print("Analyzing chromosome: " + c)
        chrSize = df_sizes.loc[c, "length"]
        
        # select only raw scores within current chromosome being analyzed
        df_rs_c = df_rs[df_rs.index.str.split(":").str[0] == c]
        df_rs_c.index = df_rs_c.index.str.split(":").str[-1].astype(int)
        
        # compute enrichment peak regions within current chromosome
        df_es = compute_enrichment_score(df_rs_c, chrSize, nullDists, \
                                         stepSize=stepSize, windowSize=windowSize, \
                                         metric=metric)
        df_es["chr"] = c
        df_es.set_index("chr", inplace=True)
        esTot.append(df_es)

    df_es_out = pd.concat(esTot)

    
    # do multiple hypothesis testing correction on enrichment significances
    print(" --->>> Performing Multiple Hypothesis Testing Correction <<<--- ")
    reject, fdrs, _, _ = multipletests(df_es_out["p"], method="fdr_bh")
    df_es_out["FDR"] = fdrs


    # Collect genes within each peak
    peakGenesTot, peakSigsTot = {}, []
    for c in mainChrs:
        df_es = df_es_out[df_es_out.index == c]
        peakZones, peakSigs = peak_region_finder(df_es, fdrThresh=fdrThresh)
        peakSigs.index = [c+":"+x for x in peakSigs.index]
        
        # get peak genes
        peakGenes = get_genes_in_zones(df_ref, peakZones, c, direction=direction)
        
        # store results
        peakGenesTot.update(peakGenes)
        peakSigsTot.append(peakSigs)

    peakSigsTot = pd.concat(peakSigsTot)

    
    # write out results (i.e. genes that are contained within peak regions)
    if writeToDisk:
        print(" --->>> Writing DEGEF Results To File <<<--- ")
    
        #if (os.path.exists(outDirName)):
        #    # if output directory exists already, then moves it to a temporary
        #    # directory before deleting it.
        #    tmp = tempfile.mktemp(dir=os.path.dirname(outDirName))
        #    shutil.move(outDirName, tmp)
        #    shutil.rmtree(tmp)

        # make output directory outDirName if it does not exist yet
        # then fill it up with output files        
        if not (os.path.exists(outDirName)):
            print("Creating output directory: " + outDirName)
            os.mkdir(outDirName)

        print("Writing parameters to "+outDirName+"/runtime_arguments.txt")
        write_params_to_file([fnDEG, fnRef, fnSize, direction, metric, p, \
                              logFC, maxWindow, windowSize, stepSize, fdrThresh],\
                             outDirName)

        print("Writing peak genes to "+outDirName+"/peaks_output.txt")
        write_peakGenes_to_file(peakGenesTot, peakSigsTot, outDirName+"/peaks_output.txt")
        
        print("Writing enrichment scores to "+outDirName+"/EnrichmentScores.txt")
        df_es_out.to_csv(outDirName + "/EnrichmentScores.txt", sep="\t")
        
        print("Writing peak FDRs to "+outDirName+"/PeakSignificances.txt")
        peakSigsTot.to_csv(outDirName + "/PeakSignificances.txt", sep="\t")
    
        print("Generating visual plots at each chromosome")
        for c in mainChrs:
            print("\t Plotting Chromosome " + c + " ... ")
            chrSize = df_sizes.loc[c, "length"]
            outPlotName = outDirName + "/chr"+c+"_"+direction+".png"
            visualize_results(df_ref, df_rs, df_es_out, peakSigsTot, c, chrSize, \
                              direction=direction, p=p, logFC=logFC, \
                              fdrThresh=fdrThresh, outname=outPlotName)
            
    print("DEGEF finished running: {0} peaks identified".format(len(peakSigsTot)))

    return peakGenesTot, peakSigsTot


def visualize_enrichment_score(df_ref, c, direction="up", metric="count", \
                               windowSize=1e6, stepSize=5e5, fdrThresh=0.05):
    """Given dataframe of dataframe of DE result on a particular chromosome, 
    compute enrichment score, summary score, and plot it out on axes."""
    # Perform the DEGEF computations for this chromosome:
    chrSize=df_sizes.loc[c, "length"]

    # compute raw scores    
    df_rs = compute_raw_score(df_ref, direction=direction, metric=metric)
    df_rs_c = df_rs[df_rs.index.str.split(":").str[0] == c]
    df_rs_c.index = df_rs_c.index.str.split(":").str[-1].astype(int)
    
    nullDists = compute_null_distributions(df_rs["RS"], metric=metric, \
                                           maxWindow=10, nboot=6400)

    # compute the enrichment-scores
    df_es = compute_enrichment_score(df_rs_c, chrSize, nullDists, windowSize, \
                                     stepSize, metric=metric)
    df_es.loc[:,"mid"] = (df_es["start"] + df_es["end"]) / 2
    df_es["chr"] = c
    df_es.set_index("chr", inplace=True)
    
    # multiple hypothesis testing correction
    reject, fdrs, _, _ = multipletests(df_es["p"], method="fdr_bh")
    df_es["FDR"] = fdrs

    # plot peaks that surpass FDR threshold
    peakZones, peakSigs = peak_region_finder(df_es, fdrThresh=fdrThresh)
    peakSigs.index = [c+":"+x for x in peakSigs.index]

    visualize_results(df_ref, df_rs, df_es, peakSigs, c, chrSize, \
                      direction=direction, fdrThresh=fdrThresh, )
    
# global variables
mainChrs = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]

if __name__ == "__main__":
    
    # parse and echo user-specified runtime parameters
    purpose = "Search for local genomic neighborhoods of differential gene \
                expression enrichment of differential gene expression"
    parser = argparse.ArgumentParser(description=purpose)

    # parse input/output file locations
    parser.add_argument("-o", "--outputDirectory", type=str, required=True, \
        help="Path to output directory where DEGEF results will be written to")
    parser.add_argument("-i", "--inputDEGfile", type=str, required=True, \
        help="Path to limma-voom with TREAT output table containing \
        significances, FDRs, and foldchanges for each transcript")
    parser.add_argument("-l", "--refTranscriptLoci", type=str, required=True, \
        default="./ref/hg38_ensemble_loci_modified.tsv", \
        help="Reference file reporting the loci of all transcripts in the genome")
    parser.add_argument("-c", "--refChromSizes", type=str, required=True, \
        default="./ref/chrom_hg38.sizes", \
        help="Reference file reporting the sizes of each chromosome in the genome")
    
    # parse actual running parameters
    parser.add_argument("-m", "--scoringMetric", type=str, \
        choices=["count", "significance", "foldchange"], default="foldchange", \
        help="Scoring metric for raw scores. Count/Significance/FoldChange")
    parser.add_argument("-d", "--scoringDirection", type=str, \
        choices=["up", "dn", "mx"], default="up", \
        help="Directionality of expression that is enriched. Up/Down/Mixed.")
    parser.add_argument("-w", "--windowSize", type=int, default=1e6, \
        help="Size of genomic window used to compute enrichment score")
    parser.add_argument("-s", "--stepSize", type=int, default=2e5, \
        help="Size the genomic window shifts by when computing enrichment scores")
    parser.add_argument("-p", "--limmaPValueThreshold", type=float, default=0.05, \
        help="Adjusted p-value threshold to classify genes as up/down/not \
        regulated.")
    parser.add_argument("-f", "--limmaLogFoldChangeThreshold", type=float, \
        default=0.58, help="Log fold-change threshold to classify genes as \
        up/down/not regulated.")
    parser.add_argument("-x", "--maxGeneInWindow", type=int, default=30, \
        help="Max # genes within a window before computing null distribution \
        of ES's using CLT rather than bootstrapping ")
    parser.add_argument("-n", "--nBootstraps", type=int, default=6400, \
        help="Number of bootstraps to use when empirically estimating ")
    parser.add_argument("-q", "--fdrPeakThreshold", type=float, default=0.05, \
        help="False Discovery Rate (FDR) threshold defining DEGEF enrichment peaks")

    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_usage(sys.stderr)
        sys.exit(1)

    # read in file paths of input and reference files
    # also read in desired output directory path
    outDirName = args.outputDirectory
    fnDEG = args.inputDEGfile
    fnRef = args.refTranscriptLoci
    fnSize = args.refChromSizes

    # set hyperparameters -- default = 100kb windows w/10kb steps
    # for testing purposes (fast) = 1Mb windows w/500kb steps
    metric = args.scoringMetric
    direction = args.scoringDirection
    windowSize = args.windowSize
    stepSize = args.stepSize
    p = args.limmaPValueThreshold
    logFC = args.limmaLogFoldChangeThreshold
    maxGeneInWindow  = args.maxGeneInWindow
    nboot = args.nBootstraps
    sigThresh = args.fdrPeakThreshold
    
    # read in chromosome sizes from fnSize
    df_sizes = pd.read_table(fnSize, index_col=0, header=None)
    df_sizes.columns = ["length"]
    if df_sizes.index.str[3:][0] == 'chr':
        df_sizes.index = df_sizes.index.str[3:]
    
    # Run DEGEF
    peakGenesTot, peakSigsTot = run_enrichment_finder(fnDEG, fnRef, fnSize, 
                                    direction=direction, metric=metric, p=p,\
                                    logFC=logFC, maxWindow=maxGeneInWindow, \
                                    nboot=nboot,windowSize=windowSize, \
                                    stepSize=stepSize, writeToDisk=True, \
                                    outDirName=outDirName)


    ### Miscellaneous code used for testing/debugging; now commented out
    # df_ref = load_data(fnDEG, fnRef)
    # df_rs = compute_enrichment_score(df_ref, direction="up", metric="count")
    # c = "5"
    # visualize_enrichment_score(df_ref, c, direction="up", metric=metric, windowSize=windowSize, stepSize=stepSize)
