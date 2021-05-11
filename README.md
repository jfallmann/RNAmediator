# RIssmed  [![GitHub](https://img.shields.io/github/tag/jfallmann/RIssmed.svg)](https://github.com/jfallmann/RIssmed) [![Build Status](https://github.com/jfallmann/RIssmed/.github/workflows/python-app.yml/badge.svg?branch=main)](https://github.com/jfallmann/RIssmed/actions) [![Bioconda](https://anaconda.org/bioconda/RIssmed/badges/version.svg)](https://anaconda.org/bioconda/RIssmed)

The RNA Interaction via secondary structure mediation (RIssmed) tool suite
analyses the change of RNA secondary structure upon binding other molecules.

It is available as suite of commandline tools, the source code of RIssmed is open source and available via GitHub (License GPL-3).

### <u>Installation via bioconda - recommended</u>

RIssmed can be installed with all tool dependencies via [conda](https://conda.io/docs/install/quick.html). Once you have conda installed simply type:

    conda create -n RIssmed -c conda-forge -c bioconda RIssmed

Activate the environment in which RIssmed was installed to use it:

    conda activate RIssmed


### <u>Usage</u>
    
RIssmed provides different tools to answer different research questions. We will discuss the target application of the tools
and the required input in the following points. 

* CalcTempDiffs.py

    Analyses how the positions-wise probability of base-pairing changes for an RNA sequence upon temperature change. 
    
    * Input: Sequence string
    
    * Output: Plot of relative position-wise probability change
    

* CollectConsResults.py
 
    Computes statistics for change of position-wise probability of being base-paired, with binding constraint for a set of provided genes.
    
    * Input:
    
    * Output:


* CollectWindowResults.py 

    Computes statistics for change of position-wise probability of being base-paired, using windows, with binding constraint for a set of provided genes.
    
    * Input:
    
    * Output:



* ConstraintFold.py
    
    Computes RNA secondary structure prediction with defined hard binding constraint for a set of sequences.
    
    * Input:
    
    * Output:


* ConstraintPLFold.py

    Computes position-wise probability of being base-paired, using windows, with hard binding constraint for a set of sequences.

    * Input:
    
    * Output:


* FoldWindows.py 

    Calculate base pairing probs of given seqs or random seqs for given window size, span and region
    
    * Input:
    
    * Output:


* risvis.py

    Visualize base pairing probabilties before and after constraining.
    
    * Input:
    
    * Output:


    
