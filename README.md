# Pol1_termination_MS
Supplemental code for Pol1 termination MS. This repository provides step-by-step guidance to reproduce data in the MS.
It is divided into two stages (cluster and desktop) according to computational demand. 
All steps require git and conda installed

# Cluster processing
## Getting started

Clone this repository
```
git clone git@github.com:tturowski/Pol1_termination_MS.git
```

Create and activate conda environment (you can use mamba instead)
```
conda env create -f envs/snakemake.yml
conda activate snakemake
```
Run SnakeMake file to process raw files (-c determines number of CPUs to use)
```
snakemake -c64 --use-conda -s SM_CRACprocessing3end_all.py
```
**NOTE:** Very first run of the SnakeMake file ```SM_CRACprocessing3end_all.py``` will initialize new conda environments. This may take a several minutes.
