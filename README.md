# **P**redictive **H**TS **E**pigenetic **B**iomarkers (PHEB)

Environment
--- 
For running the analysis, you can use [conda](https://docs.conda.io/en/latest/) for installing its requirements. In the desired location, type:
```
git clone https://github.com/aljoshoh/pheb.git
cd pheb/
source environment/make_env.sh
``` 
This will install a virtual conda environment called `pheb_r3.5.3`. It contains an installation of `Rstudio`.

Analysis
---

```
paths.R
```
Set path for relevant datasets for analysis

```
drug_preprocessing.R:
```
Imports drug HTS cell viability values from GDSC1/2 and groups into cancer types in metadata/GDSC_create_drugs_AUC_2/response*

```
meth_preprocessing.R:
```
Imports raw GDSC methylation data, preprocesses and groups into cancer types in metadata/GDSC/methyl_preprocessed*

```
cpg_analyze_dmp.R:
```
Perform drug DMP analysis on the GDSC data running for each cancer type separately


```
cpg_analyze_compp.R:
```
Run comb-p algorithm for calling drug DMRs on spatially correlated p-values from DMP analysis


```
cpg_analyze_tcga.R:
```
Download, process and store TCGA methylation and gene expression data

```
elmer.R:
```
Run correlation analysis between drug DMRs and proximal gene expression, for each GDSC and TCGA datasets

```
CTRP_validation.R:
```
Compare CTRP drug HTS to drug DMRs

```
CCLE_validation.R:
```
Process CCLE RRBS data and compare with drug DMRs

```
pool_results.R:
```
Combine metadata for the results object, finding overlaps between GDSC and TCGA for patient drug DMRs

```
hits_genomics.R:
```
Screen for associations between patient drug DMRs and somatic mutations curated by the GDSC

