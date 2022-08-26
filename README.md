# **P**redictive **H**TS **E**pigenetic **B**iomarkers (PHEB)

Environment
--- 
For running the analysis, you can use [conda](https://docs.conda.io/en/latest/) for installing its requirements. In the desired location, type:
```
git clone https://github.com/aljoshoh/pheb.git
cd pheb/
source environment/make_env.sh
``` 
This will install a virtual conda environment called `pheb_r3.5.3`. It contains an installation of `Rstudio`. Alternatively, one can use the attached Dockerfile to build a virtual Docker environment that has the requirements installed.

Analysis
---
The chunks provided in `workflow.Rmd` run different substeps of the analysis.
