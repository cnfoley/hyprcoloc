# HyPrColoc
Hypothesis Prioritisation in multi-trait Colocalization (HyPrColoc).

## Functions
* hyprcoloc - identifies clusters of colocalized traits and candidate causal SNPs using the HyPrColoc Bayesian divisive clustering algorithm.

## Installation
1. install.packages("devtools")
2. library(devtools)
3. install_github("cnfoley/hyprcoloc", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
4. library(hyprcoloc)
5. browseVignettes("hyprcoloc")

## Example
\# Regression coefficients and standard errors from ten GWAS studies (Traits 1-5, 6-8 & 9-10 colocalize)  
betas <- hyprcoloc::test.betas  
head(betas)  
ses <- hyprcoloc::test.ses  
head(ses)  
  
\# Trait names and SNP IDs  
traits <- paste0("T", 1:10)  
rsid <- rownames(betas)  

\# Colocalization analysis  
hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid)  

## Citations
* HyPrColoc: Foley CN, Staley JR, et al. A fast and efficient colocalization algorithm for identifying shared genetic risk factors across multiple traits. BioRxiv 2019. doi: https://doi.org/10.1101/592238
* HyPrColoc software: Foley, CN and Staley JR. (2020, November 27). cnfoley/hyprcoloc: First release of software (Version v1.0.0). Zenodo. http://doi.org/10.5281/zenodo.4293559
* LD blocks: Berisa T & Pickrell JK. Approximately independent linkage disequilibrium blocks in human populations. Bioinformatics 2016; 32(2):283-285
