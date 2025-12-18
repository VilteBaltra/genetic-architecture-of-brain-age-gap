### Genetic architecture of brain age gap

This repository contains scripts and workflows for estimating brain age gaps (BAG), performing genome-wide association studies (GWASs) in individual cohorts, meta-analysing GWAS summary statistics, running genetic correlations, Mendelian randomisation, and polygenic score analyses.  

---

#### Repository structure:

**0.run-cohort-level-enigma-GWAS**   
   Scripts for deriving Han’s BAG phenotype and running GWAS in individual cohorts.

**1.format-cohort-level-enigma-sumstats**  
   Formatting GWAS summary statistics.

**2.metal-GWAS-meta-analysis-enigma**  
   Scripts to perform METAL-based meta-analysis of ENIGMA GWAS (BAG Han).

**3.genomic-SEM-6-BAGs**  
   Analysis scripts for genomic structural equation modeling of six BAG GWASs.

**4.ldsc-BAGfactor**  
   Scripts for running LDSC regression between BAG factor and 32 other traits.

**5.ldsc-BAGHan**  
   LDSC analyses specifically for Han (ENIGMA) BAG.

**6.Mendelian-randomisation**  
   Scripts for forward and reverse Mendelian randomisation analyses.

**7.polygenic-scores**  
   Scripts to derive polygenic scores in independent cohorts.

**plots**  
   Plotting scripts for visualising GWAS, PGS, and MR results.




# Genetic Architecture of Brain Age Gap

This repository contains scripts and workflows for estimating **brain age gaps (BAG)**, performing **GWAS** in individual cohorts, **meta-analysing GWAS summary statistics**, running **genetic correlations**, **Mendelian randomisation**, and **polygenic score analyses**.  

---

## Overview

![Analysis Flowchart](/Users/vb506/Documents/Github_scripts_brainage_upload/Figure_1_flowchart_map.jpg)  
*Figure: High-level workflow from phenotype derivation to downstream analyses. Light brown boxes represent genome-wide association studies (GWASs). Maroon boxes represent post-GWAS analyses. Summary statistics for BAGHan have been obtained as part of the present study (see map in the top right). Summary statistics for BAGLeonardsen, BAGWen, BAGSmith, BAGJawinski, and BAGKaufmann have been obtained from previously published studies. BAG = brain age gap; PheWAS = phenome-wide association study; UKBB = UK Biobank; GenR = Generation R study.*

---

## Repository Structure

| Folder | Description |
|--------|-------------|
| `0.run-cohort-level-enigma-GWAS` | Scripts for deriving Han’s BAG phenotype and running GWAS in individual cohorts. |
| `1.format-cohort-level-enigma-sumstats` | Formatting and QC of GWAS summary statistics. |
| `2.metal-GWAS-meta-analysis-enigma` | METAL-based meta-analysis of ENIGMA GWAS (BAG Han). |
| `3.genomic-SEM-6-BAGs` | Genomic structural equation modeling of six BAG GWASs. |
| `4.ldsc-BAGfactor` | LDSC regression between BAG factor and 32 other traits. |
| `5.ldsc-BAGHan` | LDSC analyses specifically for Han (ENIGMA) BAG. |
| `6.Mendelian-randomisation` | Forward and reverse Mendelian randomisation analyses. |
| `7.polygenic-scores` | Scripts to derive polygenic scores in independent cohorts. |
| `plots` | Example plotting scripts for GWAS, PGS, and MR results. |

---

## Notes

- Large files (e.g., Singularity containers) are tracked via **Git LFS**.  
- Most scripts require pre-formatted input files; see instructions within each folder.  
- Complementary UK Biobank–specific brain age workflows: [pjawinski/enigma_brainage](https://github.com/pjawinski/enigma_brainage).
