## Plotting scripts

This folder contains scripts used to generate key figures for the manuscript. Scripts are written in **R**, **Python**, and **shell**, and generally assume that all upstream analyses (GWAS, MR, PGS, genetic correlations, etc.) have already been run and that cleaned result tables are available.

---

### Contents

#### Main manuscript figures

* **rg_among_six_BAGs_Figure-2A.py**: Plots genetic correlations among the six brain age gap (BAG) phenotypes (Figure 2A)

* **manhattan-plot-Figure-2C.R**
  Generates the Manhattan plot for the brain age GWAS used in the main manuscript (Figure 2C)

* **rg_BAGfactor_with_33traits_Figure-2D.py**
  Plots genetic correlations between the BAG factor and 33 external traits (Figure 2D)

* **MR-combined-33traits-plots-main-manuscript-Figure-3.py**
  Generates the combined forward and reverse MR plots across 33 traits for the main manuscript (Figure 3)

* **PGS-and-trio-plots-Figure-4AtoC_Figure-5.py**
  Generates polygenic score (PGS) and trio-based plots for Figures 4A–4C and Figure 5

* **PGS_associations_BAG_models_UKBB_Figure-4E.py**
  Plots associations between PGSs and BAG phenotypes in UK Biobank (Figure 4E).

---

#### Supplementary figures

* **rg_BAGHan_with_33traits_Figure-S2.sh**
  Script to generate genetic correlation plots the BAGHan and 33 traits (Supplementary Figure S2)

* **plot-forward-MR-inputGWASs-Figure-S7.R**
  Plots forward MR results for BAG factor and each of the six input GWASs (Supplementary Figure S7)

* **plot-reverse-MR-inputGWASs-Figure-S8.R**
  Plots reverse MR results for BAG factor and each of the six input GWASs (Supplementary Figure S8)

* **MRlap-IVW-df-plotting-forward-and-reverse-Figure-S9.py**
  Plots forward and reverse MR results using MRlap estimates (Supplementary Figure S9).

